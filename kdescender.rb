#! /usr/bin/ruby
# todo: use better velvet/oases params

require 'rubygems'
require 'trollop'
require 'bio'

# options
opts = Trollop::options do
  version "v0.0.1a"
  banner <<-EOS

--------------------- k-descender ---------------------

run assembler (by default, velvet-oases) with descending
k-mers in a range (default, 90:25, step 6) for paired-end
reads

Richard Smith, April 2013

-------------------------------------------------------

EOS
  opt :left, "left reads file in fastq format", :required => true, :type => String
  opt :right, "right reads file in fastq format", :required => true, :type => String
  opt :insertsize, "mean insert size",  :required => true, :type => Integer
  opt :insertsd, "insert size standard deviation", :required => true, :type => Integer
  opt :threads, "number of threads for bowtie", :required => true, :type => Integer

end

trap("SIGTERM") {
  puts "process killed!"
}

trap("SIGINT") {
  puts "process interrupted!"
}

trap("SIGQUIT") {
  puts "process quitting!"
}

puts "running velvet/oases k-descent from 90 to 25"

# set velvet thread limit
File.open("setmpthreads.sh", "w") do |f|
  f.puts "export OMP_NUM_THREADS=#{opts.threads}"
  f.puts "export OMP_THREAD_LIMIT=#{opts.threads}"
end

`sh setmpthreads.sh`

l = opts.left
r = opts.right

first = true
t0 = Time.now # time the whole run
toclean = nil # don't clean up during the first round
(25..90).step(6).reverse_each do |k|
  t1 = Time.now # time this assembly
  puts "running assembly with k = #{k}"
  longstr = first ? "" : "-fasta -long longreads.fa "
  `~/apps/velveth k#{k} #{k} -fastq -shortPaired #{l} #{r} #{longstr}> k#{k}.log`
  `~/apps/velvetg k#{k} -read_trkg yes -cov_cutoff 3 -scaffolding yes -min_contig_lgth 200 -ins_length #{opts.insertsize} -ins_length_sd #{opts.insertsd} -min_pair_count 6 >> k#{k}.log`
  `~/apps/oases k#{k} -ins_length #{opts.insertsize} -ins_length_sd #{opts.insertsd} -scaffolding yes -cov_cutoff 3 -min_trans_lgth 200 -edgeFractionCutoff 0.1 -min_pair_count 6 >> k#{k}.log`
  # skip to next k if we didn't get any contigs
  if File.size("k#{k}/transcripts.fa") == 0
    puts "no contigs assembled, descending..."
    next
  end
  # print assembly stats
  puts "assembly report for k#{k}:"
  puts `leaff --stats k#{k}/transcripts.fa`
  # filter out reads mapping to the current assembly
  puts "building bowtie index for k#{k} assembly"
  `/usr/local/bin/bowtie2-build k#{k}/transcripts.fa k#{k} >> k#{k}.log`
  puts "filtering out reads mapping to k#{k} assembly"
  `/usr/local/bin/bowtie2 --very-fast -p #{opts.threads} --un-conc k#{k} --no-discordant --phred64 -x k#{k} -1 #{l} -2 #{r} > /dev/null`
  puts "storing unmapped reads for next step of descent"
  l = "k#{k}.1"
  r = "k#{k}.2"
  # clean up reads from last assembly
  if toclean
    `rm #{toclean}*`
  end
  # prepare the current assembly to be used as long reads in the next
  puts "triplicating assembly to meet coverage requirements"
  File.open("longreads.fa", "w") do |lr|
    Bio::FastaFormat.open("k#{k}/transcripts.fa").each do |record|
      (1..3).each do |n|
        lr.puts record.naseq.to_fasta(record.entry_id, 60)
      end
    end
  end
  # clean up unneeded files generated during oases assembly
  puts "cleaning up assembly files"
  Dir.chdir("k#{k}") do
    `rm Sequences Roadmaps PreGraph Graph2 LastGraph`
  end
  toclean = "k#{k}"
  delta = Time.now - t1
  first = false
  puts "assembly with k = #{k} completed in #{delta} seconds"
end

delta = Time.now - t0
puts "all velvet/oases runs completed in #{delta} seconds"
