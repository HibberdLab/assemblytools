#! /usr/bin/ruby
# todo: use better velvet/oases params

require 'rubygems'
require 'trollop'

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

l = opts.left
r = opts.right

t0 = Time.now
(25..90).step(6).reverse_each do |k|
  t1 = Time.now
  puts "running assembly with k = #{k}"
  `~/apps/velveth k#{k} #{k} -fastq -shortPaired #{l} #{r} > k#{k}.log`
  `~/apps/velvetg k#{k} -read_trkg yes -cov_cutoff 3 -scaffolding yes -min_contig_lgth 200 -ins_length #{opts.insertsize} -ins_length_sd #{opts.insertsd} -min_pair_count 6 >> k#{k}.log`
  `~/apps/oases k#{k} -ins_length #{opts.insertsize} -ins_length_sd #{opts.insertsd} -scaffolding yes -cov_cutoff 3 -min_trans_lgth 200 -edgeFractionCutoff 0.1 -min_pair_count 6 >> k#{k}.log`
  if File.size("k#{k}/transcripts.fa") == 0
    puts "no contigs assembled, descending..."
    next
  end
  puts "building bowtie index for k#{k} assembly"
  `/usr/local/bin/bowtie2-build k#{k}/transcripts.fa k#{k} >> k#{k}.log`
  puts "filtering out reads mapping to k#{k} assembly"
  `/usr/local/bin/bowtie2 --very-fast -p 20 --un-conc k#{k} --no-discordant --phred64 -x k#{k} -1 #{l} -2 #{r} > /dev/null`
  puts "storing unmapped reads for next step of descent"
  l = "k#{k}.1"
  r = "k#{k}.2"
  delta = Time.now - t1
  puts "assembly with k = #{k} completed in #{delta} seconds"
end
delta = Time.now - t0
puts "all velvet/oases runs completed in #{delta} seconds"
