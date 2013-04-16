#! /usr/bin/ruby
# todo: use better velvet/oases params

require 'rubygems'
require 'trollop'
require 'bio'

# options
$opts = Trollop::options do
  version "v0.0.1a"
  banner <<-EOS

--------------------- k-descender ---------------------

run assembler (by default, velvet-oases) with descending
k-mers in a range (default, 90:25, step 6) for paired-end
reads

Richard Smith, April 2013

-------------------------------------------------------

EOS
  opt :assembler, "assembler to use (choice of: velvetoases or soapdt)", :required => true, :type => String
  opt :left, "left reads file in fastq format", :required => true, :type => String
  opt :right, "right reads file in fastq format", :required => true, :type => String
  opt :insertsize, "mean insert size",  :required => true, :type => Integer
  opt :insertsd, "insert size standard deviation", :required => true, :type => Integer
  opt :threads, "number of threads for bowtie", :required => true, :type => Integer
  opt :mink, "minimum K-mer size (must be odd)", :default => 25, :type => Integer
  opt :maxk, "maximum K-mer size (must be odd)", :default => 85, :type => Integer
  opt :kstep, "K-mer descent step size (must be even)", :default => 10, :type => Integer
end

# catch exit statuses (cluster sniping detector!)
trap("SIGTERM") {
  puts "process killed!"
}
trap("SIGINT") {
  puts "process interrupted!"
}
trap("SIGQUIT") {
  puts "process quitting!"
}

# ensure valid assembler is specified
validassemblers = ['velvetoases', 'soapdt']
unless validassemblers.include? $opts.assembler
  abort("Assembler option '#{$opts.assembler}' invalid, please choose from [#{validassemblers.join(', ')}]")
end

# assembler functions - should accept a k-mer size, left and right read files
# and a boolean flag for the first assembly
def runSoap(k, l, r, first)
  # make config file
  File.open("soapdt.config", "w") do |conf|
    conf.puts "max_rd_len=20000"
    conf.puts "[LIB]"
    conf.puts "avg_ins=#{$opts.insertsize}"
    conf.puts "reverse_seq=0"
    conf.puts "asm_flags=3"
    conf.puts "rank=2"
    conf.puts "q1=#{l}"
    conf.puts "q2=#{r}"
    if !first
      conf.puts "[LIB]"
      conf.puts "asm_flags=2"
      conf.puts "rank=1" # prioritise the higher-k contigs in scaffolding
      conf.puts "longreads.fa"
    end
  end
  
  # construct command
  cmd = "/applications/soapdenovo-trans/SOAPdenovo-Trans-127mer all"
  cmd += " -s soapdt.config" # config file
  cmd += " -a 60" # memory assumption
  cmd += " -o k#{k}" # output directory
  cmd += " -K #{k}" # kmer size
  cmd += " -p #{$opts.threads}" # number of threads
  cmd += " -d 3" # minimum kmer frequency
  cmd += " -F" # fill gaps in scaffold
  cmd += " -M 1" # strength of contig flattening
  cmd += " -D 1" # delete edges with coverage no greater than
  cmd += " -L 200" # minimum contig length
  # cmd += " -u" # unmask high coverage contigs before scaffolding
  cmd += " -e 2" # delete contigs with coverage no greater than
  cmd += " -t 5" # maximum number of transcripts from one locus
  # cmd += " -S" # scaffold structure exists

  # run command
  `#{cmd} > k#{k}.log`

  # cleanup unneeded files
  `mkdir k#{k}`
  `mv k#{k}.scafSeq k#{k}/transcripts.fa`
  `rm k#{k}.*`
end

# set thread number for velvet
if $opts.assembler == 'velvetoases'
  `sh setmpthreads.sh`
end

def runVelvetOases(k, l, r, first)
  longstr = first ? "" : "-fasta -long longreads.fa "
  `~/apps/velveth k#{k} #{k} -fasta -shortPaired #{l} #{r} #{longstr}> k#{k}.log`
  `~/apps/velvetg k#{k} -read_trkg yes -cov_cutoff 3 -scaffolding yes -min_contig_lgth 200 -ins_length #{$opts.insertsize} -ins_length_sd #{$opts.insertsd} -min_pair_count 6 >> k#{k}.log`
  `~/apps/oases k#{k} -ins_length #{$opts.insertsize} -ins_length_sd #{$opts.insertsd} -scaffolding yes -cov_cutoff 3 -min_trans_lgth 200 -edgeFractionCutoff 0.1 -min_pair_count 6 >> k#{k}.log`
  # cleanup unneeded files
  Dir.chdir("k#{k}") do
    `rm Sequences Roadmaps PreGraph Graph2 LastGraph`
  end
end
puts "running #{$opts.assembler} k-descent from 90 to 25"

l = $opts.left
r = $opts.right

first = true
t0 = Time.now # time the whole run
toclean = nil # don't clean up during the first round
($opts.mink..$opts.maxk).step($opts.kstep).reverse_each do |k|
  t1 = Time.now # time this assembly
  puts "running assembly with k = #{k}"
  case $opts.assembler
  when 'velvetoases'
    runVelvetOases(k, l, r, first)
  when 'soapdt'
    runSoap(k, l, r, first)
  end
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
  `/usr/local/bin/bowtie2 --very-fast -f -p #{$opts.threads} --un-conc k#{k} --no-discordant --phred64 -x k#{k} -1 #{l} -2 #{r} > /dev/null`
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
  toclean = "k#{k}"
  delta = Time.now - t1
  first = false
  puts "assembly with k = #{k} completed in #{delta} seconds"
end

delta = Time.now - t0
puts "all velvet/oases runs completed in #{delta} seconds"
