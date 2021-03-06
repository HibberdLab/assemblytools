#!/usr/bin/env ruby

require 'rubygems'
require 'trollop'

# post-process assembly

# options
opts = Trollop::options do
  version "v0.0.1a"
  banner <<-EOS

----------------- post-process assembly ------------------

post-process an assembly by sorting, flattening, dereplicating
clustering, and renaming the transcripts to a standard format

--------------------------------------------------------

EOS
  opt :input, "Assembly file in fasta format", :required => true, :type => String
  opt :idthreshold, "Redundancy ID similarity threshold (default: 0.95)", :default => 0.95, :type => Float
  opt :polyacutoff, "Number of polyA/T to use as cutoff (default: 12)", :default => 12, :type => Integer
  opt :cap3, "Run CAP3 to further scaffold the assembly (default: no)"
  opt :nocleanup, "Keep intermediate files (default: no)"
  opt :onlyrename, "Only rename the sequences (default: no)"
end

infile = opts.input.split(/\//)
indir = infile.length > 1 ? infile[0..infile.length-2].join("/") + "/" : ''
infile = infile[infile.length-1]

# sequence renaming
def renameSeqs(outfile, infile)
  File.open(outfile, "w") do |output|
    i = 0
    File.open(infile, "r").each do |line|
      if line =~ /^>/
        line = ">transcript_#{i}"
        i += 1
      end
      output << line.gsub(/\n/, '') + "\n"
    end
  end
end

if opts.onlyrename
  renameSeqs("renamed.#{infile}", infile)
  exit
end

usearch = `which usearch`

# remove polyA
puts "removing polyA/T tails of #{opts.polyacutoff} or longer"
`prinseq-lite.pl -fasta #{opts.input} -out_good #{indir}trimmed.#{infile} -trim_tail_left #{opts.polyacutoff} -trim_tail_right #{opts.polyacutoff}`
`mv #{indir}trimmed.#{infile}.fasta #{indir}trimmed.#{infile}`

# flatten the assembly
puts "sorting by length"
`#{usearch} --sortbylength #{indir}trimmed.#{infile} --output #{indir}sorted.#{infile} -quiet`
`rm #{indir}trimmed.#{infile}` unless opts.nocleanup
puts "removing completely redundant transcripts (100% sequence similarity)"
`#{usearch} --cluster_smallmem #{indir}sorted.#{infile} --id 1 --centroids #{indir}centroids.#{infile} -quiet`
`rm #{indir}sorted.#{infile}` unless opts.nocleanup

if opts.cap3
  # `cap3 merged\/cap_in.fasta >> assembly_report.txt`
  # `rm centroids.#{infile}` unless opts.nocleanup
  # system "cat merged\/cap_in.fasta.cap.contigs > merged\/cap_in.fasta";
  # system "cat merged\/cap_in.fasta.cap.singlets >> merged\/cap_in.fasta";
end

puts "standardising transcript names"
renameSeqs "#{indir}renamed.#{infile}", "#{indir}centroids.#{infile}"
`rm #{indir}centroids.#{infile}` unless opts.nocleanup

puts "re-sorting by length"
`#{usearch} --sortbylength #{indir}renamed.#{infile} --output #{indir}sorted_renamed.#{infile} -quiet`
`rm #{indir}renamed.#{infile}` unless opts.nocleanup
puts "clustering with similarity cutoff #{opts.idthreshold}"
`#{usearch} --cluster_smallmem #{indir}sorted_renamed.#{infile} --id #{opts.idthreshold} --centroids #{indir}final.#{infile} -quiet`
`rm #{indir}sorted_renamed.#{infile}` unless opts.nocleanup

puts "assembly post-processing finished! final assembly in file #{indir}final.#{infile}"
