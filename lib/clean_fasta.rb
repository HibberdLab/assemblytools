#!/usr/bin/env ruby

#
# clean up a fasta file
#
# option to replace the names of the contigs with sequential numbers
# option to split contigs that contain X numbers of Ns into 2 contigs
#


require 'rubygems'
require 'trollop'
require 'bio'

def split_around_N(sequence, threshold)
  seqs = []
  count_n = 0
  count_base = 0
  sequence.chars.each_with_index do |c,i|
    if c == "N"
      count_n += 1
      if count_n > threshold
        count_base -= (count_n-1)
        if count_base > 0
          seqs << sequence.slice(i-count_base-2, count_base)
        end
        count_base = 0
      else
        count_base += 1
      end
    else
      count_n = 0
      count_base += 1
    end
  end
  if count_base > 0
    seqs << sequence.slice(sequence.length-count_base, count_base)
  end
  return seqs
end

if __FILE__ == $0

  opts = Trollop::options do
    version "v0.0.1a"
    opt :fasta, "Fasta file to be cleaned", :required => true, :type => String
    opt :rename, "Rename contigs"
    opt :uppercase, "Make all bases uppercase"
    opt :threshold, "How many Ns in a row required to split", :default => 2, :type => :int
    opt :minimum, "Minimum size of new contigs formed", :default => 100, :type => :int
  end

  Trollop::die :fasta, "must exist" if !File.exist?(opts[:fasta]) if opts[:fasta]

  genome = Bio::FastaFormat.open(opts.fasta)
  count = 0
  genome.each do |entry|
    a = split_around_N(entry.seq, opts.threshold)
    a.each_with_index do |seq, i|
      if seq.length >= opts.minimum
        if opts.rename
          puts ">#{count}"
        else
          puts ">#{entry.definition}-#{i}"
        end
        if opts.uppercase
          puts "#{seq.upcase}"
        else
          puts "#{seq}"
        end
        count += 1
      end
    end
  end
end