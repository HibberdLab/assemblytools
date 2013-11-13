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
    opt :threshold, "How many Ns in a row required to split", :default => 2, :type => :int
  end


  Trollop::die :fasta, "must exist" if !File.exist?(opts[:fasta]) if opts[:fasta]

  genome = Bio::FastaFormat.open(opts.fasta)
  count=0
  genome.each do |entry|
    # puts "#{entry.definition}"

    a = split_around_N(entry.seq, opts.threshold)
    puts a
    # if opts.rename
    #   a.each do |seq|
    #     puts ">#{count}"
    #     puts "#{seq}"
    #     count+=1
    #   end
    # else
    #   a.each_with_index do |seq, i|
    #     puts ">#{entry.definition}#{i}"
    #     puts "#{seq}"
    #   end
    # end

  end

end