#! /usr/bin/ruby

require 'rubygems'
require 'trollop'

# options
opts = Trollop::options do
  version "v0.0.1a"
  banner <<-EOS

--------------------- merge oases --------------------

merge assemblies with velvet-oases

-------------------------------------------------------

EOS
  opt :outdir, "output directory", :required => true, :type => String
  opt :infiles, "list of input files", :required => true, :type => String
  opt :mergek, "merge k-mer size", :required => true, :type => Integer
end

puts "running velveth..."
`~/apps/velveth #{opts.outdir}/ #{opts.mergek} -long #{opts.infiles}`
puts "running velvetg..."
`~/apps/velvetg #{opts.outdir}/ -read_trkg yes -conserveLong yes -scaffolding yes`
puts "merging with oases"
`~/apps/oases #{opts.outdir}/ -merge yes -scaffolding yes -min_trans_lgth 200`
puts "finalising merged assembly"
`~/scripts/finaliseassembly.rb --input=#{opts.outdir}/transcripts.fa`
puts "done! final processed assembly is in final.transcripts.fa"
