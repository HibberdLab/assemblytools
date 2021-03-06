#!/usr/bin/env ruby

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

velveth = "~/apps/velveth"
velvetg = "~/apps/velvetg"
oases = "~/apps/oases"
finaliseassembly = File.dirname(__FILE__) + '/finaliseassembly.rb'

# check all dependencies exist
def missing(dep)
  puts "could not find #{dep} - please check the path and correct the script"
end

[velveth, velvetg, oases, finaliseassembly].each do |d|
  missing(d) if not File.exist? d
end

# run
puts "running velveth..."
`#{velveth} #{opts.outdir}/ #{opts.mergek} -long #{opts.infiles}`
puts "running velvetg..."
`#{velvetg} #{opts.outdir}/ -read_trkg yes -conserveLong yes -scaffolding yes`
puts "merging with oases"
`#{oases} #{opts.outdir}/ -merge yes -scaffolding yes -min_trans_lgth 200`
puts "finalising merged assembly"
`#{finaliseassembly} --input=#{opts.outdir}/transcripts.fa`
puts "done! final processed assembly is in final.transcripts.fa"
