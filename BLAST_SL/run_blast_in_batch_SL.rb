#! /usr/bin/env ruby

# old name: haha.rb


#############################################
require 'getoptlong'
require 'parallel'

require 'Dir'
require 'util'


#############################################
indir = nil
gap_argu = '-ungapped'
cpu = 10
outdir = nil
is_force = false


#############################################
opts = GetoptLong.new(
  ['--indir', GetoptLong::REQUIRED_ARGUMENT],
  ['--with_gap', GetoptLong::NO_ARGUMENT],
  ['--outdir', GetoptLong::REQUIRED_ARGUMENT],
  ['--force', GetoptLong::NO_ARGUMENT],
  ['--cpu', GetoptLong::REQUIRED_ARGUMENT],
)

opts.each do |opt, value|
  case opt
    when '--indir'
      indir = value
    when '--with_gap'
      gap_argu = ''
    when '--outdir'
      outdir = value
    when '--force'
      is_force = true
    when '--cpu'
      cpu = value.to_i
  end
end


#############################################
mkdir_with_force(outdir, is_force)

infiles = read_infiles(indir)


#############################################
Parallel.map(infiles, in_threads: cpu) do |infile|
  c=getCorename(infile)
  outfile = File.join(outdir, c+'.blast6')
  `blastn -query ../SL.fas -subject #{infile} -word_size 4 -evalue 100 #{gap_argu} -strand plus -outfmt 6 > #{outfile}`
end


