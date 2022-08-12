#! /bin/env ruby


require 'getoptlong'
require 'parallel'

require 'Dir'


######################################################################################
indir = nil
outdir = nil
is_force = false
cpu = 1


######################################################################################
opts = GetoptLong.new(
  ['--indir', GetoptLong::REQUIRED_ARGUMENT],
  ['--outdir', GetoptLong::REQUIRED_ARGUMENT],
  ['--force', GetoptLong::NO_ARGUMENT],
  ['--cpu', GetoptLong::REQUIRED_ARGUMENT],
)


opts.each do |opt, value|
  case opt
    when '--indir'
      indir = value
    when '--outdir'
      outdir = value
    when '--force'
      is_force = true
    when '--cpu'
      cpu = value.to_i
  end
end


######################################################################################
mkdir_with_force(outdir, is_force)
infiles = read_infiles(indir)

results = Parallel.map(infiles, in_threads: cpu) do |infile|
  a = `awk '/SL-1/{a++; if(a==1){e=$12}; if($12==e){print}}' #{infile} | shuf | head -1`.chomp
  #puts "awk '/SL-1/{a++; if(a==1){e=$12}; if($12==e){print}}' #{infile} | shuf | head -1"
end


######################################################################################
outfile = File.join(outdir, 'blast_1st')
out_fh = File.open(outfile, 'w')
results.each do |i|
  next if i !~ /\w/
  out_fh.puts i
end

out_fh.close


