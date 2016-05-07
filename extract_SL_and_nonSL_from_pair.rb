#! /bin/env ruby

require "getoptlong"


#################################################################
infile = nil
outdir = nil
sep = '-'

sl_genes = Array.new
parents = Array.new
parents_all = Array.new


#################################################################
opts = GetoptLong.new(
  ["-i", GetoptLong::REQUIRED_ARGUMENT],
  ["--outdir", GetoptLong::REQUIRED_ARGUMENT],
  ["--sep", GetoptLong::REQUIRED_ARGUMENT],
)

opts.each do |opt, value|
  case opt
    when "-i"
      infile = value
    when "--outdir"
      outdir = value
    when '--sep'
      sep = value
  end
end

`mkdir -p #{outdir}`


#################################################################
File.open(infile, "r").each_line do |line|
  line.chomp!
  line_arr = line.split("\t")
  sl_gene, num_of_SL, pair = line_arr
  num_of_SL = num_of_SL.to_i
  sl_genes << sl_gene
  parent = nil
  genes = pair.split(sep) 
  if genes[0] == sl_gene
    parent = genes[1]
  else
    parent = genes[0]
  end
  parents_all << parent
  if num_of_SL == 1
    parents << parent
  else
    parents << ''
  end
end

parents.uniq!
parents_all.uniq!
sl_genes.uniq!


outfile_sl = File.join([outdir, "genes_with_SLs.e0-0.001.AlnLength10-50.list"])
outfile_parent = File.join([outdir, "parent.e0-0.001.AlnLength10-50.list"])
outfile_parent_all = File.join([outdir, "parent_all.e0-0.001.AlnLength10-50.list"])
out_fh1 = File.open(outfile_sl, "w")
out_fh2 = File.open(outfile_parent, "w")
out_fh3 = File.open(outfile_parent_all, "w")
sl_genes.map{|i|out_fh1.puts i}
parents.map{|i|out_fh2.puts i}
parents_all.map{|i|out_fh3.puts i}


