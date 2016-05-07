#! /bin/env ruby

require "getoptlong"


##################################################
infile = nil
sep = "-"
exon_count_files = Array.new

exon_counts = Hash.new


##################################################
opts = GetoptLong.new(
  ["-i", GetoptLong::REQUIRED_ARGUMENT],
  ["--sep", GetoptLong::REQUIRED_ARGUMENT],
  ["--exon_counts", "--exon_count", GetoptLong::REQUIRED_ARGUMENT],
)

opts.each do |opt, value|
  case opt
    when "-i"
      infile = value
    when "--exon_counts", "--exon_count"
      exon_count_files << value
    when "--sep"
      sep = value
  end
end


##################################################
exon_count_files.each do |exon_count_file|
  File.open(exon_count_file, 'r').each_line do |line|
    line.chomp!
    gene, exon_count = line.split("\t").values_at(0,2)
    exon_counts[gene] = exon_count.to_i
  end
end


File.open(infile, 'r').each_line do |line|
  line.chomp!
  line_arr = line.split("\t")
  sl_gene, num_of_sl_genes, pair = line_arr
  next if num_of_sl_genes == "2"
  genes = pair.split(sep)
  parent = nil
  if genes[0] == sl_gene
    parent = genes[1]
  else
    parent = genes[0]
  end
  a = [exon_counts[sl_gene], exon_counts[parent]].map{|i|i.to_s}.join("\t")
  puts [[sl_gene, parent].join(sep), a].join("\t")
end


