#! /bin/env ruby


require "getoptlong"


##################################################################
infile = nil
list_file = nil

gene_lists = Hash.new


##################################################################
opts = GetoptLong.new(
  ["-i", "--in", GetoptLong::REQUIRED_ARGUMENT],
  ["--list", GetoptLong::REQUIRED_ARGUMENT],
)

opts.each do |opt, value|
  case opt
    when "-i", "--in"
      infile = value
    when "--list"
      list_file = value
  end
end


##################################################################
File.open(list_file, "r").each_line do |line|
  line.chomp!
  gene_lists[line] = ""
end


File.open(infile, "r").each_line do |line|
  line.chomp!
  sl_gene, num_of_sl_gene, pair = line.split("\t")
  genes = pair.split("-")
  num_of_sl_gene = genes.all?{|i|gene_lists.include?(i)} ? 2 : 1
  puts [sl_gene, num_of_sl_gene, pair].map{|i|i.to_s}.join("\t")
end


