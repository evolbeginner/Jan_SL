#! /bin/env ruby2.1

# calculate correlations for gene pairs 

#############################################################################
require 'getoptlong'
require 'basic_math'


inputs = Array.new
list_file = nil
gene_sep = "\t"
is_sort_genes_in_pair = true
fields = [0,1]
field_sep = "\t"
is_log = true
allele_RegExp = nil

gene_pairs = Array.new
value_info = Hash.new{|h,k|h[k]={}}
fpkm_added = 0.0000001


#############################################################################
def read_list_file(list_file, sep="\t", is_sort=false, allele_RegExp=nil)
  gene_pairs = Array.new
  fh = File.open(list_file, 'r')
  fh.each_line do |line|
    line.chomp!
    genes = line.split(sep)
    genes.sort! if is_sort
    genes = genes.map{|gene|del_allele(gene, allele_RegExp)} if ! allele_RegExp.nil?
    gene_pairs.push genes
  end
  fh.close
  return gene_pairs
end


def read_numbers_file(file, value_info, counter, fields=[0,1], sep="\t", allele_RegExp=nil, is_log=true, fpkm_added=0.0000001)
  fh = File.open(file, 'r')
  fh.each_line do |line|
    next if $. == 1
    line.chomp!
    gene_name, value  = line.split(sep).values_at(fields[0], fields[1])
    gene_name = del_allele(gene_name, allele_RegExp) if ! allele_RegExp.nil?
    value = value.to_f
    if not is_log
      value_info[gene_name][counter] = Math::log(value+fpkm_added)
    else
      value_info[gene_name][counter] = value
    end
  end
  fh.close
  return(value_info)
end


def del_allele(gene_name, allele_RegExp)
  new_gene_name = gene_name.sub(/#{allele_RegExp}/, "")
  return(new_gene_name)
end


#############################################################################
opts = GetoptLong.new(
  ["-i", "--in", "--input", GetoptLong::REQUIRED_ARGUMENT],
  ["--input_find", GetoptLong::REQUIRED_ARGUMENT],
  ["--list", "--gene_list", GetoptLong::REQUIRED_ARGUMENT],
  ["--gene_sep", GetoptLong::REQUIRED_ARGUMENT],
  ["--no_sort_genes_in_pair", GetoptLong::NO_ARGUMENT],
  ["-f", "--field", "--fields", GetoptLong::REQUIRED_ARGUMENT],
  ["--field_sep", GetoptLong::REQUIRED_ARGUMENT],
  ["--no_log", GetoptLong::NO_ARGUMENT],
  ["--allele_RegExp", "--allele_reg_exp", GetoptLong::REQUIRED_ARGUMENT],
)

opts.each do |opt,value|
  case opt
    when '-i', '--in', '--input'
      value.split(',').each do |i|
        inputs.push i
      end
    when '--input_find'
      find_sentence = value
      IO.popen(find_sentence).each_line do |file|
        inputs.push file.chomp!
      end
    when '--list', '--gene_list'
      list_file = value
    when '--gene_sep'
      gene_sep = value
    when '--no_sort_genes_in_pair'
      is_sort_genes_in_pair = false
    when '-f', '--field', '--fields'
      fields = value.split(',').map{|i|i.to_i-1}
    when '--field_sep'
      field_sep = value
    when '--no_log'
      is_log = false
    when "--allele_RegExp", "--allele_reg_exp"
      allele_RegExp = value
  end
end


#############################################################################
gene_pairs = read_list_file(list_file, gene_sep, is_sort_genes_in_pair, allele_RegExp)

inputs.each_with_index do |file, index|
  value_info = read_numbers_file(file, value_info, index+1, fields, field_sep, allele_RegExp, is_log, fpkm_added)
end


gene_pairs.each do |pair_genes|
  values = pair_genes.map{|gene|value_info[gene].values}
  #next if values[0] == values[1]
  next if values[0].size != values[1].size
  puts pair_genes.join("-") + "\t" + spearman_correlate(values[0], values[1]).to_s
  #puts pair_genes.join("-") + "\t" + pearson_correlate(values[0], values[1]).to_s
end


