#! /bin/env ruby


require "getoptlong"


################################################################################
infile = nil
tandem_slr_file = nil
num_of_genes_in_clstr = 0


rela = Hash.new{|h,k|h[k]={}}
slrgs = Hash.new
slr_counts = Hash.new
best_one = Hash.new
retrogenes = Hash.new
dna_based_duplicates = Hash.new


################################################################################
def read_tandem_slr_file(infile)
  slr_counts = Hash.new
  File.open(infile).each_line do |line|
    line.chomp!
    line_arr = line.split("\t")
    num_of_slrgs = (line_arr.size-1)/2
    slrg = line_arr[0]
    slr_counts[slrg] = num_of_slrgs
  end
  return(slr_counts)
end


################################################################################
opts = GetoptLong.new(
  ["-i", GetoptLong::REQUIRED_ARGUMENT],
  ["--tandem_slr", GetoptLong::REQUIRED_ARGUMENT],
  ["-n", GetoptLong::REQUIRED_ARGUMENT],
)

opts.each do |opt, value|
  case opt
    when "-i"
      infile = value
    when "--tandem_slr"
      tandem_slr_file = value
    when "-n"
      num_of_genes_in_clstr = value.to_i
  end
end


################################################################################
slr_counts = read_tandem_slr_file(tandem_slr_file)


File.open(infile, "r").each_line do |line|
  line.chomp!
  line_arr = line.split("\t")
  slrg, num_of_slrg, pair = line_arr
  genes = pair.split("-")
  parent = nil
  genes.each do |gene|
    if gene != slrg
      parent  = gene
    end
  end
  type = nil
  if num_of_slrg == "2"
    if slr_counts.include?(slrg) and slr_counts.include?(parent)
      if slr_counts[slrg] == slr_counts[parent]
        type = "DNA"
      else
        type = "RNA"
      end
    elsif slr_counts.include?(slrg) and not slr_counts.include?(parent)
      type = "RNA"
    elsif not slr_counts.include?(slrg) and slr_counts.include?(parent)
      type = "RNA"
    elsif not slr_counts.include?(slrg) and not slr_counts.include?(parent)
      type = "DNA"
    end
  else
    type = "RNA"
  end
  rela[slrg][parent] = type
  best_one[slrg] = parent
  slrgs[slrg] = ""
end


################################################################################
rela.each_pair do |gene, v|
  if v.values.include?('RNA')
    retrogenes[gene] = ""
  else
    dna_based_duplicates[gene] = ""
  end
end

puts retrogenes.size
puts dna_based_duplicates.size


