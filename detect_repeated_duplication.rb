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

clstrs = Hash.new{|h,k|h[k]=[]}
final_clstrs = Hash.new()


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


def identify_clstr(gene1, rela, already={}, genes=[])
  already[gene1] = ""
  rela[gene1].each_pair do |gene2, v2|
    next if already.include?(gene2)
    genes << gene2
    rela[gene2].each_pair do |gene3, v3|
      if gene1 == gene3
        ;
      else
        genes << identify_clstr(gene2, rela, already, genes).flatten
      end
    end
  end
  return(genes)
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
  rela[parent][slrg] = type
  best_one[slrg] = parent
  slrgs[slrg] = ""
end


################################################################################
rela.each_pair do |gene, v|
  clstr = identify_clstr(gene, rela)
  clstr << gene
  clstr.flatten!
  clstr.sort!.uniq!
  final_clstrs[clstr] = ""
end


################################################################################
final_clstrs.each_pair do |clstr, v1|
  if clstr.size < num_of_genes_in_clstr
    next
  end
  puts clstr.join("\t")
  clstr.each do |gene|
    if slrgs.include?(gene)
      rela[gene].each_pair do |parent, v2|
        best = nil
        if parent == best_one[gene]
          best = "best"
        end
        if v2 == "DNA" and [gene, parent] != [gene, parent].sort
          next
        end
        if not slr_counts.include?(gene)
          puts [gene+'-'+parent, v2, best].flatten.join("\t")
        else
          puts [gene+'-'+parent, v2, best, slr_counts[gene].to_s].flatten.join("\t")
        end
      end
    end
  end
end


