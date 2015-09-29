#! /env/bin ruby

require 'getoptlong'


###################################################
is_parent = false
is_SL = false
final_pair_yn00_file = nil
infile = nil

ks_min = 0
ks_max = 1000
gene_pairs = Hash.new


###################################################
class Gene_pair
  attr_accessor :SL, :parent
end


###################################################
opts = GetoptLong.new(
  ['--parent', GetoptLong::NO_ARGUMENT],
  ['--SL', '--sl', GetoptLong::NO_ARGUMENT],
  ['--final_pair_yn00', GetoptLong::REQUIRED_ARGUMENT],
  ['-i', '--in', GetoptLong::REQUIRED_ARGUMENT],
  ['--ks', GetoptLong::REQUIRED_ARGUMENT],
)

opts.each do |opt, value|
  case opt
    when '--parent'
      is_parent = true
    when '--SL', '--sl'
      is_SL = true
    when '--final_pair_yn00'
      final_pair_yn00_file = value
    when '-i', '--in'
      infile = value
    when '--ks'
      kss = value.split('-').map{|i| i.to_f}
      ks_min, ks_max = kss
  end
end


###################################################
File.open(final_pair_yn00_file, 'r').each_line do |line|
  #symbB.v1.2.000501.t1	1	symbB.v1.2.000501.t1-symbB.v1.2.007689.t1
  line.chomp!
  sl_gene, gene_pair = line.split("\t").values_at(0,2)
  genes = gene_pair.split("-")
  gene_pairs[$.] = Gene_pair.new
  parent_gene = genes.select{|i| i != sl_gene }
  gene_pairs[$.].SL = sl_gene
  gene_pairs[$.].parent = parent_gene
end


File.open(infile, 'r').each_line do |line|
  #symbB.v1.2.000022.t1_symbB.v1.2.012916.t1	0.8123	4.4974	0.1806
  line.chomp!
  gene_pair, ka, ks, ka_ks = line.split("\t")
  if ks =~ /\d/
    ks = ks.to_f
    if ks > ks_min and ks <= ks_max
      if gene_pairs.include?($.)
        if is_SL
          puts gene_pairs[$.].SL
        elsif is_parent
          puts gene_pairs[$.].parent
        end
      end
    end
  end
end



