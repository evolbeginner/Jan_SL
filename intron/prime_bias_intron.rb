#! /bin/env ruby


require "getoptlong"


############################################################
gff_file = nil
attr = "Parent|ID"
include_list_file = nil

genes_included = Hash.new
cds_info = Hash.new{|h,k|h[k]=[]}
cds_length_info = Hash.new{|h,k|h[k]=[]}
cds_lengths = Hash.new
intron_posi_info = Hash.new{|h,k|h[k]=[]}
intron_prop_info = Hash.new{|h,k|h[k]=[]}


############################################################
def read_gff_file(gff_file, attr)
  strand_info = Hash.new
  cds_info = Hash.new{|h,k|h[k]=[]}
  File.open(gff_file, "r").each_line do |line|
    line.chomp!
    #scaffold5.1|size591573	.	CDS	1	51	.	-	.	Parent=symbB.v1.2.000001
    line_arr = line.split("\t")
    gene = nil
    strand = line_arr[6]
    if line_arr[-1] =~ /#{attr}=([^;]+)/
      gene = $1
    end
    start, stop = line_arr.values_at(3,4).map{|i|i.to_i}
    cds_info[gene] << [start, stop]
    strand_info[gene] = strand
  end
  return([cds_info, strand_info])
end


def get_cds_length_info(cds_info)
  cds_length_info = Hash.new{|h,k|h[k]=[]}
  cds_info.each_pair do |gene, v|
    v.each do |cds_arr|
      cds_length_info[gene] << cds_arr[-1] - cds_arr[0] + 1
    end
  end
  return(cds_length_info)
end


def read_gene_list(gene_list_file)
  genes = Hash.new
  File.open(gene_list_file, "r").each_line do |line|
    line.chomp!
    genes[line] = ""
  end
  return(genes)
end


############################################################
opts = GetoptLong.new(
  ["--gff", GetoptLong::REQUIRED_ARGUMENT],
  ["--attr", GetoptLong::REQUIRED_ARGUMENT],
  ["--include_list", GetoptLong::REQUIRED_ARGUMENT],
)

opts.each do |opt, value|
  case opt
    when '--gff'
      gff_file = value
    when '--attr'
      attr = value
    when '--include_list'
      include_list_file = value
  end
end


############################################################
genes_included = read_list(include_list_file) if not include_list_file.nil?

cds_info, strand_info = read_gff_file(gff_file, attr)

cds_length_info = get_cds_length_info(cds_info)


cds_length_info.each_pair do |gene, v|
  intron_posi = 0
  v[0, v.size-1].each do |cds_length|
    intron_posi = intron_posi + cds_length 
    intron_posi_info[gene] << intron_posi
  end
  cds_lengths[gene] = v.inject(0){|sum, i| sum=sum+i}
end


intron_posi_info.each_pair do |gene, intron_posi_arr|
  intron_posi_arr.each do |intron_posi|
    prop = intron_posi/cds_lengths[gene].to_f
    prop = 1 - prop if strand_info[gene] == "-"
    intron_prop_info[gene] << prop
  end
end


intron_prop_info.each_pair do |gene, intron_prop_arr|
  if not genes_included.empty?
    next if not genes_included.include?(gene)
  end
  intron_prop_arr.each do |intron_prop|
    puts [gene, intron_prop].map{|i|i.to_s}.join("\t")
  end
end


