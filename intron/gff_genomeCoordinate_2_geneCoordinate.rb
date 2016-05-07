#! /bin/env ruby

require "getoptlong"


###############################################################
gff_file = nil
gene_name = nil
gene_list_file = nil
attribute = "ID"
features = Array.new

cds_info = Hash.new{|h,k|h[k]={}}


###############################################################
def renumber_coordinates(cds_info)
  cds_info.each_pair do |attr, v|
    strand = v["strand"]
    if strand == "+"
      v["info"].sort_by!{|i|i[0]}
      first_posi = v["info"][0][0]
      v["info"].map!{|i|i.map!{|j|j-first_posi+1}}
    else
      v["info"].sort_by!{|i|i[0]}.reverse!
      first_posi = v["info"][0][-1]
      v["info"].map!{|i|i.map!{|j|first_posi-j+1}}
      v["info"].map!{|i|i.reverse}
    end
  end
end


###############################################################
opts = GetoptLong.new(
  ["--gff", "--gff3", GetoptLong::REQUIRED_ARGUMENT],
  ["--gene_name", "--gene", GetoptLong::REQUIRED_ARGUMENT],
  ["--gene_list", GetoptLong::REQUIRED_ARGUMENT],
  ["--attribute", "--attr", GetoptLong::REQUIRED_ARGUMENT],
  ["--feature", GetoptLong::REQUIRED_ARGUMENT],
)

opts.each do |opt, value|
  case opt
    when "--gff", "--gff3"
      gff_file = value
    when "--gene_name", "--gene"
      gene_name = value
    when "--gene_list"
      gene_list_file = value
    when "--attribute", "--attr"
      attribute = value
    when "--feature"
      features << value
  end
end


###############################################################
File.open(gff_file).each_line do |line|
  line.chomp!
  line_arr = line.split("\t")
  chr, feature, start, stop, strand, attributes_str = line_arr.values_at(0,2,3,4,6,-1)
  start = start.to_i
  stop = stop.to_i
  next if not features.include?(feature)
  attr = $1 if attributes_str =~ /#{attribute}=([^;]+)/
  cds_info[attr]["strand"] = strand
  cds_info[attr]["chr"] = chr
  if ! cds_info[attr].include?("info")
    cds_info[attr]["info"] = Array.new
  end
  cds_info[attr]["info"] << [start, stop]
end


renumber_coordinates(cds_info)


cds_info.each_pair do |gene, v|
  v["info"].each do |i|
    puts [v["chr"], ".", features[0], i[0], i[1], ".", v['strand'], ".", "#{attribute}=#{gene}"].join("\t")
  end
end


