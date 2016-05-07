#! /bin/env ruby

require "getoptlong"


#################################################################
gff_file = nil
longest_rela_file = nil
features = ["CDS"]
attr = "Parent";
coordinates = Hash.new{|h,k|h[k]=[]}
chr_info = Hash.new
strand_info = Hash.new

rela = Hash.new


#################################################################
opts = GetoptLong.new(
  ["--gff", GetoptLong::REQUIRED_ARGUMENT],
  ["--longest_rela", GetoptLong::REQUIRED_ARGUMENT],
  ["--feature", GetoptLong::REQUIRED_ARGUMENT],
  ["--attr", GetoptLong::REQUIRED_ARGUMENT],
)


opts.each do |opt, value|
  case opt
    when '--gff'
      gff_file = value
    when '--longest_rela'
      longest_rela_file = value
    when '--feature'
      features << value
    when '--attr'
      attr = value
  end
end


#################################################################
File.open(longest_rela_file).each_line do |line|
  line.chomp!
  line_arr = line.split("\t")
  rela[line_arr[1]] = line_arr[0]
end


File.open(gff_file).each_line do |line|
  line.chomp!
  line_arr = line.split("\t")
  #scaffold5.1|size591573	.	CDS	1	51	.	-	.	Parent=symbB.v1.2.000001
  start, stop = line_arr.values_at(3,4).map{|i|i.to_i}
  chr = line_arr[0]
  strand = line_arr[6]
  next if not features.include?(line_arr[2])
  if line_arr[-1] =~ /#{attr}=([^;]+)/
    if rela.include?($1)
      chr_info[$1] = chr
      coordinates[$1] << start << stop
      strand_info[$1] = strand
    end
  end
end


coordinates.each_pair do |gene, v|
  start, stop = v.minmax
  puts [chr_info[gene], '.', 'gene', start.to_s, stop.to_s, '.', strand_info[gene], '.', "ID="+rela[gene]].join("\t")
end


