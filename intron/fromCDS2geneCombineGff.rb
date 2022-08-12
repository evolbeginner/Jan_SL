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
upstream = nil

rela = Hash.new


#################################################################
opts = GetoptLong.new(
  ["--gff", GetoptLong::REQUIRED_ARGUMENT],
  ["--longest_rela", GetoptLong::REQUIRED_ARGUMENT],
  ["--feature", GetoptLong::REQUIRED_ARGUMENT],
  ["--attr", GetoptLong::REQUIRED_ARGUMENT],
  ["--upstream", GetoptLong::REQUIRED_ARGUMENT],
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
    when '--upstream'
      upstream = value.to_i
  end
end


#################################################################
if not longest_rela_file.nil?
  File.open(longest_rela_file).each_line do |line|
    line.chomp!
    line_arr = line.split("\t")
    rela[line_arr[1]] = line_arr[0]
  end
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
    if not rela.empty?
      next if not rela.include?($1)
    end
    chr_info[$1] = chr
    coordinates[$1] << start << stop
    strand_info[$1] = strand
  end
end


coordinates.each_pair do |gene, v|
  start, stop = v.minmax
  out_gene = rela.empty? ? gene : rela[gene]

  if not upstream.nil?
    if strand_info[gene] == '+'
      new_start = [start - upstream, 1].max
      new_stop = start - 1
    else
      new_stop = stop + upstream
      new_start = stop + 1
    end
    puts [chr_info[gene], '.', 'upstream', new_start.to_s, new_stop.to_s, '.', strand_info[gene], '.', "ID="+out_gene].join("\t")

  else
    puts [chr_info[gene], '.', 'gene', start.to_s, stop.to_s, '.', strand_info[gene], '.', "ID="+out_gene].join("\t")

  end
end


