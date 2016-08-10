#! /bin/env ruby


require "getoptlong"


#####################################################
exon_count_file = nil
gff_file = nil
cufflinks_files = Array.new
attr = 'Parent'
features = Array.new
min_exon_num = 1

exon_numbers = Hash.new
chr_info = Hash.new


#####################################################
def read_exon_count_file(exon_count_file)
  exon_numbers = Hash.new
  File.open(exon_count_file, 'r').each_line do |line|
    line.chomp!
    line_arr = line.split("\t")
    gene, exon_num = line_arr.values_at(0,2)
    exon_num = exon_num.to_i
    exon_numbers[gene] = exon_num
  end
  return(exon_numbers)
end


def read_gff_file(gff_file, attr, features=[])
  chr_info = Hash.new{|h,k|h[k]={}}
  File.open(gff_file, 'r').each_line do |line|
    line.chomp!
    next if line =~ /^#/
    line_arr = line.split("\t")
    chr, feature, start, stop, strand, attr_str = line_arr.values_at(0,2,3,4,6,-1)
    next if not features.include?(feature)
    if attr_str =~ /#{attr}=([^;]+)/
      gene = $1
      if ! chr_info.include?(gene)
        chr_info[gene]['starts'] = Array.new
        chr_info[gene]['stops'] = Array.new
      end
      chr_info[gene]['chr'] = chr
      chr_info[gene]['strand'] = strand
      chr_info[gene]['starts'] << start.to_i
      chr_info[gene]['stops'] << stop.to_i
    end
  end
  chr_info.each_pair do |gene, v|
    v['start'] = v['starts'].min
    v['stop'] = v['stops'].max
  end
  return(chr_info)
end



#####################################################
opts = GetoptLong.new(
  ['--exon_count', GetoptLong::REQUIRED_ARGUMENT],
  ['--gff', GetoptLong::REQUIRED_ARGUMENT],
  ['--cufflinks', GetoptLong::REQUIRED_ARGUMENT],
  ['--attr', GetoptLong::REQUIRED_ARGUMENT],
  ['--feature', GetoptLong::REQUIRED_ARGUMENT],
  ['--exon_num', GetoptLong::REQUIRED_ARGUMENT],
)

opts.each do |opt, value|
  case opt
    when '--exon_count'
      exon_count_file = value
    when '--gff'
      gff_file = value
    when '--cufflinks'
      value.split(',').each do |cufflinks_file|
        cufflinks_files << cufflinks_file
      end
    when '--attr'
      attr = value
    when '--feature'
      value.split(',').each do |i|
        features << i
      end
    when '--min_exon_num'
      min_exon_num = value.to_i
  end
end


#####################################################
exon_numbers = read_exon_count_file(exon_count_file)

chr_info = read_gff_file(gff_file, attr, features)

chr_info.each_pair do |gene, v|
  next if not exon_numbers.include?(gene)
  next if exon_numbers[gene] < min_exon_num
  puts [gene, exon_numbers[gene], v['chr'], v['start'], v['stop']].map{|i|i.to_s}.join("\t")
end

exit
cufflinks_files.each do |cufflinks_file|
  File.open(cufflinks_file).each_line do |line|
    line.chomp!
  end
end


