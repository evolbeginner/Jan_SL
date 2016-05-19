#! /bin/env ruby

require 'getoptlong'
require 'String'


#################################################################################################
yn00_file=nil
exon_num_file=nil
bit_scores=Hash.new
kss=Hash.new
separator='-'
yn00_fields=Array.new
exon_num_fields=Array.new
lists_included=Array.new
items_included=Hash.new
lists_excluded=Array.new
items_excluded=Hash.new
suffix = nil
is_strict_sl = false


#################################################################################################
class String
  def numeric?
    Float(self) != nil rescue false
  end
end


#################################################################################################
opts = GetoptLong.new(
  ['-y', "--yn00_file","--yn00", "--kaks_file", GetoptLong::REQUIRED_ARGUMENT],
  ["-e", "--exon",GetoptLong::REQUIRED_ARGUMENT],
  ["--sep","--separator",GetoptLong::REQUIRED_ARGUMENT],
  ["--yn00_field", "--kaks_field", GetoptLong::REQUIRED_ARGUMENT],
  ["--exon_field", '--exon_num_field', '--exonNum_field', GetoptLong::REQUIRED_ARGUMENT],
  ["--list_included",GetoptLong::REQUIRED_ARGUMENT],
  ["--list_excluded",GetoptLong::REQUIRED_ARGUMENT],
  ["--suffix",GetoptLong::REQUIRED_ARGUMENT],
  ["--strict_sl",GetoptLong::NO_ARGUMENT],
)


opts.each do |opt,value|
  case opt
    when '-y', '--yn00_file', "--yn00", "--kaks_file"
      yn00_file=value
    when "-e", "--exon"
      exon_num_file=value
    when '--sep', '--separator'
      separator=value
    when '--yn00_field', '--kaks_field'
      yn00_fields = value.split(',').map{|i|i.to_i}
    when '--exon_field', '--exon_num_field', '--exonNum_field'
      exon_num_fields = value.split(',').map{|i|i.to_i}
    when '--list_included'
      lists_included.push value
    when '--list_excluded'
      lists_excluded.push value
    when '--suffix'
      suffix = value
    when '--strict_sl'
      is_strict_sl = true
  end
end

if not yn00_file or not exon_num_file
  raise "yn00_file and exon_num_file have to be given!"
end

exon_num_fields=[1,3] if exon_num_fields.empty?
yn00_fields=[1,3] if yn00_fields.empty?


#################################################################################################
if ! lists_included.empty?
  lists_included.each do |file|
    File.open(file,'r').each_line do |line|
      items_included[line.chomp!]=1
    end
  end
end


if ! lists_excluded.empty?
  lists_excluded.each do |file|
    File.open(file,'r').each_line do |line|
      items_excluded[line.chomp!]=1
    end
  end
end


#################################################################################################
File.open(yn00_file,'r').each_line do |line|
  next if line =~ /^#/
  line.chomp!
  line_arr = line.split("\t")
  next if line_arr.size < 4
  pair, bit_score = line_arr.values_at(yn00_fields[0]-1,yn00_fields[1]-1).map{|i|i.numeric? ? i.to_f : i}
  genes = pair.split('-')
  genes.each do |gene|
    if ! items_included.empty?
      next if ! items_included.include?(gene)
    end
    if ! items_excluded.empty?
      next if items_excluded.include?(gene)
    end

    if ! suffix.nil?
      gene.sub!(/#{suffix}/, "")
    end
    bit_scores[gene] = (! bit_scores.include? gene or bit_score>bit_scores[gene]) ? bit_score : bit_scores[gene]
  end
end


File.open(exon_num_file,'r').each_line do |line|
  line.chomp!
  line_arr = line.split("\t")
  pair, ks = line_arr.values_at(exon_num_fields[0]-1, exon_num_fields[1]-1)
  if is_strict_sl
    next if pair.split(separator).select{|i|i if bit_scores.include?(i)}.size >= 2
  end
  pair.split(separator).map{|i| kss[i]=ks if ks =~ /\d/}
end


bit_scores.keys.select{|i|kss.include? i}.each do |k|
  puts [k,bit_scores[k],kss[k]].join("\t")
end


