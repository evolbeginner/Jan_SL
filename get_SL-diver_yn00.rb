#! /bin/env ruby

require 'getoptlong'
require 'String'


#################################################################################################
blast8_file=nil
kaks_file=nil
bit_scores=Hash.new
kss=Hash.new
separator='_'
blast8_fields=Array.new
kaks_fields=Array.new
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
  ["--blast8_file","--blast_file",GetoptLong::REQUIRED_ARGUMENT],
  ["--kaks_file",GetoptLong::REQUIRED_ARGUMENT],
  ["--sep","--separator",GetoptLong::REQUIRED_ARGUMENT],
  ["--blast8_field",GetoptLong::REQUIRED_ARGUMENT],
  ["--kaks_field",GetoptLong::REQUIRED_ARGUMENT],
  ["--list_included",GetoptLong::REQUIRED_ARGUMENT],
  ["--list_excluded",GetoptLong::REQUIRED_ARGUMENT],
  ["--suffix",GetoptLong::REQUIRED_ARGUMENT],
  ["--strict_sl",GetoptLong::NO_ARGUMENT],
)


opts.each do |opt,value|
  case opt
    when '--blast8_file', "--blast_file"
      blast8_file=value
    when '--kaks_file'
      kaks_file=value
    when '--sep', '--separator'
      separator=value
    when '--blast8_field'
      blast8_fields=value.split(',').map{|i|i.to_i}
    when '-f','--kaks_field'
      value.split(',').each do |i|
        kaks_fields.push i.to_i
      end
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

if not blast8_file or not kaks_file
  raise "blast8_file and kaks_file have to be given!"
end

kaks_fields=[1,3] if kaks_fields.empty?
blast8_fields=[2,12] if blast8_fields.empty?


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
File.open(blast8_file,'r').each_line do |line|
  next if line =~ /^#/
  line.chomp!
  subject,bit_score = line.split("\t").values_at(blast8_fields[0]-1,blast8_fields[1]-1).map{|i|i.numeric? ? i.to_f : i}
  if ! items_included.empty?
    next if ! items_included.include?(subject)
  end
  if ! items_excluded.empty?
    next if items_excluded.include?(subject)
  end

  if ! suffix.nil?
    subject.sub!(/#{suffix}/, "")
  end
  bit_scores[subject] = (! bit_scores.include? subject or bit_score>bit_scores[subject]) ? bit_score : bit_scores[subject]
end


File.open(kaks_file,'r').each_line do |line|
  line.chomp!
  line_arr = line.split("\t")
  pair, ks = line_arr.values_at(kaks_fields[0]-1, kaks_fields[1]-1)
  if is_strict_sl
    next if pair.split(separator).select{|i|i if bit_scores.include?(i)}.size >= 2
  end
  pair.split(separator).map{|i| kss[i]=ks if ks =~ /\d/}
end


bit_scores.keys.select{|i|kss.include? i}.each do |k|
  puts [k,bit_scores[k],kss[k]].join("\t")
end


