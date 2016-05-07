#! /bin/env ruby

require 'getoptlong'
require 'Hash'


########################################################################
class String
  def numeric?
    Float(self) != nil rescue false
  end
end


########################################################################
blast8_file=nil
evalue_cutoff = 1
include_list = nil
suffix = nil

posi_subject=Hash.new
genes_included = Hash.new


########################################################################
opts = GetoptLong.new(
  ["--blast8", "--blast_file","--blast8_file",GetoptLong::REQUIRED_ARGUMENT],
  ["-e", "--evalue", "--e_value", GetoptLong::REQUIRED_ARGUMENT],
  ["--include_list", GetoptLong::REQUIRED_ARGUMENT],
  ["--suffix", GetoptLong::REQUIRED_ARGUMENT],
)

opts.each do |opt,value|
  case opt
    when '--blast8', '--blast_file', '--blast8_file'
      blast8_file=value
    when "-e", "--evalue", "--e_value"
      evalue_cutoff = value.to_f
    when "--include_list"
      include_list = value
    when "--suffix"
      suffix = value
  end
end


########################################################################
if ! include_list.nil?
  File.open(include_list, "r").each_line do |line|
    line.chomp!
    if ! suffix.nil?
      line.sub!(/#{suffix}$/, '')
    end
    genes_included[line] = ""
  end
end


File.open(blast8_file,'r').each_line do |line|
  next if line =~ /^#/
  line.chomp!
  subject,subject_start,subject_stop,evalue,bit_score = line.split("\t").values_at(1,8,9,10,11).map{|i|i.numeric? ? i.to_f : i}
  if ! suffix.nil?
    subject.sub!(/#{suffix}$/, '')
  end
  if not include_list.nil?
    next if not genes_included.include?(subject)
  end
  subject_start=subject_start.to_i
  subject_stop=subject_stop.to_i
  subject_range=subject_start..subject_stop
  evalue = evalue.to_f
  next if evalue > evalue_cutoff
  posi_subject[subject]=multi_D_Hash(2) if not posi_subject.include? subject
  posi_subject[subject][subject_range]['bit_score']=bit_score
end


posi_subject.each_pair do |subject,ranges|
  next if ranges.size == 1
  posi_infos = ranges.keys.sort_by{|r|r.to_a[0]}.map{|k|k}
  bit_score_infos  = ranges.keys.sort_by{|r|r.to_a[0]}.map{|k|ranges[k]['bit_score']}
  puts [subject,bit_score_infos].join("\t")+['',posi_infos].join("\t")
  #{|a|a[0]}.map{|i|posi_subject[subject][i]['bit_score'].to_s}].join("\t")
end


