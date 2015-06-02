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
posi_subject=Hash.new

opts = GetoptLong.new(
  ["--blast_file","--blast8_file",GetoptLong::REQUIRED_ARGUMENT],
)

opts.each do |opt,value|
  case opt
    when '--blast_file', '--blast8_file'
      blast8_file=value
  end
end

File.open(blast8_file,'r').each_line do |line|
  line.chomp!
  subject,subject_start,subject_stop,bit_score = line.split("\t").values_at(1,8,9,11).map{|i|i.numeric? ? i.to_f : i}
  subject_start=subject_start.to_i
  subject_stop=subject_stop.to_i
  subject_range=subject_start..subject_stop
  posi_subject[subject]=multi_D_Hash(2) if not posi_subject.include? subject
  posi_subject[subject][subject_range]['bit_score']=bit_score
end

posi_subject.each_pair do |subject,ranges|
  next if ranges.size == 1
  posi_infos = ranges.keys.sort_by{|r|r.to_a[0]}.map{|k|k}
  bit_score_infos  = ranges.keys.sort_by{|r|r.to_a[0]}.map{|k|ranges[k]['bit_score']}
  puts [subject,bit_score_infos].join("\t")
  puts ['',posi_infos].join("\t")
  #{|a|a[0]}.map{|i|posi_subject[subject][i]['bit_score'].to_s}].join("\t")
end


