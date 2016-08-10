#! /bin/env ruby

require 'getoptlong'
require 'Hash'
require 'bio'


########################################################################
class String
  def numeric?
    Float(self) != nil rescue false
  end
end


def read_seq_file(seq_file)
  seq_info = Hash.new
  Bio::FlatFile.open(seq_file).each_entry do |f|
    seq_info[f.definition] = f.seq
  end
  return(seq_info)
end


########################################################################
blast8_file=nil
evalue_cutoff = 1
include_list = nil
suffix = nil
queries = Array.new
evalue_higher_level_cutoff = 0.001
seq_file = nil
is_output_seq = false
is_output_query = false

posi_subject=Hash.new
genes_included = Hash.new
final_subjects = Hash.new


########################################################################
opts = GetoptLong.new(
  ["--blast8", "--blast_file","--blast8_file",GetoptLong::REQUIRED_ARGUMENT],
  ["-e", "--evalue", "--e_value", GetoptLong::REQUIRED_ARGUMENT],
  ["--include_list", GetoptLong::REQUIRED_ARGUMENT],
  ["--suffix", GetoptLong::REQUIRED_ARGUMENT],
  ["--query", GetoptLong::REQUIRED_ARGUMENT],
  ["--evalue_higher_level", "--e0", GetoptLong::REQUIRED_ARGUMENT],
  ["--seq", "--seq_file", GetoptLong::REQUIRED_ARGUMENT],
  ["--output_query", GetoptLong::NO_ARGUMENT],
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
    when "--query"
      value.split(",").map{|i|queries << i}
    when "--evalue_higher_level", "--e0"
      evalue_higher_level_cutoff = value.to_f
    when "--seq", "--seq_file"
      seq_file = value
      is_output_seq = true
    when "--output_query"
      is_output_query = true
  end
end


########################################################################
if not seq_file.nil?
  seq_info = read_seq_file(seq_file)
end


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
  query, subject,subject_start,subject_stop,evalue,bit_score = line.split("\t").values_at(0,1,8,9,10,11).map{|i|i.numeric? ? i.to_f : i}
  next if not queries.include?(query) if not queries.empty?
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
  posi_subject[subject][subject_range]['bit_score'] = bit_score
  posi_subject[subject][subject_range]['evalue'] = evalue.to_f
  posi_subject[subject][subject_range]["query"] = query
end


posi_subject.each_pair do |subject,ranges|
  next if ranges.size == 1
  posi_infos = ranges.keys.sort_by{|r|r.to_a[0]}.map{|k|k}
  bit_score_infos  = ranges.keys.sort_by{|r|r.to_a[0]}.map{|k|ranges[k]['bit_score']}
  evalues = ranges.keys.sort_by{|r|r.to_a[0]}.map{|k|ranges[k]['evalue']}
  tmp_queries = ranges.keys.sort_by{|r|r.to_a[0]}.map{|k|ranges[k]['query']}
  next if not evalues.any?{|i|i.to_f <= evalue_higher_level_cutoff}
  final_subjects[subject] = ""
  if is_output_query
    puts [subject,tmp_queries,bit_score_infos].join("\t")+['',posi_infos].join("\t")
  else
    puts [subject,bit_score_infos].join("\t")+['',posi_infos].join("\t")
  end
  #{|a|a[0]}.map{|i|posi_subject[subject][i]['bit_score'].to_s}].join("\t")
end


########################################################################
if is_output_seq
  final_subjects.each_key do |subject|
    ranges = posi_subject[subject].keys
    min, max = ranges.map{|i|i.to_a.first}.min, ranges.map{|i|i.to_a.last}.max
    puts ">#{subject}"
    #p [min,max]
    puts seq_info[subject][min-1, max-min]
    #ranges.each do |range|
    #  start, stop = range.to_a.values_at(0,-1)
    #  puts seq_info[subject][start-1, stop-start]
    #end
  end
end


