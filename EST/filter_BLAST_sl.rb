#! /bin/env ruby


require "getoptlong"
require "bio"


###########################################################################
gff_file = nil
blast8_file = nil
sl_seq_file = nil

gff_info = Hash.new{|h,k|h[k]={}}
sl_seq_objs = Hash.new


###########################################################################
opts = GetoptLong.new(
  ["--gff", GetoptLong::REQUIRED_ARGUMENT],
  ["--blast", "--blast8", "-b", GetoptLong::REQUIRED_ARGUMENT],
  ["--sl_seq", GetoptLong::REQUIRED_ARGUMENT],
)

opts.each do |opt, value|
  case opt
    when "--gff"
      gff_file = value
    when "--blast", "--blast8", "-b"
      blast8_file = value
    when "--sl_seq"
      sl_seq_file = value
  end
end


###########################################################################
Bio::FlatFile.open(sl_seq_file).each_entry do |f|
  sl_seq_objs[f.definition] = f
end


File.open(gff_file, "r").each_line do |line|
  line.chomp!
  line_arr = line.split("\t")
  chr, start, stop, strand = line_arr.values_at(0,3,4,6)
  gff_info[chr]["stop"] = stop
  gff_info[chr]["strand"] = strand
end


File.open(blast8_file, "r").each_line do |line|
  line.chomp!
  next if line =~ /^#/
  line_arr = line.split("\t")
  query, subject, query_start, query_stop, subject_start, subject_stop  = line_arr.values_at(0,1,6,7,8,9)
  is_pass = false
  if gff_info[subject]['strand'] == '+'
    if subject_start.to_i <= 1
      is_pass = true
    end
  elsif gff_info[subject]['strand'] == '-'
    if (subject_start.to_i - gff_info[subject]["stop"].to_i).abs <= 1
      is_pass = true
    end 
  end
  is_pass = false if query_stop.to_i < sl_seq_objs[query].seq.size-3
  puts line if is_pass
end


