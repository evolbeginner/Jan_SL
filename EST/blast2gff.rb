#! /bin/env ruby


require "getoptlong"
require "bio"

require "SSW_bio"


#######################################################
infile = nil
feature = nil
identity_cutoff = 0
coverage_cutoff = 0
seq_file = nil

seq_objs = Hash.new


#######################################################
opts = GetoptLong.new(
  ["-i", GetoptLong::REQUIRED_ARGUMENT],
  ["--feature", GetoptLong::REQUIRED_ARGUMENT],
  ["--identity", GetoptLong::REQUIRED_ARGUMENT],
  ["--cov", "--coverage", GetoptLong::REQUIRED_ARGUMENT],
  ["--seq", GetoptLong::REQUIRED_ARGUMENT],
)


opts.each do |opt, value|
  case opt
    when '-i'
      infile = value
    when '--feature'
      feature = value
    when '--identity'
      identity_cutoff = value.to_f
    when '--cov'
      coverage_cutoff = value.to_f
    when '--seq'
      seq_file = value
  end
end


if not seq_file.nil?
  seq_objs = read_seq_file(seq_file)
end


#######################################################
File.open(infile).each_line do |line|
  next if line =~ /^#/
  line.chomp!
  line_arr = line.split("\t")
	#symbB.v1.2.011972	gi|186963336|gb|FE864493.1|FE864493	85.97	563	77	2	327	888	1	562	1e-170  601
  query, subject, identity, q_start, q_stop, s_start, s_stop  = line_arr.values_at(0,1,2,6,7,8,9)
  identity = identity.to_f
  q_start = q_start.to_i
  q_stop = q_stop.to_i
  s_start = s_start.to_i
  s_stop = s_stop.to_i
  strand = "+"
  cov = nil
  if not seq_objs.empty?
    cov = (q_stop-q_start+1)/seq_objs[query].seq.length.to_f
    next if cov < coverage_cutoff
  end
  next if identity < identity_cutoff * 100
  if s_start > s_stop
    s_start, s_stop = s_stop, s_start 
    strand = "-"
  end
  puts [subject, '.', feature, s_start, s_stop, '.', strand, '.', 'ID='+query].map{|i|i.to_s}.join("\t")
end


