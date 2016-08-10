#! /bin/env ruby


require "getoptlong"
require "bio"


############################################################
def get_seq_objs(seq_file)
  seq_objs = Hash.new
  Bio::FlatFile.open(seq_file).each_entry do |f|
    seq_objs[f.definition] = f
  end
  return(seq_objs)
end


############################################################
infile = nil
blast_file_sl = nil
evalue_cutoff = 1e-100
identity_cutoff = 95.0
coverage_cutoff = 0
seq_file = nil
subject_seq_file = nil
query_start_max = 10

seq_objs = Hash.new
subject_seq_objs = Hash.new
blast_info = Hash.new{|h,k|h[k]=Hash.new{|h,k|}}
evalue_info = Hash.new{|h,k|h[k]=Hash.new{|h,k|}}


############################################################
opts = GetoptLong.new(
  ["-i", GetoptLong::REQUIRED_ARGUMENT],
  ["--blast8_sl", GetoptLong::REQUIRED_ARGUMENT],
  ["--identity", GetoptLong::REQUIRED_ARGUMENT],
  ["--cov", "--coverage", GetoptLong::REQUIRED_ARGUMENT],
  ["-e", "--evalue", GetoptLong::REQUIRED_ARGUMENT],
  ["--seq", "--seq_file", GetoptLong::REQUIRED_ARGUMENT],
  ["--subject", "--subject_file", GetoptLong::REQUIRED_ARGUMENT],
  ["--query_start_max", GetoptLong::REQUIRED_ARGUMENT],
)

opts.each do |opt, value|
  case opt
    when '-i'
      infile = value
    when '--blast8_sl'
      blast_sl_file = value
    when '--identity'
      identity_cutoff = value.to_f
      if identity_cutoff < 1
        identity_cutoff *= 100
      end
    when '--cov', '--coverage'
      coverage_cutoff = value.to_f
    when "-e", "--evalue"
      evalue_cutoff = value.to_f
    when "--seq", "--seq_file"
      seq_file = value
    when "--subject", "--subject_file"
      subject_seq_file = value
    when "--query_start_max"
      query_start_max = value.to_i
  end
end


############################################################
seq_objs = get_seq_objs(seq_file)

subject_seq_objs = get_seq_objs(subject_seq_file)


File.open(infile, "r").each_line do |line|
  line.chomp!
  #Smin_SL	Smin_SL	100.00	1550	0	0	18	1567	1	1550	0.0	 2863
  line_arr = line.split("\t")
  query, subject, identity, query_start, query_stop, subject_start, subject_stop, evalue = line_arr.values_at(0,1,2,6,7,8,9,10)
  identity = identity.to_f
  query_start = query_start.to_i
  query_stop = query_stop.to_i
  subject_start = subject_start.to_i
  subject_stop = subject_stop.to_i
  evalue = evalue.to_f

  next if evalue > evalue_cutoff
  next if identity < identity_cutoff
  next if not seq_objs.include?(query)
  coverage = (query_stop - query_start + 1).abs.to_f/seq_objs[query].seq.length
  next if coverage < coverage_cutoff
  next if query_start > query_start_max
  blast_info[subject][query] = [subject_start, subject_stop]
  evalue_info[subject][query] = evalue

  #scaffold5.1|size591573	.	gene	1325	20615	.	-	.	ID=symbB.v1.2.000001
end


############################################################
blast_info.each_pair do |est, v|
  query = v.keys.sort_by{|i|evalue_info[est][i]}[0]
  #next if v.size >= 2
  strand = "+"
  subject_start, subject_stop = v[query][0, 2].map{|i|i.to_i}
  if subject_start < subject_stop
    #new_start = [subject_start-300, 1].max; new_stop = subject_start - 1
    new_start = 1; new_stop = 50
  else
    #new_stop = [subject_start+300, subject_seq_objs[est].seq.size].min; new_start = subject_start + 1
    new_stop = subject_seq_objs[est].seq.size; new_start = new_stop - 49
    strand = "-"
  end

  puts [est, '.', 'EST', new_start, new_stop, '.', strand, '.', 'ID='+query].map{|i|i.to_s}.join("\t")
end


