#! /bin/env ruby


require "getoptlong"
require "bio"


##################################################################
infile = nil
queries = Array.new
evalue_cutoff = 1e-5
outdir = nil
sl_seq = nil
sl_truncation_length = 0
is_force = false

evalues = Array.new
q_starts = Array.new
q_ends = Array.new
q_reverse_ends = Array.new
s_starts = Array.new
s_ends = Array.new
lengths = Array.new
q_start_props = Array.new
q_end_props = Array.new
length_props = Array.new


##################################################################
def read_seq(seq_file)
  seq_objs = Hash.new
  Bio::FlatFile.open(seq_file).each_entry do |f|
    seq_objs[f.definition] = f
  end
  return(seq_objs)
end


def output_array(list, outfile)
  out_fh = File.open(outfile, 'w')
  list.each do |i|
    out_fh.puts i
  end
  out_fh.close
end


##################################################################
opts = GetoptLong.new(
  ["-i", GetoptLong::REQUIRED_ARGUMENT],
  ["-e", GetoptLong::REQUIRED_ARGUMENT],
  ["--query", GetoptLong::REQUIRED_ARGUMENT],
  ["--outdir", GetoptLong::REQUIRED_ARGUMENT],
  ["--sl", GetoptLong::REQUIRED_ARGUMENT],
  ["--sl_truncation_length", GetoptLong::REQUIRED_ARGUMENT],
  ["--force", GetoptLong::NO_ARGUMENT],
)


##################################################################
opts.each do |opt, value|
  case opt
    when "-i"
      infile = value
    when "-e"
      evalue_cutoff = value.to_f
    when "--query"
      value.split(",").each do |i|
        queries << i
      end
    when "--outdir"
      outdir = value
    when "--sl"
      sl_seq = value
    when "--sl_truncation_length"
      sl_truncation_length = value.to_i
    when "--force"
      is_force = true
  end
end


raise "infile has to be given. Exiting ......" if infile.nil?
if is_force
  `rm -rf #{outdir}`
end
`mkdir -p #{outdir}`


sl_seq_objs = read_seq(sl_seq)


##################################################################
File.open(infile, "r").each_line do |line|
  line.chomp!
  #Smin_SL	symbB.v1.2.000037	100.00	21	0	0	2	22	223	243	5e-09	39.9
  line_arr = line.split("\t")
  query, aln_length, q_start, q_end, s_start, s_end, evalue = line_arr.values_at(0,3,6,7,8,9,10)
  q_reverse_end = sl_seq_objs[query].seq.size - q_end.to_i + 1
  
  (next if not queries.include?(query)) unless queries.empty?

  if evalue.to_f > evalue_cutoff
    next
  end

  evalues << evalue
  q_starts << q_start
  q_ends << q_end
  q_reverse_ends << q_reverse_end
  s_starts << (300-s_start.to_i+1)
  s_ends << (300-s_end.to_i+1)
  lengths << (q_end.to_i - q_start.to_i + 1).to_s
  q_start_props << q_start.to_f/(sl_seq_objs[query].seq.size - sl_truncation_length).to_f
  q_end_props << q_end.to_f/sl_seq_objs[query].seq.size.to_f
  length_props << (q_end.to_i - q_start.to_i + 1).to_f/(sl_seq_objs[query].seq.size - sl_truncation_length).to_f
end


##################################################################
q_start_outfile = File.join([outdir, "q_start"])
q_end_outfile = File.join([outdir, "q_end"])
q_reverse_end_outfile = File.join([outdir, "q_reverse_end"])
s_start_outfile = File.join([outdir, "s_start"])
s_end_outfile = File.join([outdir, "s_end"])
length_outfile = File.join([outdir, "length"])
q_start_prop_outfile = File.join([outdir, "q_start_prop"])
q_end_prop_outfile = File.join([outdir, "q_end_prop"])
length_prop_outfile = File.join([outdir, "length_prop"])


output_array(q_starts, q_start_outfile)
output_array(q_ends, q_end_outfile)
output_array(q_reverse_ends, q_reverse_end_outfile)
output_array(s_starts, s_start_outfile)
output_array(s_ends, s_end_outfile)
output_array(lengths, length_outfile)
output_array(q_start_props, q_start_prop_outfile)
output_array(q_end_props, q_end_prop_outfile)
output_array(length_props, length_prop_outfile)


