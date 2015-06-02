#! /bin/env ruby

require 'String'
require 'bio'
require 'getoptlong'

seqfile=nil
query_VS_EST_file=nil
seqs={}
distance_from_query_CDS_start={}
distance_cutoff=0

########################################################
def read_SeqFile(seqfile)
  seqs_info=Hash.new{|h,k|h[k]=Hash.new}
  fh = Bio::FlatFile.open(seqfile)
  fh.each_entry do |f|
    seqs_info[f.definition]['length']=f.seq.size
  end
  return(seqs_info)
end

def read_query_VS_EST_file(query_VS_EST_file,seqs_info)
  distance_from_query_CDS_start=Hash.new()
  fh = File.open(query_VS_EST_file,'r')
  fh.each_line do |line|
    line.chomp!
    query,subject,direction,query_start,subject_start = line.split("\t").map{|i| i.numeric? ? i.to_i : i}
    pair = [query,subject].join('#')
    distance_from_SubjectStart = direction=='+' ? subject_start : seqs_info[subject]['length']-subject_start
    distance_from_query_CDS_start[pair] = distance_from_SubjectStart-(query_start-1)
  end
  return(distance_from_query_CDS_start)
end

########################################################
opts=GetoptLong.new(
  ['--seq','--seqfile',GetoptLong::REQUIRED_ARGUMENT],
  ['--in','--infile',GetoptLong::REQUIRED_ARGUMENT],
  ['--distance','--distance_cutoff',GetoptLong::REQUIRED_ARGUMENT],
)

opts.each do |opt,value|
  case opt
    when '--seq','--seqfile'
      seqfile=value
    when '--in','--infile'
      query_VS_EST_file=value
    when '--distance', 'distance_cutoff'
      distance_cutoff=value.to_i
  end
end

########################################################
seqs_info = read_SeqFile(seqfile)

distance_from_query_CDS_start = read_query_VS_EST_file(query_VS_EST_file,seqs_info)

distance_from_query_CDS_start.delete_if{|k,v| v < distance_cutoff}
distance_from_query_CDS_start.each_pair do |k,v|
  puts [k,v].join("\t")
end


