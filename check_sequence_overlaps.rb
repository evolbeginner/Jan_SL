#! /bin/env ruby

require 'bio'
require 'getoptlong'

seq_files=Array.new
seqs=Hash.new

####################################################
def get_seqs_info(seq_files)
  puts 'reading seqs ......'
  seqs=Hash.new
  seq_files.each_with_index do |seq_file,index|
    puts seq_file
    seqs[index]=Hash.new{|h,k|h[k]=Array.new}
    fh = Bio::FlatFile.open(File.expand_path(seq_file))
    fh.each_entry do |f|
      seqs[index][f.seq].push f.definition
      #seqs[index][f.seq] = 1
    end
  end
  return seqs
end

####################################################
opts = GetoptLong.new(
  ['--seq',GetoptLong::REQUIRED_ARGUMENT],
)

opts.each do |opt,value|
  case opt
    when '--seq'
      seq_files.push value
  end
end

####################################################
seqs = get_seqs_info seq_files[0,2]

# overlap_seqs = seqs[1].keys & seqs[0].keys

seqs[0].each do |seq,title|
  seqs[1].to_a.map{|k,v|k.to_s.include?(seq) and puts [title,v].join("\t")}
end

puts overlap_seqs.size

