#! /bin/env ruby


require "getoptlong"
require "bio"


#####################################################################
def parse_bl2seq(bl2seq_str, evalues)
  return(evalues) if bl2seq_str.nil? or bl2seq_str == ""
  bl2seq_arr = bl2seq_str.split("\n")
  bl2seq_arr.each do |line|
    next if line =~ /^#/
    line_arr = line.split("\t")
    query, subject = line_arr.values_at(0,1)
    evalue = line_arr[-2].to_f
    pair = [query,subject].join("\t")
    evalues[pair] << evalue
  end
  return(evalues)
end


#####################################################################
infile = nil
seq_files = Array.new
outdir = nil

seqs = Hash.new
evalues = Hash.new{|h,k|h[k]=Array.new}


#####################################################################
opts = GetoptLong.new(
  ["-i", GetoptLong::REQUIRED_ARGUMENT],
  ["--seq", GetoptLong::REQUIRED_ARGUMENT],
  ["--outdir", GetoptLong::REQUIRED_ARGUMENT],
)

opts.each do |opt, value|
  case opt
    when '-i'
      infile = value
    when '--seq'
      seq_files << value
    when '--outdir'
      outdir = value
  end
end


`mkdir -p #{outdir}`


Signal.trap("TERM") {
  `rm -rf #{outdir}`
  exit
}

Signal.trap("INT") {
  `rm -rf #{outdir}`
  exit
}


#####################################################################
seq_files.each do |seq_file|
  Bio::FlatFile.open(seq_file).each_entry do |f|
    seqs[f.definition] = f.seq
  end
end


File.open(infile, "r").each_line do |line|
  #143	shared	symbB.v1.2.003826	137	1611-2062	symbB.v1.2.040251	142	3477-4071
  line.chomp!
  next if line !~ /\tshared\t/
  genes = Array.new
  posis = Array.new
  full_names = Array.new
  outfiles = Array.new
  line_arr = line.split("\t")
  genes = line_arr.values_at(2,5)
  posis = line_arr.values_at(4,7)
  genes.each_with_index do |gene, index|
    full_names << [gene, posis[index]].join(':')
  end
  full_names.each_with_index do |full_name, index|
    outfiles << outdir + "/" + index.to_s+".fas"
    fh = File.open(outfiles[index], 'w')
    fh.puts ">" + full_name
    fh.puts seqs[full_name]
    fh.close
  end
  begin
    a = `bl2seq -i #{outfiles[0]} -j #{outfiles[1]} -p blastn -D 1 -e 0.1 -W 7 -S 1 -F F`
    evalues = parse_bl2seq(a, evalues)
  rescue
    ;
  end
end


#####################################################################
evalues.each_pair do |gene, v|
  puts [gene, v.min.to_s].join("\t")
end

`rm -rf #{outdir}`


