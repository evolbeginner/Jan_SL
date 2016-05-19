#! /bin/env ruby


###################################################################
BEGIN{
  file_name=__FILE__
  $: << File.join([File.dirname(file_name),'lib'])
}


###################################################################
require "getoptlong"
require 'bio'

require 'detect_mini_repeat'


###################################################################
def parse_bl2seq(bl2seq_str, evalues)
  return(evalues) if bl2seq_str.nil? or bl2seq_str == ""
  bl2seq_arr = bl2seq_str.split("\n")
  bl2seq_arr.each do |line|
    next if line =~ /^#/
    line_arr = line.split("\t")
    query, subject = line_arr.values_at(0,1)
    evalue = line_arr[-2].to_f
    evalues[query] << evalue
  end
  return(evalues)
end


###################################################################
gff_file = nil
gene_seq_file = nil
attr = "ID"
out_intron_seq_file = nil
out_gff_file = nil
outdir = nil
prefix = nil

gene_seq_info = Hash.new
cds_coordinates = Hash.new{|h,k|h[k]=[]}
intron_coordinates = Hash.new{|h,k|h[k]=[]}
evalues = Hash.new{|h,k|h[k]=[]}


###################################################################
def read_seq(seq_file)
  seq_info = Hash.new
  Bio::FlatFile.open(seq_file).each_entry do |f|
    seq_info[f.definition] = f.seq
  end
  return(seq_info)
end


###################################################################
opts = GetoptLong.new(
  ['-i', '--gff', GetoptLong::REQUIRED_ARGUMENT],
  ['--seq_file', '--seq', '--gene_seq', GetoptLong::REQUIRED_ARGUMENT],
  ['--attr', GetoptLong::REQUIRED_ARGUMENT],
  ['--out_seq', '--out_intron', '--out_intron_seq', GetoptLong::REQUIRED_ARGUMENT],
  ['--out_gff', GetoptLong::REQUIRED_ARGUMENT],
  ['--outdir', GetoptLong::REQUIRED_ARGUMENT],
  ['--prefix', GetoptLong::REQUIRED_ARGUMENT],
)

opts.each do |opt, value|
  case opt
    when '-i', '--gff'
      gff_file = value
    when '--attr'
      attr = value
    when '--seq_file', '--seq', '--gene_seq'
      gene_seq_file = value
    when '--out_seq', '--out_intron', '--out_intron_seq'
      out_intron_seq_file = value
    when '--out_gff'
      out_gff_file = value
    when '--outdir'
      outdir = value
    when '--prefix'
      prefix = value
  end
end


`mkdir -p #{outdir}`

tmp_outdir = File.join(["#{outdir}", "tmp"])
`mkdir -p #{tmp_outdir}`
=begin
Signal.trap("TERM") {
  `rm -rf #{tmp_outdir}`
  exit
}

Signal.trap("INT") {
  `rm -rf #{tmp_outdir}`
  exit
}
=end

if ! prefix.nil?
  out_intron_seq_file = File.join([outdir, prefix+'.'+"intron.fas"])
  out_gff_file = File.join([outdir, prefix+'.'+"gene_converted.intron.gff"])
  out_intron_NHEJ_file = File.join([outdir, prefix+'.'+"intron_repeat.list"])
end


###################################################################
gene_seq_info = read_seq(gene_seq_file)

File.open(gff_file).each_line do |line|
  line.chomp!
  gene = nil
  line_arr = line.split("\t")
  if line_arr[-1] =~ /#{attr}=([^;]+)/
    gene = $1
  end
  start, stop = line_arr.values_at(3,4)
  cds_coordinates[gene] << [start, stop].map{|i|i.to_i}
end


cds_coordinates.each_pair do |gene, v|
  v.sort_by!{|i|i[0]}
  v.each_with_index do |list, index|
    if index != v.size-1
      next_list = v[index+1]
      start = list[1] + 1
      stop = next_list[0] - 1
      intron_coordinates[gene] << [start, stop]
    end
  end
end


###################################################################
gff_out_fh = File.open(out_gff_file, 'w')
intron_seq_out_fh = File.open(out_intron_seq_file, 'w')
intron_NHEJ_out_fh = File.open(out_intron_NHEJ_file, 'w')


intron_coordinates.each_pair do |gene, v|
  v.each do |start, stop|
    intron_name = gene + ":" + [start,stop].map{|i|i.to_s}.join("-")

    gff_out_fh.puts [gene, ".", "intron", start, stop, ".", "+", ".", attr+"="+gene].map{|i|i.to_s}.join("\t")
    intron_seq_out_fh.puts ">" + intron_name
    gene_seq = gene_seq_info[gene]
    intron_start = start-1
    intron_stop = stop-1
    intron_seq_out_fh.puts gene_seq[intron_start, intron_stop-intron_start]

    prime_5_outfile = tmp_outdir+"/"+"intron_boundary_5prime.fas"
    prime_3_outfile = tmp_outdir+"/"+"intron_boundary_3prime.fas"
    intron_boundary_seq_5prime_out_fh = File.open(prime_5_outfile, 'w')
    intron_boundary_seq_3prime_out_fh = File.open(prime_3_outfile, 'w')
    intron_boundary_seq_5prime_out_fh.puts ">" + intron_name
    intron_boundary_seq_3prime_out_fh.puts ">" + intron_name
    intron_boundary_seq_5prime_out_fh.puts gene_seq[intron_start-11, 20]
    intron_boundary_seq_3prime_out_fh.puts gene_seq[stop-9, 20]
    intron_boundary_seq_5prime_out_fh.close
    intron_boundary_seq_3prime_out_fh.close
    
    a = `bl2seq -i #{prime_5_outfile} -j #{prime_3_outfile} -p blastn -D 1 -e 0.1 -S 1 -F F -W7`
    evalues = parse_bl2seq(a, evalues)
    
    #repeat_finder_obj = Repeat_finder.new
    #repeat_finder_obj.length(5)
    #if repeat_finder_obj.find(gene_seq[intron_start-6, 10], gene_seq[stop-4, 10])
    #  intron_NHEJ_out_fh.puts [gene, intron_name].join("\t")
    #end
  end
end


evalues.each do |intron_name, v|
  gene = (intron_name.split(':'))[0]
  intron_NHEJ_out_fh.puts [gene, intron_name, v.min.to_s].join("\t")
end


gff_out_fh.close
intron_seq_out_fh.close
intron_NHEJ_out_fh.close


