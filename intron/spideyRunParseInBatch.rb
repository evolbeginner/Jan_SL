#! /bin/env ruby

require "getoptlong"
require "bio"


######################################################################
infile = nil
sep = "\t"
gene_file = nil
mRNA_file = nil
outdir = nil


######################################################################
def output_seq(gene, outfile, seq_info)
  out_fh = File.open(outfile, 'w')
  out_fh.puts ">"+gene
  out_fh.puts seq_info[gene]
  out_fh.close
end


def read_seq(infile)
  seq_info = Hash.new{|h,k|h[k]={}}
  Bio::FlatFile.open(infile).each_entry do |f|
    seq_info[f.definition] = f.seq
  end
  return(seq_info)
end


######################################################################
opts = GetoptLong.new(
  ["-i", GetoptLong::REQUIRED_ARGUMENT],
  ["--sep", GetoptLong::REQUIRED_ARGUMENT],
  ["--gene", "--gene_file", GetoptLong::REQUIRED_ARGUMENT],
  ["--mRNA", "--mRNA_file",  GetoptLong::REQUIRED_ARGUMENT],
  ["--outdir",  GetoptLong::REQUIRED_ARGUMENT],
)

opts.each do |opt, value|
  case opt
    when '-i'
      infile  = value
    when '--sep'
      sep = value
    when '--gene', '--gene_file'
      gene_file = value
    when '--mRNA', '--mRNA_file'
      mRNA_file = value
    when '--outdir'
      outdir = value
  end
end


`mkdir -p #{outdir}`
mRNA_outdir = File.join(outdir, "mRNA")
`mkdir -p #{mRNA_outdir}`
gene_outdir = File.join(outdir, "gene")
`mkdir -p #{gene_outdir}`


######################################################################
gene_seq_info = read_seq(gene_file)

mRNA_seq_info = read_seq(mRNA_file)

File.open(infile, 'r').each_line do |line|
  line.chomp!
  genes = line.split(sep)
  pair = genes.join("-")
  genes.each do |gene|
    mRNA_output = File.join([mRNA_outdir, gene+".mRNA.fas"])
    output_seq(gene, mRNA_output, mRNA_seq_info)
    gene_output = File.join([gene_outdir, gene+".gene.fas"])
    output_seq(gene, gene_output, gene_seq_info)
  end
end


