#! /bin/env ruby

require "getoptlong"
require "SSW_bio"


#############################################################################
dirname = File::dirname($0)
pal2nal = "pal2nal.pl"
intron_align = File.join(dirname, "intron_align.rb")


#############################################################################
infile = nil
outdir = nil
gff = nil
aligner = "muscle"
sep = "\t"
cds_seq = nil
pep_seq = nil
pairs = Hash.new


#############################################################################
opts = GetoptLong.new(
  ['-i', GetoptLong::REQUIRED_ARGUMENT],
  ['--outdir', GetoptLong::REQUIRED_ARGUMENT],
  ['--gff', GetoptLong::REQUIRED_ARGUMENT],
  ['--cds', '--cds_seq', GetoptLong::REQUIRED_ARGUMENT],
  ['--pep', '--pep_seq', GetoptLong::REQUIRED_ARGUMENT],
  ['--sep', GetoptLong::REQUIRED_ARGUMENT],
  ['--aligner', GetoptLong::REQUIRED_ARGUMENT],
)


opts.each do |opt, value|
  case opt
    when "-i"
      infile = value
    when "--outdir"
      outdir = value
    when "--gff"
      gff = value
    when "--sep"
      sep = value
    when "--aligner"
      aligner = value
    when "--cds", "--cds_seq"
      cds_seq = value
    when "--pep", "--pep_seq"
      pep_seq = value
  end
end


raise "outdir has to be given! Exting ......" if outdir.nil?
`mkdir -p #{outdir}`


#############################################################################
cds_info = read_seq_file(cds_seq)
pep_info = read_seq_file(pep_seq)


File.open(infile, "r").each_line do |line|
  line.chomp!
  line_arr = line.split(sep).values_at(0,1)
  pairs[line_arr] = ""
end


#############################################################################
pairs.each_key do |genes|
  puts genes.join("-")
  pep_fasta = File.join(outdir, genes.join("-")+".pep.fas")
  cds_fasta = File.join(outdir, genes.join("-")+".cds.fas")
  pep_aln = File.join(outdir, genes.join("-")+".pep.aln")
  paml = File.join(outdir, genes.join("-")+".paml")

  out_fh_pep = File.open(pep_fasta, 'w')
  out_fh_cds = File.open(cds_fasta, 'w')
  genes.each do |gene|
    out_fh_pep.puts pep_info[gene]
    out_fh_cds.puts cds_info[gene]
  end
  out_fh_pep.close
  out_fh_cds.close
  `#{aligner} -in #{pep_fasta} -out #{pep_aln} -quiet 2>/dev/null`
  #puts "ruby #{intron_align} -i #{pep_aln} --gff #{gff}"
  `#{pal2nal} #{pep_aln} #{cds_fasta} -output paml > #{paml}`
end


