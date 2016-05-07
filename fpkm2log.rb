#! /bin/env ruby

require "getoptlong"
require "basic_math"


###############################################################
infiles = Array.new
outfile = nil
cufflinks_gene_field = 4
include_list_file = nil
pair_file = nil
suffix = nil

fpkms = Hash.new{|h,k|h[k]=[]}
log_transferred_fpkms = Hash.new
genes_included = Hash.new


###############################################################
opts = GetoptLong.new(
  ["-i", "--in", GetoptLong::REQUIRED_ARGUMENT],
  ["-o", "--out", GetoptLong::REQUIRED_ARGUMENT],
  ["--cufflinks_gene_field", GetoptLong::REQUIRED_ARGUMENT],
  ["--include_list", GetoptLong::REQUIRED_ARGUMENT],
  ["--pair", "--pair_file", GetoptLong::REQUIRED_ARGUMENT],
  ["--suffix", "--gene_suffix", GetoptLong::REQUIRED_ARGUMENT],
)

opts.each do |opt, value|
  case opt
    when "-i", "--in"
      value.split(",").map{|i|infiles << i}
    when "-o", "--out"
      outfile = value
    when "--cufflinks_gene_field"
      cufflinks_gene_field = value.to_i
    when "--include_list"
      include_list_file = value
    when "--pair", "--pair_file"
      pair_file = value
    when "--suffix", "--gene_suffix"
      suffix = value
  end
end

outdir = File::dirname(outfile)
`mkdir -p #{outdir}`


###############################################################
if ! include_list_file.nil?
  File.open(include_list_file, "r").each_line do |line|
    line.chomp!
    gene = line.split("\t")[0]
    if not suffix.nil?
      gene.sub!(/#{suffix}/, "")
    end
    genes_included[gene] = ""
  end
end


infiles.each do |infile|
  File.open(infile, "r").each_line do |line|
    line.chomp!
    gene, fpkm = line.split("\t").values_at(cufflinks_gene_field, -1)
    if ! genes_included.empty?
      next if ! genes_included.include?(gene)
    end
    fpkm = fpkm.to_f
    fpkms[gene] << fpkm
  end
end


fpkms.each_pair do |gene, arr|
  next if arr.size != infiles.size
  fpkm_average = arr.inject(:+)/arr.size.to_f
  log_tranferred_fpkm = Math.log2(fpkm_average+0.1)
  log_transferred_fpkms[gene] = log_tranferred_fpkm
end


###############################################################
out_fh = File.open(outfile, "w")
log_transferred_fpkms.each_pair do |gene, fpkm|
  out_fh.puts [gene, fpkm.to_s].join("\t")
end
out_fh.close


