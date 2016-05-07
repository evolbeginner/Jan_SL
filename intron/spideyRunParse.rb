#! /bin/env ruby


require "getoptlong"


###################################################################################
SPIDEY = "~/software/sequence_analysis/spidey.linux"

gene_file = nil
mRNA_file = nil
gene_id = nil


###################################################################################
opts = GetoptLong.new(
  ["--gene", GetoptLong::REQUIRED_ARGUMENT],
  ["--mRNA", GetoptLong::REQUIRED_ARGUMENT],
  ["--id", "--ID", "--gene_id", "--gene_ID", GetoptLong::REQUIRED_ARGUMENT],
)

opts.each do |opt, value|
  case opt
    when "--gene" 
      gene_file = value
    when "--mRNA"
      mRNA_file = value
    when "--id", "--ID", "--gene_id", "--gene_ID"
      gene_id = value
  end
end


###################################################################################
a = `#{SPIDEY} -i #{gene_file} -m #{mRNA_file} -p 1`
a.lines.each do |line|
  line.chomp!
  if line =~ /Exon \d+: (\d+)[-](\d+)/
    start = $1
    stop = $2
    chr = "chr"
    #1	TAIR10	CDS	1	2016	.	+	0	ID=1;
    puts [chr, ".", "CDS", start, stop, ".", "+", ".",  "ID="+gene_id].join("\t")
  end
end

<<EOF
Exon 1: 1-24 (gen)  1-24 (mRNA)  id 100.0% mismatches 0 gaps 0  splice site (d  a): 1  0
Exon 2: 61-2412 (gen)  25-2376 (mRNA)  id 100.0% mismatches 0 gaps 0  splice site (d  a): 0  0
EOF


