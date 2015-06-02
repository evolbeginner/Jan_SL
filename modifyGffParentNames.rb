#! /bin/env ruby

corenames=Array.new
coordinates=Hash.new
gene_loci=Hash.new

################################################################
class String
  def corename()
    gene = self
    return(gene.sub!(/\.[^.]+$/,''))
  end
end

################################################################
in_gff = ARGV[0]
File.open(in_gff, 'r').each_line do |line|
  line.chomp!
  next if line =~ /^#/
  chr, start, stop, attribute = line.split("\t").values_at(3,4,8)
  #coordinates[chr + '|' + [,start,stop].join("-")] = ''
  attribute =~ /(.+=)(\.[^.]+$)/
  corenames.push $2
  gene_loci[$2.corename] = Hash.new if not gene_loci.include?($2.corename)
  gene_loci[$2.corename][chr]=''
  if $1.corename != corenames[-1]
    if gene_loci.include?($1.corename)
      if ! gene_loci[$1.corename].include?(chr)
        gene_loci[$1.corename]=
      end
    end
  end
  
  if  
end

