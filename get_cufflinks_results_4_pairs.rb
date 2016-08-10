#! /bin/env ruby

require 'getoptlong'

require 'Dir'


BEGIN{
  file_name=__FILE__
  $: << [File.dirname(file_name),'lib'].join('/')
}
require 'read_infiles_file'


##############################################################################
def read_cufflinks_results(cufflinks_file, cufflinks_gene_field=0, fpkms={})
  File.open(cufflinks_file,'r').each_line do |line|
    next if $. == 1
    line.chomp!
    tracking_id, fpkm = line.split("\t").values_at(cufflinks_gene_field-1,9)
    fpkm = fpkm.to_f
    fpkms[tracking_id] << fpkm
    #fpkms[gene_id]=fpkm
  end
  return(fpkms)
end


def get_genes_from_list(gene_list, allele_RegExp)
  genes_included=Hash.new
  File.open(gene_list,'r').each_line do |line|
    line.chomp!
    line.sub!(/#{allele_RegExp}/,'')
    genes_included[line]=''
  end
  return genes_included
end


def output_fpkm_4_each_gene(fpkmsh, outfile, genes_included=Array.new)
  fh=File.open(outfile,'w')
  fpkmsh.each_pair do |gene, fpkm|
    if ! genes_included.empty?
      next if not genes_included.include?(gene)
    end
    fh.puts [gene, fpkm.to_s].join("\t")
  end
  fh.close
end


def get_genes_with_SL_info(pairs_file, genes_included, gene_sep="-", gene_pair_file_line_sep="\t", allele_RegExp)
  genes=Hash.new{|h,k|h[k]=Array.new}
  File.open(pairs_file,'r').each_line do |line|
    sl_containing_gene, num_of_SL_containing_genes, pair = nil, nil, nil
    line.chomp!
    line_array = line.split(gene_pair_file_line_sep)

    if line_array.size == 3
      sl_containing_gene, num_of_SL_containing_genes, pair  = line_array
    elsif line_array.size == 2
      sl_containing_gene, sl_non_containing_gene = line_array
      pair = line_array.join(gene_sep)
    else
      puts "Warming! Lines with more than three columns are found!"
      next
    end

    if ! allele_RegExp.nil? and ! sl_containing_gene.nil?
      sl_containing_gene.sub!(/#{allele_RegExp}/, '')
      sl_non_containing_gene.sub!(/#{allele_RegExp}/, '') if ! sl_non_containing_gene.nil?
    end

    next if num_of_SL_containing_genes == 2

    pair.split(gene_sep).each do |gene|
      gene.sub!(/#{allele_RegExp}/, '') if ! allele_RegExp.nil?
      if gene == sl_containing_gene
        if ! genes_included.empty?
          break if not genes_included.include? gene
        end
        genes['with_SL'].push gene
      else
        genes['without_SL'].push gene
      end
    end
  end
  return genes
end


def get_fpkm_results(genes,fpkms)
  fpkm_results=Hash.new{|h,k|h[k]=Array.new}
  genes['with_SL'].each_with_index do |gene_with_SL,index1|
    gene_without_SL = genes['without_SL'][index1]
    next if not fpkms.include? gene_with_SL
    next if not fpkms.include? gene_without_SL
    pair_name = [gene_with_SL, gene_without_SL].join("-")
    fpkm_results["0"].push [pair_name, fpkms[gene_with_SL],fpkms[gene_without_SL]]
  end
  return fpkm_results
end


def output_pair_fpkm(fpkm_results,pairs_fpkm_outdir)
  fpkm_results.each_pair do |key,value|
    fh = File.open(File.join(pairs_fpkm_outdir,key.to_s), 'w')
    value.each do |array|
      fh.puts array.join("\t")
    end
    fh.close
  end
end


##############################################################################
cufflinks_files=Array.new
pairs_file=nil
genes_lists_included=Array.new
force=false
outdir="./"
gene_pair_file_line_sep = "\t"
gene_sep = "-"
allele_RegExp = nil
cufflinks_gene_field = 0

fpkm_samples=Array.new
genes=Hash.new{|h,k|h[k]=Array.new}
fpkm_results=Hash.new{|h,k|h[k]=Array.new}
genes_included=Hash.new
fpkms = Hash.new{|h,k|h[k]=[]}


##############################################################################
opts=GetoptLong.new(
  ['--cufflinks','--cufflinks_results',GetoptLong::REQUIRED_ARGUMENT],
  ["--read_infiles", "--read_infile", "--read_cufflinks", GetoptLong::REQUIRED_ARGUMENT],
  ['--cufflinks_gene_field',GetoptLong::REQUIRED_ARGUMENT],
  ['--pairs_file','--pairs',GetoptLong::REQUIRED_ARGUMENT],
  ['--genes_included','--genes_list_included',GetoptLong::REQUIRED_ARGUMENT],
  ['--outdir',GetoptLong::REQUIRED_ARGUMENT],
  ['--gene_pair_file_line_sep',GetoptLong::REQUIRED_ARGUMENT],
  ['--gene_sep',GetoptLong::REQUIRED_ARGUMENT],
  ["--allele_RegExp", "--allele_reg_exp", GetoptLong::REQUIRED_ARGUMENT],
  ['--force',GetoptLong::NO_ARGUMENT],
)

opts.each do |opt,value|
  case opt
    when '--cufflinks', 'cufflinks_results'
      value.split(',').each do |file|
        cufflinks_files.push file
      end
    when '--read_infiles', "--read_infile", "--read_cufflinks"
      cufflinks_files = read_infiles_file(value, cufflinks_files)
    when '--cufflinks_gene_field'
      cufflinks_gene_field = value.to_i
    when '--pairs', '--pairs_file'
      pairs_file=value
    when '--genes_included','--gene_list_included'
      genes_lists_included.push value
    when '--outdir'
      outdir=value
    when '--gene_pair_file_line_sep'
      gene_pair_file_line_sep = value
    when '--gene_sep'
      gene_sep = value
    when '--allele_RegExp', 'allele_reg_exp'
      allele_RegExp = value
    when '--force'
      force=true
  end
end


mkdir_with_force(outdir,force)
fpkm_outdir=File.join([outdir,'fpkm'])
Dir.mkdir(fpkm_outdir)
pairs_fpkm_outdir=File.join([outdir,'pairs_fpkm'])
Dir.mkdir(pairs_fpkm_outdir)

if genes_lists_included.empty?
  raise "--genes_included needs to be specified! Exiting ......"
end


##############################################################################
cufflinks_files.each do |cufflinks_file|
  fpkms = read_cufflinks_results(cufflinks_file, cufflinks_gene_field, fpkms)
end


if ! genes_lists_included.empty?
  genes_lists_included.map{|list_file|genes_included.merge! get_genes_from_list(list_file, allele_RegExp)}
end


fpkms.each do |gene, v|
  fpkms[gene] = v.inject{|sum,i|sum+i}/v.size
end

output_fpkm_4_each_gene(fpkms, fpkm_outdir+'/0')

output_fpkm_4_each_gene(fpkms, fpkm_outdir+'/1', genes_included)

genes = get_genes_with_SL_info(pairs_file, genes_included, gene_sep, gene_pair_file_line_sep, allele_RegExp)

fpkm_results = get_fpkm_results(genes, fpkms)

output_pair_fpkm(fpkm_results, pairs_fpkm_outdir)


