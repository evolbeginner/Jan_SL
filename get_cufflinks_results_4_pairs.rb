#! /bin/env ruby

require 'getoptlong'
require 'Dir'


##############################################################################
def read_cufflinks_results(cufflinks_file)
  fpkms=Hash.new
  File.open(cufflinks_file,'r').each_line do |line|
    next if $. == 1
    line.chomp!
    tracking_id,gene_id,fpkm = line.split("\t").values_at(0,3,9)
    fpkm=fpkm.to_f
    fpkms[tracking_id]=fpkm
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


def output_fpkm_4_each_gene(fpkm_samples,output_dir="./")
  fpkm_samples.each_with_index do |fpkm,index|
    fh=File.open(File.join([output_dir,index.to_s]),'w')
    fpkm.each_pair do |gene,value|
      fh.puts [value.to_s,gene].join("\t") if value != 0
    end
    fh.close
  end
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


def get_fpkm_results(genes,fpkm_samples)
  fpkm_results=Hash.new{|h,k|h[k]=Array.new}
  genes['with_SL'].each_with_index do |gene_with_SL,index1|
    gene_without_SL = genes['without_SL'][index1]
    fpkm_samples.each_with_index do |fpkm_sample,index2|
      next if not fpkm_sample.include? gene_with_SL
      next if not fpkm_sample.include? gene_without_SL
      fpkm_results[index2].push [fpkm_sample[gene_with_SL],fpkm_sample[gene_without_SL]]
    end
  end
  return fpkm_results
end


def output_fpkm_results(fpkm_results,output_dir)
  fpkm_results.keys.each do |key|
    fh = File.open(File.join([output_dir,key.to_s]),'w')
    fpkm_results[key].each do |pair_array|
      k_0=0
      pair_array.each do |i|
        k_0+=1 if i!=0
      end
      puts pair_array.join("\t") if k_0>=1
    end
    fh.close
  end
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

fpkm_samples=Array.new
genes=Hash.new{|h,k|h[k]=Array.new}
fpkm_results=Hash.new{|h,k|h[k]=Array.new}
genes_included=Hash.new


opts=GetoptLong.new(
  ['--cufflinks','--cufflinks_results',GetoptLong::REQUIRED_ARGUMENT],
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
  fpkm_samples.push read_cufflinks_results(cufflinks_file)
end

if ! genes_lists_included.empty?
  genes_lists_included.map{|list_file|genes_included.merge! get_genes_from_list(list_file, allele_RegExp)}
end

output_fpkm_4_each_gene(fpkm_samples, fpkm_outdir)

genes = get_genes_with_SL_info(pairs_file, genes_included, gene_sep, gene_pair_file_line_sep, allele_RegExp)

fpkm_results = get_fpkm_results(genes, fpkm_samples)

output_pair_fpkm(fpkm_results, pairs_fpkm_outdir)


