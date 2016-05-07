#! /bin/env ruby2.1

# get best-reciprocal last hits partly based on the criteria defined in Plant Cell paper 'Pollen-Specific Activation of Arabidopsis Retrogenes Is Associated with Global Transcriptional Reprogramming'
# created by Sishuo Wang from the Department of Botany, the University of British Columbia
# E-mail: sishuowang@hotmail.ca
# To see the usage, please type "ruby get_duplicate_pairs_based_on_PlantCell.rb -h"
# The basic steps are as follows
# 1. OrthoMCL result where only two genes are included in a cluster
# 2. Best reciprocal BLAST hits
# 3. For genes in a given a gene list, if a they are not included in gene pairs detected in the first two steps, see whether they have 'good' BLAST hits. If they have, choose the best hits (e.g., those with the lowest e_values) and make them into a pair.


################################################################################
require 'getoptlong'


################################################################################
orthomcl_file=nil
blast_file=nil
gene_list_file=nil
exon_info_file=nil
exon_info_file_format=nil
gene_suffix=nil
gene_prefix=nil
seq_file=nil
seq_file_suffix=nil
outfile1=nil
outfile2=nil
is_corename=false

blast_args=Hash.new
gene_list=Array.new
blast_pairs=Hash.new
pairs=Hash.new{|h,k|h[k]=Array.new}

seq_info=Hash.new
is_paralogs=nil
counters_included={1=>'',2=>''}


################################################################################
class Is_paralogs
  def initialize(blast_args=Hash.new)
    @blast_args=blast_args
    @is_paralogs = lambda do |blast_args,blast_result|
      blast_args.each_pair do |k,v|
        case k
          when :coverage
            return false if blast_result[k] < v
          when :e_value
            return false if blast_result[k] > v
        end
      end
      return true
    end
  end
  def currying()
    @is_paralogs.curry.(@blast_args)
  end
end


################################################################################
def output_result(outfile, genes_InOrNotIn_pairs, counters_included)
  out_fh = outfile.nil? ? STDOUT : File.open(outfile, 'w')
  genes_InOrNotIn_pairs['IN'].each do |gene,v1|
    v1.each do |counter,v2|
      next if ! counters_included.include?(counter)
      out_fh.puts [gene,counter,v2].flatten.join("\t")
    end
  end
end


def read_seq_file(seq_file, seq_file_suffix=nil)
  require 'bio'
  seq_info=Hash.new{|h,k|h[k]={}}
  fh = Bio::FlatFile.open(seq_file)
  fh.each_entry do |f|
    title = f.definition
    title.sub!(/#{seq_file_suffix}/,'') if seq_file_suffix =~ /./
    seq_info[title]['length']=f.seq.length
  end
  return seq_info
end


def read_gene_list_file(gene_list_file=nil)
  gene_list=Hash.new
  fh=File.open(gene_list_file,'r')
  while(line=fh.gets) do
    line.chomp!
    gene_list[line]=''
  end
  return gene_list
end


def determine_gene_suffix(genes, gene_suffix='')
    gene_suffix='' if gene_suffix.nil?
    genes_without_suffix=genes.map{|i| i.sub(/#{gene_suffix}$/,'')}
    return genes_without_suffix[0] == genes_without_suffix[1] ? true : false
end


def get_pair(genes)
  return genes.sort.join('_')
end


def parse_orthomcl_file(infile, gene_prefix=nil, gene_suffix='', is_corename=false)
  orthoMCL_pairs=Array.new
  fh=File.open(infile,'r')
  while(line=fh.gets) do
    line.chomp!
    genes=line.split("\t")
    genes.map!{|i|i.sub(/^#{gene_prefix}/,'')} if gene_prefix
    next if genes.size > 2
    next if determine_gene_suffix(genes,gene_suffix)
    next if genes[0].corename(is_corename) == genes[1].corename(is_corename)
    pair=get_pair genes
    orthoMCL_pairs.push pair
  end
  return orthoMCL_pairs
end


def get_blast_pairs(blast_file, gene_prefix=nil, gene_suffix='', seq_info={}, is_paralogs=nil, gene_list=[], is_corename=false)
  all_pairs_of_a_gene=multi_D_Hash(3)
  fh=File.open(blast_file,'r')
  while(line=fh.gets)
    line.chomp!
    query, subject, aln_length, e_value, bit_score = line.split("\t").values_at(0,1,3,10,11)
    query, subject = [query,subject].map{|i|i.sub(/^#{gene_prefix}/,'')} if gene_prefix
    next if determine_gene_suffix([query,subject], gene_suffix)
    next if query.corename(is_corename) == subject.corename(is_corename)

    e_value=e_value.to_f
    if ! seq_info.empty?
      coverage=(aln_length.to_f/[query,subject].map{|i|seq_info[i]['length']}.min).to_f
    else
      coverage=1
    end
    if ! is_paralogs.nil? then
      next if not is_paralogs.currying.curry.({:e_value=>e_value,:coverage=>coverage})
    end
    all_pairs_of_a_gene[query][subject][:e_value] = e_value.to_f
    all_pairs_of_a_gene[query][subject][:bit_score] = bit_score.to_f
    all_pairs_of_a_gene[query][subject][:coverage] = (aln_length.to_f/[query,subject].map{|i|seq_info[i]['length']}.min).to_f if ! seq_info.empty?
  end
  return all_pairs_of_a_gene
end


def get_best_hits_of_a_gene(blast_pairs)
  best_blast_pairs=Hash.new{|h,k|h[k]=Array.new}
  blast_pairs.each_pair do |gene1,v|
    bit_scores = v.keys.sort_by{|gene2| v[gene2][:bit_score]}.map{|gene2|[gene2,v[gene2][:bit_score]]}.reverse
    best_blast_pairs[gene1].push bit_scores[0][0]
    bit_scores.drop(1).each do |gene, bit_score|
      if bit_score >= bit_scores[0][1] then
        best_blast_pairs[gene1].push gene
      else
        break
      end
    end
  end
  return best_blast_pairs
end


def  get_best_reciprocal_pairs(best_blast_pairs)
  # best_blast_pairs[gene1] = [best_hit1, best_hit2, ...]
  best_pairs=Array.new
  best_blast_pairs.each_pair do |gene1,best_hits|
    best_hits.each do |best_hit|
      if best_blast_pairs.include?(best_hit) and best_blast_pairs[best_hit].include? gene1 then
        pair=get_pair [gene1,best_hit]
        best_pairs.push pair
      end
    end
  end
  return best_pairs
end


def select_pairs_from_list(gene_list,pairs)
  pairs_selected=multi_D_Hash(2)
  genes_InOrNotIn_pairs=Hash.new{|h,k|h[k]=Hash.new}
  pairs.each_key do |type|
    pairs[type].each do |pair; counter|
      genes=Array.new
      counter=0
      pair.split('_').each do |gene|
        if gene_list.include? gene then
          counter+=1
          genes.push gene
          #genes_InOrNotIn_pairs['IN'][gene]=Array.new if genes_InOrNotIn_pairs['IN'][gene].nil?
          #genes_InOrNotIn_pairs['IN'][gene].push pair
        end
      end
      if counter>=1 then
        genes.each do |gene|
          genes_InOrNotIn_pairs['IN'][gene]=Hash.new if genes_InOrNotIn_pairs['IN'][gene].nil?
          genes_InOrNotIn_pairs['IN'][gene][counter]=Array.new if genes_InOrNotIn_pairs['IN'][gene][counter].nil?
          genes_InOrNotIn_pairs['IN'][gene][counter].push pair
          genes_InOrNotIn_pairs['IN'][gene][counter].uniq!
          pairs_selected[type][counter]=Array.new if pairs_selected[type][counter].nil?
          pairs_selected[type][counter].push pair
        end
      end
    end
  end
  genes_InOrNotIn_pairs['NOTIN'] = gene_list.keys - genes_InOrNotIn_pairs['IN'].keys
  return([pairs_selected,genes_InOrNotIn_pairs])
end


def select_pairs_from_multi_blast_hits(blast_pairs, genes_not_in_pairs, gene_list, is_paralogs)
  # genes_not_in_pairs: those not included in pairs
  second_best_pairs=Array.new
  genes_on_list=Hash.new{|h,k|h[k]=Hash.new}

  genes_not_in_pairs.each do |gene|
    best=Hash.new
    pair=nil
    counter=0
    if blast_pairs.include? gene then
      blast_pairs[gene].each_pair do |paralog,v|
        is_go_on=false
        if ! is_paralogs.nil? then
          next if not is_paralogs.currying.curry.(v)
        end
        if best.empty?
          is_go_on=true
        elsif v[:e_value] < best[:e_value]
          is_go_on=true
        elsif v[:e_value] == best[:e_value]
          next if v[:coverage] < best[:coverage]
          is_go_on=true
        end
        if is_go_on then
          pair=[gene,paralog].join('_')
          counter = gene_list.include?(paralog) ? 2 : 1
          best=v
        end
      end
    if ! best.empty? then
      second_best_pairs.push pair
      genes_on_list[gene][counter]=Array.new if genes_on_list[gene][counter].nil?
      genes_on_list[gene][counter].push pair
    end
    end
  end
  return ([second_best_pairs,genes_on_list])
end


def multi_D_Hash(num_of_dimensions)
  if num_of_dimensions > 1 then
    hash=Hash.new{|h,k| h[k]=multi_D_Hash(num_of_dimensions-1)}
    # multi_D_Hash(num_of_dimensions-1)
  elsif num_of_dimensions == 1 then
    Hash.new
  else
    raise "num_of_dimensions error!"
  end
end


def show_help()
  basename = File.basename($0) 
  puts "Usage:\truby #{$0}"
  puts <<EOF
  --orthomcl_file|--orthomcl
  --blast_file|--blast
  --gene_list|--gene_list_file
  --prefix
  --suffix
  --blast_args
  --seq|--seq_file
  --seq_file_suffix
  --o1|--out1|--outfile1|--output1
  --o2|--out2|--outfile2|--output2
  --corename
  -h|--help
EOF
  exit(1)
end


class String
  def corename(is_corename=false)
    gene = self
    if is_corename
      gene.gsub!(/\.[^.]+$/,'')
    end
    return(gene)
  end
end


################################################################################
show_help() if not ARGV[0]

opts = GetoptLong.new(
  ['--orthomcl_file', '--orthomcl', GetoptLong::REQUIRED_ARGUMENT],
  ['--blast_file', '--blast', GetoptLong::REQUIRED_ARGUMENT],
  ['--gene_list', '--gene_list_file', GetoptLong::REQUIRED_ARGUMENT],
  ['--prefix', GetoptLong::REQUIRED_ARGUMENT],
  ['--suffix', GetoptLong::REQUIRED_ARGUMENT],
  ['--blast_args', GetoptLong::REQUIRED_ARGUMENT],
  ['--seq', '--seq_file', GetoptLong::REQUIRED_ARGUMENT],
  ['--seq_file_suffix', GetoptLong::REQUIRED_ARGUMENT],
  ['--o1', '--output1', '--out1', '--outfile1', GetoptLong::REQUIRED_ARGUMENT],
  ['--o2', '--output2', '--out2', '--outfile2', GetoptLong::REQUIRED_ARGUMENT],
  ['--corename', GetoptLong::NO_ARGUMENT],
  ['-h', '--help', GetoptLong::NO_ARGUMENT],
)

opts.each do |opt,value|
  case opt
    when '--orthomcl_file', '--orthomcl'
      orthomcl_file=value
    when '--blast_file', '--blast'
      blast_file=value
    when '--gene_list', '--gene_list_file'
      gene_list_file=value
    when '--prefix'
      gene_prefix=value
    when '--suffix'
      gene_suffix=value 
    when '--seq_file_suffix'
      seq_file_suffix=value
    when '--blast_args'
      value.split(/,/).each do |arg| # arg e.g. e_value:0.001
        k,v = arg.split(':')
        blast_args[k.to_sym] = v.to_f
      end
    when '--seq_file', '--seq'
      seq_file=value
    when '--counters', '--counters_included'
      counters_included[value.to_i]=''
    when '--o1', '--output1', '--out1', '--outfile1'
      outfile1 = value
    when '--o2', '--output2', '--out2', '--outfile2'
      outfile2 = value
    when '--no_corename'
      is_corename=false
    when '--corename'
      is_corename=true
    when '-h', '--help'
      show_help()
  end
end


# check params
if not gene_list_file
  raise "gene_list_file or exon_info_file has to be given"
end


################################################################################
if gene_list_file
  gene_list = read_gene_list_file(gene_list_file)
end
raise "gene list cannot be obtained!" if gene_list.empty?


seq_info = read_seq_file(seq_file,seq_file_suffix) if seq_file
is_paralogs = Is_paralogs.new(blast_args) if blast_args

pairs['orthomcl'] = parse_orthomcl_file(orthomcl_file, gene_prefix, gene_suffix, is_corename)

blast_pairs = get_blast_pairs(blast_file, gene_prefix, gene_suffix, seq_info, is_paralogs, gene_list, is_corename)
best_blast_pairs = get_best_hits_of_a_gene(blast_pairs)
pairs['blast'] = get_best_reciprocal_pairs(best_blast_pairs)


################################################################################
if gene_list_file then
  pairs_selected, genes_InOrNotIn_pairs = select_pairs_from_list(gene_list,pairs)
  output_result(outfile1, genes_InOrNotIn_pairs, counters_included) # output

  pairs['blast2nd'], second_genes_on_list = select_pairs_from_multi_blast_hits(blast_pairs,genes_InOrNotIn_pairs['NOTIN'],gene_list,is_paralogs)
  genes_InOrNotIn_pairs['IN'].merge!(second_genes_on_list)
  genes_InOrNotIn_pairs['NOTIN'] = gene_list.keys - genes_InOrNotIn_pairs['IN'].keys
  output_result(outfile2, genes_InOrNotIn_pairs, counters_included) # output
end


# pairs.map{|type,v|puts v.size}

