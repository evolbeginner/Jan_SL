#! /bin/env ruby

BEGIN{
  file_name=__FILE__
  $: << File.join([File.dirname(file_name),'lib'])
}

require 'getoptlong'
require 'bio'
require 'read_seq_objs'
require 'Dir'


####################################################################
def read_gene_pairs(gene_pair_file, pair_sep='-')
  gene_pairs=Array.new
  File.open(gene_pair_file,'r').each_line do |line|
    line.chomp!
    line_array = line.split("\t")
    target_gene = line_array[0]
    paralog = line_array[2].split(pair_sep).select{ |i|i != target_gene }[0]
    gene_pairs.push([target_gene,paralog])
  end
  return gene_pairs
end


def find_single_with_multi_gene_pairs(gene_pairs, gff_info)
  return(gene_pairs.each.select{|target_gene, paralog| gff_info[paralog]['starts'].size >= 3})
end


def find_spanning_exons(single_with_multi_gene_pairs, gff_info, cds_seq_objs, outdir)
  tmp_aln_dir = File.join(outdir,'tmp_aln')
  aln_dir = File.join(outdir,'aln')
  Dir.mkdir(tmp_aln_dir)
  Dir.mkdir(aln_dir)
  single_with_multi_gene_pairs.each do |target,paralog|
    tmp_seq_fas = File.join(tmp_aln_dir,'tmp.fas')
    seq_aln = File.join(aln_dir,[target,paralog].join("-")) + '.aln'
    fh = File.open(tmp_seq_fas,'w')
    fh.puts '>' + target + "\n" + cds_seq_objs[target].seq
    fh.puts '>' + paralog+ "\n" + cds_seq_objs[paralog].seq
    fh.close
    `muscle -in #{tmp_seq_fas} -out #{seq_aln} -quiet 2>/dev/null`
    parse_aln_file(seq_aln,gff_info)
  end
end


def parse_aln_file(aln_file, gff_info)
  alns=Array.new # alns[0]:target, alsn[1]:paralog
  Bio::FlatFile.open(aln_file).each_entry do |f|
    alns.push f
  end
  # match_start corresponds to the index in the array
  match_start = alns[0].seq =~ /([^-] .+ [^-]) (?:[-]+|$)/x
  alns[0].seq =~ /(?:^|[-]+) (.+ [^-]) (?:[-]+|$)/x
  match_length=$1.size
  # for aln[1]
  aln_1_subseqs=Hash.new
  aln_1_subseqs_length=Hash.new
  aln_1_subseqs['matched_seq']      = alns[1].seq[match_start,match_length]
  aln_1_subseqs['pre_matched_seq']  = alns[1].seq[0,match_start]
  aln_1_subseqs.each_pair do |key,value|
    aln_1_subseqs_length[key] = value.gsub(/[-]/,'').size
  end
  matched_start, matched_stop = aln_1_subseqs_length['pre_matched_seq'], aln_1_subseqs_length.values.reduce(:+)
  matched_range = matched_start..matched_stop
  counter = 0
  gff_info[alns[1].definition]['starts'].each do |exon_start|
    exon_start = (exon_start-1)/3+1
    #print exon_start; print "\t"
    counter+=1 if matched_range.include?(exon_start)
    break if exon_start > matched_stop
  end
  puts [alns[0].definition, alns[1].definition, counter.to_s].join("\t") if counter >= 4
end


def get_gene_list_by_num_of_exons(exon_info_file,max_num_of_exons=1,exon_info_file_format='list',fields=[1,4],separator="\t")
  gene_list = Hash.new
  num_of_exons_hash=Hash.new
  fh=File.open(exon_info_file, 'r')
  fh.each_line do |line|
    line.chomp!
    gene,num_of_exons = line.split(separator).values_at(fields[0]-1,fields[1]-1)
    num_of_exons=num_of_exons.to_i
    if num_of_exons <= max_num_of_exons
      gene_list[gene] = ''
    end
  end
  fh.close
  return gene_list
end


def find_pairs(orthomcl_output,blast_output,prefix,suffix,seq_file,gene_list,outfile,
                seq_file_suffix='\|.+',blast_args={'e_value'=>1e-10,'coverage'=>0.3})
  dirname = File.dirname($0)
  get_duplicate_pairs_based_on_PlantCell = File.join([dirname, "get_duplicate_pairs_based_on_PlantCell.rb"])
  blast_args_content = %w[e_value coverage].map{|i|blast_args[i]}.join(",")
  execute = "time ruby2.1 #{get_duplicate_pairs_based_on_PlantCell} --orthomcl #{orthomcl_output} --blast #{blast_output} --seq #{seq_file} --prefix #{prefix} --seq_file_suffix #{seq_file_suffix} --blast_args #{blast_args_content} --gene_list #{gene_list} --o2 #{outfile}"
  puts "\nCommand:"
  puts execute
  puts
  system("#{execute} > /dev/null")
end


def get_gff_info(gff=nil, features=[], attributes=[])
  gff_info=Hash.new
  gff=gff
  posi_info=Hash.new

  fh=File.open(gff,'r')
  while(line=fh.gets) do
    line.chomp!
    # note that start and stop here are string.
    seqid, feature, start, stop, strand, attribute =line.split("\t").values_at(0,2,3,4,6,8)
    next if ! features.include? feature
    attribute.scan(%r{([^=;]+)=([^=;]+)}).each do |item|
      item =~ /([^=;]+)=([^=;]+)/
      next if ! attributes.include? item[0]
      gff_info[$2]=Hash.new if not gff_info.include? $2
      if ! gff_info[$2].include?('starts')
        gff_info[$2]['starts']=Array.new
        gff_info[$2]['lengths']=Array.new
        gff_info[$2]['strand']=strand
      end
      gff_info[$2]['starts'].push(start.to_i)
      gff_info[$2]['lengths'].push((stop.to_i-start.to_i).abs)
    end
  end
  fh.close
  return gff_info
end


def show_help()
  basename=File.basename($0)
  puts "Usage of #{basename}: ruby #{basename}"
  puts <<EOF
  <--cds_file|seq_file>
  <--gff_file>
  <--exon_info_file>
EOF
  exit
end


####################################################################
gff_file=nil
cds_file=nil
gff_info=Hash.new
exon_info_file=nil
exon_info_file_format='list'
features=Array.new
attributes=Array.new
max_num_of_exons=1
gene_prefix="''"
gene_suffix="''"
seq_file_suffix="''"
orthomcl_file=nil
blast_file=nil
blast_seq_file=nil
blast_args=Hash.new
outdir="./"
force=false
gene_pair_outfile=nil
gene_pair_infile=nil
pair_sep='-'

gene_list=Hash.new
max_num_of_exons=1
singletons=Hash.new
cds_seq_objs=Hash.new


################################
show_help if not ARGV[0]

opts=GetoptLong.new(
  ['--cds_file',GetoptLong::REQUIRED_ARGUMENT],
  ['--gff_file',GetoptLong::REQUIRED_ARGUMENT],
  ['--exon_info_file',GetoptLong::REQUIRED_ARGUMENT],
  ['--exon_info_file_format',GetoptLong::REQUIRED_ARGUMENT],
  ['--features',GetoptLong::REQUIRED_ARGUMENT],
  ['--attributes',GetoptLong::REQUIRED_ARGUMENT],
  ['--max_num_of_exons',GetoptLong::REQUIRED_ARGUMENT],
  ['--orthomcl_file','--orthomcl',GetoptLong::REQUIRED_ARGUMENT],
  ['--blast_file','--blast',GetoptLong::REQUIRED_ARGUMENT],
  ['--blast_seq_file',GetoptLong::REQUIRED_ARGUMENT],
  ['--prefix',GetoptLong::REQUIRED_ARGUMENT],
  ['--suffix',GetoptLong::REQUIRED_ARGUMENT],
  ['--blast_args',GetoptLong::REQUIRED_ARGUMENT],
  ['--seq_file_suffix',GetoptLong::REQUIRED_ARGUMENT],
  ['--outdir',GetoptLong::REQUIRED_ARGUMENT],
  ['--force',GetoptLong::NO_ARGUMENT],
  ['--gene_pair_infile',GetoptLong::REQUIRED_ARGUMENT],
  ['--pair_sep',GetoptLong::REQUIRED_ARGUMENT],
)

opts.each do |opt,value|
  case opt
    when '--cds_file', '--seq_file'
      cds_file=value
    when '--gff_file'
      gff_file=value
    when '--exon_info_file'
      exon_info_file=value
    when '--exon_info_file_format'
      exon_info_file_format=value
    when '--features'
      features.split(',') do |feature|
        features.push(feature)
      end
    when '--attributes'
      features.split(',') do |feature|
        attributes.push(feature)
      end
    when '--max_num_of_exons'
      max_num_of_exons=value.to_i
    when '--orthomcl_file', '--orthomcl'
      orthomcl_file=value
    when '--blast_file', '--blast'
      blast_file=value
    when '--blast_seq_file'
      blast_seq_file=value
    when '--prefix'
      gene_prefix = value
      gene_prefix = "'" + value + "'"
    when '--suffix'
      gene_suffix=value 
    when '--seq_file_suffix'
      seq_file_suffix=value
    when '--blast_args'
      value.split(/,/).each do |arg| # arg e.g. e_value:0.001
        k,v = arg.split(':')
        blast_args[k.to_sym] = v.to_f
      end
    when '--outdir'
      outdir=value
    when '--force'
      force=true
    when '--gene_pair_infile'
      gene_pair_infile=value
    when '--pair_sep'
      pair_sep=value
    when '-h', '--help'
      show_help()
  end
end

if features.empty?
  features=['CDS']
end
if attributes.empty?
  attributes=['Parent']
end

mkdir_with_force(outdir,force)
gene_list_outfile = File.join([outdir,"gene_list"])
gene_pair_outfile=File.join([outdir,'gene_pairs.list'])


####################################################################
if not (exon_info_file) or (not cds_file)
  raise "exon_info_file or cds_file not given!"
end

gene_list = get_gene_list_by_num_of_exons(exon_info_file,max_num_of_exons,'',[1,4],"\t")

gff_info.merge! get_gff_info(gff_file, features, attributes)

gff_info.each_pair do |gene,value|
  # p value['starts'] if gene == 'symbB.v1.2.000296.t1'
  value['starts'].reverse! if value['strand'] == '-'
  first = value['starts'][0]
  new_starts  = [1]
  new_start   = 1
  value['lengths'].each_with_index do |length|
    new_start += length
    new_starts.push new_start
  end
  value['starts'] = new_starts
  value['starts'].pop
  singletons[gene] = '' if value['starts'].size <= max_num_of_exons
end

fh=File.open(gene_list_outfile,'w')
singletons.each_key do |i|
  fh.puts i
end
fh.close

if not gene_pair_infile
  find_pairs(orthomcl_file,blast_file,gene_prefix,gene_suffix,blast_seq_file,gene_list_outfile,gene_pair_outfile,seq_file_suffix="'\\|.+'",blast_args)
  gene_pair_infile = gene_pair_outfile
end

cds_seq_objs = read_seq_objs(cds_file,nil,seq_file_suffix)

gene_pairs = read_gene_pairs(gene_pair_infile, pair_sep)

single_with_multi_gene_pairs = find_single_with_multi_gene_pairs(gene_pairs,gff_info)

find_spanning_exons(single_with_multi_gene_pairs,gff_info,cds_seq_objs,outdir)


