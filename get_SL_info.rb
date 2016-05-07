#! /bin/env ruby

require 'getoptlong'

blast_result=nil
exon_counting_file=nil
gff=nil
features=Array.new
attributes=Hash.new
is_output_subject_info=false
cutoff={
  'e_value' =>  1e-10,
  'identity'  =>  0,
  'aln_length'  =>  0,
}
subject_info=Hash.new
exon_info=Hash.new
gff_info=Hash.new
upstream=nil
downstream=nil
is_no_middle=false

###############################################
def get_subject_info(subject_info={},blast_result=nil, cutoff={})
  subject_info=subject_info
  blast_result=blast_result
  cutoff=cutoff
  fh=File.open(blast_result,'r')
  while(line=fh.gets) do
    line.chomp!
    query, subject, identity, aln_length, num_of_mismatch, num_of_gap_opening, query_start, query_end, subject_start, subject_end, e_value, score  =line.split(/\t/)
    e_value=e_value
    identity=identity
    aln_length=aln_length
    next if e_value.to_f > cutoff['e_value']
    next if identity.to_i < cutoff['identity']
    next if aln_length.to_f < cutoff['aln_length']
    subject_info[subject] ||= Hash.new{|hash,key| hash[key]=Array.new()}
    #subject_info[subject]['posi'] ||= Array.new()
    subject_info[subject]['posi'].push [(subject_start..subject_end)]
    subject_info[subject]['identity'].push identity
    subject_info[subject]['score'].push score
  end
  fh.close
  return(subject_info)
end

def output_subject_info(subject_info={}, exon_info={}, gff_info={})
  subject_info.keys.sort{|a,b| [subject_info[a]['posi'].size]<=>[subject_info[b]['posi'].size]}.each do |subject|
  #subject_info.keys.sort{|a,b| subject_info[a]['score'][0].to_f<=>subject_info[b]['score'][0].to_f}.each do |subject|
    next if subject.nil?
    print subject+"\t"+subject_info[subject]['posi'].join("\t")+"\t"
    print subject_info[subject]['score'].join("\t")+"\t"
    if exon_info.include? subject then
      print exon_info[subject]['num_of_exons'].to_s+"\t"
      print exon_info[subject]['num_of_exons']/exon_info[subject]['length'].to_f
    end
    puts
  end
end


def get_gff_info(gff_info={}, gff=nil, features=[], attributes=[], upstream=nil, downstream=nil, is_no_middle=false)
  gff_info=gff_info
  gff=gff
  features=features
  
  posi_info=Hash.new
  fh=File.open(gff,'r')
  while(line=fh.gets) do
    line.chomp!
    seqid, feature, start, stop, strand, attribute =line.split("\t").values_at(0,2,3,4,6,8)
=begin
    if upstream then
      posi_info['upstream']=Hash.new
      posi_info['upstream']['start']=start-upstream
      posi_info['upstream']['stop']=start-1
    end
    if downstream then
      posi_info['downstream']=Hash.new
      posi_info['downstream']['start']=stop
      posi_info['downstream']['stop']=stop+downstream
    end
    if ! is_no_middle then
      posi_info['middle']=Hash.new
      posi_info['middle']['start']=start
      posi_info['middle']['stop']=stop
    end
=end
    attribute.scan(/([^=;]+)=([^=;]+)/).each do |item|
      next if ! attributes.include? $1
      gff_info[$2]=Hash.new
      gff_info[$2]['posi']=posi_info
      gff_info[$2]['strand']=strand
      gff_info[$2]['length']=stop.to_i-start.to_i+1
      #puts gff_info[$1]['length']
    end
  end
  fh.close
  return gff_info
end


def get_exon_info(exon_info={}, exon_counting_file)
  exon_info=exon_info
  fh=File.open(exon_counting_file,'r')
  while(line=fh.gets) do
    line.chomp!
    name, length, num_of_exons = line.split(/\t/).values_at(0,1,3)
    exon_info[name]=Hash.new()
    exon_info[name]['num_of_exons']=num_of_exons.to_i
    exon_info[name]['length']=length.to_i
  end
  fh.close
  return exon_info
end


###############################################
opts = GetoptLong.new(
  ['--blast_result',GetoptLong::REQUIRED_ARGUMENT],
  ['--exon_counting_file',GetoptLong::REQUIRED_ARGUMENT],
  ['-e', '--e_value','--evalue',GetoptLong::REQUIRED_ARGUMENT],
  ['--identity',GetoptLong::REQUIRED_ARGUMENT],
  ['--aln_length',GetoptLong::REQUIRED_ARGUMENT],
  ['--output_subject_info',GetoptLong::NO_ARGUMENT],
  ['--gff',GetoptLong::REQUIRED_ARGUMENT],
  ['--feature',GetoptLong::REQUIRED_ARGUMENT],
  ['--attributes',GetoptLong::REQUIRED_ARGUMENT],
  ['--upstream',GetoptLong::REQUIRED_ARGUMENT],
  ['--downstream',GetoptLong::REQUIRED_ARGUMENT],
  ['--no_middle',GetoptLong::NO_ARGUMENT],
)

opts.each do |opt,value|
  case opt
    when '--blast_result'
      blast_result=value
    when '--gff'
      gff=value
    when '--feature'
      features.push value
    when '--attributes'
      attributes[value]=1
    when '--identity'
      cutoff['identity']=value.to_f
    when '-e', '--e_value', '-evalue'
      cutoff['e_value']=value.to_f
    when '--aln_length'
      cutoff['aln_length']=value.to_i
    when '--exon_counting_file'
      exon_counting_file=value
    when '--output_subject_info'
      is_output_subject_info=true
    when '--upstream'
      upstream=value
    when 'downstream'
      downstream=value
    when 'no_middle'
      is_no_middle=true
    end
end

###############################################
exon_info=get_exon_info(exon_info, exon_counting_file) if ! exon_counting_file.nil?

gff_info = get_gff_info(gff_info, gff, features, attributes, upstream, downstream, is_no_middle) if gff

subject_info = get_subject_info(subject_info, blast_result, cutoff)

output_subject_info(subject_info, exon_info, gff_info) if is_output_subject_info

