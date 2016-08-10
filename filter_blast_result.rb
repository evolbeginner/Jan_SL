#! /bin/env ruby

BEGIN{
  file_name=__FILE__
  $: << [File.dirname(file_name),'lib'].join('/')
}


###################################################################################
require 'getoptlong'
require 'filter_blast_output'
require 'read_seq_objs'


###################################################################################
def generate_regions_range(keys, posi) # keys: %w[query subject]
  regions_range=Hash.new{|h,k|h[k]=Hash.new}
  keys.each do |key|
    regions_range[key]['start'] = posi[0]
    regions_range[key]['stop'] = posi[1]
  end
  return regions_range
end


###################################################################################
indir=nil
outdir=nil
action=nil
e_value_cutoffs=[0,10000]
identity_range=[-1,100]
aligned_length_range=[0,10000]
blast_outfmt=6
force=false
regions_range = Hash.new{|h,k|h[k]=Hash.new}
donor_args=Hash.new
seq_files=Array.new
seq_objs=Hash.new
queries = Array.new


###################################################################################
opts=GetoptLong.new(
  ['--indir',GetoptLong::REQUIRED_ARGUMENT],
  ['--outdir',GetoptLong::REQUIRED_ARGUMENT],
  ['--action',GetoptLong::REQUIRED_ARGUMENT],
  ['--e_value',GetoptLong::REQUIRED_ARGUMENT],
  ['--blast_outfmt',GetoptLong::REQUIRED_ARGUMENT],
  ['--identity_range',GetoptLong::REQUIRED_ARGUMENT],
  ['--aligned_length_range',GetoptLong::REQUIRED_ARGUMENT],
  ['--regions_range',GetoptLong::REQUIRED_ARGUMENT],
  ['--regions_range_query',GetoptLong::REQUIRED_ARGUMENT],
  ['--regions_range_subject',GetoptLong::REQUIRED_ARGUMENT],
  ['--donor_args',GetoptLong::REQUIRED_ARGUMENT],
  ['--seq_files',GetoptLong::REQUIRED_ARGUMENT],
  ['--query', GetoptLong::REQUIRED_ARGUMENT],
  ['--force',GetoptLong::NO_ARGUMENT],
)

opts.each do |opt,value|
  case opt
    when '--indir'
      indir=value
    when '--outdir'
      outdir=value
    when '--action'
      action=value
    when '--e_value'
      e_value_cutoffs=value.split(/[-,]/).map{|i| i.to_f}.sort
    when '--blast_outfmt'
      blast_outfmt=value
    when '--identity_range'
      identity_range=value.split(/[-,]/).map{|i| i.to_f}.sort
    when '--aligned_length_range'
      aligned_length_range=value.split(/[-,]/).map{|i| i.to_i}.sort
    when '--regions_range'
      start,stop=value.split('-').map{|i|i.to_i}
      regions_range.merge! generate_regions_range(%w[query subject],[start,stop])
    when '--regions_range_query'
      start,stop=value.split('-').map{|i|i.to_i}
      regions_range.merge! generate_regions_range(%w[query],[start,stop])
    when '--regions_range_subject'
      start,stop=value.split('-').map{|i|i.to_i}
      regions_range.merge! generate_regions_range(%w[subject],[start,stop])
    when '--donor_args'
      length_range,evalue_range,donor_site = value.split(',').values_at(0,1,2)
      min_len,max_len = length_range.split('-').map{|i|i.to_i}
      donor_args['length_range']=min_len..max_len
      min_evalue,max_evalue = evalue_range.split('-').map{|i|i.to_f}
      donor_args['e_value_range']={'min'=>min_evalue,"max"=>max_evalue}
      donor_args['donor_site']=donor_site
    when '--seq_files'
      value.split(',').each do |file|
        seq_files.push file
      end
    when '--query'
      value.split(',').map{|i|queries << i}
    when '--force'
      force=true
  end
end

if Dir.exists?(outdir)
  if force
    `rm -rf #{outdir}`
    `mkdir -p #{outdir}`
  end
else
  `mkdir -p #{outdir}`
end

if ! seq_files.empty?
  seq_files.each do |seq_file|
    seq_objs.merge! read_seq_objs(seq_file)
  end
end


###################################################################################
Dir.foreach(indir) do |file|
  next if file =~ /^\./
  next if file =~ /^tmp/
  file_FullName=[indir, file].join("/")
  filterBlastOutput_object = FilterBlastOutput.new(file_FullName, blast_outfmt, seq_objs)
  filterBlastOutput_object.filter(e_value_cutoffs, identity_range, aligned_length_range, regions_range, action, outdir, queries, donor_args)
end


