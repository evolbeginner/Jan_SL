#! /bin/env ruby

require 'getoptlong'
require 'bio'

##############################################################
def get_genes_on_list(list_file)
  lines = File.open(list_file,'r').readlines.map{|line|line.chomp!}.reject{|line|line=~/^[#]/}
  return lines
end

def get_genes_on_lists_in_hash(lists)
  acceptor_hash = Hash.new
  lists.each do |list_file|
    get_genes_on_list(list_file).each do |gene|
      acceptor_hash[gene]=''
    end
  end
  return acceptor_hash
end

##############################################################
pair_files=Array.new
pair_files_included=Array.new
pair_files_excluded=Array.new
lists_included=Array.new
lists_excluded=Array.new
separator='_'
field=1

genes_on_lists=Hash.new{|h,k|h[k]={}}
pairs_on_lists=Hash.new{|h,k|h[k]={}}
pairs_on_lists=Hash.new{||}
#pairs_on_lists['included']

opts = GetoptLong.new(
  ['--pair', '--pair_file',GetoptLong::REQUIRED_ARGUMENT],
  ['--list_included', '--lists_included',GetoptLong::REQUIRED_ARGUMENT],
  ['--list_excluded', '--lists_excluded',GetoptLong::REQUIRED_ARGUMENT],
  ['--pair_included', '--pairs_included',GetoptLong::REQUIRED_ARGUMENT],
  ['--pair_excluded', '--pairs_excluded',GetoptLong::REQUIRED_ARGUMENT],
  ['--sep', '--separator',GetoptLong::REQUIRED_ARGUMENT],
  ['-f', '--field',GetoptLong::REQUIRED_ARGUMENT],
)

opts.each do |opt,value|
  case opt
    when '--pair', '--pair_file'
      value.split(',').map{|file|pair_files.push(file)}
    when '--lists_included', '--list_included'
      value.split(',').map{|file|lists_included.push(file)}
    when '--lists_excluded', '--list_excluded'
      value.split(',').map{|file|lists_excluded.push(file)}
    when '--pair_included', '--pairs_included'
      value.split(',').map{|file|pair_files_included.push(file)}
    when '--pair_excluded', '--pairs_excluded'
      value.split(',').map{|file|pair_files_excluded.push(file)}
    when '--sep', '--separator'
      separator=value
    when '-f', '--field'
      field=value.to_i
  end
end


##############################################################
genes_on_lists['included'] = get_genes_on_lists_in_hash(lists_included)
genes_on_lists['excluded'] = get_genes_on_lists_in_hash(lists_excluded)
pairs_on_lists['included'] = get_genes_on_lists_in_hash(pair_files_included)
pairs_on_lists['excluded'] = get_genes_on_lists_in_hash(pair_files_excluded)


pair_files.each do |pair_file|
  File.open(pair_file,'r').each_line do |line|
    line.chomp!
    no_pass=false

    pair = line
    if pairs_on_lists['excluded'].size >= 1
      no_pass=true if pairs_on_lists['included'].include? pair
    end
    if pairs_on_lists['included'].size >= 1
      no_pass=true if ! pairs_on_lists['included'].include? pair
    end
    break if no_pass

    line.split("\t")[field-1].split(separator).each do |gene|
      if genes_on_lists['excluded'].size >= 1
        no_pass=true if genes_on_lists['excluded'].include? gene
      end
      if genes_on_lists['included'].size >= 1
        no_pass=true if ! genes_on_lists['included'].include? gene
      end
    end

    puts line if not no_pass
  end
end


