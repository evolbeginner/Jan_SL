#! /bin/env ruby

require 'getoptlong'


###############################################################
pair_file=nil
gene_list_file=nil
field=1
pair_sep='-'
is_corename=false
is_output_each_containing = false

genes_list=Hash.new()

opts = GetoptLong.new(
  ['--pair','--pair_file',GetoptLong::REQUIRED_ARGUMENT],
  ['--list_file',"--gene_list_file",GetoptLong::REQUIRED_ARGUMENT],
  ['-f',"--field",GetoptLong::REQUIRED_ARGUMENT],
  ['--pair_sep','--sep',GetoptLong::REQUIRED_ARGUMENT],
  ['--corename',GetoptLong::NO_ARGUMENT],
  ['--each_containing',GetoptLong::NO_ARGUMENT],
)

opts.each do |opt,value|
  case opt
    when '--pair', '--pair_file'
      pair_file=value
    when '--list_file', '--gene_list_file'
      gene_list_file=value
    when '-f', '--field'
      field=value.to_i
    when '--sep', '--pair_sep'
      pair_sep=value
    when '--corename'
      is_corename=true
    when '--each_containing'
      is_output_each_containing = true
  end
end


###############################################################
File.open(gene_list_file,'r').each_line do |line|
  line.chomp!
  if field
    target = line.split('\t')[field-1]
  else
    target = line
  end
  genes_list[target]=''
end


File.open(pair_file,'r').each_line do |line|
  line.chomp!
  line_arr = line.split("\t")
  sl_gene, pair = line_arr.values_at(0,2)
  #pair = (line.split("\t"))[0]
  genes = pair.split(pair_sep)
  if is_output_each_containing
    genes.each do |i|
      ele = is_corename ? i.sub!(/\.[^.]+$/,'') : i
      puts ele if genes_list.include?(ele)
    end
  else
    if ! genes.select{|i| ele = is_corename ? i.sub!(/\.[^.]+$/,'') : i; genes_list.include?(ele)}.empty?
      #p genes.select{|i| ele = is_corename ? i.sub!(/\.[^.]+$/,'') : i}
      puts pair
    end
  end
end


