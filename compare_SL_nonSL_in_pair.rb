#! /bin/env ruby


require "getoptlong"


#################################################################
def read_list(list_file)
  genes = Hash.new
  File.open(list_file, "r").each_line do |line|
    line.chomp!
    genes[line] = ""
  end
  return(genes)
end



#################################################################
list1 = nil
list2 = nil

pair_file = nil


#################################################################
opts = GetoptLong.new(
  ["--i1", GetoptLong::REQUIRED_ARGUMENT],
  ["--i2", GetoptLong::REQUIRED_ARGUMENT],
  ["-p", GetoptLong::REQUIRED_ARGUMENT],
)

opts.each do |opt, value|
  case opt
    when "--i1"
      list1 = value
    when "--i2"
      list2 = value
    when "-p"
      pair_file = value
  end
end


#################################################################
genes1 = read_list(list1)
genes2 = read_list(list2)


File.open(pair_file, "r").each_line do |line|
  line.chomp!
  line_arr = line.split("\t")
  sl_gene, pair = line_arr.values_at(0,2)
  parent = nil
  pair.split('-').each do |gene|
    if gene != sl_gene
      parent = gene
    end
  end
  
  is_pass = false
  [sl_gene, parent].each do |i|
    is_pass = true if genes1.include?(i) or genes2.include?(i)
  end

  next unless is_pass
  [sl_gene, parent].each do |gene|
    output = nil
    if genes1.include?(gene)
      if genes2.include?(gene)
        output = "both"
      else
        output = "SL"
      end
    elsif genes2.include?(gene)
      output = "nonSL"
    else
      output = "-"
    end
    print [gene,output].join("\t")
    print "\t"
  end
  print line_arr[1]
  puts
end


