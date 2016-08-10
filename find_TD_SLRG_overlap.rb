#! /bin/env ruby


require "getoptlong"


#############################################################
infile = nil
pair_file = nil
field2 = 3
sep = "\t"
is_final = false

clstr_rela = Hash.new{|h,k|h[k] = Hash.new{|h,k|h[k]={}}}


#############################################################
opts = GetoptLong.new(
  ["-i", "--td", "--TD", "--td_clstr", GetoptLong::REQUIRED_ARGUMENT],
  ["-p", GetoptLong::REQUIRED_ARGUMENT],
  ["--f2", GetoptLong::REQUIRED_ARGUMENT],
  ["--sep", GetoptLong::REQUIRED_ARGUMENT],
  ["--is_final", GetoptLong::NO_ARGUMENT],
)

opts.each do |opt, value|
  case opt
    when "-i", "--td", "--TD", "--td_clstr"
      infile = value
    when "-p"
      pair_file = value
    when "--f2"
      field2 = value.to_i
    when "--sep"
      sep = value
    when "--final"
      is_final = true
  end
end


#############################################################
File.open(infile, "r").each_line do |line|
  line.chomp!
  line_arr = line.split("\t")
  line_arr.each do |i|
    line_arr.each do |j|
      clstr_rela[i][j] = "" if i != j
    end
  end
end


if field2 == 3
  is_final = true
  sep = "-"
end


File.open(pair_file, "r").each_line do |line|
  line.chomp!
  line_arr = line.split("\t")
  gene = nil
  if is_final
    pair = line_arr[field2-1]
    gene = line_arr[0]
  else
    pair = line
  end
  genes = pair.split(sep)
  if clstr_rela.include?(genes[0])
    if clstr_rela[genes[0]].include?(genes[1])
      puts pair
    end
  end
end



