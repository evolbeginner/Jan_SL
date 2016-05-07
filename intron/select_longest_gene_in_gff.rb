#! /bin/env ruby

require "getoptlong"


#################################################################
gff_file = nil
longest_rela_file = nil
features = ["CDS"]
attr = "Parent";

rela = Hash.new


#################################################################
opts = GetoptLong.new(
  ["--gff", GetoptLong::REQUIRED_ARGUMENT],
  ["--longest_rela", GetoptLong::REQUIRED_ARGUMENT],
  ["--feature", GetoptLong::REQUIRED_ARGUMENT],
  ["--attr", GetoptLong::REQUIRED_ARGUMENT],
)


opts.each do |opt, value|
  case opt
    when '--gff'
      gff_file = value
    when '--longest_rela'
      longest_rela_file = value
    when '--feature'
      features << value
    when '--attr'
      attr = value
  end
end


#################################################################
File.open(longest_rela_file).each_line do |line|
  line.chomp!
  line_arr = line.split("\t")
  rela[line_arr[1]] = line_arr[0]
end


File.open(gff_file).each_line do |line|
  line.chomp!
  line_arr = line.split("\t")
  next if not features.include?(line_arr[2])
  if line_arr[-1] =~ /#{attr}=([^;]+)/
    if rela.include?($1)
      puts line
    end
  end
end


