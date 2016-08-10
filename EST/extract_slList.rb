#! /bin/env ruby


require "getoptlong"


#################################################################
sl_gff = nil
subject_gff = nil

est_rela = Hash.new
genes = Hash.new


#################################################################
opts = GetoptLong.new(
  ["--sl_gff", GetoptLong::REQUIRED_ARGUMENT],
  ["--subject_gff", GetoptLong::REQUIRED_ARGUMENT],
)


opts.each do |opt, value|
  case opt
    when "--sl_gff"
      sl_gff = value
    when "--subject_gff"
      subject_gff = value
  end
end


#################################################################
File.open(subject_gff, "r").each_line do |line|
  line.chomp!
  chr, attr = line.split("\t").values_at(0, 8)
  attr =~ /ID=([^;]+)/
  gene = $1
  est_rela[chr] = gene
end


File.open(sl_gff, "r").each_line do |line|
  line.chomp!
  chr, attr = line.split("\t").values_at(0, 8)
  next if not est_rela.include?(chr)
  genes[est_rela[chr]] = ""
end


genes.keys.sort.each do |gene|
  puts gene
end


