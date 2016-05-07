#! /bin/env ruby

require "getoptlong"


#########################################################
infile = nil
pair_info_file = nil

pairs = Hash.new
counters = Hash.new{|h,k|h[k]=0}


#########################################################
opts = GetoptLong.new(
  ['-i', '--in', GetoptLong::REQUIRED_ARGUMENT],
  ['-p', '--pair', GetoptLong::REQUIRED_ARGUMENT],
)

opts.each do |opt, value|
  case opt
    when '-i', '--in'
      infile = value
    when '-p', '--pair'
      pair_info_file = value
  end
end


#########################################################
if not pair_info_file.nil?
  File.open(pair_info_file).each_line do |line|
    line.chomp!
    line_arr = line.split("\t")
    if line_arr[1] == "2"
      pairs[line_arr[2]] = ""
    end
  end
end


#########################################################
in_fh = File.open(infile, "r")
in_fh.each_line do |line|
  line.chomp!
  if ($.-1) %4 == 0
    if not pairs.empty?
      next if pairs.include?(line)
    end
    in_fh.each_line do |line|
      line.chomp!
      break if line =~ /^$/
      line_arr = line.split("\t")
      next if line_arr.size < 3
      counters["shared"]+=line_arr.size-1 if line =~ /^shared/
      counters["diff"]+=(line_arr.size-1)/2 if line =~ /^diff/
    end
  end
end
in_fh.close


p counters

