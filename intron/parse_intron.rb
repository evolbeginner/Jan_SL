#! /bin/env ruby

require "getoptlong"


#########################################################
infile = nil
pair_info_file = nil
yn00_file = nil
ks_cutoff = nil

pairs = Hash.new
counters = Hash.new{|h,k|h[k]=Hash.new{|h,k|h[k]=0}}
final_counter = Hash.new{|h,k|h[k]=0}
total_counter = Hash.new{|h,k|h[k]=0}
yn00_info = Hash.new


#########################################################
def read_yn00_file(yn00_file)
  yn00_info = Hash.new{|h,k|h[k]={}}
  File.open(yn00_file, "r").each_line do |line|
    line.chomp!
    line_arr = line.split("\t")
    pair, ka, ks, kaks = line_arr
    ka, ks, kaks = ka.to_f, ks.to_f, kaks.to_f
    yn00_info[pair]["ka"] = ka
    yn00_info[pair]["ks"] = ks
    yn00_info[pair]["kaks"] = kaks
  end
  return(yn00_info)
end


#########################################################
opts = GetoptLong.new(
  ['-i', '--in', GetoptLong::REQUIRED_ARGUMENT],
  ['-p', '--pair', GetoptLong::REQUIRED_ARGUMENT],
  ['--yn00', GetoptLong::REQUIRED_ARGUMENT],
  ['--ks', GetoptLong::REQUIRED_ARGUMENT],
)

opts.each do |opt, value|
  case opt
    when '-i', '--in'
      infile = value
    when '-p', '--pair'
      pair_info_file = value
    when '--yn00'
      yn00_file = value
    when '--ks'
      ks_cutoff = value.to_f
  end
end


#########################################################
yn00_info = read_yn00_file(yn00_file) if not yn00_file.nil?


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
pair = nil
in_fh = File.open(infile, "r")
in_fh.each_line do |line|
  line.chomp!
  next if line =~ /^$/
  if line.split("\t").size == 1 
    pair = line
  end

  if not pairs.empty?
    if pairs.include?(pair)
      next
    end
  end

  if not yn00_info.empty?
    next if not yn00_info.include?(pair)
    if yn00_info[pair]['ks'] > ks_cutoff
      next
    end
  end
  
  unless line.split("\t").size < 5
    line_arr = line.split("\t")
    if line_arr[1] == "shared"
      counters[pair]["shared"] += 1
      total_counter["shared"] += 1
    elsif line_arr[1] == "diff"
      counters[pair]["diff"] += 1
      total_counter["diff"] += 1
    end
  end
end
in_fh.close


counters.each_pair do |pair, v|
  if not v.include?("shared") and v.include?("diff")
    final_counter["diff"] += 1
  elsif v.include?("shared") and not v.include?("diff")
    final_counter["shared"] += 1
  elsif v.include?("shared") and v.include?("diff")
    final_counter["both"] += 1
  end
end


p final_counter
p total_counter
p counters.keys.size


