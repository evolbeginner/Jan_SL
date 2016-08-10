#! /bin/env ruby


require "getoptlong"


########################################################################################
infile = nil

freqs = Hash.new{|h,k|h[k]=0}
1.upto(10).each do |i|
  freqs[i] = 0
end
freqs[100] = 0


########################################################################################
opts = GetoptLong.new(
  ["-i", GetoptLong::REQUIRED_ARGUMENT],
)

opts.each do |opt, value|
  case opt
    when "-i"
      infile = value
  end
end


########################################################################################
File.open(infile, 'r').each_line do |line|
  line.chomp!
  num = line.to_i
  if num < 10
    freqs[num] += 1
  else
    freqs[100] += 1
  end
end


freqs.keys.sort.each do |num|
  puts [num, freqs[num]/freqs.values.inject(0){|i,j|i+j}.to_f].map{|i|i.to_s}.join("\t")
end


