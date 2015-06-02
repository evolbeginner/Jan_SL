#! /bin/env ruby

require 'getoptlong'

###########################################################################
input=nil
ranges_cutoffs=Array.new
ranges=Array.new
outdir='./'
is_force=false
results=Hash.new{|h,k|h[k]=Array.new}
separator="\t"
f1=1
f2=2
suffix=nil

opts=GetoptLong.new(
  ['--in',GetoptLong::REQUIRED_ARGUMENT],
  ['--ranges_cutoff',GetoptLong::REQUIRED_ARGUMENT],
  ['--outdir',GetoptLong::REQUIRED_ARGUMENT],
  ['--sep',GetoptLong::REQUIRED_ARGUMENT],
  ['--f1',GetoptLong::REQUIRED_ARGUMENT],
  ['--f2',GetoptLong::REQUIRED_ARGUMENT],
  ['--suffix',GetoptLong::REQUIRED_ARGUMENT],
  ['--force',GetoptLong::NO_ARGUMENT],
)

opts.each do |opt,value|
  case opt
    when '--in'
      input=value
    when '--ranges_cutoff'
      value.split(',').map{|i|ranges_cutoffs.push i.to_i}
    when '--outdir'
      outdir=File.expand_path(value)
    when '--sep'
      separator=value
    when '--f1'
      f1=value.to_i
    when '--f2'
      f2=value.to_i
    when '--suffix'
      suffix=value
    when '--force'
      is_force=true
  end
end

raise "input has to be given!" if not input
raise "ranges have to be given!" if not ranges_cutoffs
if not Dir.exists?(outdir)
  Dir.mkdir outdir
else
  if is_force
    `rm -rf #{outdir}`
    Dir.mkdir outdir
  else
    raise "outdir #{outdir} has already existed!"
  end
end

ranges_cutoffs.sort.each_with_index do |ele,index|
  next if index == ranges_cutoffs.size-1
  ranges.push ranges_cutoffs[index]..ranges_cutoffs[index+1]
end

###########################################################################
File.open(input,'r').each_line do |line|
  line.chomp!
  line_array=line.split(separator).values_at(f1-1,f2-1).map{|i|i.to_f}
  ranges.each do |range| 
    if range.include? line_array[0]
      results[range].push line_array[1]
      break
    end
  end
end

results.each_pair do |key,value|
  output_name_array=[outdir,key.to_s]
  output_name_array.push suffix if suffix
  output_name = File.join(output_name_array)
  if File.exists?(output_name)
    if not is_force
      raise "#{output_name} has already existed!"
    end
  end
  fh = File.open(output_name, 'w')
  fh.puts value
  fh.close
end


