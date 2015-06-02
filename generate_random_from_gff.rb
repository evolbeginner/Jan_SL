#! /bin/env ruby

require 'getoptlong'


###############################################################
gff = nil
num = nil
features_included = []
attributes_included = []

final_set = Hash.new


###############################################################
opts = GetoptLong.new(
  ["--gff", GetoptLong::REQUIRED_ARGUMENT],
  ["--feature", GetoptLong::REQUIRED_ARGUMENT],
  ["--attribute", GetoptLong::REQUIRED_ARGUMENT],
  ["--num", GetoptLong::REQUIRED_ARGUMENT],
)

opts.each do |opt, value|
  case opt
    when '--gff'
      gff = value
    when '--feature'
      value.split(',').each do |i|
        features_included.push i
      end
    when '--attribute'
      value.split(',').each do |i|
        attributes_included.push i
      end
    when '--num'
      num = value.to_i
  end
end


raise "gff has to be given! Exiting ......" if gff.nil?
raise "attributes have to be given! Exiting ......" if attributes_included.empty?
raise "num has to be given! Exiting ......" if num.nil?


###############################################################
in_fh = File.open(gff, 'r')
in_fh.each_line do |line|
  line.chomp!
  next if line =~ /^#/
  #scaffold5.1|size591573	AUGUSTUS	intron	1	1324	0.42	-	.	Parent=symbB.v1.2.000001.t1
  genome, feature, start, stop, strand, attributes  = line.split(/\t/).values_at(0,2,3,4,6,8)

  if not features_included.empty?
    next if not features_included.include?(feature)
  end

  attributes.split(";").each do |attribute|
    if attribute =~ /([^=]+)=([^;]+)/
      if attributes_included.include?($1)
        final_set[$2] = ''
      end
    end
  end

  #break if $. >= 10000
end

in_fh.close


final_set_keys = final_set.keys
1.upto(num) do |i|
  puts final_set_keys.sample(2).join("\t")
end


