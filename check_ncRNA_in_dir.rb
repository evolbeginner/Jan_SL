#! /bin/env ruby2.1

require "getoptlong"
require "triez"


############################################################
gff_file = nil
features = Array.new

targets = Triez.new
is_explicit = true


############################################################
opts = GetoptLong.new(
  ["-i", "--gff", GetoptLong::REQUIRED_ARGUMENT],
  ["--feature", GetoptLong::REQUIRED_ARGUMENT],
  ["--implicit", GetoptLong::NO_ARGUMENT],
)

opts.each do |opt, value|
  case opt
    when "-i", "--gff"
      gff_file = value
    when "--feature"
      value.split(",").map{|i|features << i}
    when "--implicit"
      is_explicit = false
  end
end


raise "features empty! Exiting ......" if features.empty?


############################################################
in_fh = STDIN
in_fh.each_line do |file_basename|
  next if file_basename =~ /^\./
  file_basename.chomp!
  targets[file_basename] = 1
end
in_fh.close


### ---------------------------------------------------- ###
old_feature = "!!"
old_id = "!!"

File.open(gff_file, "r").each_line do |line|
  line.chomp!
  line_arr = line.split("\t")
  feature = line_arr[2]
  id = nil
  if features.include?(old_feature)
    if line_arr[-1] =~ /ID=([^;]+)/
      id = $1
      if is_explicit
        puts id if targets.has_key?(id)
      else
        targets.search_with_prefix(id){|i|puts i}
      end
    end
  end
  old_feature = feature
  old_id = id
end



