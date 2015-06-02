#! /bin/env ruby

require 'getoptlong'

attributes=Hash.new
features=Array.new
is_first_feature=false
gff=nil
upstream,downstream=nil,nil

####################################################
opts=GetoptLong.new(
  ['--gff',GetoptLong::REQUIRED_ARGUMENT],
  ['--upstream',GetoptLong::REQUIRED_ARGUMENT],
  ['--downstream',GetoptLong::REQUIRED_ARGUMENT],
  ['--feature',GetoptLong::REQUIRED_ARGUMENT]
)

opts.each do |opt,value|
  case opt
    when '--gff'
      gff=value
    when '--upstream'
      upstream=value.to_i
    when '--downstream'
      downstream=value.to_i
    when '--feature'
      features.push value
    when '--first_feature'
      is_first_feature=true
  end
end

####################################################
fh=File.open(gff,'r')
while(line=fh.gets) do
  line.chomp!
  #exists=false
  lines=line.split
  seqid, feature, start, stop, strand, attribute = lines.values_at(0,2,3,4,6,8)
  start=start.to_i
  stop=stop.to_i
  next if not features.include? feature
<<EOF
  if ! attributes.empty? then
    attribute.scan(/([^=;]+)=([^=;]+)/).each do |item|
      next if ! attributes.include? $1
      if is_first_feature then
        if strand == '+' then
          exists=true if gff_info.include? $2
        end
      end
      gff_info[$2]=1
    end
  end
  next if exists
EOF

  if upstream then
    case strand
      when '+'
        new_start=start-upstream
        new_stop=start-1
      when '-'
        new_start=stop+1
        new_stop=stop+upstream
    end
  end
  if downstream then
    case strand
      when '+'
        new_start=stop+1
        new_stop=stop+downstream
      when '-'
        new_start=start-1
        new_stop=start-downstream
    end
  end

  next if new_start <= 0
  print lines.values_at(0,1,2).join("\t")+"\t";
  print [new_start, new_stop].join("\t")+"\t";
  print lines[5,10].join("\t")+"\n"
end

