#! /bin/env ruby
# convert output file of RepeatMasker to gff format
# convert genome.fa.out to the gff format

#######################################################################

RepeatMasker_out=ARGV[0]
File.open(RepeatMasker_out,'r').each_line do |line|
  line.chomp!
  next if line !~ /\w/
  div,chr,start,stop,attribute = line.split("\s").values_at(1,4,5,6,10)
  new_attribute = "ID=" + attribute + "|" + chr + ":" + start + "-" + stop
  puts [chr, 'RepeatMasker', 'TE', start, stop, div, '+', '.', new_attribute].join("\t")
  # scaffold10.1|size534768	RepeatMasker	similarity	355822	355866	 2.3	+	.	Target="Motif:(TGATT)n" 1 45
end

