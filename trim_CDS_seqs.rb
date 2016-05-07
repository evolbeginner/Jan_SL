#! /bin/env ruby

# to trim sequences based on the number of nucleotides of the CDS sequences
# Written by Sishuo Wang from the department of Botany, the University of British Columbia (sishuowang@hotmail.ca)

###########################################################################
require 'bio'
require 'getoptlong'

cds_file=nil
is_translate=false
codon_table=1

###########################################################################
def trim_cds_seq(seq=nil, is_translate=false, codon_table=1)
  #numberOfBpToTrim = seq.length % 3
  numberOfBpToTrim = 0
  if is_translate then
    naSeq_obj = Bio::Sequence::NA.new(seq)
    1.upto(3) do |i|
      aaSeq=naSeq_obj.translate(i)
      numberOfBpToTrim=i-1 and break if aaSeq !~ /\*(?=.)/
    end
  end
  seq.sub!(/^.{#{numberOfBpToTrim}}/,'')
  return seq
end


def usage
  script_basename=File.basename $0
  puts "Usage of #{script_basename}: ruby #{script_basename} arguments"
  print <<EOF
Arguments
Mandatory arguments:
-i|--cds_file   input file (CDS sequences)
Optional arguments:
--translate     whether to trim sequences based on translated sequences
                default: false
--codon_table   the number of the corresponding codon table
                only functional when '--translate' is specified
                default: 1
-h|--help       see usage
EOF
  puts
  print <<EOF
Codon table
1. "Standard (Eukaryote)"
2. "Vertebrate Mitochondrial"
3. "Yeast Mitochondorial"
4. "Mold, Protozoan, Coelenterate Mitochondrial and Mycoplasma/Spiroplasma"
5. "Invertebrate Mitochondrial"
6. "Ciliate Macronuclear and Dasycladacean"
9. "Echinoderm Mitochondrial"
10. "Euplotid Nuclear"
11. "Bacteria"
12. "Alternative Yeast Nuclear"
13. "Ascidian Mitochondrial"
14. "Flatworm Mitochondrial"
15. "Blepharisma Macronuclear"
16. "Chlorophycean Mitochondrial"
21. "Trematode Mitochondrial"
22. "Scenedesmus obliquus mitochondrial"
23. "Thraustochytrium Mitochondrial"
EOF
  puts "\nPlease write e-mail to sishuowang@hotmail.ca if you have any question and/or suggestion."
  puts
  exit
end


###########################################################################
opts=GetoptLong.new(
  ['-i', '--cds_file', GetoptLong::REQUIRED_ARGUMENT],
  ['--translate', GetoptLong::NO_ARGUMENT],
  ['--codon_table', GetoptLong::REQUIRED_ARGUMENT],
  ['-h','--help', GetoptLong::NO_ARGUMENT],
)

opts.each do |opt,value|
  case opt
    when '-i', '--cds_file'
      cds_file=value
    when '--translate'
      is_translate=true
    when '--codon_table'
      codon_table=value
    when '-h', '--help'
      usage
  end
end

###########################################################################
fh = Bio::FlatFile.open(cds_file)
fh.each_entry do |f|
  f.definition
  new_seq=trim_cds_seq(f.seq, is_translate, codon_table)
  puts ['>', f.definition].join('')
  puts new_seq
end

