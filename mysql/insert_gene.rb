#! /bin/env ruby2.1


require 'mysql'
require 'getoptlong'
require "bio"


##########################################################################
infiles = Array.new
orgn = nil
exon_file = nil
pair_file = nil
gff_file = nil
attr = "ID"
pep_seq_file = nil
up300_seq_file = nil
suffix = nil

gff_info = Hash.new


##########################################################################
def read_gff(gff_file, attr, suffix)
  gff_info = Hash.new{|h,k|h[k]={}}
  File.open(gff_file, 'r').each_line do |line|
    line.chomp!
    line_arr = line.split("\t")
    if line_arr[-1] =~ /#{attr}=([^;]+)/
      gene = $1
      gene.sub!(/#{suffix}/, "") unless suffix.nil?
      #scaffold39267.1	AUGUSTUS	transcript	1	175	0.16	-	.	ID=symbB.v1.2.044366.t1;Parent=symbB.v1.2.044366
      gff_info[gene]['strand'] = line_arr[6]
      gff_info[gene]['chr'] = line_arr[0]
      if ! gff_info[gene].include?('locations')
        gff_info[gene]['locations'] = Array.new
      end
      gff_info[gene]['locations'] << line_arr[3].to_i
      gff_info[gene]['locations'] << line_arr[4].to_i
    end
  end
  gff_info.each_key do |gene|
    gff_info[gene]['start'] = gff_info[gene]['locations'].min
    gff_info[gene]['end'] = gff_info[gene]['locations'].max
  end
  return(gff_info)
end


def read_gene_lists(infiles)
  genes = Hash.new
  infiles.each do |infile|
    File.open(infile, 'r').each_line do |line|
      line.chomp!
      genes[line] = ""
    end
  end
  return(genes)
end


def read_exon_file(exon_file)
  exon_num = Hash.new
  File.open(exon_file, 'r').each_line do |line|
    line.chomp!
    line_arr = line.split("\t")
    exon_num[line_arr[0]] = line_arr[2].to_i
  end
  return(exon_num)
end


def read_pair_file(pair_file)
  parent_info = Hash.new
  File.open(pair_file, 'r').each_line do |line|
    line.chomp!
    parent = nil
    #symbB.v1.2.000323	1	symbB.v1.2.000323-symbB.v1.2.012087
    line_arr = line.split("\t")
    sl_gene = line_arr[0]
    pair = line_arr[2]
    pair.split('-').each do |i|
      if sl_gene != i
        parent = i
        break
      end
    end
    parent_info[sl_gene] = parent
  end
  return(parent_info)
end


def read_seq_file(seq_file)
  seq_info = Hash.new
  Bio::FlatFile.open(seq_file).each_entry do |f|
    seq_info[f.definition] = f.seq
  end
  return(seq_info)
end


##########################################################################
opts = GetoptLong.new(
  ["-i", '--in', GetoptLong::REQUIRED_ARGUMENT],
  ["--orgn", '--ORGN', GetoptLong::REQUIRED_ARGUMENT],
  ["--exon", '--exon_file', '--num_exon', GetoptLong::REQUIRED_ARGUMENT],
  ["-p", '--pair', GetoptLong::REQUIRED_ARGUMENT],
  ["--gff", GetoptLong::REQUIRED_ARGUMENT],
  ["--attr", GetoptLong::REQUIRED_ARGUMENT],
  ["--pep", "--pep_seq", GetoptLong::REQUIRED_ARGUMENT],
  ["--up300", "--up300_seq", GetoptLong::REQUIRED_ARGUMENT],
  ["--suffix", GetoptLong::REQUIRED_ARGUMENT],
)

opts.each do |opt, value|
  case opt
    when '-i', '--in'
      infiles << value
    when '--orgn', '--ORGN'
      orgn = value
    when '--exon', '--exon_file', '--num_exon'
      exon_file = value
    when '-p', '--pair'
      pair_file = value
    when '--gff'
      gff_file = value
    when '--attr'
      attr = value
    when '--pep', '--pep_seq'
      pep_seq_file = value
    when '--up300', '--up300_seq'
      up300_seq_file = value
    when '--suffix'
      suffix = value
  end
end


##########################################################################
pep_seq_info = read_seq_file(pep_seq_file)

up300_seq_info = read_seq_file(up300_seq_file)

genes = read_gene_lists(infiles)

gff_info = read_gff(gff_file, attr, suffix)

exon_num = read_exon_file(exon_file)

parent_info = read_pair_file(pair_file)


##########################################################################
begin
  con = Mysql.new 'localhost', 'sswang', '123456', 'SL'
  #con.query("drop table IF EXISTS Genes")
  con.query("CREATE TABLE IF NOT EXISTS \
            Genes(
              pri VARCHAR(50),
              id VARCHAR(50),
              orgn VARCHAR(30),
              parent VARCHAR(50),
              exon_num INT,
              chr VARCHAR(40),
              start INT,
              end INT,
              pep_seq VARCHAR(5000),
              up300_seq VARCHAR(350)
              )"
           )

  genes.each_key do |gene|
    if not exon_num.include?(gene)
      puts "Warming: gene #{gene} not included in exon_num file."
      next
    end

    if not parent_info.include?(gene)
      parent = ''
    else
      parent = "#{parent_info[gene]}"
    end

    pri = [orgn, gene].join('-')
    con.query("INSERT INTO Genes(pri, id, orgn, parent, exon_num, chr, start, end, pep_seq, up300_seq) \
              VALUES(\"#{pri}\", \"#{gene}\", \"#{orgn}\", \"#{parent}\", #{exon_num[gene]}, \"#{gff_info[gene]['chr']}\", #{gff_info[gene]['start']}, #{gff_info[gene]['end']}, \"#{pep_seq_info[gene]}\", \"#{up300_seq_info[gene]}\")"
             )
  end

  rs = con.query("SELECT * FROM Genes")

  rs.each do |row|
    puts row.join("\t")
  end    
 
rescue Mysql::Error => e
  puts e.errno
  puts e.error
    
ensure
  con.close if con
end


