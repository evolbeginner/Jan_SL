#! /bin/env ruby


require "getoptlong"
require "Hash"


###################################################################
infile = nil
pair_info_file = nil
intron_diff_file = nil
cds_num_file = nil
kaks_file = nil
blast8_file = nil
evalue_cutoff = 0.001
is_show_gene = false

final_gene_sets = multi_D_Hash(3)
gene_info = Hash.new


###################################################################
class Gene_info
  attr_accessor :cds_num, :ks
  def initialize()
    ;
  end
end


class Infile
  attr_accessor

  def initialize(infile)
    @infile = infile
  end

  def parse_file
    intron_info = Hash.new{|h,k|h[k]={}}
    in_fh = File.open(@infile, "r")
    in_fh.each_line do |line|
      line.chomp!
      line_arr = line.split("\t")
      intron_info[line_arr[0]][line_arr[1]] = ""
    end
    in_fh.close
    return(intron_info)
  end
end


class Pair_info_file
  attr_reader :sl_gene_info, :parent_info

  def initialize(infile)
    @infile = infile
    @sl_gene_info = Hash.new{|h,k|h[k]={}}
    @parent_info = Hash.new{|h,k|h[k]={}}
  end

  def set_intron_info(intron_info)
    @intron_info = intron_info
  end

  def select_intron()
    counter = 0
    @sl_gene_info.each_pair do |sl_gene, v|
      if @intron_info.include?(sl_gene)
        counter += @intron_info[sl_gene].size
      end
    end
    #puts counter
  end

  def parse_file
    File.open(@infile).each_line do |line|
      line.chomp!
      line_arr = line.split("\t")
      sl_gene, num_of_SL, pair = line_arr
      genes = pair.split('-')
      parent = nil
      genes.each do |gene|
        parent = gene and break if gene != sl_gene
      end
      @sl_gene_info[sl_gene]['num_of_SL'] = num_of_SL.to_i
      @parent_info[parent]['num_of_SL'] = num_of_SL.to_i
    end
  end
end


class Intron_diff_file
  attr_accessor :intron_posi_info

  def initialize(infile)
    @infile = infile
    @intron_posi_info = Hash.new{|h,k|h[k]=Hash.new{|h,k|h[k]=[]}}
  end

  def get_diff_shared
    a=0
    File.open(@infile, 'r').each_line do |line|
      line.chomp!
      #26	diff	symbB.v1.2.000269	26	80-2389
      line_arr = line.split("\t")
      next if line_arr.size < 5
      a+=1
      0.upto((line_arr.size-2)/3-1) do |i|
        type, gene, intron_posi = line_arr.values_at(1,2+3*i,4+3*i)
        @intron_posi_info[type][gene] << intron_posi
      end
    end
  end
end


class Blast8_file
  attr_accessor :passed
  def initialize(blast8_file)
    @blast8_file = blast8_file
  end
  def read_file
    @evalue_info = Hash.new{|h,k|h[k]=[]}
    File.open(@blast8_file, 'r').each_line do |line|
      line.chomp!
      next if line =~ /^#/
      line_arr = line.split("\t")
      query, subject, evalue = line_arr.values_at(0,1,-2)
      @evalue_info[subject] << evalue.to_f
    end
    return(@evalue_info)
  end
  def filter_evalue(evalue_cutoff=0.001)
    @passed = Hash.new
    @evalue_info.each_pair do |gene, v|
      @passed[gene] = "" if v.any?{|i|i<=evalue_cutoff}
    end
    #p [1, @passed.size]
    #return(@passed)
  end
end


###################################################################
def read_cds_num_file(cds_num_file, gene_info)
  File.open(cds_num_file, 'r').each_line do |line|
    line.chomp!
    line_arr = line.split("\t")
    gene, cds_num = line_arr.values_at(0,2)
    cds_num = cds_num.to_i
    if not gene_info.include?(gene)
      gene_info_obj = Gene_info.new()
      gene_info[gene] = gene_info_obj
    end
    gene_info[gene].cds_num = cds_num
  end
  return(gene_info)
end


def read_kaks_file(kaks_file, gene_info)
  File.open(kaks_file, "r").each_line do |line|
    line.chomp!
    #symbB.v1.2.000022-symbB.v1.2.012916	0.8123	4.4974	0.1806
    line_arr = line.split("\t")
    pair, ka, ks, kaks = line_arr
    ka, ks, kaks = ka.to_f, ks.to_f, kaks.to_f
    genes = pair.split("-")
    genes.each do |gene|
      if not gene_info.include?(gene)
        gene_info_obj = Gene_info.new(gene)
        gene_info[gene] = gene_info_obj
      end
      gene_info[gene].ks = ks
    end
  end
  return(gene_info)
end


###################################################################
opts = GetoptLong.new(
  ['-i', '--in', GetoptLong::REQUIRED_ARGUMENT],
  ['-p', '--pair_info', GetoptLong::REQUIRED_ARGUMENT],
  ['--intron_diff', GetoptLong::REQUIRED_ARGUMENT],
  ['--cds_num', '--exon_num', GetoptLong::REQUIRED_ARGUMENT],
  ['--KaKs', '--KaKs_file', '--kaks', '--kaks_file', GetoptLong::REQUIRED_ARGUMENT],
  ['-b', '--blast8', GetoptLong::REQUIRED_ARGUMENT],
  ['-e', '--evalue', GetoptLong::REQUIRED_ARGUMENT],
  ['--show', '--show_gene', GetoptLong::NO_ARGUMENT],
)

opts.each do |opt, value|
  case opt
    when '-i', '--in'
      infile = value
    when '-p', '--pair_info'
      pair_info_file = value
    when '--intron_diff'
      intron_diff_file = value
    when '--cds_num', '--exon_num'
      cds_num_file = value
    when '--KaKs', '--KaKs_file', '--kaks', '--kaks_file'
      kaks_file = value
    when '-b', '--blast8'
      blast8_file = value
    when '-e', '--evalue'
      evalue_cutoff = value.to_f
    when '--show', '--show_gene'
      is_show_gene = true    
  end
end


###################################################################
infile_obj = Infile.new(infile)
intron_info = infile_obj.parse_file

pair_info_file_obj = Pair_info_file.new(pair_info_file)
pair_info_file_obj.parse_file
pair_info_file_obj.set_intron_info(intron_info)

intron_diff_file_obj = Intron_diff_file.new(intron_diff_file)
intron_diff_file_obj.get_diff_shared

#p intron_diff_file_obj.intron_posi_info["diff"].map{|i,v|v.size}.reduce(:+)
#p intron_diff_file_obj.intron_posi_info["shared"].map{|i,v|v.size}.reduce(:+)

if not blast8_file.nil?
  blast8_file_obj = Blast8_file.new(blast8_file)
  blast8_file_obj.read_file
  blast8_file_obj.filter_evalue(evalue_cutoff)
end


gene_info = read_cds_num_file(cds_num_file, gene_info)
gene_info = read_kaks_file(kaks_file, gene_info)


###################################################################
intron_diff_file_obj.intron_posi_info.each_pair do |type, v|
  v.each_pair do |gene, positions|
    positions.each do |position|
      intron_full_name = gene + ':' + position
      #next if not pair_info_file_obj.parent_info.include?(gene)
      next if not pair_info_file_obj.sl_gene_info.include?(gene)
      unless blast8_file_obj.nil?
        next if not blast8_file_obj.passed.include?(gene)
      end
      next if pair_info_file_obj.sl_gene_info[gene]['num_of_SL'] == 2
      final_gene_sets[type]["all"][intron_full_name] = ""

      if gene_info[gene].methods.include?('ks') and gene_info[gene].ks <= 2
        final_gene_sets[type]["ks"][intron_full_name] = ""
      end
      if gene_info[gene].cds_num <= 5
        final_gene_sets[type]["cds"][intron_full_name] = ""
      end
      if not (not intron_info.include?(gene) or not intron_info[gene].include?(intron_full_name))
        final_gene_sets[type]["on_list"][intron_full_name] = ""
      end
    end
  end
  
  #intersect_list = final_gene_sets[type]["ks"].keys & final_gene_sets[type]["cds"].keys & final_gene_sets[type]["on_list"].keys
  intersect_list = ['ks', 'cds', 'on_list'].inject(Array.new){|i,j| i.empty? ? final_gene_sets[type][j].keys : i & final_gene_sets[type][j].keys}
  intersect_list.each do |i|
    final_gene_sets[type]["intersect"][i] = ""
  end
end


final_gene_sets.each_pair do |type, v|
  puts type
  puts ["all", v["all"].size].map{|i|i.to_s}.join("\t")
  if is_show_gene
    v['all'].each_key do |intron|
      puts intron
    end
  end
  puts ["ks", v["ks"].size].map{|i|i.to_s}.join("\t")
  puts ["cds", v["cds"].size].map{|i|i.to_s}.join("\t")
  puts ["on_list", v["on_list"].size].map{|i|i.to_s}.join("\t")
  puts ["intersect", v["intersect"].size].map{|i|i.to_s}.join("\t")
end


