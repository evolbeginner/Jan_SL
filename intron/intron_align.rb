#! /bin/env ruby


require "getoptlong"
require "bio"


#################################################
def get_well_aligned_posi(aln_info, aligned_posi, prop_of_diff_sites_cutoff, prop_of_gaps_cutoff)
  diff_sites = Hash.new
  gap_sites = Hash.new
  well_aligned_posi = Hash.new
  site_info = Hash.new{|h,k|h[k]={}}

  aln_info.each_pair do |gene, sites|
    sites.each_with_index do |site, index|
      site_info[index+1][site] = ""
    end
  end

  site_info.keys.sort.each do |posi|
    diff_sites[posi] = "" if site_info[posi].size != 1
    gap_sites[posi] = "" if site_info[posi].include?('-')
  end

  aligned_posi.keys.sort.each do |aln_posi|
    diff_counter = 0
    down_counter = 0
    up_counter = 0
    gap_counter = 0
    (aln_posi-1).to_i.downto(aln_posi-5).select{|i|i>0 and i<=aligned_posi.keys.max}.each do |i|
      down_counter += 1
      diff_counter += 1 if diff_sites.include?(i) #or not aligned_posi.include?(i)
      gap_counter += 1 if not aligned_posi.include?(i)
    end
    (aln_posi+1).to_i.upto(aln_posi+5).select{|i|i>0 and i<=aligned_posi.keys.max}.each do |i|
      up_counter += 1
      diff_counter += 1 if diff_sites.include?(i) #or not aligned_posi.include?(i)
      gap_counter += 1 if not aligned_posi.include?(i)
    end
    effect_counter = up_counter + down_counter
    identity = 1 - diff_counter/effect_counter.to_f
    prop_of_gaps = gap_counter/effect_counter.to_f
    if prop_of_gaps < prop_of_gaps_cutoff and identity > 1-prop_of_diff_sites_cutoff
      #p [aln_posi, effect_counter, identity, prop_of_gaps] 
      well_aligned_posi[aln_posi] = ""
    end
  end

  return(well_aligned_posi)
end


def get_intron_in_aln(aligned_posi, real_aligned_posi, intron_info)
  intron_in_aln = Hash.new{|h,k|h[k]=[]}
  aligned_posi.keys.sort.each do |aln_posi|
    real_aligned_posi.each_pair do |seq_title, v|
      real_posi = v[aln_posi]
      if intron_info[seq_title].include?(real_posi)
        intron_in_aln[aln_posi] << seq_title
        # puts [seq_title, real_posi, intron_info[seq_title][real_posi]].map{|i|i.to_s}.join("\t")
      end
    end
  end
  return(intron_in_aln)
end


def get_aln_info(aln_info)
  aln_posi_info = Hash.new{|h,k|h[k]=[]}
  aligned_posi = Hash.new # aligned positions in all sequences in MSA
  real_aligned_posi = Hash.new{|h,k|h[k]={}}

  aln_info.each_pair do |seq_title, arr|
    num_of_spaces = 0
    arr.each_with_index do |site, index|
      aln_posi = index + 1
      if site != '-'
        aln_posi_info[aln_posi] << seq_title
        real_aln_posi = aln_posi - num_of_spaces
        real_aligned_posi[seq_title][aln_posi] = real_aln_posi
      else
        num_of_spaces += 1
      end
    end
  end

  aln_posi_info.each_pair do |aln_posi, v|
    if v.size == aln_info.keys.size
      aligned_posi[aln_posi] = ""
    end
  end
  return([aln_posi_info, aligned_posi, real_aligned_posi])
end


def read_aln_file(infile, aln_info)
  Bio::FlatFile.open(infile).each_entry do |f|
    (f.seq).split("").each_with_index do |site, index|
      aln_info[f.definition].push site
    end
  end
  return(aln_info)
end


def read_indir(indir, aln_info, extensions=[])
  Dir.foreach(indir) do |file|
    extname = File.extname(file)
    extname.sub!(/^./, "")
    next if not extensions.include?(extname)
    file_full_name = File.join([indir, file])
    aln_info = read_aln_file(file_full_name, aln_info)
  end
  return(aln_info)
end


def read_gff(gff_file, cds_info)
  File.open(gff_file, "r").each_line do |line|
    line.chomp!
    line_arr = line.split("\t")
    #chr1	AUGUSTUS	CDS	1	38	0.92	+	.	ID=2;
    type, start, stop, attr = line_arr.values_at(2,3,4,8)
    id = nil
    if type == "CDS"
      start = start.to_i
      stop = stop.to_i
      id = $1 if attr =~ /(?:Parent|ID)=([^;]+)/
      (start..stop).each do |posi|
        cds_info[id][posi] = ""
      end
    end
  end
  return(cds_info)
end


def parse_cds_info(cds_info)
  intron_info = Hash.new{|h,k|h[k]={}}
  cds_info.each_pair do |seq_title, v1|
    is_exon = true
    pre_na_posi = 0
    pre_aa_posi = 0
    aa_posi = nil
    first_na_posi = 1
    v1.each_pair do |na_posi, v2|
      if na_posi != pre_na_posi+1
        first_na_posi = na_posi
        intron_info[seq_title][pre_aa_posi] = [pre_na_posi+1, na_posi-1].map{|i|i.to_s}.join('-')
      end
      rela_aa_posi = (na_posi-first_na_posi+1)/3

      if (na_posi-first_na_posi+1)%3 != 0
        aa_posi = pre_aa_posi + 1
      else
        aa_posi = pre_aa_posi + 1
        pre_aa_posi = aa_posi
      end
      #p [na_posi, aa_posi, pre_aa_posi]
      pre_na_posi = na_posi
    end
  end
  return(intron_info)
end


#################################################
infile = nil
indir = nil
extensions = Array.new
gff_files = Array.new
prop_of_diff_sites_cutoff = 0.5
prop_of_gaps_cutoff = 0.3
intron_sliding_cutoff = 5

aln_info = Hash.new{|h,k|h[k]=Array.new}
cds_info = Hash.new{|h,k|h[k]=Hash.new}


#################################################
opts = GetoptLong.new(
  ["-i", "--in", "--aln", GetoptLong::REQUIRED_ARGUMENT],
  ["--indir", GetoptLong::REQUIRED_ARGUMENT],
  ["--extension", GetoptLong::REQUIRED_ARGUMENT],
  ["--gff", GetoptLong::REQUIRED_ARGUMENT],
  ["--prop_gaps", "--prop_gap", GetoptLong::REQUIRED_ARGUMENT],
  ["--identity", GetoptLong::REQUIRED_ARGUMENT],
  ["--prop_diff_sites", "--prop_diff_site", GetoptLong::REQUIRED_ARGUMENT],
  ["--intron_sliding", "--intron_sliding_cutoff", GetoptLong::REQUIRED_ARGUMENT],
)


opts.each do |opt, value|
  case opt
    when "-i", "--in", "--aln"
      infile = value
    when '--indir'
      indir = value
    when '--extension'
      extensions << value
    when "--gff"
      gff_files << value
    when "--identity"
      prop_of_diff_sites_cutoff = 1-value.to_f
    when "--diff_site", "--diff_sites"
      prop_of_diff_sites_cutoff = value.to_f
    when "--prop_gaps", "--prop_gap"
      prop_of_gaps_cutoff = value.to_f
    when "--intron_sliding", "--intron_sliding_cutoff"
      intron_sliding_cutoff = value.to_i
  end
end


#################################################
if not indir.nil?
  aln_info = read_indir(indir, aln_info, extensions)
else
  aln_info = read_aln_file(infile, aln_info)
end

gff_files.each do |gff_file|
  cds_info = read_gff(gff_file, cds_info)
end

intron_info = parse_cds_info(cds_info)

aln_posi_info, aligned_posi, real_aligned_posi = get_aln_info(aln_info)

well_aligned_posi = get_well_aligned_posi(aln_info, aligned_posi, prop_of_diff_sites_cutoff, prop_of_gaps_cutoff)

intron_in_aln = get_intron_in_aln(aligned_posi, real_aligned_posi, intron_info)


#################################################
intron_posi = Hash.new{|h,k|h[k]={}}
intron_in_aln.each_pair.sort_by{|aln_posi,v|aln_posi}.each do |aln_posi, v|
  if v.size == aln_info.keys.size
    intron_posi["shared"][aln_posi] = []
    v.each do |seq_title|
      real_posi = real_aligned_posi[seq_title][aln_posi]
      intron_posi["shared"][aln_posi] << [seq_title, real_posi, intron_info[seq_title][real_posi]]
    end
  else
    intron_posi["diff"][aln_posi] = Array.new if not intron_posi["diff"].include?(aln_posi)
    v.each do |seq_title|
      real_posi = real_aligned_posi[seq_title][aln_posi]
      intron_posi["diff"][aln_posi] << [seq_title, real_posi, intron_info[seq_title][real_posi]]
    end
  end
end

#p intron_posi


#################################################
intron_sliding_aln_posi = Hash.new
intron_posi['diff'].keys.sort.each do |aln_posi|
  titles = intron_posi['diff'][aln_posi].sort_by{|i|i[0].to_i}
  (aln_posi-intron_sliding_cutoff).upto(aln_posi+intron_sliding_cutoff).each do |i|
    next if not intron_posi['diff'].include?(i)
    titles2 = intron_posi['diff'][i].sort_by{|j|j[0].to_i}
    if titles != titles2
      intron_sliding_aln_posi[aln_posi] = ""
    end
  end
end


puts aln_info.keys.sort.join("-")
intron_posi["shared"].each_pair do |aln_posi, v|
  next if not well_aligned_posi.include?(aln_posi)
  print [aln_posi, "shared"].map{|i|i.to_s}.join("\t") + "\t"
  puts v.map{|i|i.join("\t")}.join("\t")
end


intron_posi['diff'].each_pair do |aln_posi, v|
  next if not well_aligned_posi.include?(aln_posi)
  next if intron_sliding_aln_posi.include?(aln_posi)
  print [aln_posi, "diff"].map{|i|i.to_s}.join("\t") + "\t"
  puts v.map{|i|i.join("\t")}.join("\t")
  #puts [aln_posi, "diff", v[0][0], v[0][1], v[0][2]].map{|i|i.to_s}.join("\t")
end

puts


