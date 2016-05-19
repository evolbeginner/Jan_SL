#! /bin/env ruby

BEGIN{
  file_name=__FILE__
  IO.popen("ls -lh #{__FILE__}", "r").readlines.each do |line|
    line.chomp!
    if line =~ /\-\> (.+)/ then
      file_name=$1
    end
  end
  $: << [File.dirname(file_name),'lib'].join('/')
}


###############################################
require 'bio'
require 'getoptlong'
require 'filter_blast_output'

input=nil
outdir=nil
sl_seq=nil
region=nil
start, stop = nil, nil
is_remove_null=false
is_quiet=false
blast_outfmt=0
e_value_cutoff=10000
add_args=nil

############################################################################
class BL2SEQ
  @@counter=nil
  attr_accessor :seq, :sl_seq, :title, :outdir, :add_args
  attr_reader :counter

  def initialize()
    @seq=nil
    @title=nil
    @outdir=nil
    @bl2seq_output=nil
    @add_args=nil
  end

  def generate_in_fasta
    in_fasta=(@outdir+"/tmp_haha")
    #bl2seq_output=@outdir+"/bl2seq_output"+@@counter.to_s
    title.gsub!('|','--')
    bl2seq_output=@outdir+"/"+title
    fh=File.open(in_fasta, 'w')
    fh.puts ">"+title+"\n"+seq
    fh.close
    return ([in_fasta,bl2seq_output])
  end

  def find_repeats(is_remove_null, blast_outfmt=0, e_value_cutoff)
    @@counter = @@counter.nil? ? (1) : (@@counter+1)
    puts @@counter if @@counter%1000==0
    in_fasta, bl2seq_output = generate_in_fasta
    sl_counter = 0
    bl2seq_outputs_lines = Array.new
    Bio::FlatFile.open(sl_seq,'r').each_entry do |f|
      sl_counter += 1
      in_tmp_sl = File.join(@outdir,"tmp_SL.fas")
      out_SL = File.open(in_tmp_sl, 'w')
      out_SL.puts ">#{f.definition}\n#{f.seq}"
      out_SL.close
      cmd="bl2seq -i #{in_tmp_sl} -j #{in_fasta} -o #{bl2seq_output}"
      cmd += " -p blastn -D #{blast_outfmt} -e #{e_value_cutoff} -W 7 -S 1 -F F"
      cmd += ' ' + add_args if add_args
      `#{cmd}`
      bl2seq_outputs_lines.push File.open(bl2seq_output, 'r').readlines
      File.delete(in_tmp_sl)
    end

    fh = File.open(bl2seq_output,'w')
    fh.puts bl2seq_outputs_lines
    fh.close

    if is_remove_null then
      case blast_outfmt
        when '1', '6'
          `if ! grep -P '^[^#]' #{bl2seq_output}; then rm #{bl2seq_output}; fi`
        when '0'
          `if ! grep "Expect" #{bl2seq_output}; then rm #{bl2seq_output}; fi`
      end
    end

    File.delete(in_fasta)
  end
end


############################################################################
opts = GetoptLong.new(
  ['-i', '--in', '--input', GetoptLong::REQUIRED_ARGUMENT],
  ['-o', '--outdir', GetoptLong::REQUIRED_ARGUMENT],
  ['--SL', '--SL_seq', GetoptLong::REQUIRED_ARGUMENT],
  ['--region', GetoptLong::REQUIRED_ARGUMENT],
  ['--quiet', GetoptLong::NO_ARGUMENT],
  ['--remove_null', GetoptLong::NO_ARGUMENT],
  ['--blast_outfmt', GetoptLong::REQUIRED_ARGUMENT],
  ['--e_value', GetoptLong::REQUIRED_ARGUMENT],
  ['--args','--add_args', GetoptLong::REQUIRED_ARGUMENT],
)

opts.each do |opt, value|
  case opt
    when '-i', '--in', '--input'
      input=value
    when '-o', '--outdir'
      outdir=value
    when '--SL', '--SL_seq'
      sl_seq=value
    when '--region'
      region=value
    when '--quiet'
      is_quiet=true
    when '--remove_null'
      is_remove_null=true
    when '--blast_outfmt'
      blast_outfmt=value
    when '--e_value'
      e_value_cutoff=value.to_f
    when '--args', 'add_args'
      add_args=value
  end
end

raise "Params are not complete! Please check." if not input or not outdir

`rm -rf #{outdir}` if File.exists?(outdir)
`mkdir -p #{outdir}`
if ! region.nil? then
  if region=~/(\d+),(\d+)/ then
    start, stop = $1.to_i, $2.to_i
  else
    raise "region format is not correct!"
  end
end

############################################################################
fh = Bio::FlatFile.open(input)
fh.each_entry do |f|
  seq_obj = Bio::Sequence::NA.new(f.seq)
  seq = region.nil? ? f.seq : seq_obj.subseq(start,stop)
  bl2seq_obj=BL2SEQ.new()
  bl2seq_obj.title=f.definition
  bl2seq_obj.seq=seq
  bl2seq_obj.sl_seq=sl_seq
  bl2seq_obj.outdir=outdir
  bl2seq_obj.add_args=add_args if add_args
  bl2seq_obj.find_repeats(is_remove_null, blast_outfmt, e_value_cutoff)
  if not is_quiet then
    puts [f.definition,bl2seq_obj.counter].join("\t")
  end
end

