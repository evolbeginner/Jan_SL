#! /bin/env ruby

require 'String'

###################################################
class FilterBlastOutput
  def initialize(blast_output,blast_outfmt,seq_objs={})
    @blast_output=blast_output
    @blast_outfmt=blast_outfmt
    @seq_objs=seq_objs if ! seq_objs.empty?
  end

  def filter(e_value_cutoffs=[0,10000], identity_range=[-1,100], aligned_length_range=[0,10000], regions_range={}, action='rm', outdir=nil, donor_args={})
    @e_value_cutoffs = e_value_cutoffs.map{|i| i.to_f}
    @identity_range = identity_range.map{|i| i.to_f}
    @aligned_length_range = aligned_length_range.map{|i| i.to_i}
    @regions_range=regions_range
    pass=false
    case @blast_outfmt
      when '1', '6'
        start, stop = Hash.new, Hash.new
        fh=File.open(@blast_output ,'r')
        while(line=fh.gets) do
          next if line =~ /^#/;
          line.chomp!
          subject, identity, aligned_length, start['query'], stop['query'],start['suject'], stop['subject'], e_value = line.split("\t").values_at(1,2,3,6,7,8,9,10).map{|i| i.numeric? ? i.to_f : i}
          regions = {
                  'query'   =>  {'start' =>  start['query'],  'stop'  =>  stop['query']} ,
                  'subject' =>  {'start' =>  start['subject'],'stop'  =>  stop['subject']}
                  }
          aligned_length=aligned_length.to_i
          next if ! check_pass?(e_value, identity, aligned_length, regions)
          (next if ! check_donor_pass?(donor_args,aligned_length,e_value,start,stop,subject)) if ! donor_args.empty?
          pass=true
        end
      when '0'
        fh=File.open(@blast_output, 'r')
        while(line=fh.gets) do
          pass_k=0
          line.chomp!
          if line =~ /Expect\s+=\s+(.+)/ then
            e_value=$1.to_f
            line=fh.gets # read one line
            iden_length,aligned_length = $1.to_i,$2.to_i if line =~ /Identities\s+=\s+ (\d+)\/(\d+)/x
            identity=iden_length.to_f/aligned_length*100
            pass=check_pass?(e_value, identity, aligned_length, regions, donor_args)
            break if pass
          end
          #pass=true and next if e_value <= e_value_cutoffs[1] and e_value > e_value_cutoffs[0]
        end
    end

    if ! pass then
      case action
        when 'rm'
          `rm #{@blast_output}`
      end
    else
      case action
        when 'mv'
          `mv #{@blast_output} #{outdir}`
        when 'cp'
          `cp #{@blast_output} #{outdir}`
      end
    end
  end


  def check_donor_pass?(donor_args={},aligned_length,e_value,start,stop,subject)
    pass=false
    raise "seq_objs not given when performing donor_checks!" if @seq_objs.nil? or @seq_objs.empty?
    if ! donor_args.empty?
      ranges_pass=true
      if donor_args.include?('length_range')
        ranges_pass=false if donor_args['length_range'].include?(aligned_length)
      end
      if donor_args.include?('e_value_range')
        if (donor_args['e_value_range']['min']<=e_value and donor_args['e_value_range']['max']>=e_value)
          ranges_pass=false
        end
      end
      return true if ranges_pass == true

      if @seq_objs.include?(subject)
        if @seq_objs[subject].seq[stop['subject']-2,2] =~ /#{donor_args['donor_site']}/i
          pass=true
        end
      else
        raise "#{subject} is not found in seq_objs!"
      end
    end
    puts start['query'].to_i
    return(pass)
  end


  def check_pass?(e_value, identity, aligned_length, regions={})
    pass=false
    pass_k=0
    (e_value <= @e_value_cutoffs[1] and e_value > @e_value_cutoffs[0])? (pass_k+=1) : (return pass)
    (identity <= @identity_range[1] and identity > @identity_range[0])? (pass_k+=1) : (return pass)
    (aligned_length <= @aligned_length_range[1] and aligned_length > @aligned_length_range[0])? (pass_k+=1) : (return pass)
    
    if ! @regions_range.empty?
      regions.each_pair do |k1,v1| # k: 'query','subject'
        next if ! @regions_range.include? k1
        region_pass=false
        regions[k1].each_pair do |k2,v2|
          if (@regions_range[k1]['start']..@regions_range[k1]['stop']).include? v2
            region_pass=true
          end
        end
        if ! region_pass
          return pass
        end
      end
    end
    
    pass=true
    return pass
  end
end

