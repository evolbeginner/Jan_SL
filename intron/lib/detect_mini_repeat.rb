#! /bin/env ruby


#################################################################
class Repeat_finder
  attr_accessor :min_length
  attr_reader :counter

  def initialize()
    ;
  end

  def length(min_length)
    @min_length = min_length
  end

  def find(sequence1, sequence2)
    is_found = false
    1.upto(sequence1.size) do |i|
      seq = sequence1[i-1, @min_length]
      if seq.size >= @min_length
        if sequence2 =~ /#{seq}/
          is_found = true
          break
        end
      end
    end
    return(is_found)
  end
end


#################################################################

