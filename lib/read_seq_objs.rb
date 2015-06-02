require 'bio'

def read_seq_objs(seq_file,prefix=nil,suffix=nil)
  seq_objs=Hash.new
  Bio::FlatFile.open(seq_file).each_entry do |f|
    definition = f.definition
    definition.sub!(/^#{prefix}/,'') if prefix
    definition.sub!(/#{suffix}$/,'') if suffix
    seq_objs[definition]=f
  end
  return seq_objs
end
