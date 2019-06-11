# lib/sequence.rb
# Includes functions for sequence operations


module ViralSeq

  # array for all amino acid one letter abbreviations
  AMINO_ACID_LIST = ["A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y", "*"]

  # sequence class
  # =USAGE
  #   # create a sequence object
  #   seq = ViralSeq::Sequence.new('my_sequence', 'ACCTAGGTTCGGAGC')
  #
  #   # print dna sequence
  #   puts seq.dna_sequence
  #
  #   # reserce complement sequence of DNA sequence, return as a string
  #   seq.rev_complement
  #
  #   # change @dna_sequence to reverse complement DNA sequence
  #   seq.rev_complement!
  #
  #   # generate amino acid sequences. either return string or array.
  #   # starting codon option 0, 1, 2 for 1st, 2nd, 3rd reading frame.
  #   # if sequence contains ambiguities, Sequence.get_aa_array will return all possible amino acids.
  #   seq.get_aa_sequence
  #   # or
  #   seq.get_aa_array
  #
  #   # print amino acid sequence
  #   puts seq.aa_sequence

  class Sequence
    def initialize (name = ">sequence",dna_sequence ="")
      @name = name
      @dna_sequence = dna_sequence.upcase
      @aa_sequence = ""
      @aa_array = []
    end

    attr_accessor :name, :dna_sequence, :aa_sequence, :aa_array

    def rev_complement
      @dna_sequence.reverse.upcase.tr('ATCG','TAGC')
    end
    def rev_complement!
      @dna_sequence = @dna_sequence.reverse.upcase.tr('ATCG','TAGC')
    end

    def get_aa_sequence(initial_position = 0)
      @aa_sequence = ""
      require_sequence = @dna_sequence[initial_position..-1]
      base_array = []
      require_sequence.each_char {|base| base_array << base}
      while (base_array.length>=3) do
        base_3= ""
        3.times {base_3 += base_array.shift}
        @aa_sequence << amino_acid(base_3)
      end
      return @aa_sequence
    end

    # get amino acid calls, return a array.keep ambiguity calls.
    def get_aa_array(initial_position = 0)
      @aa_array = []
      require_sequence = @dna_sequence[initial_position..-1].tr('-','N')
      base_array = []
      require_sequence.each_char {|base| base_array << base}
      while (base_array.length>=3) do
        base_3= ""
        3.times{base_3 += base_array.shift}
        @aa_array<< ViralSeq.amino_acid_2(base_3)
      end
      return @aa_array
    end
    def dna_length
      @dna_sequence.length
    end
    def aa_length
      @aa_sequence.length
    end
  end

  # generate amino acid abbreviations from 3 bases, ambiguity will return "#"
  def self.amino_acid (bases)
    case bases
    when /^TT[TCY]$/
      return "F"
    when /^TT[AGR]$/
      return "L"
    when /^CT.$/
      return "L"
    when /^AT[TCAHYWM]$/
      return "I"
    when "ATG"
      return "M"
    when /^GT.$/
      return "V"
    when /^TC.$/
      return "S"
    when /^CC.$/
      return "P"
    when /^AC.$/
      return "T"
    when /^GC.$/
      return "A"
    when /^TA[TCY]$/
      return "Y"
    when /^TA[AGR]$/
      return "*"
    when /^T[GR]A$/
      return "*"
    when /^CA[TCY]$/
      return "H"
    when /^CA[AGR]$/
      return "Q"
    when /^AA[TCY]$/
      return "N"
    when /^AA[AGR]$/
      return "K"
    when /^GA[TCY]$/
      return "D"
    when /^GA[AGR]$/
      return "E"
    when /^TG[TCY]$/
      return "C"
    when "TGG"
      return "W"
    when /^CG.$/
      return "R"
    when /^AG[TCY]$/
      return "S"
    when /^[AM]G[AGR]$/
      return "R"
    when /^GG.$/
      return "G"
    when /^[ATW][CGS][CTY]$/
      return "S"
    when /^[TCY]T[AGR]$/
      return "L"
    else
      return "#"
    end
  end

  # keep ambiguities, return all possible amino acids.

  def self.amino_acid_2 (bases)
    bases_to_aa = []
    aa_list = []
    base1 = to_list(bases[0])
    base2 = to_list(bases[1])
    base3 = to_list(bases[2])
    l1 = base1.size - 1
    l2 = base2.size - 1
    l3 = base3.size - 1
    (0..l1).each do |n1|
      b1 = base1[n1]
      (0..l2).each do |n2|
        b2 = base2[n2]
        (0..l3).each do |n3|
          b3 = base3[n3]
          bases_all = b1 + b2 + b3
          bases_to_aa << bases_all
        end
      end
    end

    bases_to_aa.each do |bases|
    case bases
    when /^TT[TCY]$/
      aa =  "F"
    when /^TT[AGR]$/
      aa =  "L"
    when /^CT.$/
      aa =  "L"
    when /^AT[TCAHYWM]$/
      aa =  "I"
    when "ATG"
      aa =  "M"
    when /^GT.$/
      aa =  "V"
    when /^TC.$/
      aa =  "S"
    when /^CC.$/
      aa =  "P"
    when /^AC.$/
      aa =  "T"
    when /^GC.$/
      aa =  "A"
    when /^TA[TCY]$/
      aa =  "Y"
    when /^TA[AGR]$/
      aa =  "*"
    when /^T[GR]A$/
      aa =  "*"
    when /^CA[TCY]$/
      aa =  "H"
    when /^CA[AGR]$/
      aa =  "Q"
    when /^AA[TCY]$/
      aa =  "N"
    when /^AA[AGR]$/
      aa =  "K"
    when /^GA[TCY]$/
      aa =  "D"
    when /^GA[AGR]$/
      aa =  "E"
    when /^TG[TCY]$/
      aa =  "C"
    when "TGG"
      aa =  "W"
    when /^CG.$/
      aa =  "R"
    when /^AG[TCY]$/
      aa =  "S"
    when /^[AM]G[AGR]$/
      aa =  "R"
    when /^GG.$/
      aa =  "G"
    when /^[ATW][CGS][CTY]$/
      aa =  "S"
    when /^[TCY]T[AGR]$/
      aa =  "L"
    else
      aa =  "-"
    end
    aa_list << aa
  end
    aa_out = aa_list.uniq.join('/')
    return aa_out
  end

  # parse ambiguity bases, aka %w{W S M K R Y B D H V N}

  def self.to_list(base = "")
    list = []
    case base
    when /[A|T|C|G]/
      list << base
    when "W"
      list = ['A','T']
    when "S"
      list = ['C','G']
    when "M"
      list = ['A','C']
    when 'K'
      list = ['G','C']
    when 'R'
      list = ['A','G']
    when 'Y'
      list = ['C','T']
    when 'B'
      list = ['C','G','T']
    when 'D'
      list = ['A','G','T']
    when 'H'
      list = ['A','C','T']
    when 'V'
      list = ['A','C','G']
    when 'N'
      list = ['A','T','C','G']
    end
    return list
  end

  # ViralSeq.uniq_sequence_hash(input_sequence_hash, master_sequence_tag)
  # collapse sequence hash to unique sequence hash.
  # input_sequence_hash is a sequence hash {:name => :sequence, ...}
  # master_sequence_tag is the master tag for unique sequences
  # sequences will be named as (master_sequence_tag + "_" + Integer)
  
  def self.uniq_sequence_hash(seq = {}, sequence_name = "sequence")
    uni = ViralSeq.count(seq.values)
    new_seq = {}
    n = 1
    uni.each do |s,c|
      name = ">" + sequence_name + "_" + n.to_s + "_" + c.to_s
      new_seq[name] = s
      n += 1
    end
    return new_seq
  end

end

# functions added to Class::String for direct operation on sequence if it is a String object
# String.rc
#   # reverse complement
#   # example
#   "ACAGA".rc
#   => "TCTGT"
#
# String.mutation(error_rate)
#   # mutate a nt sequence (String class) randomly
#   # must define error rate, default value 0.01, aka 1%
# =USAGE
#   # example
#   seq = "TGGAAGGGCTAATTCACTCCCAACGAAGACAAGATATCCTTGATCTGTGGATCTACCACACACAAGGCTACTTCCCTG"
#   seq.mutation(0.05)
#   => "TGGAAGGGCTAATGCACTCCCAACGAAGACACGATATCCTTGATCTGTGGATCTACGACACACAAGGCTGCTTCCCTG"
#
# String.nt_parser
#   # parse the nucleotide sequences as a String object and return a Regexp object for possible matches
# =USAGE
#   "ATRWCG".nt_parser
#   => /AT[A|G][A|T]CG/

class String
    # direct function of calling reverse complement on String class
  def rc
      self.reverse.tr("ACTG","TGAC")
  end

  def mutation(error_rate = 0.01)
    new_string = ""
    self.split("").each do |nt|
      pool = ["A","C","T","G"]
      pool.delete(nt)
      s = error_rate * 10000
      r = rand(10000)
      if r < s
        nt = pool.sample
      end
      new_string << nt
    end
    return new_string
  end

  def nt_parser
    match = ""
    self.each_char.each do |base|
      base_array = ViralSeq.to_list(base)
      if base_array.size == 1
        match += base_array[0]
      else
        pattern = "[" + base_array.join("|") + "]"
        match += pattern
      end
    end
    Regexp.new match
  end
end
