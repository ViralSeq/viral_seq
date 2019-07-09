# functions added to Class::String for direct operation on sequence as a String object

class String

  # reverse complement
  # @return [String] reverse complement sequence
  # @example Reverse complement
  #   "ACAGA".rc
  #   => "TCTGT"

  def rc
      self.reverse.tr("ACTG","TGAC")
  end

  # mutate a nt sequence (String class) randomly
  # @param error_rate [Float] define an error rate for mutation, default to `0.01`
  # @return [String] mutated sequence as String
  # @example mutate a sequence at an error rate of 0.05
  #   seq = "TGGAAGGGCTAATTCACTCCCAACGAAGACAAGATATCCTTGATCTGTGGATCTACCACACACAAGGCTACTTCCCTG"
  #   seq.mutation(0.05)
  #   => "TGGAAGGGCTAATGCACTCCCAACGAAGACACGATATCCTTGATCTGTGGATCTACGACACACAAGGCTGCTTCCCTG"

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

  # parse the nucleotide sequences as a String object
  #   and return a Regexp object for possible matches
  # @return [Regexp] as possible matches
  # @example parse a sequence with ambiguities
  #   "ATRWCG".nt_parser
  #   => /AT[A|G][A|T]CG/

  def nt_parser
    match = ""
    self.each_char.each do |base|
      base_array = base.to_list
      if base_array.size == 1
        match += base_array[0]
      else
        pattern = "[" + base_array.join("|") + "]"
        match += pattern
      end
    end
    Regexp.new match
  end

  # parse IUPAC nucleotide ambiguity codes (W S M K R Y B D H V N) as String if String.size == 1
  # @return [Array] parsed nt bases
  # @example parse IUPAC `R`
  #   'R'.to_list
  #   => ["A", "G"]

  def to_list
    list = []
    case self.upcase
    when /[A|T|C|G]/
      list << self
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

  # compare two sequences as String objects, two sequence strings need to aligned first
  # @param seq2 [String] the sequence string to compare with
  # @return [Integer] the total number of differences as integer
  # @example compare two sequence strings, without alignment and with alignment
  #   seq1 = 'AAGGCGTAGGAC'
  #   seq2 = 'AAGCTTAGGACG'
  #   seq1.compare_with(seq2) # no alignment
  #   => 8
  #   aligned_seqs = ViralSeq::Muscle.align(seq1,seq2) # align using MUSCLE
  #   aligned_seqs[0].compare_with(aligned_seqs[1])
  #   => 4

  def compare_with(seq2)
    seq1 = self
    length = seq1.size
    diff = 0
    (0..(length-1)).each do |position|
      nt1 = seq1[position]
      nt2 = seq2[position]
      diff += 1 unless nt1 == nt2
    end
    return diff
  end
end
