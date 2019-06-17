# viral_seq/nt_variation

# contain functions to cacluate shannon's entropy, pairwise diversity, and TN93 distance
#   ViralSeq::shannons_entropy
#   ViralSeq::nucleotide_pi
#   ViralSeq::TN93

# ViralSeq.shannons_entropy(sequences)
#   # calculate Shannon's entropy, Euler's number as the base of logarithm
#   # https://en.wikipedia.org/wiki/Entropy_(information_theory)
#   # input sequences alignment in Array [:seq1, :seq2, ...] or Hash [:name => :sequence] object
#   # it works on both nt and aa sequences
#   # return a Hash object for entropy score at each position in the alignment {:position => :entropy, ...}
#   # position starts at 1.
# =Usage
#   # example
#   sequence_file = 'spec/sample_files/sample_sequence_alignment_for_entropy.fasta'
#   sequence_hash = ViralSeq.fasta_to_hash(sequence_file)
#   entropy_hash = ViralSeq.shannons_entropy(sequence_hash)
#   entropy_hash[3]
#   => 0.0
#   entropy_hash[14].round(3)
#   => 0.639
#   # This example is the sample input of LANL Entropy-One
#   # https://www.hiv.lanl.gov/content/sequence/ENTROPY/entropy_one.html?sample_input=1

# ViralSeq.nucleotide_pi(sequences)
#   # Function to calculate nucleotide diversity π.
#   # Refer to https://en.wikipedia.org/wiki/Nucleotide_diversity
#   # input sequences alignment in Array [:seq1, :seq2, ...] or Hash [:name => :sequence] object
#   # nt sequence only
#   # return π as a Float object
# =Usage
#   # example
#   sequences = %w{ AAGGCCTT ATGGCCTT AAGGCGTT AAGGCCTT AACGCCTT AAGGCCAT }
#   ViralSeq.nucleotide_pi(sequences)
#   =>

module ViralSeq
  # calculate Shannon's entropy, Euler's number as the base of logarithm

  def self.shannons_entropy(sequences)
    sequences = if sequences.is_a?(Hash)
                  sequences.values
                elsif sequences.is_a?(Array)
                  sequences
                else
                  raise ArgumentError.new("Wrong type of input sequences. it has to be Hash or Array object")
                end
    entropy_hash = {}
    seq_l = sequences[0].size
    (0..(seq_l - 1)).each do |position|
      element = []
      sequences.each do |seq|
        element << seq[position]
      end
      entropy = 0
      element.delete('*')
      element_size = element.size
      ViralSeq.count(element).each do |_k,v|
        p = v/element_size.to_f
        entropy += (-p * Math.log(p))
      end
      entropy_hash[(position + 1)] = entropy
    end
    return entropy_hash
  end

  # nucleotide pairwise diversity
  def self.nucleotide_pi(sequences)
    sequences = if sequences.is_a?(Hash)
                  sequences.values
                elsif sequences.is_a?(Array)
                  sequences
                else
                  raise ArgumentError.new("Wrong type of input sequences. it has to be Hash or Array object")
                end
    seq_length = sequences[0].size - 1
    nt_position_hash = {}
    (0..seq_length).each do |n|
      nt_position_hash[n] = []
      sequences.each do |s|
        nt_position_hash[n] << s[n]
      end
    end
    diver = 0
    com = 0
    nt_position_hash.each do |_p,nt|
      nt.delete_if {|n| n =~ /[^A|^C|^G|^T]/}
      next if nt.size == 1
      nt_count = ViralSeq.count(nt)
      combination = (nt.size)*(nt.size - 1)/2
      com += combination
      a = nt_count["A"]
      c = nt_count["C"]
      t = nt_count["T"]
      g = nt_count["G"]
      div = a*c + a*t + a*g + c*t + c*g + t*g
      diver += div
    end
    pi = (diver/com.to_f).round(5)
    return pi
  end

  # TN93 distance function. Input: sequence array, output hash: diff => counts
  def self.TN93(sequence_array = [])
    diff = []
    seq_hash = ViralSeq.count(sequence_array)
    seq_hash.values.each do |v|
      comb = v * (v - 1) / 2
      comb.times {diff << 0}
    end

    seq_hash.keys.combination(2).to_a.each do |pair|
      s1 = pair[0]
      s2 = pair[1]
      diff_temp = ViralSeq.compare_two_seq(s1,s2)
      comb = seq_hash[s1] * seq_hash[s2]
      comb.times {diff << diff_temp}
    end

    count_diff = ViralSeq.count(diff)
    out_hash = Hash.new(0)
    Hash[count_diff.sort_by{|k,_v|k}].each do |k,v|
      out_hash[k] = v
    end
    return out_hash
  end
end
