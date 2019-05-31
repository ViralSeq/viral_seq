# viral_seq/shannons_entropy

# calculate Shannon's entropy, Euler's number as the base of logarithm
# https://en.wikipedia.org/wiki/Entropy_(information_theory)

module ViralSeq
  def self.shannons_entropy(sequences)
    entropy_hash = {}
    seq_l = sequences[0].size
    seq_size = sequences.size
    (0..(seq_l - 1)).each do |position|
      element = []
      sequences.each do |seq|
        element << seq[position]
      end
      entropy = 0
      ViralSeq.count(element).each do |_k,v|
        p = v/seq_size.to_f
        entropy += (-p * Math.log(p))
      end
      entropy_hash[(position + 1)] = entropy
    end
    return entropy_hash
  end
end
