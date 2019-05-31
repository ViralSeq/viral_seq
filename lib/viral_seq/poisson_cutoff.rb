# viral_seq/poisson_cutoff.rb

# define Poission cut-off for minority variants. (Ref: Zhou, et al. J Virol 2015)
# Input is a sequence array, the redisual sequencing error rate (default = 0.0001), and a fold cut-off to determine the cut-off
# example: cut-off = 2 means that mutations appear at least 2 times are very likely to be a true mutation instead of residual methods errors.
module ViralSeq

  def self.poisson_minority_cutoff(seq_array, error_rate = 0.0001, fold_cutoff = 20)
    if seq_array.size == 0
      return 0
    else
      cut_off = 1
      l = seq_array[0].size
      rate = seq_array.size * error_rate
      count_mut = ViralSeq.variant_for_poisson(seq_array)
      max_count = count_mut.keys.max
      poisson_hash = ViralSeq.poisson_distribution(rate, max_count)

      poisson_hash.each do |k,v|
        cal = l * v
        obs = count_mut[k] ? count_mut[k] : 0
        if obs >= fold_cutoff * cal
          cut_off = k
          break
        end
      end
      return cut_off
    end
  end

  # Input sequence array. output Variant distribution for Poisson cut-off
  def self.variant_for_poisson(seq)
    seq_size = seq.size
    l = seq[0].size - 1
    var = []
    (0..l).to_a.each do |pos|
      nt = []
      seq.each do |s|
        nt << s[pos]
      end
      count_nt = ViralSeq.count(nt)
      v = seq_size - count_nt.values.max
      var << v
    end
    var_count = count(var)
    var_count.sort_by{|key,_value|key}.to_h
  end

end
