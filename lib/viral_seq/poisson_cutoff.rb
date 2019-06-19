# viral_seq/poisson_cutoff.rb
# define Poission cut-off for minority variants.
# (Ref: Zhou, et al. J Virol 2015)
#   ViralSeq::poisson_minority_cutoff

# ViralSeq.poisson_minority_cutoff(sequences, error_rate, fold_cutoff)
#   # sequences: input sequences alignment in Array [:seq1, :seq2, ...] or Hash [:name => :sequence] object
#   # error_rate: the redisual sequencing error rate (default = 0.0001),
#   # fold_cutoff: a fold cut-off to determine poisson minority cut-off. default = 20. i.e. <5% mutations from randome method error.
#   # example: cut-off = 2 means that mutations appear at least 2 times are very likely to be a true mutation instead of residual methods errors.
# =Usage
#   sequence_file = 'spec/sample_files/sample_sequence_for_poisson.fasta'
#   sequences = ViralSeq.fasta_to_hash(sequence_file)
#   ViralSeq.poisson_minority_cutoff(sequences)
#   => 2


module ViralSeq

  def self.poisson_minority_cutoff(sequences, error_rate = 0.0001, fold_cutoff = 20)
    sequences = if sequences.is_a?(Hash)
                  sequences.values
                elsif sequences.is_a?(Array)
                  sequences
                else
                  raise ArgumentError.new("Wrong type of input sequences. it has to be Hash or Array object")
                end
    if sequences.size == 0
      return 0
    else
      cut_off = 1
      l = sequences[0].size
      rate = sequences.size * error_rate
      count_mut = ViralSeq.variant_for_poisson(sequences)
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
