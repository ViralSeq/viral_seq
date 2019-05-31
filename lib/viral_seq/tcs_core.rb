# viral_seq/tcs_core
# core functions for TCS and DR pipeline

module ViralSeq
  def self.calculate_pid_cut_off(m)
    n = 0
    if m <= 10
      n = 2
    elsif m <= 8500
      n = -1.24*10**-21*m**6 + 3.53*10**-17*m**5 - 3.90*10**-13*m**4 + 2.12*10**-9*m**3 - 6.06*10**-6*m**2 + 1.80*10**-2*m + 3.15
    else
      n = 0.0079 * m + 9.4869
    end
    n = n.round
    n = 2 if n < 3
    return n
  end


  # create one consensus sequence from a sequence array with an optional majority cut-off for mixed bases.
  # example:
  # position with 15% "A" and 85% "G" will be called as "G" with 20% cut-off and as "R" with 10% cut-off.
  def self.consensus(seq_array, cutoff = 0.5)
    seq_length = seq_array[0].size
    seq_size = seq_array.size
    consensus_seq = ""
    (0..(seq_length - 1)).each do |position|
      all_base = []
      seq_array.each do |seq|
        all_base << seq[position]
      end
      base_count = ViralSeq.count(all_base)
      max_base_list = []

      base_count.each do |k,v|
        if v/seq_size.to_f >= cutoff
          max_base_list << k
        end
      end
      consensus_seq += ViralSeq.call_consensus_base(max_base_list)
    end
    return consensus_seq
  end

  #call consensus nucleotide
  def self.call_consensus_base(base_array)
    if base_array.size == 1
       base_array[0]
    elsif base_array.size == 2
      case base_array.sort!
      when ["A","T"]
        "W"
      when ["C","G"]
        "S"
      when ["A","C"]
        "M"
      when ["G","T"]
         "K"
      when ["A","G"]
         "R"
      when ["C","T"]
         "Y"
      else
         "N"
      end

    elsif base_array.size == 3
      case base_array.sort!
      when ["C","G","T"]
         "B"
      when ["A","G","T"]
         "D"
      when ["A","C","T"]
         "H"
      when ["A","C","G"]
         "V"
      else
         "N"
      end
    else
       "N"
    end
  end

  

end
