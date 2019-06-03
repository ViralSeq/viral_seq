# viral_seq/tcs_core
# core functions for TCS and DR pipeline

module ViralSeq

  def self.calculate_pid_cut_off(m, error_rate = 0.02)
    if m <= 10
      return 2
    end
    n = 0
    case error_rate
    when 0...0.0075
      n = -9.59*10**-27*m**6 + 3.27*10**-21*m**5 - 3.05*10**-16*m**4 + 1.2*10**-11*m**3 - 2.19*10**-7*m**2 + 0.004044*m + 2.273
    when 0.0075...0.015
      n = 1.09*10**-26*m**6 + 7.82*10**-22*m**5 - 1.93*10**-16*m**4 + 1.01*10**-11*m**3 - 2.31*10**-7*m**2 + 0.00645*m + 2.872
    when 0.015..0.03
      if m <= 8500
        n = -1.24*10**-21*m**6 + 3.53*10**-17*m**5 - 3.90*10**-13*m**4 + 2.12*10**-9*m**3 - 6.06*10**-6*m**2 + 1.80*10**-2*m + 3.15
      else
        n = 0.0079 * m + 9.4869
      end
    else
      raise ArgumentError.new('Error_rate has be between 0 to 0.03')
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

  # call consensus nucleotide
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

  # generate all Primer ID combinations given the length of Primer ID
  def self.generate_primer_id_pool(l=8)
    nt = ['A','T','C','G']
    pid_pool = ['A','T','C','G']
    (l-1).times do
      pid_pool = pid_pool.product(nt)
      pid_pool.collect! do |v|
        v.join("")
      end
    end
    return pid_pool
  end

  # compare two primer ID sequences. If they differ in x base, return boolean value "TURE", else, return boolean value "FALSE"
  def self.similar_pid?(pid1="",pid2="", x=0)
    l = pid1.size
    m = l - x
    n = 0
    if pid1.size != pid2.size
      return false
    else
      (0..(pid1.size - 1)).each do |k|
        if pid1[k] == pid2[k]
          n += 1
        end
      end
      if n >= m
        return true
      else
        return false
      end
    end
  end

  # primer with ambiguities to match
  def self.primer_match (primer = "")
    match = ""
    primer.each_char.each do |base|
      base_array = to_list(base)
      if base_array.size == 1
        match += base_array[0]
      else
        pattern = "[" + base_array.join("|") + "]"
        match += pattern
      end
    end
    return match
  end

  # #primer with ambiguities to match
  # def primer_match (primer = "")
  #   match = ""
  #   primer.each_char.each do |base|
  #     base_array = to_list(base)
  #     if base_array.size == 1
  #       match += base_array[0]
  #     else
  #       pattern = "[" + base_array.join("|") + "]"
  #       match += pattern
  #     end
  #   end
  #   return match
  # end

  def self.filter_similar_pid(sequence_file = "", cutoff = 10)
    seq = ViralSeq.fasta_to_hash(sequence_file)
    uni_seq = seq.values.uniq
    uni_seq_pid = {}
    uni_seq.each do |k|
      seq.each do |name,s|
        name = name[1..-1]
        if k == s
          if uni_seq_pid[k]
            uni_seq_pid[k] << [name.split("_")[0],name.split("_")[1]]
          else
            uni_seq_pid[k] = []
            uni_seq_pid[k] << [name.split("_")[0],name.split("_")[1]]
          end
        end
      end
    end

    dup_pid = []
    uni_seq_pid.values.each do |v|
      next if v.size == 1
      pid_hash = Hash[v]
      list = pid_hash.keys
      list2 = Array.new(list)
      pairs = []

      list.each do |k|
        list2.delete(k)
        list2.each do |k1|
          pairs << [k,k1]
        end
      end


      pairs.each do |p|
        pid1 = p[0]
        pid2 = p[1]
        if ViralSeq.similar_pid?(pid1,pid2,1)
          n1 = pid_hash[pid1].to_i
          n2 = pid_hash[pid2].to_i
          if n1 >= cutoff * n2
            dup_pid << pid2
            #puts pid1 + "\t" + n1.to_s + "\t" + pid2 + "\t" + n2.to_s
          elsif n2 >= cutoff * n1
            dup_pid << pid1
            #puts pid2 + "\t" + n2.to_s + "\t" + pid1 + "\t" + n1.to_s
          end
        end
      end
    end


    new_seq = {}
    seq.each do |name,s|
      pid = name.split("_")[0][1..-1]
      unless dup_pid.include?(pid)
        new_seq[name] = s
      end
    end
    return new_seq
  end

  # collapse sequences with x number of nt differences. make sure sequences are aligned. The return frequency is NOT the frequency of the collasped sequences.
  def self.collapse_sequence_by_x_nt_difference(seq_array,cutoff)
      new_seq_freq = {}
      seq_freq = ViralSeq.count(seq_array)
      if seq_freq.size == 1
          new_seq_freq = seq_freq
      else
          uniq_seq = seq_freq.keys
          unique_seq_pair = uniq_seq.combination(2)
          dupli_seq = []
          unique_seq_pair.each do |pair|
              seq1 = pair[0]
              seq2 = pair[1]
              diff = ViralSeq.compare_two_seq(seq1,seq2)
              if diff <= cutoff
                  freq1 = seq_freq[seq1]
                  freq2 = seq_freq[seq2]
                  freq1 >= freq2 ? dupli_seq << seq2 : dupli_seq << seq1
              end
          end

          seq_freq.each do |seq,freq|
              unless dupli_seq.include?(seq)
                  new_seq_freq[seq] = freq
              end
          end
          return new_seq_freq
      end
  end

  # collapse sequence hash to unique sequence hash
  def self.uniq_sequence(seq = {}, sequence_name = "sequence")
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

  # compare two sequences, return the number of different positions, NO NEED alignment

  def self.compare_two_seq(seq1 = "", seq2 = "")
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
