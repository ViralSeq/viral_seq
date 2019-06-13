# viral_seq/tcs_core
# core functions for TCS and DR pipeline
# functions to manipulate sequences including:
#   ViralSeq::calculate_pid_cut_off
#   ViralSeq::consensus
#   ViralSeq::generate_primer_id_pool
#   ViralSeq::similar_pid?
#   ViralSeq::filter_similar_pid
#   ViralSeq::collapse_sequence_by_x_nt_difference
#   ViralSeq::compare_two_seq
#   ViralSeq::gap_strip
#   ViralSeq::gap_strip_ends
#   ViralSeq::paired_join1
#   ViralSeq::paired_join2

# ViralSeq.calculate_pid_cut_off(PID_abundance, estimated_error_rate)
#   # A function to calcuate cut-off for offspring primer IDs.
#   # see reference at Zhou et al. JVI 2016.
#   # https://www.ncbi.nlm.nih.gov/pubmed/26041299
#   # PID_abundance is the abundance of a certain PID
#   # estimated_error_rate is the estimated platform error rate, 0.02 (2%) as default
#   # the model supports error rate from 0.003 to 0.03.
#   # return an abundance cut-off (Integer) for offspring Primer IDs.

# ViralSeq.consensus(seq_array, majority_cutoff)
#   # Generate a consensus sequence from a given sequence array.
#   # where seq_array is an Array of input sequences (aligned) [seq1, seq2, seq3, ...]
#   # majority_cutoff is a Float of majority cut-off. default as simply majority (0.5)
# =USAGE
#   a_consensus_sequence = ViralSeq.cosensus(seq_array, majority_cutoff)

# ViralSeq.generate_primer_id_pool(n)
#   # generate all Primer ID combinations given the length of Primer ID
#   # n is the length of the Primer ID (Integer). default value of n is 8.
# =USAGE
#   primer_id_pool = ViralSeq.generate_primer_id_pool(10) # 10 is the length of Primer ID
#   puts primer_id_pool.size  #should be 4^10
#   => 1048576

# ViralSeq.similar_pid?(pid1, pid2, base_difference)
#   # compare two primer ID sequences.
#   # If they differ in certain bases, return boolean value "TURE",
#   # else, return boolean value "FALSE"
#   # where pid1 and pid2 are two Primer IDs for comparison
#   # base_difference is an Integer for difference bases that allowed
# =USAGE
#   # example
#   ViralSeq.similar_pid?("AAGGCTACGA", "AAGGATACGA", 1)
#   => true

# ViralSeq.filter_similar_pid(sequence_fasta_file, cut_off)
#   # compare PID with sequences which have identical sequences.
#   # PIDs differ by 1 base will be recognized.
#   # if PID1 is x time (cut-off) greater than PID2, PID2 will be disgarded
#   # where sequence_fasta_file is the sequence file in fasta format
#   # each sequence tag starting with ">" and the Primer ID sequence
#   # followed by the number of Primer ID appeared in the raw sequence
#   # the information sections in the tags are separated by underscore "_"
#   # example sequence tag: >AGGCGTAGA_32_sample1_RT
#   # cut_off is the fold cut-off to remove the potential residual offspring Primer IDs
#   # default value for cut_off is 10
#   # return a new sequence hash. {sequence_name => sequence, ...}

# ViralSeq.collapse_sequence_by_x_nt_difference(sequence_array, cutoff)
#   # ollapse sequences with x number of nt differences.
#   # input an Array object of sequences, make sure sequences are aligned.
#   # return a new Array object of collapsed sequences
#   # The return frequency is NOT the frequency of the collasped sequences.

# ViralSeq.compare_two_seq(seq1, seq2)
#   # compare two sequences as String object, return the number of differences as integer
#   # sequences will NOT align
#   # can use ViralSeq.muscle_align(seq1, seq2) to get the aligned sequences
# =USAGE
#   # example
#   seq1 = 'AAGGCGTAGGAC'
#   seq2 = 'AAGCTTAGGACG'
#   puts ViralSeq.compare_two_seq(seq1, seq2)
#   => 8
#   aligned_seqs = ViralSeq.muscle_align(seq1,seq2)
#   puts ViralSeq.compare_two_seq(aligned_seqs.values[0], aligned_seqs.values[1])
#   => 4

# ViralSeq.gap_strip(sequence_hash)
#   # strip positions with gaps in the sequence alignment as Hash object {:name => sequence, ...}
# =USAGE
#   # example
#   sequence_hash = {'>seq1' => 'AACCGGTT',
#                    '>seq2' => 'A-CCGGTT',
#                    '>seq3' => 'AAC-GGTT',
#                    '>seq4' => 'AACCG-TT',
#                    '>seq5' => 'AACCGGT-'}
#   ViralSeq.gap_strip(sequence_hash)
#   => {">seq1"=>"ACGT", ">seq2"=>"ACGT", ">seq3"=>"ACGT", ">seq4"=>"ACGT", ">seq5"=>"ACGT"}

# ViralSeq.gap_strip_ends(sequence_hash)
#   # similar to ViralSeq.gap_strip , but only strip the gaps at both ends of the alignment
# =USAGE
#   # example
#   sequence_hash = {'>seq1' => 'AACCGGTT',
#                    '>seq2' => 'A-CCGGTT',
#                    '>seq3' => 'AAC-GGTT',
#                    '>seq4' => 'AACCG-TT',
#                    '>seq5' => 'AACCGGT-'}
#   ViralSeq.gap_strip_ends(sequence_hash)
#   => {">seq1"=>"AACCGGT", ">seq2"=>"A-CCGGT", ">seq3"=>"AAC-GGT", ">seq4"=>"AACCG-T", ">seq5"=>"AACCGGT"}

# ViralSeq.paired_join1(sequence_pair_hash, overlap, difference_cut_off)
#   # pair-end join function for KNOW overlap size
#   # sequence_pair_hash is a Hash object for paired sequences {:seq_name => [:r1_seq, :r2_seq], ...}
#   # can use ViralSeq::pair_fasta_to_hash to load paired r1 and r2 sequences into paired sequence hash
#   # overlap is an integer that indicate how many bases are overlapped.
#   # overlap value at 0 means no overlap. R1 and R2 will be simply put together.
#   # difference_cut_off is a Float variable for the maximum mismatch rate allowed for the overlapping region
#   # default value for difference_cut_off is 0.0, i.e. no mis-match allowed
# =USAGE
#   # example
#   paired_seqs = {">pair1"=>["GGGGGGGGGGAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA",
#                      "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAATTTTTTTTTT"],
#           ">pair2"=>["GGGGGGGGGGAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA",
#                      "AAAAAAAAAAGAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAATTTTTTTTTT"],
#           ">pair3"=>["GGGGGGGGGGAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA",
#                      "AAAAAAAAAAGGAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAATTTTTTTTTT"]}
#   ViralSeq.paired_join1(paired_seqs, 100, 0.0).keys
#   => [">pair1"]
#   ViralSeq.paired_join1(paired_seqs, 100, 0.01).keys
#   => [">pair1", ">pair2"]
#   ViralSeq.paired_join1(paired_seqs, 100, 0.02)
#   => [">pair1", ">pair2", ">pair3"]

# ViralSeq.paired_join2(seq_pair_hash, model, diff)
#   # pair-end join function for UNKNOW overlap
#   # sequence_pair_hash is a Hash object for paired sequences {:seq_name => [:r1_seq, :r2_seq], ...}
#   # can use ViralSeq::pair_fasta_to_hash to load paired r1 and r2 sequences into paired sequence hash
#   # model has two options, 1 or 2 as Integer
#   # model 1: overlap is determined based on consensus, all sequence pairs are supposed to have the same overlap size
#   # model 2: overlap is determined for each sequence pair, sequence pairs can have different size of overlap
#   # minimal overlap by model 2 set to 4 positions
#   # if the sequence overlap may be smaller than 3 bases the model will consider as no overlap.
#   # difference_cut_off is a Float variable for the maximum mismatch rate allowed for the overlapping region
#   # default value for difference_cut_off is 0.0, i.e. no mis-match allowed
# =USAGE
#   # example 1
#   paired_seqs = {">pair1"=>["GGGGGGGGGGAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA",
#                      "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAATTTTTTTTTT"],
#           ">pair2"=>["GGGGGGGGGGAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA",
#                      "AAAAAAAAAAGAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAATTTTTTTTTT"],
#           ">pair3"=>["GGGGGGGGGGAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA",
#                      "AAAAAAAAAAGGAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAATTTTTTTTTT"]}
#   ViralSeq.paired_join2(paired_seqs, 1).keys
#   => [">pair1"]
#   ViralSeq.paired_join2(paired_seqs, 1, 0.01).keys
#   => [">pair1", ">pair2"]
#
#   # example 2
#   paired_seq2 = {">pair4" => ["AAAGGGGGGG", "GGGGGGGTT"],
#                  ">pair5" => ["AAAAAAGGGG", "GGGGTTTTT"],
#                  ">pair6" => ["AAACAAGGGG", "GGGGTTTTT"]
#                  }
#   ViralSeq.paired_join2(paired_seq2, 1)
#   => {">pair4"=>"AAAGGGGGGGGGGTT", ">pair5"=>"AAAAAAGGGGTTTTT", ">pair6"=>"AAACAAGGGGTTTTT"}
#   ViralSeq.paired_join2(paired_seq2, 2)
#   => {">pair4"=>"AAAGGGGGGGTT", ">pair5"=>"AAAAAAGGGGTTTTT", ">pair6"=>"AAACAAGGGGTTTTT"}


module ViralSeq

  # calculate cut-off for offspring primer IDs.
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

  # call consensus nucleotide, used by ViralSeq.consensus
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

  # compare two primer ID sequences.
  # If they differ in x base, return boolean value "TURE",
  # else, return boolean value "FALSE"
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

  # compare PID with sequences which have identical sequences.
  # PIDs differ by 1 base will be recognized.
  # if PID1 is x time greater than PID2, PID2 will be disgarded
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
          elsif n2 >= cutoff * n1
            dup_pid << pid1
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

  # collapse sequences with x number of nt differences. make sure sequences are aligned.
  # The return frequency is NOT the frequency of the collasped sequences.
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

  # gap strip from a sequence alignment

  def self.gap_strip(sequence_alignment)
    new_seq_hash = {}
    seq_size = sequence_alignment.values[0].size
    seq_matrix = {}
    (0..(seq_size - 1)).each do |p|
      seq_matrix[p] = []
      sequence_alignment.values.each do |s|
        seq_matrix[p] << s[p]
      end
    end

    seq_matrix.delete_if do |_p, list|
      list.include?("-")
    end

    sequence_alignment.each do |n,s|
      new_s = ""
      seq_matrix.keys.each {|p| new_s += s[p]}
      new_seq_hash[n] = new_s
    end
    return new_seq_hash
  end

  # gap strip from a sequence alignment, only strip the gaps at the ends of the alignment

  def self.gap_strip_ends(sequence_alignment)
    new_seq_hash = {}
    seq_size = sequence_alignment.values[0].size
    seq_matrix = {}
    (0..(seq_size - 1)).each do |p|
      seq_matrix[p] = []
      sequence_alignment.values.each do |s|
        seq_matrix[p] << s[p]
      end
    end
    n1 = 0
    n2 = 0
    seq_matrix.each do |_p, list|
      if list.include?("-")
        n1 += 1
      else
        break
      end
    end

    seq_matrix.keys.reverse.each do |p|
      list = seq_matrix[p]
      if list.include?("-")
        n2 += 1
      else
        break
      end
    end

    sequence_alignment.each do |n,s|
      new_s = s[n1..(- n2 - 1)]
      new_seq_hash[n] = new_s
    end
    return new_seq_hash
  end

  # input paired-end sequence hash format seq_name => [r1_seq, r2_seq]
  # overlap is pre-determined
  def self.paired_join1(seq_pair_hash, overlap, diff = 0.0)
    raise ArgumentError.new(":overlap has to be Integer, input #{overlap} invalid.") unless overlap.is_a? Integer
    raise ArgumentError.new(":diff has to be float or integer, input #{diff} invalid.") unless (diff.is_a? Integer or diff.is_a? Float)
    joined_seq_hash = {}
    seq_pair_hash.each do |seq_name, seq_pair|
      r1_seq = seq_pair[0]
      r2_seq = seq_pair[1]
      if overlap.zero?
        joined_seq_hash[seq_name] = r1_seq + r2_seq
      elsif ViralSeq.compare_two_seq(r1_seq[-overlap..-1], r2_seq[0,overlap]) <= (overlap * diff)
        joined_seq_hash[seq_name] = r1_seq + r2_seq[overlap..-1]
      else
        next
      end
    end
    return joined_seq_hash
  end


  # overlap is not predetermined
  # model 1: overlap is determined based on consensus, all sequence pairs are supposed to have the same overlap size
  # model 2: overlap is determined for each sequence pair, sequence pairs can have different size of overlap
  def self.paired_join2(seq_pair_hash, model = 1, diff = 0.0)
    begin
      raise ArgumentError.new(":diff has to be float or integer, input #{diff} invalid.") unless (diff.is_a? Integer or diff.is_a? Float)
      if model == 1
        overlap = ViralSeq.determine_overlap_pid_pair(seq_pair_hash, diff)
        return ViralSeq.paired_join1(seq_pair_hash, overlap, diff)
      elsif model == 2
        joined_seq_hash = {}
        seq_pair_hash.each do |seq_name, seq_pair|
          overlap_list = []
          ViralSeq.overlap_matrix(seq_pair[0], seq_pair[1]).each do |overlap1, diff_nt|
            cut_off_base = overlap1 * diff
            overlap_list << overlap1 if diff_nt <= cut_off_base
          end
          if overlap_list.empty?
            joined_seq_hash[seq_name] = seq_pair[0] + seq_pair[1]
          else
            overlap = overlap_list.max
            joined_seq_hash[seq_name] = seq_pair[0] + seq_pair[1][overlap..-1]
          end
        end
        return joined_seq_hash
      else
        raise ArgumentError.new("Error::Wrong Overlap Model Argument. Given \'#{model}\', expected '1' or '2'.")
      end
    rescue ArgumentError => e
      puts e
      return nil
    end
  end

  # determine overlap size from a paired sequence Hash object
  def self.determine_overlap_pid_pair(seq_pair_hash, diff = 0.0)
    overlaps = []
    seq_pair_hash.each do |_seq_name, seq_pair|
      overlap_list = []
      matrix = ViralSeq.overlap_matrix(seq_pair[0], seq_pair[1])
      matrix.each do |overlap, diff_nt|
        cut_off_base = overlap * diff
        overlap_list << overlap if diff_nt <= cut_off_base
      end
      if overlap_list.empty?
        overlaps << 0
      else
        overlaps << overlap_list.max
      end
    end
    count_overlaps = ViralSeq.count(overlaps)
    max_value = count_overlaps.values.max
    max_overlap_list = []
    count_overlaps.each {|overlap, counts| max_overlap_list << overlap if counts == max_value}
    max_overlap_list.max
  end

  # input a pair of sequences as String, return a Hash object of overlapping Hash object
  # {:overlap_size => number_of_differnt_positions, ...}
  # {minimal overlap set to 4. }
  def self.overlap_matrix(sequence1, sequence2)
    min_overlap = 4
    max_overlap = [sequence1.size, sequence2.size].max
    matrix_hash = {}
    (min_overlap..max_overlap).each do |overlap|
      matrix_hash[overlap] = ViralSeq.compare_two_seq(sequence1[-overlap..-1], sequence2[0, overlap])
    end
    return matrix_hash
  end

end
