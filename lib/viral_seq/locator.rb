# viral_seq/locator.rb

# HIV sequence locator function
#   resembling HIV Sequence Locator from LANL
#   https://www.hiv.lanl.gov/content/sequence/LOCATE/locate.html
#   require MUSCLE (http://www.drive5.com/muscle) installed
#   current version only supports nucleotide sequence, not for amino acid sequence.
#
# =USAGE1
#   # Find the location of a sequence
#   ViralSeq.sequence_locator(input_sequence, reference_options, path_to_muscle)
#   # input_sequence: String of nucleotide sequence
#   # reference_options: choose a reference genome from :HXB2 (default), :NL43, or :MAC239
#   # path_to_muscle: path to the muscle executable. Default as 'muscle', add MUSCLE excutable path to $PATH to use the default setting or provide the path
#   # function returns an array of
#   #   start_location (Integer)
#   #   end_location (Integer)
#   #   percentage_of_similarity_to_reference_sequence (Float)
#   #   containing_indel? (Boolean)
#   #   aligned_input_sequence (String)
#   #   aligned_reference_sequence (String)
#   # example code
#   sequence = 'AGCAGATGATACAGTATTAGAAGAAATAAATTTGCCAGGAAGATGGAAACCAAAAATGATAGGGGGAATTGGAGGTTTTATCAAAGTAAGACAATATGATC'
#   p ViralSeq.sequence_locator(sequence, :NL43, 'muscle')
#   => [2333, 2433, 98.0, false, "AGCAGATGATACAGTATTAGAAGAAATAAATTTGCCAGGAAGATGGAAACCAAAAATGATAGGGGGAATTGGAGGTTTTATCAAAGTAAGACAATATGATC", "AGCAGATGATACAGTATTAGAAGAAATGAATTTGCCAGGAAGATGGAAACCAAAAATGATAGGGGGAATTGGAGGTTTTATCAAAGTAAGACAGTATGATC"]
#
# =USAGE2
#   # Given a pair of specific start and end positions, and an input sequence, return a sub-sequence of that range
#   # return nil if the input sequence is not in the range
#   ViralSeq.sequence_clip(input_sequence, start_position, end_position, reference_options, path_to_muscle)
#   # input_sequence: String of nucleotide sequence
#   # start_position and end_position: Integer of the start and end reference number of the sub-sequence
#   # reference_options and path_to_muscle are same as in ViralSeq.sequence_locator
#   # example code
#   seq = "CCTCAGATCACTCTTTGGCAACGACCCCTAGTTACAATAAGGGTAGGGGGGCAACTAAAGGAAGCCCTATTAGATACAGGAGCAGATGATACAGTATTAGAAGAAATAAATTTGCCAGGAAGATGGAAACCAAAAATGATAGGGGGAATTGGAGGTTTTATCAAAGTAAGACAATATGATCAGATACCCATAGAAATTTGTGGACATGAAGCTATAGGTACAGTATTAGTGGGACCTACACCTGTCAACATAATTGGGAGAAATCTGTTGACTCAGATTGGTTGCACTCTAAATTTT"
#   p ViralSeq.sequence_clip(seq, 2333, 2433, :HXB2, 'muscle')
#   => "AGCAGATGATACAGTATTAGAAGAAATAAATTTGCCAGGAAGATGGAAACCAAAAATGATAGGGGGAATTGGAGGTTTTATCAAAGTAAGACAATATGATC"


module ViralSeq

  def self.sequence_locator(seq='', ref_option = :HXB2, path_to_muscle = 'muscle')

    # ViralSeq.check_muscle(path_to_muscle)
    ori_ref = ViralSeq.locator_ref(ref_option)

    begin
      ori_ref_l = ori_ref.size
      l1 = 0
      l2 = 0

      aln_seq = ViralSeq.muscle_align(ori_ref, seq, path_to_muscle)
      aln_test = aln_seq[1]
      aln_test =~ /^(\-*)(\w.*\w)(\-*)$/
      gap_begin = $1.size
      gap_end = $3.size
      aln_test2 = $2
      ref = aln_seq[0]
      ref = ref[gap_begin..(-gap_end-1)]
      ref_size = ref.size
      if ref_size > 1.3*(seq.size)
        l1 = l1 + gap_begin
        l2 = l2 + gap_end
        max_seq = aln_test2.scan(/[ACGT]+/).max_by(&:length)
        aln_test2 =~ /#{max_seq}/
        before_aln_seq = $`
        before_aln = $`.size
        post_aln_seq = $'
        post_aln = $'.size
        before_aln_seq_size = before_aln_seq.scan(/[ACGT]+/).join("").size
        b1 = (1.3 * before_aln_seq_size).to_i
        post_aln_seq_size = post_aln_seq.scan(/[ACGT]+/).join("").size
        b2 = (1.3 * post_aln_seq_size).to_i
        if (before_aln > seq.size) and (post_aln <= seq.size)
          ref = ref[(before_aln - b1)..(ref_size - post_aln - 1)]
          l1 = l1 + (before_aln - b1)
        elsif (post_aln > seq.size) and (before_aln <= seq.size)
          ref = ref[before_aln..(ref_size - post_aln - 1 + b2)]
          l2 = l2 + post_aln - b2
        elsif (post_aln > seq.size) and (before_aln > seq.size)
          ref = ref[(before_aln - b1)..(ref_size - post_aln - 1 + b2)]
          l1 = l1 + (before_aln - b1)
          l2 = l2 + (post_aln - b2)
        end

        aln_seq = ViralSeq.muscle_align(ref, seq, path_to_muscle)
        aln_test = aln_seq[1]
        aln_test =~ /^(\-*)(\w.*\w)(\-*)$/
        gap_begin = $1.size
        gap_end = $3.size
        ref = aln_seq[0]
        ref = ref[gap_begin..(-gap_end-1)]
      end

      aln_test = aln_seq[1]
      aln_test =~ /^(\-*)(\w.*\w)(\-*)$/
      gap_begin = $1.size
      gap_end = $3.size
      aln_test = $2
      aln_test =~ /^(\w+)(\-*)\w/
      s1 = $1.size
      g1 = $2.size
      aln_test =~ /\w(\-*)(\w+)$/
      s2 = $2.size
      g2 = $1.size

      l1 = l1 + gap_begin
      l2 = l2 + gap_end
      repeat = 0

      if g1 == g2 and (s1 + g1 + s2) == ref.size
        if s1 > s2 and g2 > 2*s2
          ref = ref[0..(-g2-1)]
          repeat = 1
          l2 = l2 + g2
        elsif s1 < s2 and g1 > 2*s1
          ref = ref[g1..-1]
          repeat = 1
          l1 = l1 + g1
        end
      else
        if g1 > 2*s1
          ref = ref[g1..-1]
          repeat = 1
          l1 = l1 + g1
        end
        if g2 > 2*s2
          ref = ref[0..(-g2 - 1)]
          repeat = 1
          l2 = l2 + g2
        end
      end

      while repeat == 1
        aln_seq = ViralSeq.muscle_align(ref, seq, path_to_muscle)
        aln_test = aln_seq[1]
        aln_test =~ /^(\-*)(\w.*\w)(\-*)$/
        gap_begin = $1.size
        gap_end = $3.size
        aln_test = $2
        aln_test =~ /^(\w+)(\-*)\w/
        s1 = $1.size
        g1 = $2.size
        aln_test =~ /\w(\-*)(\w+)$/
        s2 = $2.size
        g2 = $1.size
        ref = aln_seq[0]
        ref = ref[gap_begin..(-gap_end-1)]
        l1 = l1 + gap_begin
        l2 = l2 + gap_end
        repeat = 0
        if g1 > 2*s1
          ref = ref[g1..-1]
          repeat = 1
          l1 = l1 + g1
        end
        if g2 > 2*s2
          ref = ref[0..(-g2 - 1)]
          repeat = 1
          l2 = l2 + g2
        end
      end
      ref = ori_ref[l1..(ori_ref_l - l2 - 1)]


      aln_seq = ViralSeq.muscle_align(ref, seq, path_to_muscle)
      aln_test = aln_seq[1]
      ref = aln_seq[0]

      #refine alignment

      if ref =~ /^(\-+)/
        l1 = l1 - $1.size
      elsif ref =~ /(\-+)$/
        l2 = l2 + $1.size
      end

      if (ori_ref_l - l2 - 1) >= l1
        ref = ori_ref[l1..(ori_ref_l - l2 - 1)]
        aln_seq = ViralSeq.muscle_align(ref, seq, path_to_muscle)
        aln_test = aln_seq[1]
        ref = aln_seq[0]

        ref_size = ref.size
        sim_count = 0
        (0..(ref_size-1)).each do |n|
          ref_base = ref[n]
          test_base = aln_test[n]
          sim_count += 1 if ref_base == test_base
        end
        similarity = (sim_count/ref_size.to_f*100).round(1)

        loc_p1 = l1 + 1
        loc_p2 = ori_ref_l - l2
        if seq.size != (loc_p2 - loc_p1 + 1)
            indel = true
        elsif aln_test.include?("-")
            indel = true
        else
            indel = false
        end
        return [loc_p1,loc_p2,similarity,indel,aln_test,ref]
      else
        return [0,0,0,0,0,0,0]
      end
    rescue => e
      puts "Unexpected error occured."
      puts "Exception Class: #{ e.class.name }"
      puts "Exception Message: #{ e.message }"
      puts "Exception Backtrace: #{ e.backtrace[0] }"
      puts "ViralSeq.sequence_locator returns nil"
      return nil
    end
  end

  # sequence clip function
  def self.sequence_clip(seq='', p1 = 0, p2 = 0, ref_option = :HXB2, path_to_muscle = 'muscle')
    loc = ViralSeq.sequence_locator(seq, ref_option, path_to_muscle)
    l1 = loc[0]
    l2 = loc[1]
    if (p1 >= l1) & (p2 <= l2)
        seq = loc[4]
        ref = loc[5]
        g1 = 0
        ref.each_char do |char|
            break if l1 == p1
            g1 += 1
            l1 += 1 unless char == "-"
        end
        g2 = 1
        ref.reverse.each_char do |char|
            break if l2 == p2
            g2 += 1
            l2 -= 1 unless char == "-"
        end
        return seq[g1..(-g2)].tr("-","")
    else
        return nil
    end
  end

end
