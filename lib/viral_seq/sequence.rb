
module ViralSeq

  # ViralSeq::Sequence class for sequence operation
  #
  # @example create a sequence object
  #   seq = ViralSeq::Sequence.new('my_sequence', 'ACCTAGGTTCGGAGC')
  #   => #<ViralSeq::Sequence:0x00007fd03c8c10b8 @name="my_sequence", @dna="ACCTAGGTTCGGAGC", @aa_string="", @aa_array=[]>
  #
  # @example return dna sequence as String
  #   seq.dna
  #   => "ACCTAGGTTCGGAGC"
  #
  # @example reverse complement sequence of DNA sequence
  #   seq.rc
  #   => "GCTCCGAACCTAGGT"
  #
  # @example change @dna to reverse complement DNA sequence
  #   seq.rc!
  #
  # @example translate the DNA sequence, return values for @aa_string and @aa_array
  #   seq = ViralSeq::Sequence.new('my_sequence', 'AWTCGRAGAG')
  #   seq.translate(1)
  #   seq.aa_string
  #   => "##E"
  #   seq.aa_array
  #   => ["IF", "EG", "E"]

  class Sequence
    # initialize a ViralSeq::Sequence class with sequence name (default as '>sequence')
    # and DNA sequence as String object
    def initialize (name = ">sequence",dna_sequence ="")
      @name = name
      @dna = dna_sequence.upcase
      @aa_string = ""
      @aa_array = []
    end

    # @return [String] sequence tag name
    attr_accessor :name

    # @return [String] DNA sequence
    attr_accessor :dna

    # @return [String] amino acid sequence
    attr_accessor :aa_string

    # @return [Array] amino acid sequence as an Array object,
    # ambiguity dna sequence will be translated in all possible amino acid sequence at the position
    attr_accessor :aa_array

    # @return [String] reverse compliment sequence of the @dna.
    def rev_complement
      @dna.rc
    end

    # replace the @dna with reverse complement DNA sequence.
    def rev_complement!
      @dna = @dna.rc
    end

    alias_method :rc, :rev_complement
    alias_method :rc!, :rev_complement!

    # translate @dna to amino acid sequence.
    # generate values for @aa_string and @aa_array
    # @param initial_position [Integer] option `0`, `1` or `2`, indicating 1st, 2nd, 3rd reading frames

    def translate(initial_position = 0)
      @aa_string = ""
      require_sequence = @dna[initial_position..-1]
      base_array = []
      require_sequence.each_char {|base| base_array << base}
      while (base_array.length>=3) do
        base_3= ""
        3.times {base_3 += base_array.shift}
        @aa_string << amino_acid(base_3)
      end

      @aa_array = []
      require_sequence = @dna[initial_position..-1].tr('-','N')
      base_array = []
      require_sequence.each_char {|base| base_array << base}
      while (base_array.length>=3) do
        base_3= ""
        3.times{base_3 += base_array.shift}
        @aa_array<< amino_acid_2(base_3)
      end
    end

    # @return [Integer] length of DNA sequence
    def dna_length
      @dna.length
    end

    # @return [Integer] length of amino acid sequence
    def aa_length
      @aa_string.length
    end

    # resistant mutation interpretation for a chosen region from a translated ViralSeq::Sequence object
    # @param option [Symbol] option of region to interpret, `:hcv_ns5a`, `:hiv_pr`, `:nrti`, `:nnrti`, `hiv_in`
    # @param start_aa [Integer] the starting aa number of the input sequence
    # @return [Hash] return a Hash object for SDRMs identified. :posiiton => [:wildtype_codon, :mutation_codon]
    # @example examine an HIV PR region sequence for drug resistance mutations
    #   my_seq_name = 'a_pr_seq'
    #   my_seq = 'CCTCAGATCACTCTTTGGCAACGACCCCTCGTCACAGTAAAAATAGGAGGGCAATTAAAGGAAGCTCTATTAGATACAGGAGCAGATAATACAGTATTAGAAGACATGGAGTTACCAGGAAGATGGAAACCAAAAATGATAGGGGGAATTGGAGGTTTTATCAAAGTAAGACAATATGATCAGATACCCATAGAAATCTGTGGGCATAAAACTACAGGTACAGTGTTAATAGGACCTACACCCGTCAACATAATTGGAAGAGATCTGTTGACTCAGCTTGGTTGCACTTTAAATTTT'
    #   s = ViralSeq::Sequence.new(my_seq_name, my_seq)
    #   s.translate
    #   s.sdrm(:hiv_pr)
    #   => {30=>["D", "N"], 88=>["N", "D"]}

    def sdrm(option, start_aa = 1)
      aa_array = self.aa_array
      out_hash = {}
      sdrm = sdrm_hash(option)
      aa_length = aa_array.size
      end_aa = start_aa + aa_length - 1
      (start_aa..end_aa).each do |position|
        array_position = position - start_aa
        if sdrm.keys.include?(position)
          wt_aa = sdrm[position][0]
          test_aa = aa_array[array_position]
          if test_aa.size == 1
            unless wt_aa == test_aa
              if sdrm[position][1].include?(test_aa)
                out_hash[position] = [wt_aa,test_aa]
              end
            end
          else
            test_aa_array = test_aa.split("")
            if (test_aa_array & sdrm[position][1])
              out_hash[position] = [wt_aa,test_aa]
            end
          end
        end
      end
      return out_hash
    end # end of #hcv_ns5a

    # HIV sequence locator function, resembling HIV Sequence Locator from LANL
    #   # current version only supports nucleotide sequence, not for amino acid sequence.
    # @param ref_option [Symbol], name of reference genomes, options are `:HXB2`, `:NL43`, `:MAC239`
    # @param path_to_muscle [String], path to the muscle executable, if not provided, use MuscleBio to run Muscle
    # @return [Array] an array of the following info:
    #
    #   start_location (Integer)
    #
    #   end_location (Integer)
    #
    #   percentage_of_similarity_to_reference_sequence (Float)
    #
    #   containing_indel? (Boolean)
    #
    #   aligned_input_sequence (String)
    #
    #   aligned_reference_sequence (String)
    #
    # @example identify the location of the input sequence on the NL43 genome
    #   sequence = 'AGCAGATGATACAGTATTAGAAGAAATAAATTTGCCAGGAAGATGGAAACCAAAAATGATAGGGGGAATTGGAGGTTTTATCAAAGTAAGACAATATGATC'
    #   s = ViralSeq::Sequence.new('my_sequence', sequence)
    #   loc = s.locator(:NL43)
    #   h = ViralSeq::SeqHash.new; h.dna_hash['NL43'] = loc[5]; h.dna_hash[s.name] = loc[4]
    #   rs_string = h.to_rsphylip.split("\n")[1..-1].join("\n") # get a relaxed phylip format string for display of alignment.
    #   puts "The input sequence \"#{s.name}\" is located on the NL43 nt sequence from #{loc[0].to_s} to #{loc[1].to_s}.\nIt is #{loc[2].to_s}% similar to the reference.\nIt #{loc[3]? "does" : "does not"} have indels.\nThe alignment is\n#{rs_string}"
    #   => The input sequence "my_sequence" is located on the NL43 nt sequence from 2333 to 2433.
    #   => It is 98.0% similar to the reference.
    #   => It does not have indels.
    #   => The alignment is
    #   => NL43         AGCAGATGAT ACAGTATTAG AAGAAATGAA TTTGCCAGGA AGATGGAAAC CAAAAATGAT AGGGGGAATT GGAGGTTTTA TCAAAGTAAG ACAGTATGAT C
    #   => my_sequence  AGCAGATGAT ACAGTATTAG AAGAAATAAA TTTGCCAGGA AGATGGAAAC CAAAAATGAT AGGGGGAATT GGAGGTTTTA TCAAAGTAAG ACAATATGAT C
    # @see https://www.hiv.lanl.gov/content/sequence/LOCATE/locate.html LANL Sequence Locator

    def locator(ref_option = :HXB2, path_to_muscle = false)
      seq = self.dna
      ori_ref = ViralSeq::RefSeq.get(ref_option)

      begin
        ori_ref_l = ori_ref.size
        l1 = 0
        l2 = 0

        aln_seq = ViralSeq::Muscle.align(ori_ref, seq, path_to_muscle)
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

          aln_seq = ViralSeq::Muscle.align(ref, seq, path_to_muscle)
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
          aln_seq = ViralSeq::Muscle.align(ref, seq, path_to_muscle)
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


        aln_seq = ViralSeq::Muscle.align(ref, seq, path_to_muscle)
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
          aln_seq = ViralSeq::Muscle.align(ref, seq, path_to_muscle)
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
    end # end of locator

    # Given start and end positions on the reference genome, return a sub-sequence of the target sequence in that range
    # @param p1 [Integer] start position number on the reference genome
    # @param p2 [Integer] end position number on the reference genome
    # @param ref_option [Symbol], name of reference genomes, options are `:HXB2`, `:NL43`, `:MAC239`
    # @param path_to_muscle [String], path to the muscle executable, if not provided, use MuscleBio to run Muscle
    # @return [ViralSeq::Sequence, nil] a new ViralSeq::Sequence object that of input range on the reference genome or nil
    #   if either the start or end position is beyond the range of the target sequence.
    # @example trim a sequence to fit in the range of [2333, 2433] on the HXB2 nt reference
    #   seq = "CCTCAGATCACTCTTTGGCAACGACCCCTAGTTACAATAAGGGTAGGGGGGCAACTAAAGGAAGCCCTATTAGATACAGGAGCAGATGATACAGTATTAGAAGAAATAAATTTGCCAGGAAGATGGAAACCAAAAATGATAGGGGGAATTGGAGGTTTTATCAAAGTAAGACAATATGATCAGATACCCATAGAAATTTGTGGACATGAAGCTATAGGTACAGTATTAGTGGGACCTACACCTGTCAACATAATTGGGAGAAATCTGTTGACTCAGATTGGTTGCACTCTAAATTTT"
    #   s = ViralSeq::Sequence.new('my_seq', seq)
    #   s.sequence_clip(2333, 2433, :HXB2).dna
    #   => "AGCAGATGATACAGTATTAGAAGAAATAAATTTGCCAGGAAGATGGAAACCAAAAATGATAGGGGGAATTGGAGGTTTTATCAAAGTAAGACAATATGATC"

    def sequence_clip(p1 = 0, p2 = 0, ref_option = :HXB2, path_to_muscle = false)
      loc = self.locator(ref_option, path_to_muscle)
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
          return ViralSeq::Sequence.new(self.name,seq[g1..(-g2)].tr("-",""))
      else
          return nil
      end
    end

    # start of private functions
    private

    # generate amino acid abbreviations from 3 bases, ambiguity will return "#"
    def amino_acid (bases)
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
    end # end of amino_acid

    # keep ambiguities, return all possible amino acids.

    def amino_acid_2 (bases)
      bases_to_aa = []
      aa_list = []
      base1 = bases[0].to_list
      base2 = bases[1].to_list
      base3 = bases[2].to_list
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

      bases_to_aa.each do |base|
      case base
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
      aa_out = aa_list.uniq.join
      return aa_out
    end # end of #amino_acid_2

    # sdrm position hash
    def sdrm_hash(options)
      sdrm = {}
      case options
      when :hcv_ns5a
        sdrm[28] = ['M',['T']]
        sdrm[30] = ['L',['H','K','R','Q','A','S','D']]
        sdrm[31] = ['L',['M','V','F']]
        sdrm[32] = ['P',['L']]
        sdrm[44] = ['K',['R']]
        sdrm[58] = ['H',['D','P','S']]
        sdrm[64] = ['T',['A','S']]
        sdrm[77] = ['P',['A','S']]
        sdrm[78] = ['R',['K']]
        sdrm[79] = ['T',['A']]
        sdrm[83] = ['T',['M']]
        sdrm[85] = ['S',['N','H','Y']]
        sdrm[92] = ['A',['P','T','K','E']]
        sdrm[93] = ['Y',['C','F','H','N']]
        sdrm[107] = ['K',['T','S']]
        sdrm[121] = ['I',['V']]
        sdrm[135] = ['T',['A']]
      when :nrti
        sdrm[41] = ['M',['L']]
        sdrm[65] = ['K',['R']]
        sdrm[67] = ['D',['N','G','E']]
        sdrm[69] = ['T',['D']]
        sdrm[70] = ['K',['R','E']]
        sdrm[74] = ['L',['V','I']]
        sdrm[75] = ['V',['M','T','A','S']]
        sdrm[77] = ['F',['L']]
        sdrm[115] = ['Y',['F']]
        sdrm[116] = ['F',['Y']]
        sdrm[151] = ['Q',['M']]
        sdrm[184] = ['M',['V','I']]
        sdrm[210] = ['L',['W']]
        sdrm[215] = ["T",["Y","F","I","C","D","V","E"]]
        sdrm[219] = ["K",["Q","E","N","R"]]
      when :nnrti
        sdrm[100] = ['L',['I']]
        sdrm[101] = ['K',['E','P']]
        sdrm[103] = ['K',['N','S']]
        sdrm[106] = ['V',['M','A']]
        sdrm[179] = ['V',['F','D']]
        sdrm[181] = ['Y',['C','I','V']]
        sdrm[188] = ['Y',['L','H','C']]
        sdrm[190] = ['G',['A','S','E']]
        sdrm[225] = ['P',['H']]
        sdrm[230] = ['M',['L']]
      when :hiv_pr
        sdrm[23] = ['L',['I']]
        sdrm[24] = ['L',['I']]
        sdrm[30] = ['D',['N']]
        sdrm[32] = ['V',['I']]
        sdrm[46] = ['M',['I','L']]
        sdrm[47] = ['I',['V','A']]
        sdrm[48] = ['G',['V','M']]
        sdrm[50] = ['I',['V','L']]
        sdrm[53] = ['F',['L']]
        sdrm[54] = ['I',['V','L','M','T','A','S']]
        sdrm[73] = ['G',['S','T','C','A']]
        sdrm[76] = ['L',['V']]
        sdrm[82] = ['V',['A','T','S','F','L','C','M']]
        sdrm[83] = ['N',['D']]
        sdrm[84] = ['I',['V','A','C']]
        sdrm[88] = ['N',['D','S']]
        sdrm[90] = ['L',['M']]
      when :hiv_in
        sdrm[66] = ['T',['A','I','K']]
        sdrm[74] = ['L',['M']]
        sdrm[92] = ['E',['Q']]
        sdrm[95] = ['Q',['K']]
        sdrm[97] = ['T',['A']]
        sdrm[121] = ['F',['Y']]
        sdrm[140] = ['G',['A','S','C']]
        sdrm[143] = ["Y",["C","H","R"]]
        sdrm[147] = ['S',['G']]
        sdrm[148] = ['Q',['H','K','R']]
        sdrm[155] = ['N',['S','H']]
      else raise "Input option `#{options}` for ViralSeq::Sequence.sdrm not supported"
      end
      return sdrm
    end
  end # end of ViralSeq::Sequence
end # end of ViralSeq
