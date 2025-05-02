
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
      sdrm = ViralSeq::DRMs.sdrm_hash(option)
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
    end # end of #sdrm

    # Similar to #sdrm but use a DRM list as a param

    def check_drm(drm_list_single_type)
      aa_array = self.aa_array
      out_hash = {}

      drm_list_single_type.each do |position, mut|
        wt_aa = mut[0]
        mut_aas = mut[1]
        test_aa = aa_array[position - 1]
        if test_aa.size == 1 and mut_aas.include?(test_aa)
          out_hash[position] = [wt_aa, test_aa]
        elsif test_aa.size > 1
          test_aa_array = test_aa.split("")
          mut_detected = test_aa_array & mut_aas

          if !mut_detected.empty?
            out_hash[position] = [wt_aa, mut_detected.join]
          end

        end
      end
      return out_hash
    end

    # HIV sequence locator function, resembling HIV Sequence Locator from LANL
    #   # current version only supports nucleotide sequence, not for amino acid sequence.
    # @param ref_option [Symbol], name of reference genomes, options are `:HXB2`, `:SIVmm239`
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
    # @example identify the location of the input sequence on the HXB2 genome
    #   sequence = 'AGCAGATGATACAGTATTAGAAGAAATAAATTTGCCAGGAAGATGGAAACCAAAAATGATAGGGGGAATTGGAGGTTTTATCAAAGTAAGACAATATGATC'
    #   s = ViralSeq::Sequence.new('my_sequence', sequence)
    #   loc = s.locator(:HXB2)
    #   h = ViralSeq::SeqHash.new; h.dna_hash['HXB2'] = loc[5]; h.dna_hash[s.name] = loc[4]
    #   rs_string = h.to_rsphylip.split("\n")[1..-1].join("\n") # get a relaxed phylip format string for display of alignment.
    #   puts "The input sequence \"#{s.name}\" is located on the HXB2 nt sequence from #{loc[0].to_s} to #{loc[1].to_s}.\nIt is #{loc[2].round(1).to_s}% similar to the reference.\nIt #{loc[3]? "does" : "does not"} have indels.\nThe alignment is\n#{rs_string}"
    #   => The input sequence "my_sequence" is located on the HXB2 nt sequence from 2333 to 2433.
    #   => It is 97.0% similar to the reference.
    #   => It does not have indels.
    #   => The alignment is
    #   => HXB2         AGCAGATGAT ACAGTATTAG AAGAAATGAA TTTGCCAGGA AGATGGAAAC CAAAAATGAT AGGGGGAATT GGAGGTTTTA TCAAAGTAAG ACAGTATGAT C
    #   => my_sequence  AGCAGATGAT ACAGTATTAG AAGAAATAAA TTTGCCAGGA AGATGGAAAC CAAAAATGAT AGGGGGAATT GGAGGTTTTA TCAAAGTAAG ACAATATGAT C
    # @see https://www.hiv.lanl.gov/content/sequence/LOCATE/locate.html LANL Sequence Locator
    def locator(ref_option = :HXB2, algorithm = 1)
      seq = self.dna
      ref = ref_option.to_s
      begin
        loc = VirustLocator::Locator.exec(seq, "nt", algorithm, ref).split("\t")
        loc[0] = loc[0].to_i
        loc[1] = loc[1].to_i
        loc[2] = loc[2].to_f.round(1)
        if loc[3].to_s.downcase == "true"
          loc[3] = true
        else
          loc[3] = false
        end
      rescue => e
        puts "Unexpected error occured."
        puts "Exception Class: #{ e.class.name }"
        puts "Exception Message: #{ e.message }"
        puts "Exception Backtrace: #{ e.backtrace[0] }"
        puts "ViralSeq.sequence_locator returns nil"
        return nil
      end
      return loc
    end #end of locator

    # Given start and end positions on the reference genome, return a sub-sequence of the target sequence in that range
    # @param p1 [Integer] start position number on the reference genome
    # @param p2 [Integer] end position number on the reference genome
    # @param ref_option [Symbol], name of reference genomes, options are `:HXB2`, `:SIVmm239`
    # @param path_to_muscle [String], path to the muscle executable, if not provided, use MuscleBio to run Muscle
    # @return [ViralSeq::Sequence, nil] a new ViralSeq::Sequence object that of input range on the reference genome or nil
    #   if either the start or end position is beyond the range of the target sequence.
    # @example trim a sequence to fit in the range of [2333, 2433] on the HXB2 nt reference
    #   seq = "CCTCAGATCACTCTTTGGCAACGACCCCTAGTTACAATAAGGGTAGGGGGGCAACTAAAGGAAGCCCTATTAGATACAGGAGCAGATGATACAGTATTAGAAGAAATAAATTTGCCAGGAAGATGGAAACCAAAAATGATAGGGGGAATTGGAGGTTTTATCAAAGTAAGACAATATGATCAGATACCCATAGAAATTTGTGGACATGAAGCTATAGGTACAGTATTAGTGGGACCTACACCTGTCAACATAATTGGGAGAAATCTGTTGACTCAGATTGGTTGCACTCTAAATTTT"
    #   s = ViralSeq::Sequence.new('my_seq', seq)
    #   s.sequence_clip(2333, 2433, :HXB2).dna
    #   => "AGCAGATGATACAGTATTAGAAGAAATAAATTTGCCAGGAAGATGGAAACCAAAAATGATAGGGGGAATTGGAGGTTTTATCAAAGTAAGACAATATGATC"

    def sequence_clip(p1 = 0, p2 = 0, ref_option = :HXB2)
      loc = self.locator(ref_option)
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

  end # end of ViralSeq::Sequence
end # end of ViralSeq
