
module ViralSeq

  # Class for paired-end sequences.
  # @example initialize a new SeqHashPair object from a directory containing paired-end sequences
  #   my_seqhashpair = ViralSeq::SeqHashPair.fa('my_seq_directory')
  # @example join the paired-end sequences with an overlap of 100 bp
  #   my_seqhashpair.join1(100)
  # @example join the paired-end sequences with unknown overlap, each pair of sequences has its own overlap size
  #   my_seqhashpair.join2(model: :indiv)

  class SeqHashPair

    # initialize SeqHashPair object with @dna_hash, @title and @file

    def initialize (dna_hash = {}, title = "", file = [])
      @dna_hash = dna_hash
      @title = title
      @file = file
    end

    # @return [Hash] Hash object for :name => [:r1_sequence_string, :r2_sequence_string]

    attr_accessor :dna_hash

    # @return [String] the title of the SeqHash object.
    # default as the directory basename if SeqHash object is initialized using ::fa

    attr_accessor :title

    # @return [String] the r1 and r2 files that are used to initialize SeqHash object, if they exist

    attr_accessor :file

    # initialize a new ViralSeq::SeqHashPair object from a directory containing paired sequence files in the FASTA format
    # @param indir [String] directory containing paired sequence files in the FASTA format,
    #
    #     Paired sequence files need to have "r1" and "r2" in their file names
    #
    #     Example for the file structure
    #       ├───lib1
    #           │     lib1_r1.txt
    #           │     lib1_r2.txt
    #     The sequence taxa should only differ by last 3 characters to distinguish r1 and r2 sequence.
    # @return [ViralSeq::SeqHashPair] new SeqHashPair object from the paired FASTA sequence files
    # @example initialize a new SeqHashPair object from a directory containing paired-end sequences
    #   my_seqhashpair = ViralSeq::SeqHashPair.fa('spec/sample_paired_seq')

    def self.new_from_fasta(indir)
      files = Dir[indir + "/*"]
      r1_file = ""
      r2_file = ""
      files.each do |f|
        if File.basename(f) =~ /r1/i
          r1_file = f
        elsif File.basename(f) =~ /r2/i
          r2_file = f
        end
      end

      seq1 = ViralSeq::SeqHash.fa(r1_file).dna_hash
      seq2 = ViralSeq::SeqHash.fa(r2_file).dna_hash

      new_seq1 = seq1.each_with_object({}) {|(k, v), h| h[k[0..-4]] = v}
      new_seq2 = seq2.each_with_object({}) {|(k, v), h| h[k[0..-4]] = v}

      seq_pair_hash = {}

      new_seq1.each do |seq_name,seq|
        seq_pair_hash[seq_name] = [seq, new_seq2[seq_name]]
      end
      seq_hash = ViralSeq::SeqHashPair.new
      seq_hash.dna_hash = seq_pair_hash
      seq_hash.title = File.basename(indir,".*")
      seq_hash.file = [r1_file, r2_file]
      return seq_hash
    end # end of .new_from_fasta

    class << self
      alias_method :fa, :new_from_fasta
    end

    # the size of nt sequence hash of the SeqHashPair object
    # @return [Integer] size of nt sequence hash of the SeqHash object
    def size
      self.dna_hash.size
    end

    # Pair-end join function for KNOWN overlap size.
    # @param overlap [Integer] simple overlap value indicating how many bases are overlapped. `0` means no overlap, R1 and R2 will be simply put together.
    # overlap can also be an explicit [Hash] object for :overlap_size, :r1_overlap, :r2_overlap, :before_overlap, :after_overlap
    # @param diff [Integer, Float] the maximum mismatch rate allowed for the overlapping region. default at 0.0, i.e. no mis-match allowed.
    # @return [ViralSeq::SeqHash] a SeqHash object of joined sequences.
    # @example join paired-end sequences with different :diff cut-offs, overlap provided.
    #   paired_seqs = {">pair1"=>["GGGGGGGGGGAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA",
    #                             "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAATTTTTTTTTT"],
    #                  ">pair2"=>["GGGGGGGGGGAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA",
    #                             "AAAAAAAAAAGAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAATTTTTTTTTT"],
    #                  ">pair3"=>["GGGGGGGGGGAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA",
    #                             "AAAAAAAAAAGGAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAATTTTTTTTTT"]}
    #   my_seqhashpair = ViralSeq::SeqHashPair.new(paired_seqs)
    #   my_seqhashpair.join1(100).dna_hash.keys
    #   => [">pair1"]
    #   my_seqhashpair.join1(100,0.01).dna_hash.keys
    #   => [">pair1", ">pair2"]
    #   my_seqhashpair.join1(100,0.02).dna_hash.keys
    #   => [">pair1", ">pair2", ">pair3"]

    def join1(overlap = 0, diff = 0.0)
      raise ArgumentError.new(":diff has to be float or integer, input #{diff} invalid.") unless (diff.is_a? Integer or diff.is_a? Float)

      if overlap.is_a? Integer and overlap.zero?
        overlap = {
          overlap_size: 0,
          r1_overlapped: 0...0,
          r2_overlapped: 0...0,
          before_overlap: {
            region: :r1,
            range: 0..-1,
          } ,
          after_overlap: {
            region: :r2,
            range: 0..-1
          }
        }
      elsif overlap.is_a? Integer
        overlap = {
          overlap_size: overlap,
          r1_overlapped: -overlap..-1,
          r2_overlapped: 0..(overlap - 1),
          before_overlap: {
            region: :r1,
            range: 0..(-overlap - 1),
          } ,
          after_overlap: {
            region: :r2,
            range: overlap..-1
          }
        }
      end

      seq_pair_hash = self.dna_hash
      joined_seq = {}
      seq_pair_hash.each do |seq_name,seq_pair|
        r1_seq = seq_pair[0]
        r2_seq = seq_pair[1]

        r1_overlap = r1_seq[overlap[:r1_overlapped]]
        r2_overlap = r2_seq[overlap[:r2_overlapped]]

        overlap_size = overlap[:overlap_size]

        if (diff.zero? and r1_overlap == r2_overlap) or (!diff.zero? and r1_overlap.compare_with(r2_overlap) <= (overlap_size.abs * diff))
          if overlap[:before_overlap][:region] == :r1
            before_overlap_seq = r1_seq[overlap[:before_overlap][:range]]
          elsif overlap[:before_overlap][:region] == :r2
            before_overlap_seq = r2_seq[overlap[:before_overlap][:range]]
          end

          if overlap[:after_overlap][:region] == :r1
            after_overlap_seq = r1_seq[overlap[:after_overlap][:range]]
          elsif overlap[:after_overlap][:region] == :r2
            after_overlap_seq = r2_seq[overlap[:after_overlap][:range]]
          end
          joined_sequence = before_overlap_seq + r1_overlap + after_overlap_seq
        end

        joined_seq[seq_name] = joined_sequence if joined_sequence
      end

      joined_seq_hash = ViralSeq::SeqHash.new
      joined_seq_hash.dna_hash = joined_seq
      joined_seq_hash.title = self.title + "_joined"
      joined_seq_hash.file = File.dirname(self.file[0]) if self.file.size > 0
      return joined_seq_hash
    end # end of join1


    # Pair-end join function for UNKNOWN overlap.
    # @param model [Symbol] models used to determine the overlap, `:con`, `:indiv`
    #
    #   model `:con`: overlap is determined based on consensus, all sequence pairs are supposed to have the same overlap size
    #
    #     note: minimal overlap as 4 bases.
    #   model `:indiv`: overlap is determined for each sequence pair, sequence pairs can have different size of overlap
    # @param diff (see #join1)
    # @return (see #join1)
    # @example join paired-end sequences, overlap NOT provided
    #   paired_seq2 = {">pair4" => ["AAAGGGGGGG", "GGGGGGGTT"],
    #                  ">pair5" => ["AAAAAAGGGG", "GGGGTTTTT"],
    #                  ">pair6" => ["AAACAAGGGG", "GGGGTTTTT"] }
    #   my_seqhashpair = ViralSeq::SeqHashPair.new(paired_seq2)
    #   my_seqhashpair.join2.dna_hash
    #   => {">pair4"=>"AAAGGGGGGGGGGTT", ">pair5"=>"AAAAAAGGGGTTTTT", ">pair6"=>"AAACAAGGGGTTTTT"}
    #   my_seqhashpair.join2(model: :indiv).dna_hash
    #   => {">pair4"=>"AAAGGGGGGGTT", ">pair5"=>"AAAAAAGGGGTTTTT", ">pair6"=>"AAACAAGGGGTTTTT"}

    def join2(model: :con, diff: 0.0)
      seq_pair_hash = self.dna_hash
      begin
        raise ArgumentError.new(":diff has to be float or integer, input #{diff} invalid.") unless (diff.is_a? Integer or diff.is_a? Float)
        if model == :con
          overlap = determine_overlap_pid_pair(seq_pair_hash, diff)
          return self.join1(overlap, diff)
        elsif model == :indiv
          joined_seq = {}
          seq_pair_hash.each do |seq_name, seq_pair|
            r1_seq = seq_pair[0]
            r2_seq = seq_pair[1]
            overlap_list = []

            overlap_matrix(r1_seq, r2_seq).each do |overlap1, diff_nt|
              cut_off_base = overlap1[:overlap_size] * diff
              overlap_list << overlap1 if diff_nt <= cut_off_base
            end

            if overlap_list.empty?
              joined_seq[seq_name]  = seq_pair[0] + seq_pair[1]
            else
              overlap_to_use = overlap_list.sort_by{|k| k[:overlap_size].abs}.reverse[0]

              if overlap_to_use[:before_overlap][:region] == :r1
                before_overlap_seq = r1_seq[overlap_to_use[:before_overlap][:range]]
              elsif overlap_to_use[:before_overlap][:region] == :r2
                before_overlap_seq = r2_seq[overlap_to_use[:before_overlap][:range]]
              end

              if overlap_to_use[:after_overlap][:region] == :r1
                after_overlap_seq = r1_seq[overlap_to_use[:after_overlap][:range]]
              elsif overlap_to_use[:after_overlap][:region] == :r2
                after_overlap_seq = r2_seq[overlap_to_use[:after_overlap][:range]]
              end
              joined_seq[seq_name] = before_overlap_seq + r1_seq[overlap_to_use[:r1_overlapped]] + after_overlap_seq
            end
          end

          joined_seq_hash = ViralSeq::SeqHash.new
          joined_seq_hash.dna_hash = joined_seq
          joined_seq_hash.title = self.title + "_joined"
          joined_seq_hash.file = File.dirname(self.file[0]) if self.file.size > 0
          return joined_seq_hash
        else
          raise ArgumentError.new("Error::Wrong Overlap Model Argument. Given \`#{model}\`, expected `:con` or `:indiv`.")
        end
      rescue ArgumentError => e
        puts e
        return nil
      end
    end # end of join2

    private
    # determine overlap size from @dna_hash
    def determine_overlap_pid_pair(seq_pair_hash, diff = 0.02)
      overlaps = []
      seq_pair_hash.each do |_seq_name, seq_pair|
        overlap_list = []
        matrix = overlap_matrix(seq_pair[0], seq_pair[1])
        matrix.each do |overlap_positions, diff_nt|
          overlap = overlap_positions[:overlap_size].abs
          cut_off_base = overlap * diff
          overlap_list << overlap_positions if diff_nt <= cut_off_base
        end

        if overlap_list.empty?
          overlaps <<    {
            overlap_size: 0,
            r1_overlapped: 0...0,
            r2_overlapped: 0...0,
            before_overlap: {
              region: :r1,
              range: 0..-1,
            } ,
            after_overlap: {
              region: :r2,
              range: 0..-1
            }
          }
        else
          overlaps << overlap_list.sort_by{|k| k[:overlap_size].abs}.reverse[0]
        end

      end
      count_overlaps = overlaps.count_freq
      max_value = count_overlaps.values.max
      max_overlap_list = []
      count_overlaps.each {|overlap, counts| max_overlap_list << overlap if counts == max_value}
      max_overlap_list.sort_by{|k| k[:overlap_size].abs}.reverse[0]
    end # end pf determine_overlap_pid_pair

    # input a pair of sequences as String, return a Hash object of overlapping Hash object
    # {:overlap_size => number_of_differnt_positions, ...}
    # {minimal overlap set to 4. }
    def overlap_matrix(sequence1, sequence2)
      list = overlap_list(sequence1.size, sequence2.size)
      matrix_hash = {}
      list.each do |l|
        range1 = l[:r1_overlapped]
        range2 = l[:r2_overlapped]
        matrix_hash[l] = sequence1[range1].compare_with(sequence2[range2])
      end
      matrix_hash
    end

    # given two [Integer], return all possible overlaping ranges in an [Array]
    def overlap_list(l1, l2)
      return_list = []
      min_overlap = 4
      max_overlap = [l1, l2].min
      diff = (l1 - l2).abs
      max_reverse = l1/2

      (min_overlap..max_overlap).each do |overlap|
        return_list<< {
          overlap_size: overlap,
          r1_overlapped: (l1-overlap)..(l1-1),
          r2_overlapped: 0..(overlap -1),
          before_overlap: {region: :r1, range: 0..(l1 - overlap - 1)},
          after_overlap: {region: :r2, range: overlap..(l2-1)}
        }
      end

      if l1 >= l2
        (1..diff).each do |overlap|
          return_list << {
            overlap_size: max_overlap,
            r1_overlapped: (diff - overlap)..(l1-1-overlap),
            r2_overlapped: 0..(l2-1),
            before_overlap: {region: :r1, range: 0...(diff - overlap)},
            after_overlap: {region: :r1, range: (l1-overlap)...l1},
        }
        end
      else
        (1..diff).each do |overlap|
          return_list << {
            overlap_size: max_overlap,
            r1_overlapped: 0..(l1-1),
            r2_overlapped: overlap..(max_overlap + overlap - 1),
            before_overlap: {region: :r2, range: 0...overlap},
            after_overlap: {region: :r2, range: (max_overlap + overlap)...l2},
        }
        end
      end

      (max_reverse..(max_overlap-1)).reverse_each do |overlap|
        return_list << {
          overlap_size: overlap,
          r1_overlapped: 0..(overlap -1),
          r2_overlapped: (l2-overlap)..(l2-1),
          before_overlap: {region: :r2, range: 0..(l2-overlap-1)},
          after_overlap: {region: :r1, range: overlap..(l1-1)},
        }
      end

      return_list
    end # end of overlap_list

  end # end of SeqHashPair
end # end of ViralSeq
