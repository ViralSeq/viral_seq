
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
    # @param overlap [Integer] how many bases are overlapped. `0` means no overlap, R1 and R2 will be simply put together.
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
      seq_pair_hash = self.dna_hash
      raise ArgumentError.new(":overlap has to be Integer, input #{overlap} invalid.") unless overlap.is_a? Integer
      raise ArgumentError.new(":diff has to be float or integer, input #{diff} invalid.") unless (diff.is_a? Integer or diff.is_a? Float)
      joined_seq = {}
      seq_pair_hash.each do |seq_name,seq_pair|
        r1_seq = seq_pair[0]
        r2_seq = seq_pair[1]
        if overlap.zero?
          joined_sequence = r1_seq + r2_seq
        elsif diff.zero?
          if r1_seq[-overlap..-1] == r2_seq[0,overlap]
            joined_sequence= r1_seq + r2_seq[overlap..-1]
          end
        elsif r1_seq[-overlap..-1].compare_with(r2_seq[0,overlap]) <= (overlap * diff)
          joined_sequence= r1_seq + r2_seq[overlap..-1]
        else
          next
        end
        joined_seq[seq_name] = joined_sequence
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
            overlap_list = []
            overlap_matrix(seq_pair[0], seq_pair[1]).each do |overlap1, diff_nt|
              cut_off_base = overlap1 * diff
              overlap_list << overlap1 if diff_nt <= cut_off_base
            end
            if overlap_list.empty?
              joined_seq[seq_name] = seq_pair[0] + seq_pair[1]
            else
              overlap = overlap_list.max
              joined_seq[seq_name] = seq_pair[0] + seq_pair[1][overlap..-1]
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
    def determine_overlap_pid_pair(seq_pair_hash, diff = 0.0)
      overlaps = []
      seq_pair_hash.each do |_seq_name, seq_pair|
        overlap_list = []
        matrix = overlap_matrix(seq_pair[0], seq_pair[1])
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
      count_overlaps = overlaps.count_freq
      max_value = count_overlaps.values.max
      max_overlap_list = []
      count_overlaps.each {|overlap, counts| max_overlap_list << overlap if counts == max_value}
      max_overlap_list.max
    end # end pf determine_overlap_pid_pair

    # input a pair of sequences as String, return a Hash object of overlapping Hash object
    # {:overlap_size => number_of_differnt_positions, ...}
    # {minimal overlap set to 4. }
    def overlap_matrix(sequence1, sequence2)
      min_overlap = 4
      max_overlap = [sequence1.size, sequence2.size].min
      matrix_hash = {}
      (min_overlap..max_overlap).each do |overlap|
        matrix_hash[overlap] = sequence1[-overlap..-1].compare_with(sequence2[0, overlap])
      end
      return matrix_hash
    end # end of overlap_matrix

  end # end of SeqHashPair
end # end of ViralSeq
