
module ViralSeq
  # alignment using MUSCLE alignment program
  # @see http://www.drive5.com/muscle MUSCLE download link

  module Muscle
    # check if path_to_muscle is correct, prompt error messages if MUSCLE is not found.
    # @param path_to_muscle [String] path to muscle excutable
    # @return [boolean]

    def self.check_muscle?(path_to_muscle)
      begin
        `#{path_to_muscle} -version`
        return true
      rescue Errno::ENOENT
        puts "
              Error: MUSCLE is not found for at the provided {path_to_muscle}!!
              MUSLCE can be download at http://www.drive5.com/muscle
              Add MUSCLE excutable path to $PATH using
              $  export PATH=$PATH:/path/to/muscle
              or
              provide path_to_MUSCLE in the function arguments\n
              "
        return false
      end
    end # end of .check_muscle?

    # align a sequence with reference sequence Strings
    # @param ref_seq [String] reference sequence
    # @param test_seq [String] test sequence
    # @param path_to_muscle [String], path to MUSCLE excutable. if not provided (as default), it will use RubyGem::MuscleBio
    # @return [Array] a pair of [:ref_seq_aligned, :test_seq_aligned] or nil
    #   if the cannot find MUSCLE excutable
    # @example
    #   seq1 = 'AAGGCGTAGGAC'
    #   seq2 = 'AAGCTTAGGACG'
    #   aligned_seqs = ViralSeq::Muscle.align(seq1,seq2)
    #   => ["AAGGCGTAGGAC-", "-AAGCTTAGGACG"]

    def self.align(ref_seq = "", test_seq = "", path_to_muscle = false)
      temp_dir = Dir.home
      temp_name = "_"  + SecureRandom.alphanumeric
      temp_file = File.join(temp_dir, temp_name)
      temp_aln = File.join(temp_dir, (temp_name + "_aln"))
      name = ">test"
      temp_in = File.open(temp_file,"w")
      temp_in.puts ">ref"
      temp_in.puts ref_seq
      temp_in.puts name
      temp_in.puts test_seq
      temp_in.close
      if path_to_muscle
        unless ViralSeq::Muscle.check_muscle?(path_to_muscle)
          File.unlink(temp_file)
          return nil;
        end
        print `#{path_to_muscle} -in #{temp_file} -out #{temp_aln} -quiet`
      else
        MuscleBio.run("muscle -in #{temp_file} -out #{temp_aln} -quiet")
      end
      aln_seq_hash = ViralSeq::SeqHash.fa(temp_aln).dna_hash
      File.unlink(temp_file)
      File.unlink(temp_aln)
      return [aln_seq_hash[">ref"], aln_seq_hash[">test"]]
    end # end of .align
  end # end of ViralSeq::Muscle

end # end of ViralSeq
