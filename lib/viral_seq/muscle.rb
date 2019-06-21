# viral_seq/muscle.rb
# wrapper for MUSCLE (http://www.drive5.com/muscle)
# Including Methods as:
#   ViralSeq::check_muscle
#   ViralSeq::muscle_align
#   ViralSeq::muscle_align_multi

# ViralSeq.check_muscle?(path_to_muscle)
#   # check if the path_to_muscle provided is valid,
#   # prompt error messages if MUSCLE is not found.

# ViralSeq.muscle_align(reference_seq, test_sequence, path_to_muscle)
#   # takes a reference sequence and a test sequence as String object
#   # without specification on path_to_muscle, MuscleBio will be called to run Muscle
#   # specify path_to_muscle if other source of muscle needed
#   # returns aligned reference sequence and test sequences

# ViralSeq.muscle_align_multi(sequence_hash, path_to_muscle)
#   # input a sequence_hash object {:name=>:sequence,...}
#   # without specification on path_to_muscle, MuscleBio will be called to run Muscle
#   # specify path_to_muscle if other source of muscle needed
#   # return aligned sequences an hash

module ViralSeq

  # check if path_to_muscle is correct
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
  end

  def self.muscle_align(ref_seq = "", test_seq = "", path_to_muscle = false)
    temp_dir=File.dirname($0)
    temp_file = temp_dir + "/_temp_muscle_in"
    temp_aln = temp_dir + "/_temp_muscle_aln"
    name = ">test"
    temp_in = File.open(temp_file,"w")
    temp_in.puts ">ref"
    temp_in.puts ref_seq
    temp_in.puts name
    temp_in.puts test_seq
    temp_in.close
    if path_to_muscle
      unless ViralSeq.check_muscle?(path_to_muscle)
        File.unlink(temp_file)
        return nil;
      end
      print `#{path_to_muscle} -in #{temp_file} -out #{temp_aln} -quiet`
    else
      MuscleBio.run("muscle -in #{temp_file} -out #{temp_aln} -quiet")
    end
    aln_seq_hash = ViralSeq.fasta_to_hash(temp_aln)
    File.unlink(temp_file)
    File.unlink(temp_aln)
    return [aln_seq_hash[">ref"], aln_seq_hash[">test"]]
  end

  def self.muscle_align_multi(seq_hash = {}, path_to_muscle = false)
    temp_dir=File.dirname($0)
    temp_file = temp_dir + "/_temp_muscle_in"
    temp_aln = temp_dir + "/_temp_muscle_aln"
    File.open(temp_file, 'w'){|f| seq_hash.each {|k,v| f.puts k; f.puts v}}
    if path_to_muscle
      unless ViralSeq.check_muscle?(path_to_muscle)
        File.unlink(temp_file)
        return nil
      end
      print `#{path_to_muscle} -in #{temp_file} -out #{temp_aln} -quiet`
    else
      MuscleBio.run("muscle -in #{temp_file} -out #{temp_aln} -quiet")
    end
    out_seq_hash = ViralSeq.fasta_to_hash(temp_aln)
    File.unlink(temp_file)
    File.unlink(temp_aln)
    return out_seq_hash
  end
end
