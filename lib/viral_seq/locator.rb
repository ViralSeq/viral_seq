# viral_seq/locator.rb


module ViralSeq

  def self.sequence_locator(seq='', ref_option = 'HXB2', path_to_muscle = 'muscle')

    # check if path_to_muscle is correct
    begin
      `#{path_to_muscle} -version`
    rescue Errno::ENOENT
      puts "
            Error: MUSCLE is not found!!\n
            MUSLCE can be download at http://www.drive5.com/muscle\n
            Add MUSCLE excutable path to $PATH using\n
            $  export PATH=$PATH:/path/to/muscle\n
            or\n
            provide path_to_MUSCLE in the ViralSeq::sequence_locator arguments\n
            "
      return nil
    end

    # check if reference option is correct
    begin
      case ref_option
      when 'HXB2'
        ori_ref = HXB2.dup
      when 'NL43'
        ori_ref = NL43.dup
      when 'MAC239'
        ori_ref = MAC239.dup
      else
        raise StandardError.new("reference sequence not recognized, choose from 'HXB2' (default), 'NL43', or 'MAC239'.")
      end
    rescue StandardError => e
      puts e.message
      return nil
    end

    begin
      ori_ref_l = ori_ref.size

      temp_dir = File.dirname($0)
      temp_file = temp_dir + "/_temp_muscle_in"
      temp_aln = temp_dir + "/_temp_muscle_aln"

      l1 = 0
      l2 = 0
      name = ">test"
      temp_in = File.open(temp_file,"w")
      temp_in.puts ">ref"
      temp_in.puts ori_ref
      temp_in.puts name
      temp_in.puts seq
      temp_in.close

      print `#{path_to_muscle} -in #{temp_file} -out #{temp_aln} -quiet`
      aln_seq = ViralSeq.fasta_to_hash(temp_aln)
      aln_test = aln_seq[name]
      aln_test =~ /^(\-*)(\w.*\w)(\-*)$/
      gap_begin = $1.size
      gap_end = $3.size
      aln_test2 = $2
      ref = aln_seq[">ref"]
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
        temp_in = File.open(temp_file,"w")
        temp_in.puts ">ref"
        temp_in.puts ref
        temp_in.puts name
        temp_in.puts seq
        temp_in.close
        print `#{path_to_muscle} -in #{temp_file} -out #{temp_aln} -quiet`
        aln_seq = ViralSeq.fasta_to_hash(temp_aln)
        aln_test = aln_seq[name]
        aln_test =~ /^(\-*)(\w.*\w)(\-*)$/
        gap_begin = $1.size
        gap_end = $3.size
        ref = aln_seq[">ref"]
        ref = ref[gap_begin..(-gap_end-1)]
      end
      aln_seq = ViralSeq.fasta_to_hash(temp_aln)
      aln_test = aln_seq[name]
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
      ref = aln_seq[">ref"]
      ref = ref[gap_begin..(-gap_end-1)]

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
        temp_in = File.open(temp_file,"w")
        temp_in.puts ">ref"
        temp_in.puts ref
        temp_in.puts name
        temp_in.puts seq
        temp_in.close
        print `#{path_to_muscle} -in #{temp_file} -out #{temp_aln} -quiet`
        aln_seq = ViralSeq.fasta_to_hash(temp_aln)
        aln_test = aln_seq[name]
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
        ref = aln_seq[">ref"]
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

      temp_in = File.open(temp_file,"w")
      temp_in.puts ">ref"
      temp_in.puts ref
      temp_in.puts name
      temp_in.puts seq
      temp_in.close
      print `#{path_to_muscle} -in #{temp_file} -out #{temp_aln} -quiet`
      aln_seq = ViralSeq.fasta_to_hash(temp_aln)
      aln_test = aln_seq[name]
      ref = aln_seq[">ref"]

      #refine alignment

      if ref =~ /^(\-+)/
        l1 = l1 - $1.size
      elsif ref =~ /(\-+)$/
        l2 = l2 + $1.size
      end

      if (ori_ref_l - l2 - 1) >= l1
        ref = ori_ref[l1..(ori_ref_l - l2 - 1)]
        temp_in = File.open(temp_file,"w")
        temp_in.puts ">ref"
        temp_in.puts ref
        temp_in.puts name
        temp_in.puts seq
        temp_in.close
        print `#{path_to_muscle} -in #{temp_file} -out #{temp_aln} -quiet`
        aln_seq = ViralSeq.fasta_to_hash(temp_aln)
        aln_test = aln_seq[name]
        ref = aln_seq[">ref"]

        ref_size = ref.size
        sim_count = 0
        (0..(ref_size-1)).each do |n|
          ref_base = ref[n]
          test_base = aln_test[n]
          sim_count += 1 if ref_base == test_base
        end
        similarity = (sim_count/ref_size.to_f*100).round(1)
        print `rm -f #{temp_file}`
        print `rm -f #{temp_aln}`
        loc_p1 = l1 + 1
        loc_p2 = ori_ref_l - l2
        if seq.size != (loc_p2 - loc_p1 + 1)
            indel = true
        elsif aln_test.include?("-")
            indel = true
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

end
