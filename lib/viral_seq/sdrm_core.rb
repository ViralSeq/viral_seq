# viral_seq/sdrm_core.rb
# core functions for HIV SDRM analysis using MPID-DR protocol. 


module ViralSeq

  # drug resistant mutation summary. input: amino acid array and starting codon, output, hash of summary
  def self.sdrm_nrti(aa_array,start_aa=1)
    out_hash = {}
    sdrm = {}
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
          test_aa_array = test_aa.split("/")
          if (test_aa_array & sdrm[position][1])
            out_hash[position] = [wt_aa,test_aa]
          end
        end

      end
    end
    return out_hash
  end

  def self.sdrm_nnrti(aa_array,start_aa=1)
    out_hash = {}
    sdrm = {}
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
          test_aa_array = test_aa.split("/")
          if (test_aa_array & sdrm[position][1])
            out_hash[position] = [wt_aa,test_aa]
          end
        end

      end
    end
    return out_hash
  end

  #HIV protease surveillance mutations

  def self.hiv_protease(aa_array,start_aa=1)
    out_hash = {}
    sdrm = {}
    sdrm[23] = ['L',['I']]
    sdrm[24] = ['L',['I']]
    sdrm[30] = ['D',['N']]
    sdrm[32] = ['V',['I']]
    sdrm[46] = ['M',['I','L','V']] # M46V not on the SDRM list but we still include it.
    sdrm[47] = ['I',['V','A']]
    sdrm[48] = ['G',['V','M']]
    sdrm[50] = ['I',['V','L']]
    sdrm[53] = ['F',['Y']]
    sdrm[54] = ['I',['V','L','M','T','A','S']]
    sdrm[73] = ['G',['S','T','C','A']]
    sdrm[76] = ['L',['V']]
    sdrm[82] = ['V',['A','T','S','F','L','C','M']]
    sdrm[83] = ['N',['D']]
    sdrm[84] = ['I',['V','A','C']]
    sdrm[85] = ['I',['V']]
    sdrm[88] = ['N',['D','S']]
    sdrm[90] = ['L',['M']]
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
          test_aa_array = test_aa.split("/")
          if (test_aa_array & sdrm[position][1])
            out_hash[position] = [wt_aa,test_aa]
          end
        end
      end
    end
    return out_hash
  end

  #HIV integrase drug resistance mutations

  def self.sdrm_int(aa_array,start_aa=1)
    out_hash = {}
    sdrm = {}
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
          test_aa_array = test_aa.split("/")
          if (test_aa_array & sdrm[position][1])
            out_hash[position] = [wt_aa,test_aa]
          end
        end

      end
    end
    return out_hash
  end

  # input sequence hash, and Poisson cutoff for minority variants.
  # HIV-1 PR region SDRM based on HIVDB.stanford.edu
  # only for MPID-DR MiSeq sequences, PR codon 1-99
  # return [substitution rate with 95% CI, halpotype abundance with 95% CI, amino acid sequence report spreadsheet]
  def self.sdrm_pr_bulk(sequences, cutoff = 0, temp_r_dir = File.dirname($0))
    region = "PR"
    rf_label = 0
    start_codon_number = 1
    n_seq = sequences.size
    mut = {}
    mut_com = []
    aa = {}
    point_mutation_list = []
    sequences.each do |name,seq|
      s = ViralSeq::Sequence.new(name,seq)
      s.get_aa_array(rf_label)
      aa_seq = s.aa_array
      aa[name] = aa_seq.join("")
      record = ViralSeq.hiv_protease(aa_seq)
      mut_com << record
      record.each do |position,mutation|
        if mut[position]
          mut[position][1] << mutation[1]
        else
          mut[position] = [mutation[0],[]]
          mut[position][1] << mutation[1]
        end
      end
    end
    mut.each do |position,mutation|
      wt = mutation[0]
      mut_list = mutation[1]
      count_mut_list = ViralSeq.count(mut_list)
      count_mut_list.each do |m,number|
        ci = ViralSeq.r_binom_CI(number, n_seq, temp_r_dir)
        label = number < cutoff ? "*" : ""
        point_mutation_list << [region, n_seq, position, wt, m, number, (number/n_seq.to_f).round(5), ci[0], ci[1], label]
      end
    end
    point_mutation_list.sort_by! {|record| record[2]}

    link = ViralSeq.count(mut_com)
    link2 = {}
    link.each do |k,v|
      pattern = []
      if k.size == 0
        pattern = ['WT']
      else
        k.each do |p,m|
          pattern << (m[0] + p.to_s + m[1])
        end
      end
      link2[pattern.join("+")] = v
    end
    linkage_list = []
    link2.sort_by{|_key,value|value}.reverse.to_h.each do |k,v|
      ci = ViralSeq.r_binom_CI(v, n_seq, temp_r_dir)
      label = v < cutoff ? "*" : ""
      linkage_list << [region, n_seq, k, v, (v/n_seq.to_f).round(5), ci[0], ci[1], label]
    end

    report_list = []

    div_aa = {}
    aa_start = start_codon_number

    aa_size = aa.values[0].size - 1

    (0..aa_size).to_a.each do |p|
      aas = []
      aa.values.each do |r1|
        aas << r1[p]
      end
      count_aas = ViralSeq.count(aas)
      div_aa[aa_start] = count_aas.sort_by{|_k,v|v}.reverse.to_h
      aa_start += 1
    end

    div_aa.each do |k,v|
      record = [region, k, n_seq]
      $amino_acid_list.each do |amino_acid|
        aa_count = v[amino_acid]
        record << (aa_count.to_f/n_seq*100).round(4)
      end
      report_list << record
    end

    return [point_mutation_list, linkage_list, report_list]
  end


  #input sequence hash, and Poisson cutoff for minority variants.
  #HIV-1 RT region SDRM based on HIVDB.stanford.edu
  #only for MPID-DR MiSeq sequences
  #RT codon 34-122, 152-236 two regions are linked.
  #return [substitution rate with 95% CI, halpotype abundance with 95% CI, amino acid sequence report spreadsheet]
  def self.sdrm_rt_bulk(sequences, cutoff = 0, temp_r_dir = File.dirname($0))
    region = "RT"
    rf_label = 1
    start_codon_number = 34
    gap = "AGACTTCAGGAAGTATACTGCATTTACCATACCTAGTATAAACAATGAGACACCAGGGATTAGATATCAGTACAATGTGCTTCCAC"

    n_seq = sequences.size
    mut_nrti = {}
    mut_nnrti = {}
    mut_com = []
    r1_aa = {}
    r2_aa = {}
    point_mutation_list = []
    sequences.each do |name,seq|
      r1 = seq[0,267]
      r2 = seq[267..-1]
      seq = r1 + gap + r2
      s = ViralSeq::Sequence.new(name,seq)
      s.get_aa_array(rf_label)
      aa_seq = s.aa_array

      r1_aa[name] = aa_seq[0,89].join("")
      r2_aa[name] = aa_seq[85..-1].join("")
      nrti = ViralSeq.sdrm_nrti(aa_seq,start_codon_number)
      nnrti = ViralSeq.sdrm_nnrti(aa_seq,start_codon_number)
      mut_com << (nrti.merge(nnrti))

      nrti.each do |position,mutation|
        if mut_nrti[position]
          mut_nrti[position][1] << mutation[1]
        else
          mut_nrti[position] = [mutation[0],[]]
          mut_nrti[position][1] << mutation[1]
        end
      end
      nnrti.each do |position,mutation|
        if mut_nnrti[position]
          mut_nnrti[position][1] << mutation[1]
        else
          mut_nnrti[position] = [mutation[0],[]]
          mut_nnrti[position][1] << mutation[1]
        end
      end
    end

    mut_nrti.each do |position,mutation|
      wt = mutation[0]
      mut_list = mutation[1]
      count_mut_list = ViralSeq.count(mut_list)
      count_mut_list.each do |m,number|
        ci = ViralSeq.r_binom_CI(number, n_seq, temp_r_dir)
        label = number < cutoff ? "*" : ""
        point_mutation_list << ["NRTI", n_seq, position, wt, m, number, (number/n_seq.to_f).round(5), ci[0], ci[1], label]
      end
    end

    mut_nnrti.each do |position,mutation|
      wt = mutation[0]
      mut_list = mutation[1]
      count_mut_list = ViralSeq.count(mut_list)
      count_mut_list.each do |m,number|
        ci = ViralSeq.r_binom_CI(number, n_seq, temp_r_dir)
        label = number < cutoff ? "*" : ""
        point_mutation_list << ["NNRTI", n_seq, position, wt, m, number, (number/n_seq.to_f).round(5), ci[0], ci[1], label]
      end
    end
    point_mutation_list.sort_by! {|record| record[2]}

    link = ViralSeq.count(mut_com)
    link2 = {}
    link.each do |k,v|
      pattern = []
      if k.size == 0
        pattern = ['WT']
      else
        k.each do |p,m|
          pattern << (m[0] + p.to_s + m[1])
        end
      end
      link2[pattern.join("+")] = v
    end
    linkage_list = []
    link2.sort_by{|_key,value|value}.reverse.to_h.each do |k,v|
      ci = ViralSeq.r_binom_CI(v, n_seq, temp_r_dir)
      label = v < cutoff ? "*" : ""
      linkage_list << [region, n_seq, k, v, (v/n_seq.to_f).round(5), ci[0], ci[1], label]
    end

    report_list = []

    div_aa = {}
    r1_aa_start = 34
    r2_aa_start = 152

    r1_aa_size = r1_aa.values[0].size - 1
    r2_aa_size = r2_aa.values[0].size - 1

    (0..r1_aa_size).to_a.each do |p|
      aas = []
      r1_aa.values.each do |r1|
        aas << r1[p]
      end
      count_aas = ViralSeq.count(aas)
      div_aa[r1_aa_start] = count_aas.sort_by{|_k,v|v}.reverse.to_h
      r1_aa_start += 1
    end

    (0..r2_aa_size).to_a.each do |p|
      aas = []
      r2_aa.values.each do |r1|
        aas << r1[p]
      end
      count_aas = ViralSeq.count(aas)
      div_aa[r2_aa_start] = count_aas.sort_by{|_k,v|v}.reverse.to_h
      r2_aa_start += 1
    end

    div_aa.each do |k,v|
      record = [region, k, n_seq]
      $amino_acid_list.each do |amino_acid|
        aa_count = v[amino_acid]
        record << (aa_count.to_f/n_seq*100).round(4)
      end
      report_list << record
    end

    return [point_mutation_list, linkage_list, report_list]
  end

  #input sequence hash, and Poisson cutoff for minority variants.
  #HIV-1 IN region SDRM based on HIVDB.stanford.edu
  #only for MPID-DR MiSeq sequences
  #IN codon 53-174
  #return [substitution rate with 95% CI, halpotype abundance with 95% CI, amino acid sequence report spreadsheet]
  def self.sdrm_in_bulk(sequences, cutoff = 0, temp_r_dir = File.dirname($0))
    region = "IN"
    rf_label = 2
    start_codon_number = 53
    n_seq = sequences.size
    mut = {}
    mut_com = []
    aa = {}
    point_mutation_list = []
    sequences.each do |name,seq|
      s = ViralSeq::Sequence.new(name,seq)
      s.get_aa_array(rf_label)
      aa_seq = s.aa_array
      aa[name] = aa_seq.join("")
      record = ViralSeq.sdrm_int(aa_seq, start_codon_number)
      mut_com << record
      record.each do |position,mutation|
        if mut[position]
          mut[position][1] << mutation[1]
        else
          mut[position] = [mutation[0],[]]
          mut[position][1] << mutation[1]
        end
      end
    end
    mut.each do |position,mutation|
      wt = mutation[0]
      mut_list = mutation[1]
      count_mut_list = ViralSeq.count(mut_list)
      count_mut_list.each do |m,number|
        ci = ViralSeq.r_binom_CI(number, n_seq, temp_r_dir)
        label = number < cutoff ? "*" : ""
        point_mutation_list << [region, n_seq, position, wt, m, number, (number/n_seq.to_f).round(5), ci[0], ci[1], label]
      end
    end
    point_mutation_list.sort_by! {|record| record[2]}

    link = ViralSeq.count(mut_com)
    link2 = {}
    link.each do |k,v|
      pattern = []
      if k.size == 0
        pattern = ['WT']
      else
        k.each do |p,m|
          pattern << (m[0] + p.to_s + m[1])
        end
      end
      link2[pattern.join("+")] = v
    end
    linkage_list = []
    link2.sort_by{|_key,value|value}.reverse.to_h.each do |k,v|
      ci = ViralSeq.r_binom_CI(v, n_seq, temp_r_dir)
      label = v < cutoff ? "*" : ""
      linkage_list << [region, n_seq, k, v, (v/n_seq.to_f).round(5), ci[0], ci[1], label]
    end

    report_list = []

    div_aa = {}
    aa_start = start_codon_number

    aa_size = aa.values[0].size - 1

    (0..aa_size).to_a.each do |p|
      aas = []
      aa.values.each do |r1|
        aas << r1[p]
      end
      count_aas = ViralSeq.count(aas)
      div_aa[aa_start] = count_aas.sort_by{|_k,v|v}.reverse.to_h
      aa_start += 1
    end

    div_aa.each do |k,v|
      record = [region, k, n_seq]
      $amino_acid_list.each do |amino_acid|
        aa_count = v[amino_acid]
        record << (aa_count.to_f/n_seq*100).round(4)
      end
      report_list << record
    end

    return [point_mutation_list, linkage_list, report_list]
  end

  # input a sequence hash, return a sequence hash with stop codons.
  def self.stop_codon_seq_hash(seq_hash, rf = 0)
    out_seq_hash = {}
    seq_hash.each do |k,v|
      sequence = Sequence.new(k,v)
      sequence.get_aa_array(rf)
      if sequence.aa_array.include?("*")
        out_seq_hash[k] = v
      end
    end
    return out_seq_hash
  end

end
