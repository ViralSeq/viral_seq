# viral_seq/hcv_dr
# HCV resistant mutation interpretation
# ViralSeq::hcv_ns5a

# ViralSeq.hcv_ns5a(amino_acid_sequence_array, start_aa_position)
#   # amino_acid_sequence_array is Array object of the amino acid sequence.
#   # can use ViralSeq::Sequence#aa_array to obtain the aa array sequence
#   # start_aa_position is the starting aa number of the input sequence as Integer

module ViralSeq
  def self.hcv_ns5a(aa_array,start_aa=1)
    out_hash = {}
    sdrm = {}
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
end
