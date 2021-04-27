module ViralSeq
  class DRMs
    class << self

      # function to retrieve sdrm positions as a hash
      # @param ref_option [Symbol], name of reference genomes, options are `:hiv_pr`, `:hiv_rt`, `:hiv_in`, `hcv_ns5a`
      # @return [Hash] Hash of :position_number => [ 'wildtype_codon', ['mutation_codons']]
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
      end # end of #sdrm_hash

      # function to export SDRM positions as json object
      # @param (see #sdrm_hash)
      # @return [String] json String of SDRM positions

      def sdrm_json(options)
        sdrm = ViralSeq::DRMs.sdrm_hash(options)
        json_array = []
        sdrm.each do |pos, muts|
          mutation = {}
          mutation[:position] = pos
          mutation[:wildtypeCodon] = muts[0]
          mutation[:mutationCodons] = muts[1]
          json_array << mutation
        end
        JSON.pretty_generate(json_array)
      end
    end
  end
end
