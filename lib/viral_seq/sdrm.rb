module ViralSeq
  class DRMs

    # function to retrieve sdrm positions as a hash, DRM list are stored at `lib/viral_seq/util/drm_list.json`
    # @param ref_option [Symbol], name of reference genomes, options are `:hiv_pr`, `:hiv_rt`, `:hiv_in`, `hcv_ns5a`
    # @return [Hash] Hash of :position_number => [ 'wildtype_codon', ['mutation_codons']]

    def self.sdrm_hash(options)
      options = options.to_s.upcase
      drm_data = JSON.parse(
        File.read(
          File.join(ViralSeq.root, 'viral_seq', 'util', 'drm_list.json')
        )
      )
      if drm_data[options]
        sdrm = {}
        drm_data[options].each do |record|
          sdrm[record["position"]] = [record["wild-type"], record["mutations"]]
        end

      else
        abort "Input option `#{options}` for ViralSeq::DRMs.sdrm_hash not supported. Program aborted.\nSupported type of mutations for '#{drm_data.keys.join(", ")}' only."
      end
      return sdrm
    end # end of #sdrm_hash

    # function to export SDRM positions as json object
    # @param (see #sdrm_hash)
    # @return [Array] json Array of SDRM positions

    def self.sdrm_json(options)
      sdrm = ViralSeq::DRMs.sdrm_hash(options)
      json_array = []
      sdrm.each do |pos, muts|
        mutation = {}
        mutation[:position] = pos
        mutation[:wildtypeCodon] = muts[0]
        mutation[:mutationCodons] = muts[1]
        json_array << mutation
      end
      return json_array
    end #end of #sdrm_json
  end
end
