
module ViralSeq

  # DRM configuration for each region

  class DrmRegionConfig

    # initialize DRM region configuration
    # @param drm_version [String] version of the instance of DrmVersion
    # @param region [String] name of the region
    # @param drm_class [Array] classes of DRMs at this region
    # @param drm_range [Hash] DRM range for each class of DRMs at this region
    # @param drm_list [Hash] List of detailed DRM mutations for each DRM classes at this region
    # @param seq_drm_corrlation [Hash] correlation of sequenced region and DRM class
    # @param ref_info [Hash] information of the reference genome, including sequence coordinates on HXB2
    def initialize(drm_version, region, drm_class, drm_range, drm_list, seq_coord, ref_info)
      @drm_version = drm_version
      @region = region
      @drm_class = drm_class
      @drm_range = drm_range
      @drm_list = drm_list
      @seq_coord = seq_coord
      @ref_info = ref_info
    end

    attr_accessor :drm_version, :region, :drm_class, :drm_range, :drm_list, :seq_coord, :ref_info

    # summarize the DRM information for the output as JSON for the specific version
    # @return [Hash] json has for DRM inforation of each position
    def drm_json
      sdrm = self.drm_list
      json_hash = {}
      sdrm.each do |drm_class, drms|
        json_hash[drm_class] = []
        drms.each do |pos, muts|
          mutation = {}
          mutation[:position] = pos
          mutation[:wildtypeCodon] = muts[0]
          mutation[:mutationCodons] = muts[1]
          json_hash[drm_class] << mutation
        end
      end
      return json_hash
    end

    # calculate the length of R1 and R2 based on the sequence coordinates
    # @return [Hash] {r1_length: [Integer], r2_length: [Integer]}
    def r1_r2_length
      seq_coord = self.seq_coord
      return nil unless seq_coord["gap"]

      r1_length = seq_coord["gap"]["minimum"] - seq_coord["minimum"]
      r2_length = seq_coord["maximum"] - seq_coord["gap"]["maximum"]

      return {r1_length: r1_length, r2_length: r2_length}
    end #end of #r1_r2_length


    # determine the reading frame number based on the sequence coordinates
    # @return [Integer] reading frame of 0, 1 or 2
    def get_reading_frame_number
      m1 = (self.seq_coord["minimum"] - self.ref_info["ref_coord"][0]) % 3
      if m1.zero?
        n1 = 0
      else
        n1 = 3 - m1
      end

      if seq_coord["gap"]
        m2 = (self.seq_coord["gap"]["maximum"] + 1 - self.ref_info["ref_coord"][0]) % 3
        if m2.zero?
          n2 = 0
        else
          n2 = 3 - m2
        end
        return [n1, n2]
      else
        return [n1]
      end
    end #end get_reading_frame_number

  end
end
