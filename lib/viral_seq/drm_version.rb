
module ViralSeq

  # DRM version configuration.
  # Configuration files are located at `lib/viral_seq/drm_versions_config.json`

  class DrmVersion

    # initialize a ViralSeq::DrmVersion instance
    # @param drm_version [String] version of the instance of DrmVersion
    # @param drm_range [Hash] region/class of DRM and the range of amino acid positions included in this version.
    # @param seq_coord [Hash] region and its amplicon positions on HXB2 reference
    # @param seq_drm_corrlation [Hash] correlation of sequenced region and DRM class
    # @param ref_info [Hash] information of the reference genome, including sequence coordinates on HXB2
    def initialize(drm_version, drm_range, seq_coord, seq_drm_correlation, ref_info)
      @drm_version = drm_version
      @drm_range = drm_range
      @seq_coord = seq_coord
      @seq_drm_correlation = seq_drm_correlation
      @ref_info = ref_info
    end

    attr_accessor :drm_version, :drm_range, :seq_coord, :seq_drm_correlation, :ref_info

    # construct an instance of ViralSeq::DrmVersion
    # @param version_config_hash [Hash] json hash of stored version configurations.
    # @return [ViralSeq::DrmVersion] an instance of constructed DrmVersion

    def self.construct(version_config_hash)
      drm_version = version_config_hash["version"]
      drm_range = version_config_hash["DRM_range"]
      seq_coord = version_config_hash["seq_coord"]
      seq_drm_correlation = version_config_hash["seq_drm_correlation"]
      ref_info = version_config_hash["ref_info"]
      ViralSeq::DrmVersion.new(drm_version, drm_range, seq_coord, seq_drm_correlation, ref_info)
    end

    # construct a specific version of ViralSeq::DrmVersion
    # @param v [String] version string
    # @return [ViralSeq::DrmVersion] an instance of constructed DrmVersion

    def self.config_version(v="v1")
      v = v.downcase
      v = "v1" if v == "v2"

      drm_config = JSON.parse(
        File.read(
          File.join( ViralSeq.root, 'viral_seq', 'util', 'drm_versions_config.json')
          )
      )

      drm_versions = {}

      drm_config.each do |config|
        drm_versions[config["version"]] = ViralSeq::DrmVersion.construct(config)
      end

      if drm_versions[v]
        drm_versions[v]
      else
        abort (
        "Version '#{v}' config not found. Program aborted. \nCurrent supported versions '#{drm_versions.keys.sort.join(", ")}'\nCheck documentations for details".red
        )
      end
    end

    # construct a ViralSeq::DrmRegionConfig instance from a specific version
    # @param region [String] name of the region
    # @return [ViralSeq::DrmRegionConfig] an instance of DrmRegionConfig

    def query_region(region)
      region = region.to_s.upcase
      drm_classes = self.seq_drm_correlation[region]

      if drm_classes.nil?
        abort "Region not recognized by the specific DRM config version. Program aborted."
      end

      drm_range = {}
      drm_list = {}

      drm_classes.each do |drm_class|
        drm_range[drm_class] = self.drm_range[drm_class]
        drm_list_single_class = ViralSeq::DRMs.sdrm_hash(drm_class)

        drm_list[drm_class] = drm_list_single_class.select { |k, _v| drm_range[drm_class].include? k }

      end


      seq_coord = self.seq_coord[region]

      ref_info = {}
      ref_info["ref_type"] = self.ref_info["ref_type"]
      ref_info["ref_coord"] = self.ref_info["ref_coord"][region]


      ViralSeq::DrmRegionConfig.new(
        self.drm_version, region, drm_classes, drm_range, drm_list, seq_coord, ref_info
      )
    end

    # summarize the DRM information for the output as JSON
    # @return [Hash] json has for DRM inforation of each position

    def pull_drm_json

      summary_json_hash = {}

      self.seq_drm_correlation.keys.each do |region|
        summary_json_hash = summary_json_hash.merge query_region(region).drm_json
      end

      summary_json_hash

    end

  end # end of class

end
