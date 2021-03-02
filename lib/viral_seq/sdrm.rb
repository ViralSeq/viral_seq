module ViralSeq
  class DRMs
    def initialize (mutation_list = {})
      @mutation_list = mutation_list
    end

    attr_accessor :mutation_list
  end

  def self.sdrm_hiv_pr(seq_hash)
  end

  def self.sdrm_hiv_rt(seq_hash)
  end

  def self.sdrm_hiv_in(seq_hash)
  end

  def self.list_from_json(file)
  end

  def self.list_from_csv(file)
  end

  def self.export_list_hiv_pr(file, format = :json)
    if foramt == :json

    end
  end

  def self.export_list_hiv_rt(file, format = :json)

  end

  def self.export_list_hiv_in(file, format = :json)

  end

  def drm_analysis(seq_hash)
    mutation_list = self.mutation_list

  end
end
