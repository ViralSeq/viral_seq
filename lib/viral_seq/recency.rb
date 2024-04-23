module ViralSeq

  # recency prediction function based on HIV MPID-NGS
  # @see https://pubmed.ncbi.nlm.nih.gov/32663847 Ref: Zhou et al. J Infect Dis. 2021

  module Recency

    # @params tcs_RT [Integer] number of TCS at the RT region
    # @params tcs_V1V3 [Integer] number of TCS at the V1V3 region
    # @params pi_RT [Float] pairwise diversity at the RT region
    # @params pi_V1V3 [Float] pairwise diversity at the V1V3 region
    # @params dist20_RT [Float] dist20 at the RT region
    # @params dist20_V1V3 [Float] dist20 at the V1V3 region
    # @return [String] determination of the recency

    def self.define(tcs_RT: nil,
                     tcs_V1V3: nil,
                     pi_RT: nil,
                     dist20_RT: nil,
                     pi_V1V3: nil,
                     dist20_V1V3: nil)
      tcs_RT ||= 0
      tcs_V1V3 ||= 0
      if (tcs_RT >= 3 && pi_RT) and (tcs_V1V3 >= 3 && pi_V1V3)
        if (pi_RT + pi_V1V3) < 0.0103
            recency = "recent"
        elsif (pi_RT + pi_V1V3) >= 0.0103 and (dist20_RT + dist20_V1V3) >= 0.006
            recency = "chronic"
        else
            recency = "indeterminant"
        end
      elsif (tcs_RT >= 3 && pi_RT) and tcs_V1V3 < 3
        if pi_RT < 0.0021
          recency = "recent"
        elsif pi_RT >= 0.0021 and dist20_RT >= 0.001
          recency = "chronic"
        else
          recency = "indeterminant"
        end
      elsif (tcs_V1V3 >= 3 && pi_V1V3)
        if pi_V1V3 >= 0.0103 and dist20_V1V3 >= 0.006
          recency = "chronic"
        else
          recency = "insufficient data"
        end
      else
        recency = "insufficient data"
      end
      return recency
    end


    def self.dpi(pi_rt, pi_v1v3)

      if pi_rt.is_a? Numeric and pi_v1v3.is_a? Numeric
        pi = pi_rt*100 + pi_v1v3*100
        model_file = "rt_v1v3_fit.Rdata"
        var = "combined_perc"
      elsif pi_rt.is_a? Numeric
        pi = pi_rt*100
        model_file = "rt_only_fit.Rdata"
        var = "rtperc"
      elsif pi_v1v3.is_a? Numeric
        pi = pi_v1v3*100
        model_file = "v1v3_only_fit.Rdata"
        var = "v1v3perc"
      else
        return ["NA", "NA", "NA"]
      end
      path = File.join("lib", "viral_seq", "util", "recency_model", model_file)
      data_str = `Rscript -e 'fit = readRDS("#{path}"); test = data.frame(#{var} = #{pi}); pre= predict(fit, test, interval = "prediction", level = 0.9); cat(pre)'`
      dpi_array = data_str.split("\s")
      dpi = dpi_array[0].to_f
      lwr = dpi_array[1].to_f
      upr = dpi_array[2].to_f
      if lwr < 0
        return [dpi, 0.0, upr]
      else
        return [dpi, lwr, upr]
      end
    end
  end
end
