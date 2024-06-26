module ViralSeq

  # class to generate recency report

  class RecencyReport

    # to generate the recency report in .pdf format.
    # @param log [Hash] Hash from the json summary string of the SDRM report
    # @param outfile [String] path to the output file
    # @return [NilClass] .pdf file generated by the method. Return nil.

    def self.generate(log, outfile)

      recency_color = {
        "recent" => "d42828",
        "chronic" => "0666bf",
        "indeterminant"=> "f78914",
        "insufficient data" => "7d7b79"
      }

      dual_infection_color = {
        "Yes" => "ffcc00",
        "No" => "339900",
        "insufficient data" => "7d7b79"
      }

      Prawn::Document.generate(outfile, margin: 75) do

        def text_format(text1, text2)
          [
            { text: text1 + "\s" * (30 - text1.size), styles: [:bold], size: 14, font: "Courier"},
            { text: text2, size: 14, styles: [:underline]}
          ]
        end

        def text_format2(text1, text2, text3, text4)
            text1 = text1.to_s
            text2 = text2.to_s
            text3 = text3.to_s
            text4 = text4.to_s

            [
              { text: "\s\s\s" + text1 + "\s"*(11-text1.size) +
            text2 + "\s"*(19-text2.size) +
            text3 + "\s"*(11-text3.size) + text4,
              size: 14,
              font: "Courier"
              }
            ]
        end

        text("Quantitative Recency Report by MPID-NGS",
        size: 18,
        align: :center,
        style: :bold
        )

        move_down 20

        formatted_text(
          text_format("Library ID:", log[:sample_id])
        )

        move_down 10

        formatted_text(
          text_format("ViralSeq Version:", ViralSeq::VERSION.to_s)

        )

        formatted_text(
          text_format("TCS Version:", ViralSeq::TCS_VERSION.to_s)
        )

        formatted_text(
          text_format("Processed Date", Time.now.strftime("%Y-%b-%d %H:%M"))
        )

        move_down 30

        text("Summary of parameters",
        size: 16,
        style: :bold
        )

        move_down 20

        formatted_text(
          [
            { text: "REGION" + "\s"*5 + "AVG. DIVERSITY" + "\s"*5 + "DIST20" + "\s"*5 + "DEPTH",
              styles: [:bold],
              size: 14,
              font: "Courier"
            },
          ]
        )

        move_down 5

        formatted_text(
          text_format2("RT", log[:pi_RT], log[:dist20_RT], log[:tcs_RT])
        )

        formatted_text(
          text_format2("V1V3", log[:pi_V1V3], log[:dist20_V1V3], log[:tcs_V1V3])
        )

        formatted_text(
          text_format2("P17", log[:pi_P17], log[:dist20_P17], log[:tcs_P17])
        )

        move_down 30

        formatted_text(
          [
            { text: "Prediction: ",
              styles: [:bold],
              size: 16,
            },

            { text: log[:recency].capitalize + " Infection",
              styles: [:bold],
              size: 16,
              color: recency_color[log[:recency]]
            },

            { text: " (9-month cutoff)",
              size: 14,
            },
          ]
        )

        move_down 20

        formatted_text(
          [
            {
              text: "Estimated Day Post Infection: ",
              styles: [:bold],
              size: 16
            },

            {
              text: log[:dpi].to_s +
              " (" + log[:dpi_lwr].to_s + "-" + log[:dpi_upr].to_s + ") Days",
              styles: [:bold],
              size: 16,
              color: recency_color[log[:recency]]
            }
          ]
        )

        move_down 20

        formatted_text(
          [
            {
              text: "Possible multivariant Infection: ",
              styles: [:bold],
              size: 16,
            },

            {
              text: log[:possible_dual_infection],
              styles: [:bold],
              size: 16,
              color: dual_infection_color[log[:possible_dual_infection]]
            }
          ]
        )

        move_down 10

        if log[:possible_dual_infection] == "Yes"

          formatted_text(
            [
              {
                text: "Warning: Days Post Infection prediction not reliable!",
                styles: [:bold],
                size: 14,
                color: "ffcc00"
              }
            ]
          )
        end
      end

    end

  end

end
