module ViralSeq
  class TcsJson
    class << self

      def generate
        puts '-'*58
        puts '| JSON Parameter Generator for ' + "TCS #{ViralSeq::TCS_VERSION}".red.bold + " by " + "Shuntai Zhou".blue.bold + ' |'
        puts '-'*58 + "\n"

        param = {}

        puts 'Enter the path to the directory that contains the MiSeq pair-end R1 and R2 .fastq or .fastq.gz file'
        print '> '
        param[:raw_sequence_dir] = gets.chomp.rstrip

        puts "Choose MiSeq Platform (1-3):\n1. 150x7x150\n2. 250x7x250\n3. 300x7x300 (default)"
        print "> "
        pf_option = gets.chomp.rstrip
        # while ![1,2,3].include?(pf_option.to_i)
        #   print "Entered MiSeq Platform #{pf_option.red.bold} not valid (choose 1-3), try again\n> "
        #   pf_option = gets.chomp.rstrip
        # end
        case pf_option.to_i
        when 1
          param[:platform_format] = 150
        when 2
          param[:platform_format] = 250
        else
          param[:platform_format] = 300
        end

        puts 'Enter the estimated platform error rate (for TCS cut-off calculation), default as ' + '0.02'.red.bold
        print '> '
        input_error = gets.chomp.rstrip.to_f
        if input_error == 0.0
          param[:platform_error_rate] = 0.02
        else
          param[:platform_error_rate] = input_error
        end

        param[:primer_pairs] = []

        loop do
          data = {}
          puts "Enter the name for the sequenced region: "
          print '> '
          data[:region] = gets.chomp.rstrip

          puts "Enter the #{"cDNA".red.bold} primer sequence: "
          print '> '
          data[:cdna] = gets.chomp.rstrip

          puts "Enter the #{"forward".blue.bold} primer sequence: "
          print '> '
          data[:forward] = gets.chomp.rstrip

          puts "Enter supermajority cut-off (0.5 - 1.0). Default Simple Majority"
          print '> '
          mj = gets.chomp.rstrip.to_f
          if (0.5..1.0).include?(mj)
            data[:majority] = mj
          else
            data[:majority] = 0
          end

          print "Need end-join? Y/N \n> "
          ej = gets.chomp.rstrip
          if ej =~ /y|yes/i
            data[:end_join] = true

            puts "End-join option? Choose from (1-4):"
            puts "1: simple join, no overlap"
            puts "2: known overlap"
            puts "3: unknow overlap, use sample consensus to determine overlap, all sequence pairs have same overlap"
            puts "4: unknow overlap, determine overlap by individual sequence pairs, sequence pairs can have different overlap"
            print "> "
            ej_option = gets.chomp.rstrip
            while ![1,2,3,4].include?(ej_option.to_i)
              puts "Entered end-join option #{ej_option.red.bold} not valid (choose 1-4), try again"
              ej_option = gets.chomp.rstrip.to_i
            end
            case ej_option.to_i
            when 1
              data[:end_join_option] = 1
              data[:overlap] = 0
            when 2
              data[:end_join_option] = 1
              print "overlap bases: \n> "
              ol = gets.chomp.rstrip.to_i
              data[:overlap] = ol
            when 3
              data[:end_join_option] = 3
            when 4
              data[:end_join_option] = 4
            end

            print "Need QC for TCS? (support for HIV-1 and SIV)? Y/N \n> "
            qc = gets.chomp.rstrip
            if qc =~ /y|yes/i
              data[:TCS_QC] = true

              data[:ref_genome] = get_ref

              print "reference 5'end ref position or position range, 0 if no need to match this end \n> "
              data[:ref_start] = gets.chomp.rstrip.to_i

              print "reference 3'end ref position or position range: 0 if no need to match this end \n> "
              data[:ref_end] = gets.chomp.rstrip.to_i

              print "allow indels? (default as yes) Y/N \n> "
              indel = gets.chomp.rstrip
              if indel =~ /n|no/i
                data[:indel] = false
              else
                data[:indel] = true
              end
            else
              data[:TCS_QC] = false
            end

            print "Need trimming to a reference genome? Y/N \n> "
            trim_option = gets.chomp.rstrip
            if trim_option =~ /y|yes/i
              data[:trim] = true
              data[:trim_ref] = get_ref

              print "reference 5'end ref position \n> "
              data[:trim_ref_start] = gets.chomp.rstrip.to_i

              print "reference 3'end ref position \n> "
              data[:trim_ref_end] = gets.chomp.rstrip.to_i

            else
              data[:trim] = false
            end

          else
            data[:end_join] = false
          end

          param[:primer_pairs] << data
          print "Do you wish to conintue? Y/N \n> "
          continue_sig = gets.chomp.rstrip
          break unless continue_sig =~ /y|yes/i

        end

        puts "\nYour JSON string is:"
        puts JSON.pretty_generate(param)

        print "\nDo you wish to save it as a file? Y/N \n> "
        save_option = gets.chomp.rstrip

        if save_option =~ /y|yes/i
          print "Path to save JSON file:\n> "
          path = gets.chomp.rstrip
          while !validate_path_name(path)
            print "Entered path no valid, try again.\n".red.bold
            print "Path to save JSON file:\n> "
            path = gets.chomp.rstrip
          end
          File.open(validate_path_name(path), 'w') {|f| f.puts JSON.pretty_generate(param)}
        end

        print "\nDo you wish to execute tcs pipeline with the input params now? Y/N \n> "

        rsp = gets.chomp.rstrip
        if rsp =~ /y/i
          return param
        else
          abort "Params json file generated. You can execute tcs pipeline using `tcs -p [params.json]`".blue
        end

      end

      private
      def get_ref
        puts "Choose reference genome (1-3):"
        puts "1. HIV-1 HXB2".red.bold
        puts "2. HIV-1 NL4-3".blue.bold
        puts "3. SIV MAC239".magenta.bold
        print "> "
        ref_option = gets.chomp.rstrip
        while ![1,2,3].include?(ref_option.to_i)
          print "Entered end-join option #{ref_option.to_s.red.bold} not valid (choose 1-3), try again\n> "
          ref_option = gets.chomp.rstrip.to_i
        end
        ref = case ref_option.to_i
              when 1
                :HXB2
              when 2
                :NL43
              when 3
                :MAC239
              end
      end # end of get_ref

      def validate_path_name(path)
        if path.empty?
          return false
        elsif File.directory? path
          return File.join(path, 'params.json')
        elsif File.directory?(File.dirname(path))
          return path
        end
      end # end of validate_path_name
    end # end of class << self
  end # end TcsJson
end # end main module
