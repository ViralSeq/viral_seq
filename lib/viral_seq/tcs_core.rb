module ViralSeq

  # Core functions for `tcs` pipeline

  class TcsCore
    class << self

      # methods to calculate TCS consensus cut-off based on the maximum numbers of PIDs and platform error rate.
      # @see https://www.ncbi.nlm.nih.gov/pubmed/26041299 reference at Zhou et al. JVI 2016.
      # @param m [Integer] PID abundance
      # @param error_rate [Float] estimated platform error rate.
      # @return [Integer] an abundance cut-off (Integer) for offspring Primer IDs.

      def calculate_cut_off(m, error_rate = 0.02)
        n = 0
        case error_rate
        when 0.005...0.015
          if m <= 10
            n = 2
          else
            n = 1.09*10**-26*m**6 + 7.82*10**-22*m**5 - 1.93*10**-16*m**4 + 1.01*10**-11*m**3 - 2.31*10**-7*m**2 + 0.00645*m + 2.872
          end

        when 0...0.005
          if m <= 10
            n = 2
          else
            n = -9.59*10**-27*m**6 + 3.27*10**-21*m**5 - 3.05*10**-16*m**4 + 1.2*10**-11*m**3 - 2.19*10**-7*m**2 + 0.004044*m + 2.273
          end

        else
          if m <= 10
            n = 2
          elsif m <= 8500
            n = -1.24*10**-21*m**6 + 3.53*10**-17*m**5 - 3.90*10**-13*m**4 + 2.12*10**-9*m**3 - 6.06*10**-6*m**2 + 1.80*10**-2*m + 3.15
          else
            n = 0.0079 * m + 9.4869
          end
        end

        n = n.round
        n = 2 if n < 3
        return n
      end

      # identify which file in the directory is R1 file, and which is R2 file based on file names
      # input as directory (Dir object or a string of path)
      # by default, .gz files will be unzipped.
      # return as an hash of {r1_file: file1, r1_file: file2}
      def r1r2(directory, unzip = true)
        files = []
        Dir.chdir(directory) { files = Dir.glob "*" }
        r1_file = ""
        r2_file = ""
        files.each do |f|
          tag = parser_file_name(f)[:tag]

          if tag.include? "R1"
            unzip ? r1_file = unzip_r(directory, f) : r1_file = File.join(directory, f)
          elsif tag.include? "R2"
            unzip ? r2_file = unzip_r(directory, f) : r2_file = File.join(directory, f)
          end
        end
        return { r1_file: r1_file, r2_file: r2_file }
      end # end of ViralSeq:TcsCore.r1r2

      # sort directories containing mulitple r1 and r2 files.
      # use the library name (first string before "_") to seperate libraries
      # out_dir is the Dir object or string of the output directory, by default named as directory + "_sorted"
      # return a hash as { with_both_r1_r2: [lib1, lib2, ...], missing_r1: [lib1, lib2, ...], missing_r2: [lib1, lib2, ...], error: [lib1, lib2, ...]}

      def sort_by_lib(directory, out_dir = directory + "_sorted")
        Dir.mkdir(out_dir) unless File.directory?(out_dir)
        files = []
        Dir.chdir(directory) {files = Dir.glob("*")}

        files.each do |file|
          path = File.join(directory,file)
          index = file.split("_")[0]
          index_dir = File.join(out_dir, index)
          Dir.mkdir(index_dir) unless File.directory?(index_dir)
          File.rename(path, File.join(index_dir, file))
        end

        return_obj = { with_both_r1_r2: [],
                       missing_r1: [],
                       missing_r2: [],
                       error: []
                      }

        libs = []
        Dir.chdir(out_dir) { libs = Dir.glob('*') }
        libs.each do |lib|
          file_check = ViralSeq::TcsCore.r1r2(File.join(out_dir, lib))
          if !file_check[:r1_file].empty? and !file_check[:r2_file].empty?
            return_obj[:with_both_r1_r2] << lib
          elsif file_check[:r1_file].empty? and !file_check[:r2_file].empty?
            return_obj[:missing_r1] << lib
          elsif file_check[:r2_file].empty? and !file_check[:r1_file].empty?
            return_obj[:missing_r2] << lib
          else
            return_obj[:error] << lib
          end
        end
        return return_obj
      end

      # sort array of file names to determine if there is potential errors
      # @param name_array [Array] array of file names
      # @return [hash] name check results

      def validate_file_name(name_array)
        errors = {
                   file_type_error: [] ,
                   missing_r1_file: [] ,
                   missing_r2_file: [] ,
                   extra_r1_r2_file: [],
                   no_region_tag: [] ,
                   multiple_region_tag: []
                 }

        passed_libs = {}

        name_with_r1_r2 = []

        name_array.each do |name|
          tag = parser_file_name(name)[:tag]
          if name !~ /\.fastq\Z|\.fastq\.gz\Z/
            errors[:file_type_error] << name
          elsif tag.count("R1") == 0 and tag.count("R2") == 0
            errors[:no_region_tag] << name
          elsif tag.count("R1") > 0 and tag.count("R2") > 0
            errors[:multiple_region_tag] << name
          elsif tag.count("R1") > 1 or tag.count("R2") > 1
            errors[:multiple_region_tag] << name
          else
            name_with_r1_r2 << name
          end
        end

        libs = {}

        name_with_r1_r2.map do |name|
          libname = parser_file_name(name)[:libname]
          libs[libname] ||= []
          libs[libname] << name
        end

        libs.each do |libname, files|
          count_r1_file = 0
          count_r2_file = 0
          files.each do |name|
            tag = parser_file_name(name)[:tag]
            if tag.include? "R1"
              count_r1_file += 1
            elsif tag.include? "R2"
              count_r2_file += 1
            end
          end

          if count_r1_file > 1 or count_r2_file > 1
            errors[:extra_r1_r2_file] += files
          elsif count_r1_file.zero?
            errors[:missing_r1_file] += files
          elsif count_r2_file.zero?
            errors[:missing_r2_file] += files
          else
            passed_libs[libname] = files
          end
        end

        file_name_with_lib_name = {}
        passed_libs.each do |lib_name, files|
          files.each do |f|
            file_name_with_lib_name[f] = lib_name
          end
        end

        passed_names = []

        passed_libs.values.each { |names| passed_names += names}

        if passed_names.size < name_array.size
          pass = false
        else
          pass = true
        end

        file_name_with_error_type = {}

        errors.each do |type, files|
          files.each do |f|
            file_name_with_error_type[f] ||= []
            file_name_with_error_type[f] << type.to_s.tr("_", "\s")
          end
        end

        file_check = []

        name_array.each do |name|
          file_check_hash = {}
          file_check_hash[:fileName] = name
          file_check_hash[:errors] = file_name_with_error_type[name]
          file_check_hash[:libName] = file_name_with_lib_name[name]

          file_check << file_check_hash
        end

        return { allPass: pass, files: file_check }
      end

      # filter r1 raw sequences for non-specific primers.
      # input r1_sh, SeqHash obj.
      # return filtered Hash of sequence name and seq pair, in the object { r1_filtered_seq: r1_filtered_seq_pair }

      def filter_r1(r1_sh, forward_primer)
        if forward_primer.match(/(N+)(\w+)$/)
          forward_n = $1.size
          forward_bio_primer = $2
        else
          forward_n = 0
          forward_bio_primer = forward_primer
        end
        forward_bio_primer_size = forward_bio_primer.size
        forward_starting_number = forward_n + forward_bio_primer_size
        forward_primer_ref = forward_bio_primer.nt_parser

        r1_passed_seq = {}
        r1_raw = r1_sh.dna_hash

        proc_filter = proc do |name|
          seq = r1_raw[name]
          next unless general_filter seq
          primer_region_seq = seq[forward_n, forward_bio_primer_size]
          if primer_region_seq =~ forward_primer_ref
            new_name = remove_tag name
            r1_passed_seq[new_name] = seq
          end
        end

        r1_raw.keys.map do |name|
          proc_filter.call name
        end

        return { r1_passed_seq: r1_passed_seq, forward_starting_number: forward_starting_number }
      end # end of filter_r1

      # filter r2 raw sequences for non-specific primers.
      # input r2_sh, SeqHash obj.
      # return filtered Hash of sequence name and seq pair, as well as the length of PID.
      def filter_r2(r2_sh, cdna_primer)
        r2_raw = r2_sh.dna_hash
        cdna_primer.match(/(N+)(\w+)$/)
        pid_length = $1.size
        cdna_bio_primer = $2
        cdna_bio_primer_size = cdna_bio_primer.size
        reverse_starting_number = pid_length + cdna_bio_primer_size
        cdna_primer_ref = cdna_bio_primer.nt_parser
        r2_passed_seq = {}
        proc_filter = proc do |name|
          seq = r2_raw[name]
          next unless general_filter seq
          primer_region_seq = seq[pid_length, cdna_bio_primer_size]
          if primer_region_seq =~ cdna_primer_ref
            new_name = remove_tag name
            r2_passed_seq[new_name] = seq
          end
        end

        r2_raw.keys.map do |name|
          proc_filter.call name
        end

        return { r2_passed_seq: r2_passed_seq, pid_length: pid_length, reverse_starting_number: reverse_starting_number }
      end # end of filter_r2



      # puts error message in the log file handler, and abort with the same infor

      def log_and_abort(log, infor)
        log.puts Time.now.to_s + "\t" + infor
        log.close
        abort infor.red.bold
      end

      # lower detection sensitivity for minority mutations given the number of TCS, calculated based on binomial distribution.
      # R required.
      # @param tcs_number [Integer] number of TCS
      # @return [Float] lower detection limit
      # @example calculate lower detection limit
      #   ViralSeq::TcsCore.detection_limit(100)
      #   => 0.0362

      def detection_limit(tcs_number)
        if ViralSeq::DETECT_SEN[tcs_number]
          return ViralSeq::DETECT_SEN[tcs_number]
        else
          dl =  `Rscript -e "library(dplyr); ifelse(#{tcs_number} > 2, (binom.test(0,#{tcs_number})['conf.int'] %>% unlist %>% unname)[2] %>% round(4) %>% cat, 0)" 2>/dev/null`
          dl.to_f
        end
      end

      private

      def unzip_r(indir, f)
        r_file = File.join(indir, f)
        if f =~ /.gz/
          `gzip -d #{r_file}`
          new_f = f.sub ".gz", ""
          r_file = File.join(indir, new_f)
        end
        return r_file
      end

      def parser_file_name(file_name)
        t = file_name.split(".")[0].split("_")
        if t.size == 1
          libname = "lib"
          tag = [ t[0].upcase ]
        else
          libname = t[0]
          tag = t[1..-1].map(&:upcase)
        end
        return {libname: libname, tag: tag}
      end

      def general_filter(seq)
        return false unless seq
        if seq.size < ($platform_sequencing_length - 10)
          return false
        elsif seq[1..-2] =~ /N/ # sequences with ambiguities except the 1st and last position removed
          return false
        elsif seq =~ /A{11}/ # a string of poly-A indicates adaptor sequence
          return false
        elsif seq =~ /T{11}/ # a string of poly-T indicates adaptor sequence
          return false
        else
          return true
        end
      end

      # remove region info tags from the raw MiSeq sequences.
      def remove_tag(seq_name)
        if seq_name =~ /\s/
          new_tag = $`
        else
          new_tag = seq_name[0..-3]
        end
      end

    end # end of class << self

  end # end of TcsCore module

end # end of main module
