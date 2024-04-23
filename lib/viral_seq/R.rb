module ViralSeq
  module R

    # check if R is installed. if R is installed, return the version number of R.

    def self.check_R
      begin
        r_version = `R --version`.split("\n")[0]
      rescue Errno::ENOENT
        abort '"R" is not installed. Install R at https://www.r-project.org/' +
              "\n`tcs_sdrm` pipeline aborted."
      end
    end # end check_R

    # check if required R packages is installed.
    def self.check_R_packages
      if system "Rscript #{File.join("lib", "viral_seq", "util", "check_env.r")}"
        return 0
      else
        raise "Non-zeor exit code. Error happens when checking required R packages."
      end
    end # end check_R_packages.

    # read sdrm rscript as a string.

    def self.get_sdrm_rscript
      File.read(File.join("lib", "viral_seq", "util", "sdrm_r.r"))
    end

  end
end
