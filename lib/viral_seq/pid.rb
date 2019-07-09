
module ViralSeq

  module PID

    # generate all Primer ID combinations given the length of Primer ID
    # @param l [Integer] the length of the Primer ID.
    # @example generate a pool of Primer IDs with length of 10
    #   primer_id_pool = ViralSeq::PID.generate_pool(10) # 10 is the length of Primer ID
    #   puts primer_id_pool.size  #should be 4^10
    #   => 1048576

    def self.generate_pool(l=8)
      nt = ['A','T','C','G']
      pid_pool = ['A','T','C','G']
      (l-1).times do
        pid_pool = pid_pool.product(nt)
        pid_pool.collect! do |v|
          v.join("")
        end
      end
      return pid_pool
    end # end of .generate_primer_id_pool

  end # end of Pid
end # end of ViralSeq
