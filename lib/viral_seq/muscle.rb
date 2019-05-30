# viral_seq/muscle.rb

# check if path_to_muscle is correct

module ViralSeq
  def self.check_muscle(path_to_muscle)
    begin
      `#{path_to_muscle} -version`
    rescue Errno::ENOENT
      puts "
            Error: MUSCLE is not found!!\n
            MUSLCE can be download at http://www.drive5.com/muscle\n
            Add MUSCLE excutable path to $PATH using\n
            $  export PATH=$PATH:/path/to/muscle\n
            or\n
            provide path_to_MUSCLE in the ViralSeq::sequence_locator arguments\n
            "
      return nil
    end
  end
end
