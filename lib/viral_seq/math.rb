
module ViralSeq

  # math functions reqruied for ViralSeq

  module Math

    # Generate values from the standard normal distribution with given mean and standard deviation
    # @see http://en.wikipedia.org/wiki/Box-Muller_transform Wikipedia explanation

    class RandomGaussian

      # generate RandomGaussian instance with given mean and standard deviation
      # @param mean [Float] mean value.
      # @param sd [Float] standard deviation value.

      def initialize(mean = 0.0, sd = 1.0, rng = lambda { Kernel.rand })
        @mean, @sd, @rng = mean, sd, rng
        @compute_next_pair = false
      end

      # generate a random number that falls in the pre-defined gaussian distribution
      # @return [Float]
      # @example generate 10 random number that falls in the a gaussian distribution with mean at 0 and standard deviation at 1.0
      #   a = RandomGaussian.new
      #   numbers = []
      #   10.times {numbers << a.rand.round(5)}
      #   numbers
      #   => [-1.83457, 1.24439, -0.30109, 0.13977, 0.61556, 1.3548, 1.72878, 2.46171, 0.97031, -0.29496]

      def rand
        if (@compute_next_pair = !@compute_next_pair)
          theta = 2 * ::Math::PI * @rng.call
          scale = @sd * ::Math.sqrt(-2 * Math.log(1 - @rng.call))
          @g1 = @mean + scale * ::Math.sin(theta)
          @g0 = @mean + scale * ::Math.cos(theta)
        else
          @g1
        end
      end

    end

    # class for poisson distribution.
    #   An event can occur 0, 1, 2, … times in an interval.
    #   The average number of events in an interval is designated λ (lambda).
    #   λ is the event rate, also called the rate parameter.
    #   The probability of observing k events in an interval is given by the equation
    #
    #   P(k events in interval) = e^(-λ) * λ^k / k!
    #
    #   λ is the average number of events per interval
    #   e is the number 2.71828... (Euler's number) the base of the natural logarithms
    #   k takes values 0, 1, 2, …
    #   k! = k × (k − 1) × (k − 2) × … × 2 × 1 is the factorial of k.
    # @see https://en.wikipedia.org/wiki/Poisson_distribution Poisson Distribution (Wikipedia).
    # @example given the mutation rate at 0.01 and sequence length of 1000 bp,
    # calculate the probablity of 3 mutations on one sequence
    #   new_poisson_dist = ViralSeq::Math::PoissonDist.new(0.01)
    #   prob_hash = new_poisson_dist.poisson_hash
    #   1000 * prob_hash[3].round(5)
    #   => 0.00017
    class PoissonDist
      # initialize with given event rate λ, default events upper limit set to 5
      def initialize(rate,k = 5)
        @rate = rate
        @k = k
        @poisson_hash = {}
        (0..k).each do |n|
          p = (rate**n * ::Math::E**(-rate))/!n
          @poisson_hash[n] = p
        end
      end

      # @return [Float] event rate λ
      attr_accessor :rate
      # @return [Integer] maxinum events number shows in @poisson_hash
      attr_accessor :k
      # @return [Hash] probablity hash of :event_number => :probablity
      attr_reader :poisson_hash
    end # end of PoissonDist

    # Use R to calculate binomial 95% confidence intervals. Require R function binom.test.
    # @example mutation M184V found in 3 out of 923 sequences, calculate 95% confidence interval
    #   freq = ViralSeq::Math::BinomCI.new(3,923)
    #   freq.mean.round(5)
    #   => 0.00325
    #   freq.lower.round(5)
    #   => 0.00067
    #   freq.upper.round(5)
    #   => 0.00947

    class BinomCI
      # initialize with numerator @n1 and denominator @n2 as Integer
      def initialize(n1, n2)
        @n1 = n1
        @n2 = n2
        @mean = n1/n2.to_f
        r_output = `Rscript -e 'binom.test(#{n1},#{n2})$conf.int[1];binom.test(#{n1},#{n2})$conf.int[2]'`
        lines = r_output.split "\n"
        @lower = lines[0].chomp[4..-1].to_f
        @upper = lines[1].chomp[4..-1].to_f
      end

      # @return [Integer] number of observations
      attr_accessor :n1
      # @return [Integer] total numbers
      attr_accessor :n2
      # @return [Float] mean
      attr_reader :mean
      # @return [Float] lower limit of 95% CI
      attr_reader :lower
      # @return [Float] upper limit of 95% CI
      attr_reader :upper

    end # end of BinomCI


    # A function to calcuate cut-off for offspring primer IDs.
    # @see https://www.ncbi.nlm.nih.gov/pubmed/26041299 reference at Zhou et al. JVI 2016.
    # @param m [Integer] PID abundance
    # @param error_rate [Float] estimated platform error rate, the model supports error rate from 0.003 to 0.03.
    # @return [Integer] an abundance cut-off (Integer) for offspring Primer IDs.

    def self.calculate_pid_cut_off(m, error_rate = 0.02)
      if m <= 10
        return 2
      end
      n = 0
      case error_rate
      when 0...0.0075
        n = -9.59*10**-27*m**6 + 3.27*10**-21*m**5 - 3.05*10**-16*m**4 + 1.2*10**-11*m**3 - 2.19*10**-7*m**2 + 0.004044*m + 2.273
      when 0.0075...0.015
        n = 1.09*10**-26*m**6 + 7.82*10**-22*m**5 - 1.93*10**-16*m**4 + 1.01*10**-11*m**3 - 2.31*10**-7*m**2 + 0.00645*m + 2.872
      when 0.015..0.03
        if m <= 8500
          n = -1.24*10**-21*m**6 + 3.53*10**-17*m**5 - 3.90*10**-13*m**4 + 2.12*10**-9*m**3 - 6.06*10**-6*m**2 + 1.80*10**-2*m + 3.15
        else
          n = 0.0079 * m + 9.4869
        end
      else
        raise ArgumentError.new('Error_rate has be between 0 to 0.03')
      end
      n = n.round
      n = 2 if n < 3
      return n
    end # end of .calculate_pid_cut_off
  end # end of Math
end # end of ViralSeq
