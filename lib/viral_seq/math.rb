# lib/math.rb

# math and statistic functions


module ViralSeq

  # count elements in a array, return a hash of {:element1 => number1, :element2 => number2, ...}
  # =Usage
  #   array = %w{cat dog monkey cat cat cat monkey}
  #   ViralSeq.count(array)
  #   => {"cat"=>4, "dog"=>1, "monkey"=>2}

  def self.count(array)
    hash = Hash.new(0)
    array.each do |element|
      hash[element] +=1
    end
    return hash
  end

  # count elements in a array, return a hash of {:element1 => frequency1, :element2 => frequency2, ...}
  # default decimal as 2
  # =Usage
  #   array = %w{cat dog monkey cat cat cat monkey}
  #   ViralSeq.count_percentage(array)
  #   => {"cat"=>0.57, "dog"=>0.14, "monkey"=>0.29}

  def self.count_percentage(array,decimal = 2)
    hash1 = Hash.new(0)
    array.each do |element|
      hash1[element] += 1
    end
    total_elements = array.size
    hash2 = Hash.new(0)
    hash1.each do |key,value|
      hash2[key] = (value/total_elements.to_f).round(decimal)
    end
    return hash2
  end

  # poisson distribution. input lambda and maximum k, return a hash with keys as k
  # default k value is 5, meaning calculate up to 5 events.
  #
  # Poisson Distribution (https://en.wikipedia.org/wiki/Poisson_distribution)
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
  #
  # =USAGE
  #   # We assume the mutaiton rate is 0.005 (event rate λ),
  #   # we would like to calculate the probablity of 3 mutations on one sequence
  #   prob_hash = ViralSeq::poisson_distribution(0.005)
  #   => {0=>0.9950124791926823, 1=>0.004975062395963412, 2=>1.243765598990853e-05, 3=>2.072942664984755e-08, 4=>2.5911783312309436e-11, 5=>2.591178331230944e-14}
  #   prob_hash[3]
  #   => 2.072942664984755e-08

  def self.poisson_distribution(rate,k = 5)
    out_hash = {}
    (0..k).each do |n|
      p = (rate**n * Math::E**(-rate))/!n
      out_hash[n] = p
    end
    return out_hash
  end


  # require R pre-installed
  # calculate binomial 95% confidence intervals by R. refer to R function binom.test
  # input number x and n, return an array as [lower_interval, upper_interval]
  #
  # =USAGE
  #   # mutation M184V found in 3 out of 923 sequences, the 95% confidence interval is
  #   ViralSeq.r_binom_CI(3, 923)
  #   => [0.02223, 0.19234]
  #
  def self.r_binom_CI(x= 0, n= 0)
    r_output = `Rscript -e 'binom.test(#{x},#{n})$conf.int[1];binom.test(#{x},#{n})$conf.int[2]'`
    lines = r_output.split "\n"
    low = lines[0].chomp[4..-1].to_f
    high = lines[1].chomp[4..-1].to_f
    return [low.round(5), high.round(5)]
  end

end

# statistic methods
# :median :sum :mean :sample_variance :stdev :upper_quartile :lower_quartile
# =USAGE
#   array = [1,2,3,4,5,6,7,8,9,10]
#   array.median
#   => 5.5
#   array.sum
#   => 55
#   array.mean
#   => 5.5
#   array.sample_variance
#   => 9.166666666666666
#   array.stdev
#   => 3.0276503540974917
#   array.upper_quartile
#   => 7.5
#   array.lower_quartile
#   => 3.5

module Enumerable
  def median
    len = self.length
    sorted = self.sort
    len % 2 == 1 ? sorted[len/2] : (sorted[len/2 - 1] + sorted[len/2]).to_f / 2
  end

  def sum
     self.inject(0){|accum, i| accum + i }
  end

  def mean
    self.sum/self.length.to_f
  end

  def sample_variance
    m = self.mean
    sum = self.inject(0){|accum, i| accum + (i-m)**2 }
    sum/(self.length - 1).to_f
  end

  def stdev
    return Math.sqrt(self.sample_variance)
  end

  def upper_quartile
    return nil if self.empty?
    sorted_array = self.sort
    u = (0.25*(3*sorted_array.length))
    if (u-u.truncate).is_a?(Integer)
      return sorted_array[(u-u.truncate)-1]
    else
      sample = sorted_array[u.truncate.abs-1]
      sample1 = sorted_array[(u.truncate.abs)]
      return sample+((sample1-sample)*(u-u.truncate))
    end
  end

  def lower_quartile
    return nil if self.empty?
    sorted_array = self.sort
    u = 0.25*sorted_array.length + 1
    if (u-u.truncate).is_a?(Integer)
      return sorted_array[(u-u.truncate)-1]
    else
      sample = sorted_array[u.truncate.abs-1]
      sample1 = sorted_array[(u.truncate.abs)]
      return sample+((sample1-sample)*(u-u.truncate))
    end
  end
end

class Integer
  def !
    if self == 0
      return 1
    else
      (1..self).inject(:*)
    end
  end
end


# Fisher's Exact Test Function Library
#
# Based on JavaScript version created by: Oyvind Langsrud
# Ported to Ruby by Bryan Donovan

module Rubystats
  class FishersExactTest

    def initialize
      @sn11    = 0.0
      @sn1_    = 0.0
      @sn_1    = 0.0
      @sn      = 0.0
      @sprob   = 0.0

      @sleft   = 0.0
      @sright  = 0.0
      @sless   = 0.0
      @slarg   = 0.0

      @left    = 0.0
      @right   = 0.0
      @twotail = 0.0
    end

    # Reference: "Lanczos, C. 'A precision approximation
    # of the gamma function', J. SIAM Numer. Anal., B, 1, 86-96, 1964."
    # Translation of  Alan Miller's FORTRAN-implementation
    # See http://lib.stat.cmu.edu/apstat/245
    def lngamm(z)
      x = 0
      x += 0.0000001659470187408462 / (z+7)
      x += 0.000009934937113930748  / (z+6)
      x -= 0.1385710331296526       / (z+5)
      x += 12.50734324009056        / (z+4)
      x -= 176.6150291498386        / (z+3)
      x += 771.3234287757674        / (z+2)
      x -= 1259.139216722289        / (z+1)
      x += 676.5203681218835        / (z)
      x += 0.9999999999995183

      return(Math.log(x)-5.58106146679532777-z+(z-0.5) * Math.log(z+6.5))
    end

    def lnfact(n)
      if n <= 1
        return 0
      else
        return lngamm(n+1)
      end
    end

    def lnbico(n,k)
      return lnfact(n) - lnfact(k) - lnfact(n-k)
    end

    def hyper_323(n11, n1_, n_1, n)
      return Math.exp(lnbico(n1_, n11) + lnbico(n-n1_, n_1-n11) - lnbico(n, n_1))
    end

    def hyper(n11)
      return hyper0(n11, 0, 0, 0)
    end

    def hyper0(n11i,n1_i,n_1i,ni)
      if n1_i == 0 and n_1i ==0 and ni == 0
        unless n11i % 10 == 0
          if n11i == @sn11+1
            @sprob *= ((@sn1_ - @sn11)/(n11i.to_f))*((@sn_1 - @sn11)/(n11i.to_f + @sn - @sn1_ - @sn_1))
            @sn11 = n11i
            return @sprob
          end
          if n11i == @sn11-1
            @sprob *= ((@sn11)/(@sn1_-n11i.to_f))*((@sn11+@sn-@sn1_-@sn_1)/(@sn_1-n11i.to_f))
            @sn11 = n11i
            return @sprob
          end
        end
        @sn11 = n11i
      else
        @sn11 = n11i
        @sn1_ = n1_i
        @sn_1 = n_1i
        @sn   = ni
      end
      @sprob = hyper_323(@sn11,@sn1_,@sn_1,@sn)
      return @sprob
    end

    def exact(n11,n1_,n_1,n)

      p = i = j = prob = 0.0

      max = n1_
      max = n_1 if n_1 < max
      min = n1_ + n_1 - n
      min = 0 if min < 0

      if min == max
        @sless  = 1
        @sright = 1
        @sleft  = 1
        @slarg  = 1
        return 1
      end

      prob = hyper0(n11,n1_,n_1,n)
      @sleft = 0

      p = hyper(min)
      i = min + 1
      while p < (0.99999999 * prob)
        @sleft += p
        p = hyper(i)
        i += 1
      end

      i -= 1

      if p < (1.00000001*prob)
        @sleft += p
      else
        i -= 1
      end

      @sright = 0

      p = hyper(max)
      j = max - 1
      while p < (0.99999999 * prob)
        @sright += p
        p = hyper(j)
        j -= 1
      end
      j += 1

      if p < (1.00000001*prob)
        @sright += p
      else
        j += 1
      end

      if (i - n11).abs < (j - n11).abs
        @sless = @sleft
        @slarg = 1 - @sleft + prob
      else
        @sless = 1 - @sright + prob
        @slarg = @sright
      end
      return prob
    end

    def calculate(n11_,n12_,n21_,n22_)
      n11_ *= -1 if n11_ < 0
      n12_ *= -1 if n12_ < 0
      n21_ *= -1 if n21_ < 0
      n22_ *= -1 if n22_ < 0
      n1_     = n11_ + n12_
      n_1     = n11_ + n21_
      n       = n11_ + n12_ + n21_ + n22_
      exact(n11_,n1_,n_1,n)
      left    = @sless
      right   = @slarg
      twotail = @sleft + @sright
      twotail = 1 if twotail > 1
      values_hash = { :left =>left, :right =>right, :twotail =>twotail }
      return values_hash
    end
  end
end
