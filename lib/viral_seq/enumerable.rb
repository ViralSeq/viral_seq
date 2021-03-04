# additional statistic/math functions to Module::Enumerable
# @example median number
#   array = [1,2,3,4,5,6,7,8,9,10]
#   array.median
#   => 5.5
# @example average number (mean)
#   array = [1,2,3,4,5,6,7,8,9,10]
#   array.mean
#   => 5.5
# @example sample variance
#   array = [1,2,3,4,5,6,7,8,9,10]
#   array.sample_variance
#   => 9.166666666666666
# @example standard deviation
#   array = [1,2,3,4,5,6,7,8,9,10]
#   array.stdev
#   => 3.0276503540974917
# @example upper quartile
#   array = [1,2,3,4,5,6,7,8,9,10]
#   array.upper_quartile
#   => 7.5
# @example lower_quartile
#   array = [1,2,3,4,5,6,7,8,9,10]
#   array.lower_quartile
#   => 3.5
# @example count frequency of elements in an array
#   array = %w{cat dog monkey cat cat cat monkey}
#   array.count_freq
#   => {"cat"=>4, "dog"=>1, "monkey"=>2}
# @example count frequency as percentage of elements in an array
#   array = %w{cat dog monkey cat cat cat monkey}
#   array.count_freq2
#   => {"cat"=>0.57, "dog"=>0.14, "monkey"=>0.29}
module Enumerable

  # generate median number
  # @return [Numeric] median number
  def median
    len = self.length
    sorted = self.sort
    len % 2 == 1 ? sorted[len/2] : (sorted[len/2 - 1] + sorted[len/2]).to_f / 2
  end

  # generate mean number
  # @return [Float] mean value
  def mean
    self.sum/self.length.to_f
  end

  # generate sample variance
  # @return [Float] sample variance
  def sample_variance
    m = self.mean
    sum = self.inject(0){|accum, i| accum + (i-m)**2 }
    sum/(self.length - 1).to_f
  end

  # generate standard deviation
  # @return [Float] standard deviation
  def stdev
    return Math.sqrt(self.sample_variance)
  end

  # generate upper quartile value
  # @return [Numeric] upper quartile value
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

  # generate lower quartile value
  # @return [Numeric] lower quartile value
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

  # tabulate elements and frequencies of an Enumerable
  # return [Hash] return a hash of :element => :freq_count

  def count_freq
    hash = Hash.new(0)
    self.each do |element|
      hash[element] +=1
    end
    return hash
  end

  # tabulate elements and frequencies (as percentage) of an Enumerable {
  # @param decimal [Integer] decimals of frequency
  # return [Hash] return a hash of :element => :percentage

  def count_freq2(decimal = 2)
    hash1 = Hash.new(0)
    self.each do |element|
      hash1[element] += 1
    end
    total_elements = self.size
    hash2 = Hash.new(0)
    hash1.each do |key,value|
      hash2[key] = (value/total_elements.to_f).round(decimal)
    end
    return hash2
  end

end
