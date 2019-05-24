# viral_seq/misc.rb

# miscellaneous methods


# copy a hash

class Hash
  def copyhash
    h = Hash.new
    self.each do |pair|
      h.store(pair[0], pair[1])
    end
    return h
  end
end

# generate values from the standard normal distribution with given mean and standard deviation
# See http://en.wikipedia.org/wiki/Box-Muller_transform
class RandomGaussian
  def initialize(mean = 0.0, sd = 1.0, rng = lambda { Kernel.rand })
    @mean, @sd, @rng = mean, sd, rng
    @compute_next_pair = false
  end

  def rand
    if (@compute_next_pair = !@compute_next_pair)
      theta = 2 * Math::PI * @rng.call
      scale = @sd * Math.sqrt(-2 * Math.log(1 - @rng.call))
      @g1 = @mean + scale * Math.sin(theta)
      @g0 = @mean + scale * Math.cos(theta)
    else
      @g1
    end
  end
end
