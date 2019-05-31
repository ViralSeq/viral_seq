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


# Tail function for file as 'tail' in bash.
def tail(path, n)
  file = File.open(path, "r")
  buffer_s = 512
  line_count = 0
  file.seek(0, IO::SEEK_END)

  offset = file.pos # we start at the end

  while line_count <= n && offset > 0
    to_read = if (offset - buffer_s) < 0
                offset
              else
                buffer_s
              end

    file.seek(offset-to_read)
    data = file.read(to_read)

    data.reverse.each_char do |c|
      if line_count > n
        offset += 1
        break
      end
      offset -= 1
      if c == "\n"
        line_count += 1
      end
    end
  end

  file.seek(offset)
  data = file.read
end


# subtract one hash (h2) from the other (h1) if the keys are identical
# example:
# h1 = {"Cat" => 100, "Dog" => 5, "Bird" => 2, "Snake" => 10}
# h2 = {"Cat" => 100, "Dog" => 5, "Bison" => 30}
# h1.difference(h2) = {"Bird" => 2, "Snake" => 10}
class Hash
  def difference(other)
    reject do |k,_v|
      other.has_key? k
    end
  end
end


# input hash A, return hash B with the unique values of hash A as keys,
# and the keys of the unique values of hash A as values of hash B
#   # example
#   hash = {1=>"A", 2=>"A", 3=>"C", 4=>"C", 5=>"T"}
#   p hash.uniq_hash
#   => {"A"=>[1, 2], "C"=>[3, 4], "T"=>[5]}

class Hash
  def uniq_hash
    uniq_values = self.values.uniq
    out_hash = {}
    uniq_values.each do |uniq_va|
      self.each do |k,v|
        if v == uniq_va
          if out_hash[uniq_va]
            out_hash[uniq_va] << k
          else
            out_hash[uniq_va] = []
            out_hash[uniq_va] << k
          end
        end
      end
    end
    return out_hash
  end
end
