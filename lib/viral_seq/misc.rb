# viral_seq/misc.rb

# miscellaneous methods
# including
#   Hash#copyhash
#   Hash#difference
#   Hash#uniq_hash
#   ViralSeq::tail

class Hash

  # Hash#copyhash
  # copy a hash
  # different from "="
  #   # example
  #   h1 = {1=>'a'}
  #   h2 = h1
  #   h3 = h1.copyhash
  #   h1.object_id == h2.object_id
  #   => true
  #   h1.object_id == h3.object_id
  #   => false

  def copyhash
    h = Hash.new
    self.each do |pair|
      h.store(pair[0], pair[1])
    end
    return h
  end

  # subtract one hash (h2) from the other (h1) if the keys are identical
  # example:
  # h1 = {"Cat" => 100, "Dog" => 5, "Bird" => 2, "Snake" => 10}
  # h2 = {"Cat" => 100, "Dog" => 5, "Bison" => 30}
  # h1.difference(h2) = {"Bird" => 2, "Snake" => 10}

  def difference(other)
    reject do |k,_v|
      other.has_key? k
    end
  end

  # input hash A, return hash B with the unique values of hash A as keys,
  # and the keys of the unique values of hash A as values of hash B
  #   # example
  #   hash = {1=>"A", 2=>"A", 3=>"C", 4=>"C", 5=>"T"}
  #   p hash.uniq_hash
  #   => {"A"=>[1, 2], "C"=>[3, 4], "T"=>[5]}

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

# Tail function for file as 'tail' in bash.
def ViralSeq.tail(path, n)
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
  file.read
end
