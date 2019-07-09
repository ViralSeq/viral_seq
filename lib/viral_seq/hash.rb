# addition methods for Class::Hash required for ViralSeq

class Hash

  # subtract one hash (h2) from the other (h1) if the keys are identical
  # @param other_hash [Hash] the hash that needs to substracted from the hash before the method
  # @return [Hash] hash after substraction
  # @example substract h2 from h1 if the keys match
  #   h1 = {"Cat" => 100, "Dog" => 5, "Bird" => 2, "Snake" => 10}
  #   h2 = {"Cat" => 100, "Dog" => 5, "Bison" => 30}
  #   h1.difference(h2)
  #   => {"Bird" => 2, "Snake" => 10}

  def difference(other_hash)
    reject do |k,_v|
      other_hash.has_key? k
    end
  end

  # return a new hash with the unique values of input hash as keys,
  # and the keys of the unique values of input hash in an array as values of the new hash
  # @return [Hash] a new hash of :uniq_value_of_orginial_hash => :array_of_keys
  # @example
  #   hash = {1=>"A", 2=>"A", 3=>"C", 4=>"C", 5=>"T"}
  #   hash.uniq_hash
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
