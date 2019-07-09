# additional functions for Class::Integer

class Integer
  # factorial method for an Integer
  # @return [Integer] factorial for given Integer
  # @example factorial for 5
  #   !5
  #   => 120
  def !
    if self == 0
      return 1
    else
      (1..self).inject(:*)
    end
  end
end
