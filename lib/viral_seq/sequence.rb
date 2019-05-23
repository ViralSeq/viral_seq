
$amino_acid_list = ["A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y", "*"]

class Sequence
  def initialize (name = ">sequence",dna_sequence ="")
    @name = name
    @dna_sequence = dna_sequence.upcase
    @aa_sequence = ""
    @aa_array = []
  end
  attr_accessor :name, :dna_sequence, :aa_sequence, :aa_array
  def rev_complement
    @dna_sequence.reverse.upcase.tr('ATCG','TAGC')
  end
  def rev_complement!
    @dna_sequence = @dna_sequence.reverse.upcase.tr('ATCG','TAGC')
  end
  def get_aa_sequence(initial_position = 0)
    @aa_sequence = ""
    require_sequence = @dna_sequence[initial_position..-1]
    base_array = []
    require_sequence.each_char {|base| base_array << base}
    while (base_array.length>=3) do
      base_3= ""
      3.times{base_3 += base_array.shift}
      @aa_sequence<< amino_acid(base_3)
    end
  end

  #get amino acid calls, return a array.keep ambiguity calls.
  def get_aa_array(initial_position = 0)
    @aa_array = []
    require_sequence = @dna_sequence[initial_position..-1].tr('-','N')
    base_array = []
    require_sequence.each_char {|base| base_array << base}
    while (base_array.length>=3) do
      base_3= ""
      3.times{base_3 += base_array.shift}
      @aa_array<< amino_acid_2(base_3)
    end
  end
  def dna_length
    @dna_sequence.length
  end
  def aa_length
    @aa_sequence.length
  end
end

def amino_acid (bases)
  case bases
  when /^TT[TCY]$/
    return "F"
  when /^TT[AGR]$/
    return "L"
  when /^CT.$/
    return "L"
  when /^AT[TCAHYWM]$/
    return "I"
  when "ATG"
    return "M"
  when /^GT.$/
    return "V"
  when /^TC.$/
    return "S"
  when /^CC.$/
    return "P"
  when /^AC.$/
    return "T"
  when /^GC.$/
    return "A"
  when /^TA[TCY]$/
    return "Y"
  when /^TA[AGR]$/
    return "*"
  when /^T[GR]A$/
    return "*"
  when /^CA[TCY]$/
    return "H"
  when /^CA[AGR]$/
    return "Q"
  when /^AA[TCY]$/
    return "N"
  when /^AA[AGR]$/
    return "K"
  when /^GA[TCY]$/
    return "D"
  when /^GA[AGR]$/
    return "E"
  when /^TG[TCY]$/
    return "C"
  when "TGG"
    return "W"
  when /^CG.$/
    return "R"
  when /^AG[TCY]$/
    return "S"
  when /^[AM]G[AGR]$/
    return "R"
  when /^GG.$/
    return "G"
  when /^[ATW][CGS][CTY]$/
    return "S"
  when /^[TCY]T[AGR]$/
    return "L"
  else
    return "#"
  end
end

#keep ambiguities

def amino_acid_2 (bases)
  bases_to_aa = []
  aa_list = []
  base1 = to_list(bases[0])
  base2 = to_list(bases[1])
  base3 = to_list(bases[2])
  l1 = base1.size - 1
  l2 = base2.size - 1
  l3 = base3.size - 1
  (0..l1).each do |n1|
    b1 = base1[n1]
    (0..l2).each do |n2|
      b2 = base2[n2]
      (0..l3).each do |n3|
        b3 = base3[n3]
        bases_all = b1 + b2 + b3
        bases_to_aa << bases_all
      end
    end
  end

  bases_to_aa.each do |bases|
  case bases
  when /^TT[TCY]$/
    aa =  "F"
  when /^TT[AGR]$/
    aa =  "L"
  when /^CT.$/
    aa =  "L"
  when /^AT[TCAHYWM]$/
    aa =  "I"
  when "ATG"
    aa =  "M"
  when /^GT.$/
    aa =  "V"
  when /^TC.$/
    aa =  "S"
  when /^CC.$/
    aa =  "P"
  when /^AC.$/
    aa =  "T"
  when /^GC.$/
    aa =  "A"
  when /^TA[TCY]$/
    aa =  "Y"
  when /^TA[AGR]$/
    aa =  "*"
  when /^T[GR]A$/
    aa =  "*"
  when /^CA[TCY]$/
    aa =  "H"
  when /^CA[AGR]$/
    aa =  "Q"
  when /^AA[TCY]$/
    aa =  "N"
  when /^AA[AGR]$/
    aa =  "K"
  when /^GA[TCY]$/
    aa =  "D"
  when /^GA[AGR]$/
    aa =  "E"
  when /^TG[TCY]$/
    aa =  "C"
  when "TGG"
    aa =  "W"
  when /^CG.$/
    aa =  "R"
  when /^AG[TCY]$/
    aa =  "S"
  when /^[AM]G[AGR]$/
    aa =  "R"
  when /^GG.$/
    aa =  "G"
  when /^[ATW][CGS][CTY]$/
    aa =  "S"
  when /^[TCY]T[AGR]$/
    aa =  "L"
  else
    aa =  "-"
  end
  aa_list << aa
end
  aa_out = aa_list.uniq.join('/')
  return aa_out
end

def to_list(base = "")
  list = []
  case base
  when /[A|T|C|G]/
    list << base
  when "W"
    list = ['A','T']
  when "S"
    list = ['C','G']
  when "M"
    list = ['A','C']
  when 'K'
    list = ['G','C']
  when 'R'
    list = ['A','G']
  when 'Y'
    list = ['C','T']
  when 'B'
    list = ['C','G','T']
  when 'D'
    list = ['A','G','T']
  when 'H'
    list = ['A','C','T']
  when 'V'
    list = ['A','C','G']
  when 'N'
    list = ['A','T','C','G']
  end
  return list
end


def array_to_hash(array)
  count = 0
  hash = Hash.new
  (array.length / 2).times do
    hash[array[count]] = array[count+1]
    count += 2
  end
  return hash
end


def count(array)
  hash = Hash.new(0)
  array.each do |element|
    hash[element] +=1
  end
  return hash
end

#count percentage of each element#
def count_percentage(array,decimal = 2)
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

def clustal_score(ref_sequence,test_sequence,open = 15, ext = 6.66)
  score = 0
  gapopen = open
  gapext = ext
  temp_dir = File.dirname($0)

  temp_file_in = temp_dir + "/temp_sequence"

  ref_name = ref_sequence.name
  ref_seq = ref_sequence.dna_sequence
  test_name = test_sequence.name
  test_seq = test_sequence.dna_sequence

  File.open(temp_file_in,'w') do |line|
    line.print "#{ref_name}\n#{ref_seq}\n#{test_name}\n#{test_seq}\n"
  end

  temp_file_out = temp_dir + "/temp_out"
  temp_screen_out = temp_dir + "/temp_screen"
  print `/applications/clustalw2 -infile=#{temp_file_in} -case=upper -outorder=input -output=gde -outfile=#{temp_file_out} >#{temp_screen_out} -gapopen=#{gapopen} -gapext=#{gapext}`
  File.open(temp_screen_out,"r") do |file|
    file.readlines.each do |line|
      if ((line =~ /^Sequences/)&&(line =~ /Score/))
        score = line.match(/\w+\s\(.+\)\s\w+\.\s\w+\:\s+(\d+)/)[1]
      end
    end
  end
  File.unlink temp_file_in
  File.unlink temp_file_out
  File.unlink temp_screen_out
  Dir.chdir(temp_dir) do
    Dir.glob("*.dnd") do |dnd|
      File.unlink(dnd)
    end
  end
  return (score.to_i)
end

def clustal_seq(ref_sequence,test_sequence,open = 15, ext = 6.66)
  score = 0
  gapopen = open
  gapext = ext
  temp_dir = File.dirname($0)

  temp_file_in = temp_dir + "/temp_sequence"

  ref_name = ref_sequence.name
  ref_seq = ref_sequence.dna_sequence
  test_name = test_sequence.name
  test_seq = test_sequence.dna_sequence

  File.open(temp_file_in,'w') do |line|
    line.print "#{ref_name}\n#{ref_seq}\n#{test_name}\n#{test_seq}\n"
  end

  temp_file_out = temp_dir + "/temp_out"
  temp_screen_out = temp_dir + "/temp_screen"
  print `/applications/clustalw2 -infile=#{temp_file_in} -case=upper -outorder=input -output=gde -outfile=#{temp_file_out} >#{temp_screen_out} -gapopen=#{gapopen} -gapext=#{gapext}`
  h = {}
  File.open(temp_file_out,"r") do |file|
    n = 0
    file.readlines.each do |line|
      if line =~ /^[\#|\%]/
        n += 1
        h[n] = ""
      else
        h[n] += line.chomp
      end
    end
  end
  ref_aligned = h[1]
  test_aligned = h[2]
  length = ref_aligned.length
  test_trimmed = ""
  (0..(length - 1)).each do |n|
    base_ref = ref_aligned[n]
    base_test = test_aligned[n]
    unless base_ref == "-"
      test_trimmed << base_test
    end
  end

  File.unlink temp_file_in
  File.unlink temp_file_out
  File.unlink temp_screen_out
  Dir.chdir(temp_dir) do
    Dir.glob("*.dnd") do |dnd|
      File.unlink(dnd)
    end
  end
  return test_trimmed
end


def substitution(base)
  bases = ['A','C','T','G']
  bases.delete(base)
  new_base = bases.sample
end

def addition(base)
  bases = ['A','C','T','G']
  new_base = base + bases.sample
end

#creat consensus with alignment. Need sequence hash as input. gap_treatment option: 1: treat gap as normal base call, 2: only return a consensus gap call when gap over 75% #
def clustal_consensus_multi(seq_hash,open = 15, ext = 6.66, gap_treatment = 1)
  gapopen = open
  gapext = ext
  temp_dir = File.dirname($0)
  temp_file_in = temp_dir + "/temp_sequence"
  f = File.open(temp_file_in,'w')
  f.puts seq_hash.flatten
  f.close

  temp_file_out = temp_dir + "/temp_out"
  temp_screen_out = temp_dir + "/temp_screen"
  print `/applications/clustalw2 -infile=#{temp_file_in} -case=upper -outorder=input -output=gde -outfile=#{temp_file_out} >#{temp_screen_out} -gapopen=#{gapopen} -gapext=#{gapext}`
  h = {}
  File.open(temp_file_out,"r") do |file|
    n = 0
    file.readlines.each do |line|
      if line =~ /^\#/
        n += 1
        h[n] = ""
      else
        h[n] += line.chomp
      end
    end
  end
  length = h[1].size
  consensus_bases = []
  (0..(length-1)).each do |n|
    bases = []
    h.values.each do |seq|
      bases << seq[n]
    end
    if gap_treatment == 1
      consensus_bases << creat_consensus_base_non_gap(bases)
    else
      consensus_bases << creat_consensus_base_gap(bases)
    end
  end
  File.unlink temp_file_in
  File.unlink temp_file_out
  File.unlink temp_screen_out
  Dir.chdir(temp_dir) do
    Dir.glob("*.dnd") do |dnd|
      File.unlink(dnd)
    end
  end
  consensus_seq = consensus_bases.join('')
end

#for 454 sequencing which have problem undercalling homopolymers and creates "gaps".#
def creat_consensus_base_gap(base_array_input)
  base_array = Array.new(base_array_input)
  consensus_base = '-'
  number_of_bases = base_array.size
  hash_temp = Hash.new(0)
  base_array.each do |base|
    hash_temp[base] += 1
  end
  number_of_gap = hash_temp["-"]
  gap_percentage = (number_of_gap.to_f/number_of_bases.to_f)
  base_array.delete("-")
  h = Hash.new(0)
  if base_array.size >0
    base_array.each do |base|
      h[base] += 1
    end
    max_number = h.values.max
    max_list = []
    h.each do |k,v|
      if v == max_number
        max_list << k
      end
    end
    maxi_list_size = max_list.size
    if maxi_list_size == 1
      consensus_base = max_list.shift
    elsif maxi_list_size >= 3
      consensus_base = "N"
    elsif maxi_list_size == 2
      if max_list.include?("A") and max_list.include?("T")
        consensus_base = "W"
      elsif max_list.include?("A") and max_list.include?("C")
        consensus_base = "M"
      elsif max_list.include?("A") and max_list.include?("G")
        consensus_base = "R"
      elsif max_list.include?("T") and max_list.include?("C")
        consensus_base = "Y"
      elsif max_list.include?("G") and max_list.include?("C")
        consensus_base = "S"
      elsif max_list.include?("T") and max_list.include?("G")
        consensus_base = "M"
      end
    end
  end
  if gap_percentage >= 0.75
    consensus_base = "-"
  end
  return consensus_base
end

def creat_consensus_base_non_gap(base_array_input)
  base_array = Array.new(base_array_input)
  consensus_base = '-'
  number_of_bases = base_array.size
  h = Hash.new(0)
  if base_array.size >0
    base_array.each do |base|
      h[base] += 1
    end
    max_number = h.values.max
    max_list = []
    h.each do |k,v|
      if v == max_number
        max_list << k
      end
    end
    maxi_list_size = max_list.size
    if maxi_list_size == 1
      consensus_base = max_list.shift
    elsif maxi_list_size >= 3
      consensus_base = "N"
    elsif maxi_list_size == 2
      if max_list.include?("A") and max_list.include?("T")
        consensus_base = "W"
      elsif max_list.include?("A") and max_list.include?("C")
        consensus_base = "M"
      elsif max_list.include?("A") and max_list.include?("G")
        consensus_base = "R"
      elsif max_list.include?("T") and max_list.include?("C")
        consensus_base = "Y"
      elsif max_list.include?("G") and max_list.include?("C")
        consensus_base = "S"
      elsif max_list.include?("T") and max_list.include?("G")
        consensus_base = "K"
      elsif max_list.include?('-')
        max_list.delete('-')
        consensus_base = max_list.shift
      end
    end
  end
  return consensus_base
end

def creat_consensus_base_for_primer_design(base_array_input)
  base_array = Array.new(base_array_input)
  consensus_base = '-'
  h = count_percentage(base_array)
  if h.values.max > 0.85
    consensus_base = h.invert[h.invert.keys.max]
  else
    max_list = []
    h.each do |k,v|
      if v > 0.1
        max_list << k
      end
    end
    maxi_list_size = max_list.size
    if maxi_list_size == 1
      consensus_base = max_list.shift
    elsif maxi_list_size >= 3
      consensus_base = "N"
    elsif maxi_list_size == 2
      if max_list.include?("A") and max_list.include?("T")
        consensus_base = "W"
      elsif max_list.include?("A") and max_list.include?("C")
        consensus_base = "M"
      elsif max_list.include?("A") and max_list.include?("G")
        consensus_base = "R"
      elsif max_list.include?("T") and max_list.include?("C")
        consensus_base = "Y"
      elsif max_list.include?("G") and max_list.include?("C")
        consensus_base = "S"
      elsif max_list.include?("T") and max_list.include?("G")
        consensus_base = "M"
      elsif max_list.include?('-')
        consensus_base = "N"
      end
    end
  end
  return consensus_base
end


def consensus_without_alignment(seq_array,gap_treatment = 1)
  length = seq_array[0].size
  consensus_bases = []
  (0..(length-1)).each do |n|
    bases = []
    seq_array.each do |seq|
      bases << seq[n]
    end
    if gap_treatment == 1
      consensus_bases << creat_consensus_base_non_gap(bases)
    else
      consensus_bases << creat_consensus_base_gap(bases)
    end
  end
  consensus_seq = consensus_bases.join('')
end

#reverse complement
class String
    def rc
        self.reverse.tr("ACTG","TGAC")
    end
end

#math defs#
module Enumerable
  def median
    len = self.length
    sorted = self.sort
    median = len % 2 == 1 ? sorted[len/2] : (sorted[len/2 - 1] + sorted[len/2]).to_f / 2
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

end

def excel_upper_quartile(array)
  return nil if array.empty?
  sorted_array = array.sort
  u = (0.25*(3*sorted_array.length+1))
  if (u-u.truncate).is_a?(Integer)
    return sorted_array[(u-u.truncate)-1]
  else
    sample = sorted_array[u.truncate.abs-1]
    sample1 = sorted_array[(u.truncate.abs)]
    return sample+((sample1-sample)*(u-u.truncate))
  end
end

def excel_lower_quartile(array)
  return nil if array.empty?
  sorted_array = array.sort
  u = (0.25*(sorted_array.length+3))
  if (u-u.truncate).is_a?(Integer)
    return sorted_array[(u-u.truncate)-1]
  else
    sample = sorted_array[u.truncate.abs-1]
    sample1 = sorted_array[(u.truncate.abs)]
    return sample+((sample1-sample)*(u-u.truncate))
  end
end

def excel_upper_90(array)
   return nil if array.empty?
  sorted_array = array.sort
  u = (0.1*(9*sorted_array.length+1))
  if (u-u.truncate).is_a?(Integer)
    return sorted_array[(u-u.truncate)-1]
  else
    sample = sorted_array[u.truncate.abs-1]
    sample1 = sorted_array[(u.truncate.abs)]
    return sample+((sample1-sample)*(u-u.truncate))
  end
end

def factorial(n)
  sum = 1
  until n == 0
    sum *= n
    n -= 1
  end
 sum
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


#poisson distribution. input lambda and maximum k, return a hash with keys as k

def poisson_distribution(rate,k = 5)
  e = 2.718281828459045
  out_hash = {}
  (0..k).each do |n|
    p = (rate**n * e**(-rate))/factorial(n)
    out_hash[n] = p
  end
  return out_hash
end

#math defs end#

#copy hash#

def copyhash(inputhash)
  h = Hash.new
  inputhash.each do |pair|
    h.store(pair[0], pair[1])
  end
  return h
end

#copy hash end#

#open file in fasta format, turn into a sequence hash

def fasta_to_hash(infile)
  f=File.open(infile,"r")
  return_hash = {}
  name = ""
  while line = f.gets do
    line.tr!("\u0000","")
    next if line == "\n"
    next if line =~ /^\=/
    if line =~ /^\>/
      name = line.chomp
      return_hash[name] = ""
    else
      return_hash[name] += line.chomp.upcase
    end
  end
  f.close
  return return_hash
end

# convert fasta file into relaxed sequencial phylip format
def fasta_hash_to_rsphylip(seqs)
  outline = ""
  names = seqs.keys
  max_name_l = (names.max.size - 1)
  max_name_l > 10 ? name_block_l = max_name_l : name_block_l = 10
  seqs.each do |k,v|
    outline += k[1..-1] + "\s" * (name_block_l - k.size + 2) + v.scan(/.{1,10}/).join("\s") + "\n"
  end
  return outline
end


#fastq file to fasta, discard quality, return a sequence hash
def fastq_to_fasta(fastq_file)
    count = 0
    sequence_a = []
    count_seq = 0

    File.open(fastq_file,'r') do |file|
      file.readlines.collect do |line|
        count +=1
        count_m = count % 4
        if count_m == 1
          line.tr!('@','>')
          sequence_a << line.chomp
          count_seq += 1
        elsif count_m == 2
          sequence_a << line.chomp
        end
      end
    end
    sequence_hash = Hash[*sequence_a]
end

#fastq file to hash, including quality. seq_name=>[seq,quality]
def fastq_to_hash(fastq_file)
    count = 0
    sequence_a = []
    quality_a = []
    count_seq = 0

    File.open(fastq_file,'r') do |file|
      file.readlines.collect do |line|
        count +=1
        count_m = count % 4
        if count_m == 1
          line.tr!('@','>')
          sequence_a << line.chomp
          quality_a << line.chomp
          count_seq += 1
        elsif count_m == 2
          sequence_a << line.chomp
        elsif count_m == 0
          quality_a << line.chomp
        end
      end
    end
    sequence_hash = Hash[*sequence_a]
    quality_hash = Hash[*quality_a]
    return_hash = {}
    sequence_hash.each do |k,v|
      return_hash[k] = [v, quality_hash[k]]
    end
    return return_hash
end

#mutation of a sequence string
class String
  def mutation(error_rate = 0.01)
    new_string = ""
    self.split("").each do |nt|
      pool = ["A","C","T","G"]
      pool.delete(nt)
      s = error_rate * 10000
      r = rand(10000)
      if r < s
        nt = pool.sample
      end
      new_string << nt
    end
    return new_string
  end
end


#calculate consensus cut-off


def calculate_cut_off(m)
  n = 0
  if m <= 10
    n = 2
  elsif m <= 8500
    n = -1.24*10**-21*m**6 + 3.53*10**-17*m**5 - 3.90*10**-13*m**4 + 2.12*10**-9*m**3 - 6.06*10**-6*m**2 + 1.80*10**-2*m + 3.15
  else
    n = 0.0079 * m + 9.4869
  end
  n = n.round
  n = 2 if n < 3
  return n
end

#drug resistant mutation summary. input: amino acid array and starting codon, output, hash of summary

def sdrm_nrti(aa_array,start_aa=1)
  out_hash = {}
  sdrm = {}
  sdrm[41] = ['M',['L']]
  sdrm[65] = ['K',['R']]
  sdrm[67] = ['D',['N','G','E']]
  sdrm[69] = ['T',['D']]
  sdrm[70] = ['K',['R','E']]
  sdrm[74] = ['L',['V','I']]
  sdrm[75] = ['V',['M','T','A','S']]
  sdrm[77] = ['F',['L']]
  sdrm[115] = ['Y',['F']]
  sdrm[116] = ['F',['Y']]
  sdrm[151] = ['Q',['M']]
  sdrm[184] = ['M',['V','I']]
  sdrm[210] = ['L',['W']]
  sdrm[215] = ["T",["Y","F","I","C","D","V","E"]]
  sdrm[219] = ["K",["Q","E","N","R"]]
  aa_length = aa_array.size
  end_aa = start_aa + aa_length - 1
  (start_aa..end_aa).each do |position|
    array_position = position - start_aa
    if sdrm.keys.include?(position)
      wt_aa = sdrm[position][0]
      test_aa = aa_array[array_position]
      if test_aa.size == 1
        unless wt_aa == test_aa
          if sdrm[position][1].include?(test_aa)
            out_hash[position] = [wt_aa,test_aa]
          end
        end
      else
        test_aa_array = test_aa.split("/")
        if (test_aa_array & sdrm[position][1])
          out_hash[position] = [wt_aa,test_aa]
        end
      end

    end
  end
  return out_hash
end

def sdrm_nnrti(aa_array,start_aa=1)
  out_hash = {}
  sdrm = {}
  sdrm[100] = ['L',['I']]
  sdrm[101] = ['K',['E','P']]
  sdrm[103] = ['K',['N','S']]
  sdrm[106] = ['V',['M','A']]
  sdrm[179] = ['V',['F','D']]
  sdrm[181] = ['Y',['C','I','V']]
  sdrm[188] = ['Y',['L','H','C']]
  sdrm[190] = ['G',['A','S','E']]
  sdrm[225] = ['P',['H']]
  sdrm[230] = ['M',['L']]
  aa_length = aa_array.size
  end_aa = start_aa + aa_length - 1
  (start_aa..end_aa).each do |position|
    array_position = position - start_aa
    if sdrm.keys.include?(position)
      wt_aa = sdrm[position][0]
      test_aa = aa_array[array_position]
      if test_aa.size == 1
        unless wt_aa == test_aa
          if sdrm[position][1].include?(test_aa)
            out_hash[position] = [wt_aa,test_aa]
          end
        end
      else
        test_aa_array = test_aa.split("/")
        if (test_aa_array & sdrm[position][1])
          out_hash[position] = [wt_aa,test_aa]
        end
      end

    end
  end
  return out_hash
end

#APOBEC sentinel mutation on RT
def sentinel_rt (aa_array,start_aa=1)
  out_hash = {}
  sen = {15=>"KRS", 16=>"I", 17=>"N", 18=>"KRS", 24=>"*", 41=>"I", 42=>"KN", 45=>"KR", 51=>"KR", 53=>"K", 71=>"*", 76=>"K", 79=>"K", 86=>"K", 88=>"*", 89=>"KN", 93=>"EKR", 99=>"EKRS", 110=>"KN", 112=>"EKRS", 113=>"K", 141=>"EKRS", 152=>"EKRS", 153=>"*", 155=>"EKRS", 185=>"KN", 186=>"KN", 190=>"R", 192=>"K", 212=>"*", 213=>"KR", 229=>"*", 231=>"KRS", 233=>"K", 239=>"*", 252=>"*", 256=>"KN", 262=>"KR", 266=>"*", 273=>"EKS", 285=>"EKR", 298=>"N", 302=>"KN", 305=>"KN", 316=>"EKR", 324=>"K", 337=>"*", 344=>"K", 352=>"EKRS", 364=>"KN", 378=>"N", 383=>"*", 384=>"EKRS", 396=>"K", 398=>"*", 401=>"*", 402=>"*", 404=>"K", 406=>"*", 410=>"*", 413=>"KN", 414=>"*", 415=>"KN", 426=>"*", 430=>"KN", 436=>"R", 438=>"K", 443=>"N", 444=>"EKRS", 449=>"K", 453=>"KR", 456=>"EKRS", 462=>"ERS", 471=>"K", 478=>"KN", 488=>"K", 490=>"R", 498=>"KN", 504=>"KRS", 511=>"K", 514=>"K", 516=>"N", 523=>"KN", 535=>"*", 541=>"EKRS", 543=>"EKRS", 544=>"ERS", 546=>"N", 549=>"KN", 555=>"EKRS"}
  aa_length = aa_array.size
  end_aa = start_aa + aa_length - 1
  (start_aa..end_aa).each do |position|
    array_position = position - start_aa
    if sen.keys.include?(position)
      test_aa = aa_array[array_position]
      if sen[position].include?(test_aa)
        out_hash[position] = test_aa
      end
    end
  end
  return out_hash
end

#HCV NS5A resistant mutation

def hcv_ns5a(aa_array,start_aa=1)
  out_hash = {}
  sdrm = {}
  sdrm[28] = ['M',['T']]
  sdrm[30] = ['L',['H','K','R','Q','A','S','D']]
  sdrm[31] = ['L',['M','V','F']]
  sdrm[32] = ['P',['L']]
  sdrm[44] = ['K',['R']]
  sdrm[58] = ['H',['D','P','S']]
  sdrm[64] = ['T',['A','S']]
  sdrm[77] = ['P',['A','S']]
  sdrm[78] = ['R',['K']]
  sdrm[79] = ['T',['A']]
  sdrm[83] = ['T',['M']]
  sdrm[85] = ['S',['N','H','Y']]
  sdrm[92] = ['A',['P','T','K','E']]
  sdrm[93] = ['Y',['C','F','H','N']]
  sdrm[107] = ['K',['T','S']]
  sdrm[121] = ['I',['V']]
  sdrm[135] = ['T',['A']]
  aa_length = aa_array.size
  end_aa = start_aa + aa_length - 1
  (start_aa..end_aa).each do |position|
    array_position = position - start_aa
    if sdrm.keys.include?(position)
      wt_aa = sdrm[position][0]
      test_aa = aa_array[array_position]
      if test_aa.size == 1
        unless wt_aa == test_aa
          if sdrm[position][1].include?(test_aa)
            out_hash[position] = [wt_aa,test_aa]
          end
        end
      else
        test_aa_array = test_aa.split("/")
        if (test_aa_array & sdrm[position][1])
          out_hash[position] = [wt_aa,test_aa]
        end
      end
    end
  end
  return out_hash
end

#HIV protease surveillance mutations

def hiv_protease(aa_array,start_aa=1)
  out_hash = {}
  sdrm = {}
  sdrm[23] = ['L',['I']]
  sdrm[24] = ['L',['I']]
  sdrm[30] = ['D',['N']]
  sdrm[32] = ['V',['I']]
  sdrm[46] = ['M',['I','L','V']]
  sdrm[47] = ['I',['V','A']]
  sdrm[48] = ['G',['V','M']]
  sdrm[50] = ['I',['V','L']]
  sdrm[53] = ['F',['Y']]
  sdrm[54] = ['I',['V','L','M','T','A','S']]
  sdrm[73] = ['G',['S','T','C','A']]
  sdrm[76] = ['L',['V']]
  sdrm[82] = ['V',['A','T','S','F','L','C','M']]
  sdrm[83] = ['N',['D']]
  sdrm[84] = ['I',['V','A','C']]
  sdrm[85] = ['I',['V']]
  sdrm[88] = ['N',['D','S']]
  sdrm[90] = ['L',['M']]
  aa_length = aa_array.size
  end_aa = start_aa + aa_length - 1
  (start_aa..end_aa).each do |position|
    array_position = position - start_aa
    if sdrm.keys.include?(position)
      wt_aa = sdrm[position][0]
      test_aa = aa_array[array_position]
      if test_aa.size == 1
        unless wt_aa == test_aa
          if sdrm[position][1].include?(test_aa)
            out_hash[position] = [wt_aa,test_aa]
          end
        end
      else
        test_aa_array = test_aa.split("/")
        if (test_aa_array & sdrm[position][1])
          out_hash[position] = [wt_aa,test_aa]
        end
      end
    end
  end
  return out_hash
end

#HIV integrase drug resistance mutations

def sdrm_int(aa_array,start_aa=1)
  out_hash = {}
  sdrm = {}
  sdrm[66] = ['T',['A','I','K']]
  sdrm[68] = ['L',['V']]
  sdrm[74] = ['L',['M']]
  sdrm[92] = ['E',['Q','G','V']]
  sdrm[95] = ['Q',['K']]
  sdrm[97] = ['T',['A']]
  sdrm[114] = ['H',['Y']]
  sdrm[118] = ['G',['R']]
  sdrm[121] = ['F',['Y']]
  sdrm[128] = ['A',['T']]
  sdrm[138] = ['E',['D','K','A']]
  sdrm[140] = ['G',['A','S','C']]
  sdrm[143] = ["Y",["H","C","R","K","G","S","A"]]
  sdrm[145] = ['P',['S']]
  sdrm[146] = ['Q',['P']]
  sdrm[147] = ['S',['G']]
  sdrm[148] = ['Q',['H','K','R']]
  sdrm[151] = ['V',['A','L','I']]
  sdrm[153] = ['S',['Y','F']]
  sdrm[155] = ['N',['S','T','H']]
  sdrm[157] = ['E',['Q']]
  sdrm[163] = ['G',['P','K']]
  aa_length = aa_array.size
  end_aa = start_aa + aa_length - 1
  (start_aa..end_aa).each do |position|
    array_position = position - start_aa
    if sdrm.keys.include?(position)
      wt_aa = sdrm[position][0]
      test_aa = aa_array[array_position]
      if test_aa.size == 1
        unless wt_aa == test_aa
          if sdrm[position][1].include?(test_aa)
            out_hash[position] = [wt_aa,test_aa]
          end
        end
      else
        test_aa_array = test_aa.split("/")
        if (test_aa_array & sdrm[position][1])
          out_hash[position] = [wt_aa,test_aa]
        end
      end

    end
  end
  return out_hash
end



class RandomGaussian
  def initialize(mean = 0.0, sd = 1.0, rng = lambda { Kernel.rand })
    @mean, @sd, @rng = mean, sd, rng
    @compute_next_pair = false
  end

  def rand
    if (@compute_next_pair = !@compute_next_pair)
      # Compute a pair of random values with normal distribution.
      # See http://en.wikipedia.org/wiki/Box-Muller_transform
      theta = 2 * Math::PI * @rng.call
      scale = @sd * Math.sqrt(-2 * Math.log(1 - @rng.call))
      @g1 = @mean + scale * Math.sin(theta)
      @g0 = @mean + scale * Math.cos(theta)
    else
      @g1
    end
  end
end

def count_list(a)
  b = a.split"\n"
  c = count(b)
  (0..7).to_a.each do |n|
    puts c[n.to_s]
  end
end

def pid_to_n(a)
  a.tr("ATCG","").split("\n").collect{|n|n.to_i}
end



def pos(s,n)
  p = []
  s.values.each do |v|
    p << v[n]
  end
  count(p)
end

#compare two primer ID sequences. If they differ in x base, return boolean value "TURE", else, return boolean value "FALSE"
def two_pid_x_base_different(pid1="",pid2="", x=0)
  l = pid1.size
  m = l - x
  n = 0
  if pid1.size != pid2.size
    return false
  else
    (0..(pid1.size - 1)).each do |k|
      if pid1[k] == pid2[k]
        n += 1
      end
    end
    if n >= m
      return true
    else
      return false
    end
  end
end


#generate all Primer ID combinations given the length of Primer ID
def generate_primer_id_pool(l=8)
  nt = ['A','T','C','G']
  pid_pool = ['A','T','C','G']
  (l-1).times do
    pid_pool = pid_pool.product(nt)
    pid_pool.collect! do |v|
      v.join("")
    end
  end
  return pid_pool
end

#irb only
#copy two columns from excel (a and b) of primer ID distribution, return mean and stdev of raw/primer ID
def return_mean_stdev_from_excel(a="",b="")
  a = a.split"\n"
  b = b.split"\n"
  h = {}
  (0..(a.size-1)).each do |n|
    h[a[n].to_i] = b[n].to_i
  end
  aa = []
  h.each {|k,v| v.times {aa << k}}
  r1 = aa.mean
  r2 = aa.stdev
  return [r1,r2]
end


#irb only
#turn strings to array with single characters.

class Array
  def single
    new_array = []
    self.each do |v|
      new_array << (v.split("")).join("\s")
    end
    return new_array
  end
end

#for kure bsub. muscle to build tree
def bsub_muscle(names)
  names.split("\n").each do |name|
    puts "bsub muscle -in #{name} -out #{name + ".fa"} -tree2 #{name + ".tre"} -maxiters 2 -cluster neighborjoining"
  end
end

#sequence locator:
#align with HXB2, return location, percentage of similarity
def sequence_locator(seq="",temp_dir=File.dirname($0))
  hxb2_ref = "TGGAAGGGCTAATTCACTCCCAACGAAGACAAGATATCCTTGATCTGTGGATCTACCACACACAAGGCTACTTCCCTGATTAGCAGAACTACACACCAGGGCCAGGGATCAGATATCCACTGACCTTTGGATGGTGCTACAAGCTAGTACCAGTTGAGCCAGAGAAGTTAGAAGAAGCCAACAAAGGAGAGAACACCAGCTTGTTACACCCTGTGAGCCTGCATGGAATGGATGACCCGGAGAGAGAAGTGTTAGAGTGGAGGTTTGACAGCCGCCTAGCATTTCATCACATGGCCCGAGAGCTGCATCCGGAGTACTTCAAGAACTGCTGACATCGAGCTTGCTACAAGGGACTTTCCGCTGGGGACTTTCCAGGGAGGCGTGGCCTGGGCGGGACTGGGGAGTGGCGAGCCCTCAGATCCTGCATATAAGCAGCTGCTTTTTGCCTGTACTGGGTCTCTCTGGTTAGACCAGATCTGAGCCTGGGAGCTCTCTGGCTAACTAGGGAACCCACTGCTTAAGCCTCAATAAAGCTTGCCTTGAGTGCTTCAAGTAGTGTGTGCCCGTCTGTTGTGTGACTCTGGTAACTAGAGATCCCTCAGACCCTTTTAGTCAGTGTGGAAAATCTCTAGCAGTGGCGCCCGAACAGGGACCTGAAAGCGAAAGGGAAACCAGAGGAGCTCTCTCGACGCAGGACTCGGCTTGCTGAAGCGCGCACGGCAAGAGGCGAGGGGCGGCGACTGGTGAGTACGCCAAAAATTTTGACTAGCGGAGGCTAGAAGGAGAGAGATGGGTGCGAGAGCGTCAGTATTAAGCGGGGGAGAATTAGATCGATGGGAAAAAATTCGGTTAAGGCCAGGGGGAAAGAAAAAATATAAATTAAAACATATAGTATGGGCAAGCAGGGAGCTAGAACGATTCGCAGTTAATCCTGGCCTGTTAGAAACATCAGAAGGCTGTAGACAAATACTGGGACAGCTACAACCATCCCTTCAGACAGGATCAGAAGAACTTAGATCATTATATAATACAGTAGCAACCCTCTATTGTGTGCATCAAAGGATAGAGATAAAAGACACCAAGGAAGCTTTAGACAAGATAGAGGAAGAGCAAAACAAAAGTAAGAAAAAAGCACAGCAAGCAGCAGCTGACACAGGACACAGCAATCAGGTCAGCCAAAATTACCCTATAGTGCAGAACATCCAGGGGCAAATGGTACATCAGGCCATATCACCTAGAACTTTAAATGCATGGGTAAAAGTAGTAGAAGAGAAGGCTTTCAGCCCAGAAGTGATACCCATGTTTTCAGCATTATCAGAAGGAGCCACCCCACAAGATTTAAACACCATGCTAAACACAGTGGGGGGACATCAAGCAGCCATGCAAATGTTAAAAGAGACCATCAATGAGGAAGCTGCAGAATGGGATAGAGTGCATCCAGTGCATGCAGGGCCTATTGCACCAGGCCAGATGAGAGAACCAAGGGGAAGTGACATAGCAGGAACTACTAGTACCCTTCAGGAACAAATAGGATGGATGACAAATAATCCACCTATCCCAGTAGGAGAAATTTATAAAAGATGGATAATCCTGGGATTAAATAAAATAGTAAGAATGTATAGCCCTACCAGCATTCTGGACATAAGACAAGGACCAAAGGAACCCTTTAGAGACTATGTAGACCGGTTCTATAAAACTCTAAGAGCCGAGCAAGCTTCACAGGAGGTAAAAAATTGGATGACAGAAACCTTGTTGGTCCAAAATGCGAACCCAGATTGTAAGACTATTTTAAAAGCATTGGGACCAGCGGCTACACTAGAAGAAATGATGACAGCATGTCAGGGAGTAGGAGGACCCGGCCATAAGGCAAGAGTTTTGGCTGAAGCAATGAGCCAAGTAACAAATTCAGCTACCATAATGATGCAGAGAGGCAATTTTAGGAACCAAAGAAAGATTGTTAAGTGTTTCAATTGTGGCAAAGAAGGGCACACAGCCAGAAATTGCAGGGCCCCTAGGAAAAAGGGCTGTTGGAAATGTGGAAAGGAAGGACACCAAATGAAAGATTGTACTGAGAGACAGGCTAATTTTTTAGGGAAGATCTGGCCTTCCTACAAGGGAAGGCCAGGGAATTTTCTTCAGAGCAGACCAGAGCCAACAGCCCCACCAGAAGAGAGCTTCAGGTCTGGGGTAGAGACAACAACTCCCCCTCAGAAGCAGGAGCCGATAGACAAGGAACTGTATCCTTTAACTTCCCTCAGGTCACTCTTTGGCAACGACCCCTCGTCACAATAAAGATAGGGGGGCAACTAAAGGAAGCTCTATTAGATACAGGAGCAGATGATACAGTATTAGAAGAAATGAGTTTGCCAGGAAGATGGAAACCAAAAATGATAGGGGGAATTGGAGGTTTTATCAAAGTAAGACAGTATGATCAGATACTCATAGAAATCTGTGGACATAAAGCTATAGGTACAGTATTAGTAGGACCTACACCTGTCAACATAATTGGAAGAAATCTGTTGACTCAGATTGGTTGCACTTTAAATTTTCCCATTAGCCCTATTGAGACTGTACCAGTAAAATTAAAGCCAGGAATGGATGGCCCAAAAGTTAAACAATGGCCATTGACAGAAGAAAAAATAAAAGCATTAGTAGAAATTTGTACAGAGATGGAAAAGGAAGGGAAAATTTCAAAAATTGGGCCTGAAAATCCATACAATACTCCAGTATTTGCCATAAAGAAAAAAGACAGTACTAAATGGAGAAAATTAGTAGATTTCAGAGAACTTAATAAGAGAACTCAAGACTTCTGGGAAGTTCAATTAGGAATACCACATCCCGCAGGGTTAAAAAAGAAAAAATCAGTAACAGTACTGGATGTGGGTGATGCATATTTTTCAGTTCCCTTAGATGAAGACTTCAGGAAGTATACTGCATTTACCATACCTAGTATAAACAATGAGACACCAGGGATTAGATATCAGTACAATGTGCTTCCACAGGGATGGAAAGGATCACCAGCAATATTCCAAAGTAGCATGACAAAAATCTTAGAGCCTTTTAGAAAACAAAATCCAGACATAGTTATCTATCAATACATGGATGATTTGTATGTAGGATCTGACTTAGAAATAGGGCAGCATAGAACAAAAATAGAGGAGCTGAGACAACATCTGTTGAGGTGGGGACTTACCACACCAGACAAAAAACATCAGAAAGAACCTCCATTCCTTTGGATGGGTTATGAACTCCATCCTGATAAATGGACAGTACAGCCTATAGTGCTGCCAGAAAAAGACAGCTGGACTGTCAATGACATACAGAAGTTAGTGGGGAAATTGAATTGGGCAAGTCAGATTTACCCAGGGATTAAAGTAAGGCAATTATGTAAACTCCTTAGAGGAACCAAAGCACTAACAGAAGTAATACCACTAACAGAAGAAGCAGAGCTAGAACTGGCAGAAAACAGAGAGATTCTAAAAGAACCAGTACATGGAGTGTATTATGACCCATCAAAAGACTTAATAGCAGAAATACAGAAGCAGGGGCAAGGCCAATGGACATATCAAATTTATCAAGAGCCATTTAAAAATCTGAAAACAGGAAAATATGCAAGAATGAGGGGTGCCCACACTAATGATGTAAAACAATTAACAGAGGCAGTGCAAAAAATAACCACAGAAAGCATAGTAATATGGGGAAAGACTCCTAAATTTAAACTGCCCATACAAAAGGAAACATGGGAAACATGGTGGACAGAGTATTGGCAAGCCACCTGGATTCCTGAGTGGGAGTTTGTTAATACCCCTCCCTTAGTGAAATTATGGTACCAGTTAGAGAAAGAACCCATAGTAGGAGCAGAAACCTTCTATGTAGATGGGGCAGCTAACAGGGAGACTAAATTAGGAAAAGCAGGATATGTTACTAATAGAGGAAGACAAAAAGTTGTCACCCTAACTGACACAACAAATCAGAAGACTGAGTTACAAGCAATTTATCTAGCTTTGCAGGATTCGGGATTAGAAGTAAACATAGTAACAGACTCACAATATGCATTAGGAATCATTCAAGCACAACCAGATCAAAGTGAATCAGAGTTAGTCAATCAAATAATAGAGCAGTTAATAAAAAAGGAAAAGGTCTATCTGGCATGGGTACCAGCACACAAAGGAATTGGAGGAAATGAACAAGTAGATAAATTAGTCAGTGCTGGAATCAGGAAAGTACTATTTTTAGATGGAATAGATAAGGCCCAAGATGAACATGAGAAATATCACAGTAATTGGAGAGCAATGGCTAGTGATTTTAACCTGCCACCTGTAGTAGCAAAAGAAATAGTAGCCAGCTGTGATAAATGTCAGCTAAAAGGAGAAGCCATGCATGGACAAGTAGACTGTAGTCCAGGAATATGGCAACTAGATTGTACACATTTAGAAGGAAAAGTTATCCTGGTAGCAGTTCATGTAGCCAGTGGATATATAGAAGCAGAAGTTATTCCAGCAGAAACAGGGCAGGAAACAGCATATTTTCTTTTAAAATTAGCAGGAAGATGGCCAGTAAAAACAATACATACTGACAATGGCAGCAATTTCACCGGTGCTACGGTTAGGGCCGCCTGTTGGTGGGCGGGAATCAAGCAGGAATTTGGAATTCCCTACAATCCCCAAAGTCAAGGAGTAGTAGAATCTATGAATAAAGAATTAAAGAAAATTATAGGACAGGTAAGAGATCAGGCTGAACATCTTAAGACAGCAGTACAAATGGCAGTATTCATCCACAATTTTAAAAGAAAAGGGGGGATTGGGGGGTACAGTGCAGGGGAAAGAATAGTAGACATAATAGCAACAGACATACAAACTAAAGAATTACAAAAACAAATTACAAAAATTCAAAATTTTCGGGTTTATTACAGGGACAGCAGAAATCCACTTTGGAAAGGACCAGCAAAGCTCCTCTGGAAAGGTGAAGGGGCAGTAGTAATACAAGATAATAGTGACATAAAAGTAGTGCCAAGAAGAAAAGCAAAGATCATTAGGGATTATGGAAAACAGATGGCAGGTGATGATTGTGTGGCAAGTAGACAGGATGAGGATTAGAACATGGAAAAGTTTAGTAAAACACCATATGTATGTTTCAGGGAAAGCTAGGGGATGGTTTTATAGACATCACTATGAAAGCCCTCATCCAAGAATAAGTTCAGAAGTACACATCCCACTAGGGGATGCTAGATTGGTAATAACAACATATTGGGGTCTGCATACAGGAGAAAGAGACTGGCATTTGGGTCAGGGAGTCTCCATAGAATGGAGGAAAAAGAGATATAGCACACAAGTAGACCCTGAACTAGCAGACCAACTAATTCATCTGTATTACTTTGACTGTTTTTCAGACTCTGCTATAAGAAAGGCCTTATTAGGACACATAGTTAGCCCTAGGTGTGAATATCAAGCAGGACATAACAAGGTAGGATCTCTACAATACTTGGCACTAGCAGCATTAATAACACCAAAAAAGATAAAGCCACCTTTGCCTAGTGTTACGAAACTGACAGAGGATAGATGGAACAAGCCCCAGAAGACCAAGGGCCACAGAGGGAGCCACACAATGAATGGACACTAGAGCTTTTAGAGGAGCTTAAGAATGAAGCTGTTAGACATTTTCCTAGGATTTGGCTCCATGGCTTAGGGCAACATATCTATGAAACTTATGGGGATACTTGGGCAGGAGTGGAAGCCATAATAAGAATTCTGCAACAACTGCTGTTTATCCATTTTCAGAATTGGGTGTCGACATAGCAGAATAGGCGTTACTCGACAGAGGAGAGCAAGAAATGGAGCCAGTAGATCCTAGACTAGAGCCCTGGAAGCATCCAGGAAGTCAGCCTAAAACTGCTTGTACCAATTGCTATTGTAAAAAGTGTTGCTTTCATTGCCAAGTTTGTTTCATAACAAAAGCCTTAGGCATCTCCTATGGCAGGAAGAAGCGGAGACAGCGACGAAGAGCTCATCAGAACAGTCAGACTCATCAAGCTTCTCTATCAAAGCAGTAAGTAGTACATGTAACGCAACCTATACCAATAGTAGCAATAGTAGCATTAGTAGTAGCAATAATAATAGCAATAGTTGTGTGGTCCATAGTAATCATAGAATATAGGAAAATATTAAGACAAAGAAAAATAGACAGGTTAATTGATAGACTAATAGAAAGAGCAGAAGACAGTGGCAATGAGAGTGAAGGAGAAATATCAGCACTTGTGGAGATGGGGGTGGAGATGGGGCACCATGCTCCTTGGGATGTTGATGATCTGTAGTGCTACAGAAAAATTGTGGGTCACAGTCTATTATGGGGTACCTGTGTGGAAGGAAGCAACCACCACTCTATTTTGTGCATCAGATGCTAAAGCATATGATACAGAGGTACATAATGTTTGGGCCACACATGCCTGTGTACCCACAGACCCCAACCCACAAGAAGTAGTATTGGTAAATGTGACAGAAAATTTTAACATGTGGAAAAATGACATGGTAGAACAGATGCATGAGGATATAATCAGTTTATGGGATCAAAGCCTAAAGCCATGTGTAAAATTAACCCCACTCTGTGTTAGTTTAAAGTGCACTGATTTGAAGAATGATACTAATACCAATAGTAGTAGCGGGAGAATGATAATGGAGAAAGGAGAGATAAAAAACTGCTCTTTCAATATCAGCACAAGCATAAGAGGTAAGGTGCAGAAAGAATATGCATTTTTTTATAAACTTGATATAATACCAATAGATAATGATACTACCAGCTATAAGTTGACAAGTTGTAACACCTCAGTCATTACACAGGCCTGTCCAAAGGTATCCTTTGAGCCAATTCCCATACATTATTGTGCCCCGGCTGGTTTTGCGATTCTAAAATGTAATAATAAGACGTTCAATGGAACAGGACCATGTACAAATGTCAGCACAGTACAATGTACACATGGAATTAGGCCAGTAGTATCAACTCAACTGCTGTTAAATGGCAGTCTAGCAGAAGAAGAGGTAGTAATTAGATCTGTCAATTTCACGGACAATGCTAAAACCATAATAGTACAGCTGAACACATCTGTAGAAATTAATTGTACAAGACCCAACAACAATACAAGAAAAAGAATCCGTATCCAGAGAGGACCAGGGAGAGCATTTGTTACAATAGGAAAAATAGGAAATATGAGACAAGCACATTGTAACATTAGTAGAGCAAAATGGAATAACACTTTAAAACAGATAGCTAGCAAATTAAGAGAACAATTTGGAAATAATAAAACAATAATCTTTAAGCAATCCTCAGGAGGGGACCCAGAAATTGTAACGCACAGTTTTAATTGTGGAGGGGAATTTTTCTACTGTAATTCAACACAACTGTTTAATAGTACTTGGTTTAATAGTACTTGGAGTACTGAAGGGTCAAATAACACTGAAGGAAGTGACACAATCACCCTCCCATGCAGAATAAAACAAATTATAAACATGTGGCAGAAAGTAGGAAAAGCAATGTATGCCCCTCCCATCAGTGGACAAATTAGATGTTCATCAAATATTACAGGGCTGCTATTAACAAGAGATGGTGGTAATAGCAACAATGAGTCCGAGATCTTCAGACCTGGAGGAGGAGATATGAGGGACAATTGGAGAAGTGAATTATATAAATATAAAGTAGTAAAAATTGAACCATTAGGAGTAGCACCCACCAAGGCAAAGAGAAGAGTGGTGCAGAGAGAAAAAAGAGCAGTGGGAATAGGAGCTTTGTTCCTTGGGTTCTTGGGAGCAGCAGGAAGCACTATGGGCGCAGCCTCAATGACGCTGACGGTACAGGCCAGACAATTATTGTCTGGTATAGTGCAGCAGCAGAACAATTTGCTGAGGGCTATTGAGGCGCAACAGCATCTGTTGCAACTCACAGTCTGGGGCATCAAGCAGCTCCAGGCAAGAATCCTGGCTGTGGAAAGATACCTAAAGGATCAACAGCTCCTGGGGATTTGGGGTTGCTCTGGAAAACTCATTTGCACCACTGCTGTGCCTTGGAATGCTAGTTGGAGTAATAAATCTCTGGAACAGATTTGGAATCACACGACCTGGATGGAGTGGGACAGAGAAATTAACAATTACACAAGCTTAATACACTCCTTAATTGAAGAATCGCAAAACCAGCAAGAAAAGAATGAACAAGAATTATTGGAATTAGATAAATGGGCAAGTTTGTGGAATTGGTTTAACATAACAAATTGGCTGTGGTATATAAAATTATTCATAATGATAGTAGGAGGCTTGGTAGGTTTAAGAATAGTTTTTGCTGTACTTTCTATAGTGAATAGAGTTAGGCAGGGATATTCACCATTATCGTTTCAGACCCACCTCCCAACCCCGAGGGGACCCGACAGGCCCGAAGGAATAGAAGAAGAAGGTGGAGAGAGAGACAGAGACAGATCCATTCGATTAGTGAACGGATCCTTGGCACTTATCTGGGACGATCTGCGGAGCCTGTGCCTCTTCAGCTACCACCGCTTGAGAGACTTACTCTTGATTGTAACGAGGATTGTGGAACTTCTGGGACGCAGGGGGTGGGAAGCCCTCAAATATTGGTGGAATCTCCTACAGTATTGGAGTCAGGAACTAAAGAATAGTGCTGTTAGCTTGCTCAATGCCACAGCCATAGCAGTAGCTGAGGGGACAGATAGGGTTATAGAAGTAGTACAAGGAGCTTGTAGAGCTATTCGCCACATACCTAGAAGAATAAGACAGGGCTTGGAAAGGATTTTGCTATAAGATGGGTGGCAAGTGGTCAAAAAGTAGTGTGATTGGATGGCCTACTGTAAGGGAAAGAATGAGACGAGCTGAGCCAGCAGCAGATAGGGTGGGAGCAGCATCTCGAGACCTGGAAAAACATGGAGCAATCACAAGTAGCAATACAGCAGCTACCAATGCTGCTTGTGCCTGGCTAGAAGCACAAGAGGAGGAGGAGGTGGGTTTTCCAGTCACACCTCAGGTACCTTTAAGACCAATGACTTACAAGGCAGCTGTAGATCTTAGCCACTTTTTAAAAGAAAAGGGGGGACTGGAAGGGCTAATTCACTCCCAAAGAAGACAAGATATCCTTGATCTGTGGATCTACCACACACAAGGCTACTTCCCTGATTAGCAGAACTACACACCAGGGCCAGGGGTCAGATATCCACTGACCTTTGGATGGTGCTACAAGCTAGTACCAGTTGAGCCAGATAAGATAGAAGAGGCCAATAAAGGAGAGAACACCAGCTTGTTACACCCTGTGAGCCTGCATGGGATGGATGACCCGGAGAGAGAAGTGTTAGAGTGGAGGTTTGACAGCCGCCTAGCATTTCATCACGTGGCCCGAGAGCTGCATCCGGAGTACTTCAAGAACTGCTGACATCGAGCTTGCTACAAGGGACTTTCCGCTGGGGACTTTCCAGGGAGGCGTGGCCTGGGCGGGACTGGGGAGTGGCGAGCCCTCAGATCCTGCATATAAGCAGCTGCTTTTTGCCTGTACTGGGTCTCTCTGGTTAGACCAGATCTGAGCCTGGGAGCTCTCTGGCTAACTAGGGAACCCACTGCTTAAGCCTCAATAAAGCTTGCCTTGAGTGCTTCAAGTAGTGTGTGCCCGTCTGTTGTGTGACTCTGGTAACTAGAGATCCCTCAGACCCTTTTAGTCAGTGTGGAAAATCTCTAGCA"
  hxb2_l = hxb2_ref.size
  head = ""
  8.times {head  << (65 + rand(25)).chr}
  temp_file = temp_dir + "/" + head + "_temp"
  temp_aln = temp_dir + "/" + head + "_temp_aln"

  l1 = 0
  l2 = 0
  name = ">test"
  temp_in = File.open(temp_file,"w")
  temp_in.puts ">ref"
  temp_in.puts hxb2_ref
  temp_in.puts name
  temp_in.puts seq
  temp_in.close

  print `muscle -in #{temp_file} -out #{temp_aln} -quiet`
  aln_seq = fasta_to_hash(temp_aln)
  aln_test = aln_seq[name]
  aln_test =~ /^(\-*)(\w.*\w)(\-*)$/
  gap_begin = $1.size
  gap_end = $3.size
  aln_test2 = $2
  ref = aln_seq[">ref"]
  ref = ref[gap_begin..(-gap_end-1)]
  ref_size = ref.size
  if ref_size > 1.3*(seq.size)
    l1 = l1 + gap_begin
    l2 = l2 + gap_end
    max_seq = aln_test2.scan(/[ACGT]+/).max_by(&:length)
    aln_test2 =~ /#{max_seq}/
    before_aln_seq = $`
    before_aln = $`.size
    post_aln_seq = $'
    post_aln = $'.size
    before_aln_seq_size = before_aln_seq.scan(/[ACGT]+/).join("").size
    b1 = (1.3 * before_aln_seq_size).to_i
    post_aln_seq_size = post_aln_seq.scan(/[ACGT]+/).join("").size
    b2 = (1.3 * post_aln_seq_size).to_i
    if (before_aln > seq.size) and (post_aln <= seq.size)
      ref = ref[(before_aln - b1)..(ref_size - post_aln - 1)]
      l1 = l1 + (before_aln - b1)
    elsif (post_aln > seq.size) and (before_aln <= seq.size)
      ref = ref[before_aln..(ref_size - post_aln - 1 + b2)]
      l2 = l2 + post_aln - b2
    elsif (post_aln > seq.size) and (before_aln > seq.size)
      ref = ref[(before_aln - b1)..(ref_size - post_aln - 1 + b2)]
      l1 = l1 + (before_aln - b1)
      l2 = l2 + (post_aln - b2)
    end
    temp_in = File.open(temp_file,"w")
    temp_in.puts ">ref"
    temp_in.puts ref
    temp_in.puts name
    temp_in.puts seq
    temp_in.close
    print `muscle -in #{temp_file} -out #{temp_aln} -quiet`
    aln_seq = fasta_to_hash(temp_aln)
    aln_test = aln_seq[name]
    aln_test =~ /^(\-*)(\w.*\w)(\-*)$/
    gap_begin = $1.size
    gap_end = $3.size
    aln_test2 = $2
    ref = aln_seq[">ref"]
    ref = ref[gap_begin..(-gap_end-1)]
    ref_size = ref.size
  end
  aln_seq = fasta_to_hash(temp_aln)
  aln_test = aln_seq[name]
  aln_test =~ /^(\-*)(\w.*\w)(\-*)$/
  gap_begin = $1.size
  gap_end = $3.size
  aln_test = $2
  aln_test =~ /^(\w+)(\-*)\w/
  s1 = $1.size
  g1 = $2.size
  aln_test =~ /\w(\-*)(\w+)$/
  s2 = $2.size
  g2 = $1.size
  ref = aln_seq[">ref"]
  ref = ref[gap_begin..(-gap_end-1)]

  l1 = l1 + gap_begin
  l2 = l2 + gap_end
  repeat = 0

  if g1 == g2 and (s1 + g1 + s2) == ref.size
    if s1 > s2 and g2 > 2*s2
      ref = ref[0..(-g2-1)]
      repeat = 1
      l2 = l2 + g2
    elsif s1 < s2 and g1 > 2*s1
      ref = ref[g1..-1]
      repeat = 1
      l1 = l1 + g1
    end
  else
    if g1 > 2*s1
      ref = ref[g1..-1]
      repeat = 1
      l1 = l1 + g1
    end
    if g2 > 2*s2
      ref = ref[0..(-g2 - 1)]
      repeat = 1
      l2 = l2 + g2
    end
  end

  while repeat == 1
    temp_in = File.open(temp_file,"w")
    temp_in.puts ">ref"
    temp_in.puts ref
    temp_in.puts name
    temp_in.puts seq
    temp_in.close
    print `muscle -in #{temp_file} -out #{temp_aln} -quiet`
    aln_seq = fasta_to_hash(temp_aln)
    aln_test = aln_seq[name]
    aln_test =~ /^(\-*)(\w.*\w)(\-*)$/
    gap_begin = $1.size
    gap_end = $3.size
    aln_test = $2
    aln_test =~ /^(\w+)(\-*)\w/
    s1 = $1.size
    g1 = $2.size
    aln_test =~ /\w(\-*)(\w+)$/
    s2 = $2.size
    g2 = $1.size
    ref = aln_seq[">ref"]
    ref = ref[gap_begin..(-gap_end-1)]
    l1 = l1 + gap_begin
    l2 = l2 + gap_end
    repeat = 0
    if g1 > 2*s1
      ref = ref[g1..-1]
      repeat = 1
      l1 = l1 + g1
    end
    if g2 > 2*s2
      ref = ref[0..(-g2 - 1)]
      repeat = 1
      l2 = l2 + g2
    end
  end
  ref = hxb2_ref[l1..(hxb2_l - l2 - 1)]

  temp_in = File.open(temp_file,"w")
  temp_in.puts ">ref"
  temp_in.puts ref
  temp_in.puts name
  temp_in.puts seq
  temp_in.close
  print `muscle -in #{temp_file} -out #{temp_aln} -quiet`
  aln_seq = fasta_to_hash(temp_aln)
  aln_test = aln_seq[name]
  ref = aln_seq[">ref"]

  #refine alignment

  if ref =~ /^(\-+)/
    l1 = l1 - $1.size
  elsif ref =~ /(\-+)$/
    l2 = l2 + $1.size
  end

  if (hxb2_l - l2 - 1) >= l1
    ref = hxb2_ref[l1..(hxb2_l - l2 - 1)]
    temp_in = File.open(temp_file,"w")
    temp_in.puts ">ref"
    temp_in.puts ref
    temp_in.puts name
    temp_in.puts seq
    temp_in.close
    print `muscle -in #{temp_file} -out #{temp_aln} -quiet`
    aln_seq = fasta_to_hash(temp_aln)
    aln_test = aln_seq[name]
    ref = aln_seq[">ref"]

    ref_size = ref.size
    sim_count = 0
    (0..(ref_size-1)).each do |n|
      ref_base = ref[n]
      test_base = aln_test[n]
      sim_count += 1 if ref_base == test_base
    end
    similarity = (sim_count/ref_size.to_f*100).round(1)
    print `rm -f #{temp_file}`
    print `rm -f #{temp_aln}`
    loc_p1 = l1 + 1
    loc_p2 = hxb2_l - l2
    if seq.size != (loc_p2 - loc_p1 + 1)
        indel = true
    elsif aln_test.include?("-")
        indel = true
    end
    return [loc_p1,loc_p2,similarity,indel,aln_test,ref]
  else
    return [0,0,0,0,0,0,0]
  end
end

#Use NL43 as reference
#align with NL43, return location, percentage of similarity
def NL43_locator(seq="",temp_dir=File.dirname($0))
  hxb2_ref = "TGGAAGGGCTAATTTGGTCCCAAAAAAGACAAGAGATCCTTGATCTGTGGATCTACCACACACAAGGCTACTTCCCTGATTGGCAGAACTACACACCAGGGCCAGGGATCAGATATCCACTGACCTTTGGATGGTGCTTCAAGTTAGTACCAGTTGAACCAGAGCAAGTAGAAGAGGCCAAATAAGGAGAGAAGAACAGCTTGTTACACCCTATGAGCCAGCATGGGATGGAGGACCCGGAGGGAGAAGTATTAGTGTGGAAGTTTGACAGCCTCCTAGCATTTCGTCACATGGCCCGAGAGCTGCATCCGGAGTACTACAAAGACTGCTGACATCGAGCTTTCTACAAGGGACTTTCCGCTGGGGACTTTCCAGGGAGGTGTGGCCTGGGCGGGACTGGGGAGTGGCGAGCCCTCAGATGCTACATATAAGCAGCTGCTTTTTGCCTGTACTGGGTCTCTCTGGTTAGACCAGATCTGAGCCTGGGAGCTCTCTGGCTAACTAGGGAACCCACTGCTTAAGCCTCAATAAAGCTTGCCTTGAGTGCTCAAAGTAGTGTGTGCCCGTCTGTTGTGTGACTCTGGTAACTAGAGATCCCTCAGACCCTTTTAGTCAGTGTGGAAAATCTCTAGCAGTGGCGCCCGAACAGGGACTTGAAAGCGAAAGTAAAGCCAGAGGAGATCTCTCGACGCAGGACTCGGCTTGCTGAAGCGCGCACGGCAAGAGGCGAGGGGCGGCGACTGGTGAGTACGCCAAAAATTTTGACTAGCGGAGGCTAGAAGGAGAGAGATGGGTGCGAGAGCGTCGGTATTAAGCGGGGGAGAATTAGATAAATGGGAAAAAATTCGGTTAAGGCCAGGGGGAAAGAAACAATATAAACTAAAACATATAGTATGGGCAAGCAGGGAGCTAGAACGATTCGCAGTTAATCCTGGCCTTTTAGAGACATCAGAAGGCTGTAGACAAATACTGGGACAGCTACAACCATCCCTTCAGACAGGATCAGAAGAACTTAGATCATTATATAATACAATAGCAGTCCTCTATTGTGTGCATCAAAGGATAGATGTAAAAGACACCAAGGAAGCCTTAGATAAGATAGAGGAAGAGCAAAACAAAAGTAAGAAAAAGGCACAGCAAGCAGCAGCTGACACAGGAAACAACAGCCAGGTCAGCCAAAATTACCCTATAGTGCAGAACCTCCAGGGGCAAATGGTACATCAGGCCATATCACCTAGAACTTTAAATGCATGGGTAAAAGTAGTAGAAGAGAAGGCTTTCAGCCCAGAAGTAATACCCATGTTTTCAGCATTATCAGAAGGAGCCACCCCACAAGATTTAAATACCATGCTAAACACAGTGGGGGGACATCAAGCAGCCATGCAAATGTTAAAAGAGACCATCAATGAGGAAGCTGCAGAATGGGATAGATTGCATCCAGTGCATGCAGGGCCTATTGCACCAGGCCAGATGAGAGAACCAAGGGGAAGTGACATAGCAGGAACTACTAGTACCCTTCAGGAACAAATAGGATGGATGACACATAATCCACCTATCCCAGTAGGAGAAATCTATAAAAGATGGATAATCCTGGGATTAAATAAAATAGTAAGAATGTATAGCCCTACCAGCATTCTGGACATAAGACAAGGACCAAAGGAACCCTTTAGAGACTATGTAGACCGATTCTATAAAACTCTAAGAGCCGAGCAAGCTTCACAAGAGGTAAAAAATTGGATGACAGAAACCTTGTTGGTCCAAAATGCGAACCCAGATTGTAAGACTATTTTAAAAGCATTGGGACCAGGAGCGACACTAGAAGAAATGATGACAGCATGTCAGGGAGTGGGGGGACCCGGCCATAAAGCAAGAGTTTTGGCTGAAGCAATGAGCCAAGTAACAAATCCAGCTACCATAATGATACAGAAAGGCAATTTTAGGAACCAAAGAAAGACTGTTAAGTGTTTCAATTGTGGCAAAGAAGGGCACATAGCCAAAAATTGCAGGGCCCCTAGGAAAAAGGGCTGTTGGAAATGTGGAAAGGAAGGACACCAAATGAAAGATTGTACTGAGAGACAGGCTAATTTTTTAGGGAAGATCTGGCCTTCCCACAAGGGAAGGCCAGGGAATTTTCTTCAGAGCAGACCAGAGCCAACAGCCCCACCAGAAGAGAGCTTCAGGTTTGGGGAAGAGACAACAACTCCCTCTCAGAAGCAGGAGCCGATAGACAAGGAACTGTATCCTTTAGCTTCCCTCAGATCACTCTTTGGCAGCGACCCCTCGTCACAATAAAGATAGGGGGGCAATTAAAGGAAGCTCTATTAGATACAGGAGCAGATGATACAGTATTAGAAGAAATGAATTTGCCAGGAAGATGGAAACCAAAAATGATAGGGGGAATTGGAGGTTTTATCAAAGTAAGACAGTATGATCAGATACTCATAGAAATCTGCGGACATAAAGCTATAGGTACAGTATTAGTAGGACCTACACCTGTCAACATAATTGGAAGAAATCTGTTGACTCAGATTGGCTGCACTTTAAATTTTCCCATTAGTCCTATTGAGACTGTACCAGTAAAATTAAAGCCAGGAATGGATGGCCCAAAAGTTAAACAATGGCCATTGACAGAAGAAAAAATAAAAGCATTAGTAGAAATTTGTACAGAAATGGAAAAGGAAGGAAAAATTTCAAAAATTGGGCCTGAAAATCCATACAATACTCCAGTATTTGCCATAAAGAAAAAAGACAGTACTAAATGGAGAAAATTAGTAGATTTCAGAGAACTTAATAAGAGAACTCAAGATTTCTGGGAAGTTCAATTAGGAATACCACATCCTGCAGGGTTAAAACAGAAAAAATCAGTAACAGTACTGGATGTGGGCGATGCATATTTTTCAGTTCCCTTAGATAAAGACTTCAGGAAGTATACTGCATTTACCATACCTAGTATAAACAATGAGACACCAGGGATTAGATATCAGTACAATGTGCTTCCACAGGGATGGAAAGGATCACCAGCAATATTCCAGTGTAGCATGACAAAAATCTTAGAGCCTTTTAGAAAACAAAATCCAGACATAGTCATCTATCAATACATGGATGATTTGTATGTAGGATCTGACTTAGAAATAGGGCAGCATAGAACAAAAATAGAGGAACTGAGACAACATCTGTTGAGGTGGGGATTTACCACACCAGACAAAAAACATCAGAAAGAACCTCCATTCCTTTGGATGGGTTATGAACTCCATCCTGATAAATGGACAGTACAGCCTATAGTGCTGCCAGAAAAGGACAGCTGGACTGTCAATGACATACAGAAATTAGTGGGAAAATTGAATTGGGCAAGTCAGATTTATGCAGGGATTAAAGTAAGGCAATTATGTAAACTTCTTAGGGGAACCAAAGCACTAACAGAAGTAGTACCACTAACAGAAGAAGCAGAGCTAGAACTGGCAGAAAACAGGGAGATTCTAAAAGAACCGGTACATGGAGTGTATTATGACCCATCAAAAGACTTAATAGCAGAAATACAGAAGCAGGGGCAAGGCCAATGGACATATCAAATTTATCAAGAGCCATTTAAAAATCTGAAAACAGGAAAATATGCAAGAATGAAGGGTGCCCACACTAATGATGTGAAACAATTAACAGAGGCAGTACAAAAAATAGCCACAGAAAGCATAGTAATATGGGGAAAGACTCCTAAATTTAAATTACCCATACAAAAGGAAACATGGGAAGCATGGTGGACAGAGTATTGGCAAGCCACCTGGATTCCTGAGTGGGAGTTTGTCAATACCCCTCCCTTAGTGAAGTTATGGTACCAGTTAGAGAAAGAACCCATAATAGGAGCAGAAACTTTCTATGTAGATGGGGCAGCCAATAGGGAAACTAAATTAGGAAAAGCAGGATATGTAACTGACAGAGGAAGACAAAAAGTTGTCCCCCTAACGGACACAACAAATCAGAAGACTGAGTTACAAGCAATTCATCTAGCTTTGCAGGATTCGGGATTAGAAGTAAACATAGTGACAGACTCACAATATGCATTGGGAATCATTCAAGCACAACCAGATAAGAGTGAATCAGAGTTAGTCAGTCAAATAATAGAGCAGTTAATAAAAAAGGAAAAAGTCTACCTGGCATGGGTACCAGCACACAAAGGAATTGGAGGAAATGAACAAGTAGATGGGTTGGTCAGTGCTGGAATCAGGAAAGTACTATTTTTAGATGGAATAGATAAGGCCCAAGAAGAACATGAGAAATATCACAGTAATTGGAGAGCAATGGCTAGTGATTTTAACCTACCACCTGTAGTAGCAAAAGAAATAGTAGCCAGCTGTGATAAATGTCAGCTAAAAGGGGAAGCCATGCATGGACAAGTAGACTGTAGCCCAGGAATATGGCAGCTAGATTGTACACATTTAGAAGGAAAAGTTATCTTGGTAGCAGTTCATGTAGCCAGTGGATATATAGAAGCAGAAGTAATTCCAGCAGAGACAGGGCAAGAAACAGCATACTTCCTCTTAAAATTAGCAGGAAGATGGCCAGTAAAAACAGTACATACAGACAATGGCAGCAATTTCACCAGTACTACAGTTAAGGCCGCCTGTTGGTGGGCGGGGATCAAGCAGGAATTTGGCATTCCCTACAATCCCCAAAGTCAAGGAGTAATAGAATCTATGAATAAAGAATTAAAGAAAATTATAGGACAGGTAAGAGATCAGGCTGAACATCTTAAGACAGCAGTACAAATGGCAGTATTCATCCACAATTTTAAAAGAAAAGGGGGGATTGGGGGGTACAGTGCAGGGGAAAGAATAGTAGACATAATAGCAACAGACATACAAACTAAAGAATTACAAAAACAAATTACAAAAATTCAAAATTTTCGGGTTTATTACAGGGACAGCAGAGATCCAGTTTGGAAAGGACCAGCAAAGCTCCTCTGGAAAGGTGAAGGGGCAGTAGTAATACAAGATAATAGTGACATAAAAGTAGTGCCAAGAAGAAAAGCAAAGATCATCAGGGATTATGGAAAACAGATGGCAGGTGATGATTGTGTGGCAAGTAGACAGGATGAGGATTAACACATGGAAAAGATTAGTAAAACACCATATGTATATTTCAAGGAAAGCTAAGGACTGGTTTTATAGACATCACTATGAAAGTACTAATCCAAAAATAAGTTCAGAAGTACACATCCCACTAGGGGATGCTAAATTAGTAATAACAACATATTGGGGTCTGCATACAGGAGAAAGAGACTGGCATTTGGGTCAGGGAGTCTCCATAGAATGGAGGAAAAAGAGATATAGCACACAAGTAGACCCTGACCTAGCAGACCAACTAATTCATCTGCACTATTTTGATTGTTTTTCAGAATCTGCTATAAGAAATACCATATTAGGACGTATAGTTAGTCCTAGGTGTGAATATCAAGCAGGACATAACAAGGTAGGATCTCTACAGTACTTGGCACTAGCAGCATTAATAAAACCAAAACAGATAAAGCCACCTTTGCCTAGTGTTAGGAAACTGACAGAGGACAGATGGAACAAGCCCCAGAAGACCAAGGGCCACAGAGGGAGCCATACAATGAATGGACACTAGAGCTTTTAGAGGAACTTAAGAGTGAAGCTGTTAGACATTTTCCTAGGATATGGCTCCATAACTTAGGACAACATATCTATGAAACTTACGGGGATACTTGGGCAGGAGTGGAAGCCATAATAAGAATTCTGCAACAACTGCTGTTTATCCATTTCAGAATTGGGTGTCGACATAGCAGAATAGGCGTTACTCGACAGAGGAGAGCAAGAAATGGAGCCAGTAGATCCTAGACTAGAGCCCTGGAAGCATCCAGGAAGTCAGCCTAAAACTGCTTGTACCAATTGCTATTGTAAAAAGTGTTGCTTTCATTGCCAAGTTTGTTTCATGACAAAAGCCTTAGGCATCTCCTATGGCAGGAAGAAGCGGAGACAGCGACGAAGAGCTCATCAGAACAGTCAGACTCATCAAGCTTCTCTATCAAAGCAGTAAGTAGTACATGTAATGCAACCTATAATAGTAGCAATAGTAGCATTAGTAGTAGCAATAATAATAGCAATAGTTGTGTGGTCCATAGTAATCATAGAATATAGGAAAATATTAAGACAAAGAAAAATAGACAGGTTAATTGATAGACTAATAGAAAGAGCAGAAGACAGTGGCAATGAGAGTGAAGGAGAAGTATCAGCACTTGTGGAGATGGGGGTGGAAATGGGGCACCATGCTCCTTGGGATATTGATGATCTGTAGTGCTACAGAAAAATTGTGGGTCACAGTCTATTATGGGGTACCTGTGTGGAAGGAAGCAACCACCACTCTATTTTGTGCATCAGATGCTAAAGCATATGATACAGAGGTACATAATGTTTGGGCCACACATGCCTGTGTACCCACAGACCCCAACCCACAAGAAGTAGTATTGGTAAATGTGACAGAAAATTTTAACATGTGGAAAAATGACATGGTAGAACAGATGCATGAGGATATAATCAGTTTATGGGATCAAAGCCTAAAGCCATGTGTAAAATTAACCCCACTCTGTGTTAGTTTAAAGTGCACTGATTTGAAGAATGATACTAATACCAATAGTAGTAGCGGGAGAATGATAATGGAGAAAGGAGAGATAAAAAACTGCTCTTTCAATATCAGCACAAGCATAAGAGATAAGGTGCAGAAAGAATATGCATTCTTTTATAAACTTGATATAGTACCAATAGATAATACCAGCTATAGGTTGATAAGTTGTAACACCTCAGTCATTACACAGGCCTGTCCAAAGGTATCCTTTGAGCCAATTCCCATACATTATTGTGCCCCGGCTGGTTTTGCGATTCTAAAATGTAATAATAAGACGTTCAATGGAACAGGACCATGTACAAATGTCAGCACAGTACAATGTACACATGGAATCAGGCCAGTAGTATCAACTCAACTGCTGTTAAATGGCAGTCTAGCAGAAGAAGATGTAGTAATTAGATCTGCCAATTTCACAGACAATGCTAAAACCATAATAGTACAGCTGAACACATCTGTAGAAATTAATTGTACAAGACCCAACAACAATACAAGAAAAAGTATCCGTATCCAGAGGGGACCAGGGAGAGCATTTGTTACAATAGGAAAAATAGGAAATATGAGACAAGCACATTGTAACATTAGTAGAGCAAAATGGAATGCCACTTTAAAACAGATAGCTAGCAAATTAAGAGAACAATTTGGAAATAATAAAACAATAATCTTTAAGCAATCCTCAGGAGGGGACCCAGAAATTGTAACGCACAGTTTTAATTGTGGAGGGGAATTTTTCTACTGTAATTCAACACAACTGTTTAATAGTACTTGGTTTAATAGTACTTGGAGTACTGAAGGGTCAAATAACACTGAAGGAAGTGACACAATCACACTCCCATGCAGAATAAAACAATTTATAAACATGTGGCAGGAAGTAGGAAAAGCAATGTATGCCCCTCCCATCAGTGGACAAATTAGATGTTCATCAAATATTACTGGGCTGCTATTAACAAGAGATGGTGGTAATAACAACAATGGGTCCGAGATCTTCAGACCTGGAGGAGGCGATATGAGGGACAATTGGAGAAGTGAATTATATAAATATAAAGTAGTAAAAATTGAACCATTAGGAGTAGCACCCACCAAGGCAAAGAGAAGAGTGGTGCAGAGAGAAAAAAGAGCAGTGGGAATAGGAGCTTTGTTCCTTGGGTTCTTGGGAGCAGCAGGAAGCACTATGGGCTGCACGTCAATGACGCTGACGGTACAGGCCAGACAATTATTGTCTGATATAGTGCAGCAGCAGAACAATTTGCTGAGGGCTATTGAGGCGCAACAGCATCTGTTGCAACTCACAGTCTGGGGCATCAAACAGCTCCAGGCAAGAATCCTGGCTGTGGAAAGATACCTAAAGGATCAACAGCTCCTGGGGATTTGGGGTTGCTCTGGAAAACTCATTTGCACCACTGCTGTGCCTTGGAATGCTAGTTGGAGTAATAAATCTCTGGAACAGATTTGGAATAACATGACCTGGATGGAGTGGGACAGAGAAATTAACAATTACACAAGCTTAATACACTCCTTAATTGAAGAATCGCAAAACCAGCAAGAAAAGAATGAACAAGAATTATTGGAATTAGATAAATGGGCAAGTTTGTGGAATTGGTTTAACATAACAAATTGGCTGTGGTATATAAAATTATTCATAATGATAGTAGGAGGCTTGGTAGGTTTAAGAATAGTTTTTGCTGTACTTTCTATAGTGAATAGAGTTAGGCAGGGATATTCACCATTATCGTTTCAGACCCACCTCCCAATCCCGAGGGGACCCGACAGGCCCGAAGGAATAGAAGAAGAAGGTGGAGAGAGAGACAGAGACAGATCCATTCGATTAGTGAACGGATCCTTAGCACTTATCTGGGACGATCTGCGGAGCCTGTGCCTCTTCAGCTACCACCGCTTGAGAGACTTACTCTTGATTGTAACGAGGATTGTGGAACTTCTGGGACGCAGGGGGTGGGAAGCCCTCAAATATTGGTGGAATCTCCTACAGTATTGGAGTCAGGAACTAAAGAATAGTGCTGTTAACTTGCTCAATGCCACAGCCATAGCAGTAGCTGAGGGGACAGATAGGGTTATAGAAGTATTACAAGCAGCTTATAGAGCTATTCGCCACATACCTAGAAGAATAAGACAGGGCTTGGAAAGGATTTTGCTATAAGATGGGTGGCAAGTGGTCAAAAAGTAGTGTGATTGGATGGCCTGCTGTAAGGGAAAGAATGAGACGAGCTGAGCCAGCAGCAGATGGGGTGGGAGCAGTATCTCGAGACCTAGAAAAACATGGAGCAATCACAAGTAGCAATACAGCAGCTAACAATGCTGCTTGTGCCTGGCTAGAAGCACAAGAGGAGGAAGAGGTGGGTTTTCCAGTCACACCTCAGGTACCTTTAAGACCAATGACTTACAAGGCAGCTGTAGATCTTAGCCACTTTTTAAAAGAAAAGGGGGGACTGGAAGGGCTAATTCACTCCCAAAGAAGACAAGATATCCTTGATCTGTGGATCTACCACACACAAGGCTACTTCCCTGATTGGCAGAACTACACACCAGGGCCAGGGGTCAGATATCCACTGACCTTTGGATGGTGCTACAAGCTAGTACCAGTTGAGCCAGATAAGGTAGAAGAGGCCAATAAAGGAGAGAACACCAGCTTGTTACACCCTGTGAGCCTGCATGGAATGGATGACCCTGAGAGAGAAGTGTTAGAGTGGAGGTTTGACAGCCGCCTAGCATTTCATCACGTGGCCCGAGAGCTGCATCCGGAGTACTTCAAGAACTGCTGACATCGAGCTTGCTACAAGGGACTTTCCGCTGGGGACTTTCCAGGGAGGCGTGGCCTGGGCGGGACTGGGGAGTGGCGAGCCCTCAGATGCTGCATATAAGCAGCTGCTTTTTGCCTGTACTGGGTCTCTCTGGTTAGACCAGATCTGAGCCTGGGAGCTCTCTGGCTAACTAGGGAACCCACTGCTTAAGCCTCAATAAAGCTTGCCTTGAGTGCTTCAAGTAGTGTGTGCCCGTCTGTTGTGTGACTCTGGTAACTAGAGATCCCTCAGACCCTTTTAGTCAGTGTGGAAAATCTCTAGCACCCAGGAGGTAGAGGTTGCAGTGAGCCAAGATCGCGCCACTGCATTCCAGCCTGGGCAAGAAAACAAGACTGTCTAAAATAATAATAATAAGTTAAGGGTATTAAATATATTTATACATGGAGGTCATAAAAATATATATATTTGGGCTGGGCGCAGTGGCTCACACCTGCGCCCGGCCCTTTGGGAGGCCGAGGCAGGTGGATCACCTGAGTTTGGGAGTTCCAGACCAGCCTGACCAACATGGAGAAACCCCTTCTCTGTGTATTTTTAGTAGATTTTATTTTATGTGTATTTTATTCACAGGTATTTCTGGAAAACTGAAACTGTTTTTCCTCTACTCTGATACCACAAGAATCATCAGCACAGAGGAAGACTTCTGTGATCAAATGTGGTGGGAGAGGGAGGTTTTCACCAGCACATGAGCAGTCAGTTCTGCCGCAGACTCGGCGGGTGTCCTTCGGTTCAGTTCCAACACCGCCTGCCTGGAGAGAGGTCAGACCACAGGGTGAGGGCTCAGTCCCCAAGACATAAACACCCAAGACATAAACACCCAACAGGTCCACCCCGCCTGCTGCCCAGGCAGAGCCGATTCACCAAGACGGGAATTAGGATAGAGAAAGAGTAAGTCACACAGAGCCGGCTGTGCGGGAGAACGGAGTTCTATTATGACTCAAATCAGTCTCCCCAAGCATTCGGGGATCAGAGTTTTTAAGGATAACTTAGTGTGTAGGGGGCCAGTGAGTTGGAGATGAAAGCGTAGGGAGTCGAAGGTGTCCTTTTGCGCCGAGTCAGTTCCTGGGTGGGGGCCACAAGATCGGATGAGCCAGTTTATCAATCCGGGGGTGCCAGCTGATCCATGGAGTGCAGGGTCTGCAAAATATCTCAAGCACTGATTGATCTTAGGTTTTACAATAGTGATGTTACCCCAGGAACAATTTGGGGAAGGTCAGAATCTTGTAGCCTGTAGCTGCATGACTCCTAAACCATAATTTCTTTTTTGTTTTTTTTTTTTTATTTTTGAGACAGGGTCTCACTCTGTCACCTAGGCTGGAGTGCAGTGGTGCAATCACAGCTCACTGCAGCCTCAACGTCGTAAGCTCAAGCGATCCTCCCACCTCAGCCTGCCTGGTAGCTGAGACTACAAGCGACGCCCCAGTTAATTTTTGTATTTTTGGTAGAGGCAGCGTTTTGCCGTGTGGCCCTGGCTGGTCTCGAACTCCTGGGCTCAAGTGATCCAGCCTCAGCCTCCCAAAGTGCTGGGACAACCGGGCCCAGTCACTGCACCTGGCCCTAAACCATAATTTCTAATCTTTTGGCTAATTTGTTAGTCCTACAAAGGCAGTCTAGTCCCCAGCAAAAAGGGGGTTTGTTTCGGGAAAGGGCTGTTACTGTCTTTGTTTCAAACTATAAACTAAGTTCCTCCTAAACTTAGTTCGGCCTACACCCAGGAATGAACAAGGAGAGCTTGGAGGTTAGAAGCACGATGGAATTGGTTAGGTCAGATCTCTTTCACTGTCTGAGTTATAATTTTGCAATGGTGGTTCAAAGACTGCCCGCTTCTGACACCAGTCGCTGCATTAATGAATCGGCCAACGCGCGGGGAGAGGCGGTTTGCGTATTGGGCGCTCTTCCGCTTCCTCGCTCACTGACTCGCTGCGCTCGGTCGTTCGGCTGCGGCGAGCGGTATCAGCTCACTCAAAGGCGGTAATACGGTTATCCACAGAATCAGGGGATAACGCAGGAAAGAACATGTGAGCAAAAGGCCAGCAAAAGGCCAGGAACCGTAAAAAGGCCGCGTTGCTGGCGTTTTTCCATAGGCTCCGCCCCCCTGACGAGCATCACAAAAATCGACGCTCAAGTCAGAGGTGGCGAAACCCGACAGGACTATAAAGATACCAGGCGTTTCCCCCTGGAAGCTCCCTCGTGCGCTCTCCTGTTCCGACCCTGCCGCTTACCGGATACCTGTCCGCCTTTCTCCCTTCGGGAAGCGTGGCGCTTTCTCATAGCTCACGCTGTAGGTATCTCAGTTCGGTGTAGGTCGTTCGCTCCAAGCTGGGCTGTGTGCACGAACCCCCCGTTCAGCCCGACCGCTGCGCCTTATCCGGTAACTATCGTCTTGAGTCCAACCCGGTAAGACACGACTTATCGCCACTGGCAGCAGCCACTGGTAACAGGATTAGCAGAGCGAGGTATGTAGGCGGTGCTACAGAGTTCTTGAAGTGGTGGCCTAACTACGGCTACACTAGAAGGACAGTATTTGGTATCTGCGCTCTGCTGAAGCCAGTTACCTTCGGAAAAAGAGTTGGTAGCTCTTGATCCGGCAAACAAACCACCGCTGGTAGCGGTGGTTTTTTTGTTTGCAAGCAGCAGATTACGCGCAGAAAAAAAGGATCTCAAGAAGATCCTTTGATCTTTTCTACGGGGTCTGACGCTCAGTGGAACGAAAACTCACGTTAAGGGATTTTGGTCATGAGATTATCAAAAAGGATCTTCACCTAGATCCTTTTAAATTAAAAATGAAGTTTTAAATCAATCTAAAGTATATATGAGTAAACTTGGTCTGACAGTTACCAATGCTTAATCAGTGAGGCACCTATCTCAGCGATCTGTCTATTTCGTTCATCCATAGTTGCCTGACTCCCCGTCGTGTAGATAACTACGATACGGGAGGGCTTACCATCTGGCCCCAGTGCTGCAATGATACCGCGAGACCCACGCTCACCGGCTCCAGATTTATCAGCAATAAACCAGCCAGCCGGAAGGGCCGAGCGCAGAAGTGGTCCTGCAACTTTATCCGCCTCCATCCAGTCTATTAATTGTTGCCGGGAAGCTAGAGTAAGTAGTTCGCCAGTTAATAGTTTGCGCAACGTTGTTGCCATTGCTACAGGCATCGTGGTGTCACGCTCGTCGTTTGGTATGGCTTCATTCAGCTCCGGTTCCCAACGATCAAGGCGAGTTACATGATCCCCCATGTTGTGCAAAAAAGCGGTTAGCTCCTTCGGTCCTCCGATCGTTGTCAGAAGTAAGTTGGCCGCAGTGTTATCACTCATGGTTATGGCAGCACTGCATAATTCTCTTACTGTCATGCCATCCGTAAGATGCTTTTCTGTGACTGGTGAGTACTCAACCAAGTCATTCTGAGAATAGTGTATGCGGCGACCGAGTTGCTCTTGCCCGGCGTCAATACGGGATAATACCGCGCCACATAGCAGAACTTTAAAAGTGCTCATCATTGGAAAACGTTCTTCGGGGCGAAAACTCTCAAGGATCTTACCGCTGTTGAGATCCAGTTCGATGTAACCCACTCGTGCACCCAACTGATCTTCAGCATCTTTTACTTTCACCAGCGTTTCTGGGTGAGCAAAAACAGGAAGGCAAAATGCCGCAAAAAAGGGAATAAGGGCGACACGGAAATGTTGAATACTCATACTCTTCCTTTTTCAATATTATTGAAGCATTTATCAGGGTTATTGTCTCATGAGCGGATACATATTTGAATGTATTTAGAAAAATAAACAAATAGGGGTTCCGCGCACATTTCCCCGAAAAGTGCCACCTGACGTCTAAGAAACCATTATTATCATGACATTAACCTATAAAAATAGGCGTATCACGAGGCCCTTTCGTCTCGCGCGTTTCGGTGATGACGGTGAAAACCTCTGACACATGCAGCTCCCGGAGACGGTCACAGCTTGTCTGTAAGCGGATGCCGGGAGCAGACAAGCCCGTCAGGGCGCGTCAGCGGGTGTTGGCGGGTGTCGGGGCTGGCTTAACTATGCGGCATCAGAGCAGATTGTACTGAGAGTGCACCATATGCGGTGTGAAATACCGCACAGATGCGTAAGGAGAAAATACCGCATCAGGCGCCATTCGCCATTCAGGCTGCGCAACTGTTGGGAAGGGCGATCGGTGCGGGCCTCTTCGCTATTACGCCAGGGGAGGCAGAGATTGCAGTAAGCTGAGATCGCAGCACTGCACTCCAGCCTGGGCGACAGAGTAAGACTCTGTCTCAAAAATAAAATAAATAAATCAATCAGATATTCCAATCTTTTCCTTTATTTATTTATTTATTTTCTATTTTGGAAACACAGTCCTTCCTTATTCCAGAATTACACATATATTCTATTTTTCTTTATATGCTCCAGTTTTTTTTAGACCTTCACCTGAAATGTGTGTATACAAAATCTAGGCCAGTCCAGCAGAGCCTAAAGGTAAAAAATAAAATAATAAAAAATAAATAAAATCTAGCTCACTCCTTCACATCAAAATGGAGATACAGCTGTTAGCATTAAATACCAAATAACCCATCTTGTCCTCAATAATTTTAAGCGCCTCTCTCCACCACATCTAACTCCTGTCAAAGGCATGTGCCCCTTCCGGGCGCTCTGCTGTGCTGCCAACCAACTGGCATGTGGACTCTGCAGGGTCCCTAACTGCCAAGCCCCACAGTGTGCCCTGAGGCTGCCCCTTCCTTCTAGCGGCTGCCCCCACTCGGCTTTGCTTTCCCTAGTTTCAGTTACTTGCGTTCAGCCAAGGTCTGAAACTAGGTGCGCACAGAGCGGTAAGACTGCGAGAGAAAGAGACCAGCTTTACAGGGGGTTTATCACAGTGCACCCTGACAGTCGTCAGCCTCACAGGGGGTTTATCACATTGCACCCTGACAGTCGTCAGCCTCACAGGGGGTTTATCACAGTGCACCCTTACAATCATTCCATTTGATTCACAATTTTTTTAGTCTCTACTGTGCCTAACTTGTAAGTTAAATTTGATCAGAGGTGTGTTCCCAGAGGGGAAAACAGTATATACAGGGTTCAGTACTATCGCATTTCAGGCCTCCACCTGGGTCTTGGAATGTGTCCCCCGAGGGGTGATGACTACCTCAGTTGGATCTCCACAGGTCACAGTGACACAAGATAACCAAGACACCTCCCAAGGCTACCACAATGGGCCGCCCTCCACGTGCACATGGCCGGAGGAACTGCCATGTCGGAGGTGCAAGCACACCTGCGCATCAGAGTCCTTGGTGTGGAGGGAGGGACCAGCGCAGCTTCCAGCCATCCACCTGATGAACAGAACCTAGGGAAAGCCCCAGTTCTACTTACACCAGGAAAGGC"
  hxb2_l = hxb2_ref.size
  head = ""
  8.times {head  << (65 + rand(25)).chr}
  temp_file = temp_dir + "/temp"
  temp_aln = temp_dir + "/temp_aln"

  l1 = 0
  l2 = 0
  name = ">test"
  temp_in = File.open(temp_file,"w")
  temp_in.puts ">ref"
  temp_in.puts hxb2_ref
  temp_in.puts name
  temp_in.puts seq
  temp_in.close

  print `muscle -in #{temp_file} -out #{temp_aln} -quiet`
  aln_seq = fasta_to_hash(temp_aln)
  aln_test = aln_seq[name]
  aln_test =~ /^(\-*)(\w.*\w)(\-*)$/
  gap_begin = $1.size
  gap_end = $3.size
  aln_test2 = $2
  ref = aln_seq[">ref"]
  ref = ref[gap_begin..(-gap_end-1)]
  ref_size = ref.size
  if ref_size > 1.3*(seq.size)
    l1 = l1 + gap_begin
    l2 = l2 + gap_end
    max_seq = aln_test2.scan(/[ACGT]+/).max_by(&:length)
    aln_test2 =~ /#{max_seq}/
    before_aln_seq = $`
    before_aln = $`.size
    post_aln_seq = $'
    post_aln = $'.size
    before_aln_seq_size = before_aln_seq.scan(/[ACGT]+/).join("").size
    b1 = (1.3 * before_aln_seq_size).to_i
    post_aln_seq_size = post_aln_seq.scan(/[ACGT]+/).join("").size
    b2 = (1.3 * post_aln_seq_size).to_i
    if (before_aln > seq.size) and (post_aln <= seq.size)
      ref = ref[(before_aln - b1)..(ref_size - post_aln - 1)]
      l1 = l1 + (before_aln - b1)
    elsif (post_aln > seq.size) and (before_aln <= seq.size)
      ref = ref[before_aln..(ref_size - post_aln - 1 + b2)]
      l2 = l2 + post_aln - b2
    elsif (post_aln > seq.size) and (before_aln > seq.size)
      ref = ref[(before_aln - b1)..(ref_size - post_aln - 1 + b2)]
      l1 = l1 + (before_aln - b1)
      l2 = l2 + (post_aln - b2)
    end
    temp_in = File.open(temp_file,"w")
    temp_in.puts ">ref"
    temp_in.puts ref
    temp_in.puts name
    temp_in.puts seq
    temp_in.close
    print `muscle -in #{temp_file} -out #{temp_aln} -quiet`
    aln_seq = fasta_to_hash(temp_aln)
    aln_test = aln_seq[name]
    aln_test =~ /^(\-*)(\w.*\w)(\-*)$/
    gap_begin = $1.size
    gap_end = $3.size
    aln_test2 = $2
    ref = aln_seq[">ref"]
    ref = ref[gap_begin..(-gap_end-1)]
    ref_size = ref.size
  end
  aln_seq = fasta_to_hash(temp_aln)
  aln_test = aln_seq[name]
  aln_test =~ /^(\-*)(\w.*\w)(\-*)$/
  gap_begin = $1.size
  gap_end = $3.size
  aln_test = $2
  aln_test =~ /^(\w+)(\-*)\w/
  s1 = $1.size
  g1 = $2.size
  aln_test =~ /\w(\-*)(\w+)$/
  s2 = $2.size
  g2 = $1.size
  ref = aln_seq[">ref"]
  ref = ref[gap_begin..(-gap_end-1)]

  l1 = l1 + gap_begin
  l2 = l2 + gap_end
  repeat = 0

  if g1 == g2 and (s1 + g1 + s2) == ref.size
    if s1 > s2 and g2 > 2*s2
      ref = ref[0..(-g2-1)]
      repeat = 1
      l2 = l2 + g2
    elsif s1 < s2 and g1 > 2*s1
      ref = ref[g1..-1]
      repeat = 1
      l1 = l1 + g1
    end
  else
    if g1 > 2*s1
      ref = ref[g1..-1]
      repeat = 1
      l1 = l1 + g1
    end
    if g2 > 2*s2
      ref = ref[0..(-g2 - 1)]
      repeat = 1
      l2 = l2 + g2
    end
  end

  while repeat == 1
    temp_in = File.open(temp_file,"w")
    temp_in.puts ">ref"
    temp_in.puts ref
    temp_in.puts name
    temp_in.puts seq
    temp_in.close
    print `muscle -in #{temp_file} -out #{temp_aln} -quiet`
    aln_seq = fasta_to_hash(temp_aln)
    aln_test = aln_seq[name]
    aln_test =~ /^(\-*)(\w.*\w)(\-*)$/
    gap_begin = $1.size
    gap_end = $3.size
    aln_test = $2
    aln_test =~ /^(\w+)(\-*)\w/
    s1 = $1.size
    g1 = $2.size
    aln_test =~ /\w(\-*)(\w+)$/
    s2 = $2.size
    g2 = $1.size
    ref = aln_seq[">ref"]
    ref = ref[gap_begin..(-gap_end-1)]
    l1 = l1 + gap_begin
    l2 = l2 + gap_end
    repeat = 0
    if g1 > 2*s1
      ref = ref[g1..-1]
      repeat = 1
      l1 = l1 + g1
    end
    if g2 > 2*s2
      ref = ref[0..(-g2 - 1)]
      repeat = 1
      l2 = l2 + g2
    end
  end
  ref = hxb2_ref[l1..(hxb2_l - l2 - 1)]

  temp_in = File.open(temp_file,"w")
  temp_in.puts ">ref"
  temp_in.puts ref
  temp_in.puts name
  temp_in.puts seq
  temp_in.close
  print `muscle -in #{temp_file} -out #{temp_aln} -quiet`
  aln_seq = fasta_to_hash(temp_aln)
  aln_test = aln_seq[name]
  ref = aln_seq[">ref"]

  #refine alignment

  if ref =~ /^(\-+)/
    l1 = l1 - $1.size
  elsif ref =~ /(\-+)$/
    l2 = l2 + $1.size
  end
  l1 = 0 if l1 < 0
  if (hxb2_l - l2 - 1) >= l1
    ref = hxb2_ref[l1..(hxb2_l - l2 - 1)]
    temp_in = File.open(temp_file,"w")
    temp_in.puts ">ref"
    temp_in.puts ref
    temp_in.puts name
    temp_in.puts seq
    temp_in.close
    print `muscle -in #{temp_file} -out #{temp_aln} -quiet`
    aln_seq = fasta_to_hash(temp_aln)
    aln_test = aln_seq[name]
    ref = aln_seq[">ref"]

    ref_size = ref.size
    sim_count = 0
    (0..(ref_size-1)).each do |n|
      ref_base = ref[n]
      test_base = aln_test[n]
      sim_count += 1 if ref_base == test_base
    end
    similarity = (sim_count/ref_size.to_f*100).round(1)
    print `rm -f #{temp_file}`
    print `rm -f #{temp_aln}`
    loc_p1 = l1 + 1
    loc_p2 = hxb2_l - l2
    if seq.size != (loc_p2 - loc_p1 + 1)
        indel = true
    elsif aln_test.include?("-")
        indel = true
    end
    return [loc_p1,loc_p2,similarity,indel,aln_test,ref]
  else
    return [0,0,0,0,0,0,0]
  end
rescue
  return [0,0,0,0,"N","N"]
end

#gene cutter. given a specific HXB2 position array [first_position,last_position], return a sequence within the range
#return nil if seqeuence not in the range.
#require sequence_locator.
def sequence_clip(loc={},positions=[])
    p1 = positions[0]
    p2 = positions[1]
    l1 = loc[0]
    l2 = loc[1]
    if (p1 >= l1) & (p2 <= l2)
        seq = loc[4]
        ref = loc[5]
        g1 = 0
        ref.each_char do |char|
            break if l1 == p1
            g1 += 1
            l1 += 1 unless char == "-"
        end
        g2 = 1
        ref.reverse.each_char do |char|
            break if l2 == p2
            g2 += 1
            l2 -= 1 unless char == "-"
        end
        return seq[g1..(-g2)].tr("-","")
    else
        return nil
    end
end

#Tail function for file.
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

#Input sequence array. output Variant distribution for Poisson cut-off
def variant_for_poisson(seq)
  seq_size = seq.size
  l = seq[0].size - 1
  var = []
  (0..l).to_a.each do |pos|
    nt = []
    seq.each do |s|
      nt << s[pos]
    end
    count_nt = count(nt)
    v = seq_size - count_nt.values.max
    var << v
  end
  var_count = count(var)
  var_count.sort_by{|key,value|key}.to_h
end

#align test sequence with a reference sequence, return aligned test sequence
def muscle_sequence(ref_seq = "", test_seq = "", temp_dir=File.dirname($0))
  temp_file = temp_dir + "/temp"
  temp_aln = temp_dir + "/temp_aln"
  name = ">test"
  temp_in = File.open(temp_file,"w")
  temp_in.puts ">ref"
  temp_in.puts ref_seq
  temp_in.puts name
  temp_in.puts test_seq
  temp_in.close
  print `muscle -in #{temp_file} -out #{temp_aln} -quiet`
  aln_seq = fasta_to_hash(temp_aln)[">test"]
  File.unlink(temp_file)
  File.unlink(temp_aln)
  return aln_seq
end

#align test sequence with a reference sequence, return aligned both sequences
def muscle_sequence2(ref_seq = "", test_seq = "", temp_dir=File.dirname($0))
  temp_file = temp_dir + "/temp"
  temp_aln = temp_dir + "/temp_aln"
  name = ">test"
  temp_in = File.open(temp_file,"w")
  temp_in.puts ">ref"
  temp_in.puts ref_seq
  temp_in.puts name
  temp_in.puts test_seq
  temp_in.close
  print `muscle -in #{temp_file} -out #{temp_aln} -quiet`
  aln_seq = fasta_to_hash(temp_aln)[">test"]
  aln_ref = fasta_to_hash(temp_aln)[">ref"]
  File.unlink(temp_file)
  File.unlink(temp_aln)
  return [aln_ref, aln_seq]
end

#collapse sequences with x number of nt differences. make sure sequences are aligned. The return frequency is NOT the frequency of the collasped sequences.
def collapse_sequence_by_x_nt_difference(seq_array,cutoff)
    new_seq_freq = {}
    seq_freq = count(seq_array)
    if seq_freq.size == 1
        new_seq_freq = seq_freq
    else
        uniq_seq = seq_freq.keys
        unique_seq_pair = uniq_seq.combination(2)
        dupli_seq = []
        unique_seq_pair.each do |pair|
            seq1 = pair[0]
            seq2 = pair[1]
            diff = compare_two_seq(seq1,seq2)
            if diff <= cutoff
                freq1 = seq_freq[seq1]
                freq2 = seq_freq[seq2]
                freq1 >= freq2 ? dupli_seq << seq2 : dupli_seq << seq1
            end
        end

        seq_freq.each do |seq,freq|
            unless dupli_seq.include?(seq)
                new_seq_freq[seq] = freq
            end
        end
        return new_seq_freq
    end
end

#collapse sequence hash to unique sequence hash
def uniq_sequence(seq = {}, sequence_name = "sequence")
  uni = count(seq.values)
  new_seq = {}
  n = 1
  uni.each do |s,c|
    name = ">" + sequence_name + "_" + n.to_s + "_" + c.to_s
    new_seq[name] = s
    n += 1
  end
  return new_seq
end

#calculate pairwise diverisity pi from a sequence hash
def nucleotide_pi(seq = {})
  seq_length = seq.values[0].size - 1
  nt_position_hash = {}
  (0..seq_length).each do |n|
    nt_position_hash[n] = []
    seq.values.each do |s|
      nt_position_hash[n] << s[n]
    end
  end
  diver = 0
  com = 0
  nt_position_hash.each do |p,nt|
    nt.delete("-")
    next if nt.size == 1
    nt_count = count(nt)
    combination = (nt.size)*(nt.size - 1)/2
    com += combination
    a = nt_count["A"]
    c = nt_count["C"]
    t = nt_count["T"]
    g = nt_count["G"]
    div = a*c + a*t + a*g + c*t + c*g + t*g
    diver += div
  end
  pi = (diver/com.to_f).round(5)
  return pi
end

def nucleotide_array_pi(seq = [])
  seq_length = seq[0].size - 1
  nt_position_hash = {}
  (0..seq_length).each do |n|
    nt_position_hash[n] = []
    seq.each do |s|
      nt_position_hash[n] << s[n]
    end
  end
  diver = 0
  com = 0
  nt_position_hash.each do |p,nt|
    nt.delete("-")
    next if nt.size == 1
    nt_count = count(nt)
    combination = (nt.size)*(nt.size - 1)/2
    com += combination
    a = nt_count["A"]
    c = nt_count["C"]
    t = nt_count["T"]
    g = nt_count["G"]
    div = a*c + a*t + a*g + c*t + c*g + t*g
    diver += div
  end
  pi = (diver/com.to_f).round(5)
  return pi
end

def a1(n = 2)
  p = 0
  q = 1
  while q < n
    p += 1/q.to_f
    q += 1
  end
  return p
end


#primer with ambiguities to match
def primer_match (primer = "")
  match = ""
  primer.each_char.each do |base|
    base_array = to_list(base)
    if base_array.size == 1
      match += base_array[0]
    else
      pattern = "[" + base_array.join("|") + "]"
      match += pattern
    end
  end
  return match
end

def filter_similar_pid(sequence_file = "", cutoff = 10)
  seq = fasta_to_hash(sequence_file)
  uni_seq = seq.values.uniq
  uni_seq_pid = {}
  uni_seq.each do |k|
    seq.each do |name,s|
      name = name[1..-1]
      if k == s
        if uni_seq_pid[k]
          uni_seq_pid[k] << [name.split("_")[0],name.split("_")[1]]
        else
          uni_seq_pid[k] = []
          uni_seq_pid[k] << [name.split("_")[0],name.split("_")[1]]
        end
      end
    end
  end

  dup_pid = []
  uni_seq_pid.values.each do |v|
    next if v.size == 1
    pid_hash = Hash[v]
    list = pid_hash.keys
    list2 = Array.new(list)
    pairs = []

    list.each do |k|
      list2.delete(k)
      list2.each do |k1|
        pairs << [k,k1]
      end
    end


    pairs.each do |p|
      pid1 = p[0]
      pid2 = p[1]
      if two_pid_x_base_different(pid1,pid2,1)
        n1 = pid_hash[pid1].to_i
        n2 = pid_hash[pid2].to_i
        if n1 >= cutoff * n2
          dup_pid << pid2
          #puts pid1 + "\t" + n1.to_s + "\t" + pid2 + "\t" + n2.to_s
        elsif n2 >= cutoff * n1
          dup_pid << pid1
          #puts pid2 + "\t" + n2.to_s + "\t" + pid1 + "\t" + n1.to_s
        end
      end
    end
  end


  new_seq = {}
  seq.each do |name,s|
    pid = name.split("_")[0][1..-1]
    unless dup_pid.include?(pid)
      new_seq[name] = s
    end
  end
  return new_seq
end

# APOBEC3G/F mutation position identification
# APOBEC3G/F pattern: GRD -> ARD
# control pattern: G[YN|RC] -> A[YN|RC]

def apobec3gf(seq = "")
  seq_length = seq.size
  apobec_position = []
  control_position = []
  (0..(seq_length - 3)).each do |n|
    tri_base = seq[n,3]
    if tri_base =~ /G[A|G][A|G|T]/
      apobec_position << n
    elsif seq[n] == "G"
      control_position << n
    end
  end
  return [apobec_position,control_position]
end

#compare two sequences, return the number of different positions, NO NEED alignment

def compare_two_seq(seq1 = "", seq2 = "")
  length = seq1.size
  diff = 0
  (0..(length-1)).each do |position|
    nt1 = seq1[position]
    nt2 = seq2[position]
    diff += 1 unless nt1 == nt2
  end
  return diff
end


#compare two sequences, return the number of different positions, need alignment

def compare_two_seq2(sequence1 = "", sequence2 = "")
  aln_seq = muscle_sequence2(sequence1,sequence2)
  seq1 = aln_seq[0]
  seq2 = aln_seq[1]
  length = seq1.size
  diff = 0
  (0..(length-1)).each do |position|
    nt1 = seq1[position]
    nt2 = seq2[position]
    diff += 1 unless nt1 == nt2
  end
  return diff
end

#input hash A, return hash B with the unique values of hash A as keys, and the keys of the unique values of hash A as values of hash B
def uniq_hash(in_hash = {})
  uniq_values = in_hash.values.uniq
  out_hash = {}
  uniq_values.each do |uniq_va|
    in_hash.each do |k,v|
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

#call consensus nucleotide
def call_consensus_base(base_array)
  if base_array.size == 1
    return base_array[0]
  elsif base_array.size == 2
    case base_array.sort!
    when ["A","T"]
      return "W"
    when ["C","G"]
      return "S"
    when ["A","C"]
      return "M"
    when ["G","T"]
      return "K"
    when ["A","G"]
      return "R"
    when ["C","T"]
      return "Y"
    else
      return "N"
    end
  elsif base_array.size == 3
    case base_array.sort!
    when ["C","G","T"]
      return "B"
    when ["A","G","T"]
      return "D"
    when ["A","C","T"]
      return "H"
    when ["A","C","G"]
      return "V"
    else
      return "N"
    end
  else
    return "N"
  end
end


#create one consensus sequence from a sequence array with an optional majority cut-off for mixed bases.
#example:
#position with 15% "A" and 85% "G" will be called as "G" with 20% cut-off and as "R" with 10% cut-off.
def consensus(seq_array, cutoff = 0.5)
  seq_length = seq_array[0].size
  seq_size = seq_array.size
  consensus_seq = ""
  (0..(seq_length - 1)).each do |position|
    all_base = []
    seq_array.each do |seq|
      all_base << seq[position]
    end
    base_count = count(all_base)
    max_base_list = []

    base_count.each do |k,v|
      if v/seq_size.to_f >= cutoff
        max_base_list << k
      end
    end
    consensus_seq += call_consensus_base(max_base_list)
  end
  return consensus_seq
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

=begin
APOBEC3g/f G to A hypermutation
APOBEC3G/F pattern: GRD -> ARD
control pattern: G[YN|RC] -> A[YN|RC]

Two ways to identify hypermutation
1. Fisher's exact test on the frequencies of G to A mutation at A3G positons vs. non-A3G positions
2. Poisson distribution of G to A mutations at A3G positions, outliers sequences

Input sequences hash
Output [hypermutation sequence hash, lines for .csv file ("Sequence,Muts,Out of,Controls,Out of,Rate Ratio,Fishers Exact P-value")
=end

def a3g_hypermut_seq_hash(seq_hash)
  #mut_hash number of apobec3g/f mutations per sequence
  mut_hash = {}
  hm_hash = {}
  out_hash = {}

  #total G->A mutations at apobec3g/f positions.
  total = 0

  #make specimen consensus
  ref = consensus_without_alignment(seq_hash.values)

  #obtain apobec3g positions and control positions
  apobec = apobec3gf(ref)
  mut = apobec[0]
  control = apobec[1]

  seq_hash.each do |k,v|
    a = 0 #muts
    b = 0 #potential mut sites
    c = 0 #control muts
    d = 0 #potenrial controls
    mut.each do |n|
      next if v[n] == "-"
      if v[n] == "A"
        a += 1
        b += 1
      else
        b += 1
      end
    end
    mut_hash[k] = a
    total += a

    control.each do |n|
      next if v[n] == "-"
      if v[n] == "A"
        c += 1
        d += 1
      else
        d += 1
      end
    end
    rr = (a/b.to_f)/(c/d.to_f)

    t1 = b - a
    t2 = d - c

    fet = Rubystats::FishersExactTest.new
    fisher = fet.calculate(t1,t2,a,c)
    perc = fisher[:twotail]
    info = k + "," + a.to_s + "," + b.to_s + "," + c.to_s + "," + d.to_s + "," + rr.round(2).to_s + "," + perc.to_s
    out_hash[k] = info
    if perc < 0.05
      hm_hash[k] = info
    end
  end

  rate = total.to_f/(seq_hash.size)

  count_mut = count(mut_hash.values)
  maxi_count = count_mut.values.max

  poisson_hash = poisson_distribution(rate,maxi_count)

  cut_off = 0
  poisson_hash.each do |k,v|
    cal = seq_hash.size * v
    obs = count_mut[k]
    if obs >= 20 * cal
      cut_off = k
      break
    elsif k == maxi_count
      cut_off = maxi_count
    end
  end

  mut_hash.each do |k,v|
    if v > cut_off
      hm_hash[k] = out_hash[k]
    end
  end

  hm_seq_hash = {}
  hm_hash.keys.each do |k|
    hm_seq_hash[k] = seq_hash[k]
  end
  return [hm_seq_hash,hm_hash]
end

#input a sequence hash, return a sequence hash with stop codons.
def stop_codon_seq_hash(seq_hash, rf = 0)
  out_seq_hash = {}
  seq_hash.each do |k,v|
    sequence = Sequence.new(k,v)
    sequence.get_aa_array(rf)
    if sequence.aa_array.include?("*")
      out_seq_hash[k] = v
    end
  end
  return out_seq_hash
end

#define Poission cut-off for minority variants. (Ref: Zhou, et al. J Virol 2015)
#Input is a sequence array, the redisual sequencing error rate (default = 0.0001), and a fold cut-off to determine the cut-off
#example: cut-off = 2 means that mutations appear at least 2 times are very likely to be a true mutation instead of residual methods errors.
def poisson_minority_cutoff(seq_array, error_rate = 0.0001, fold_cutoff = 20)
  if seq_array.size == 0
    return 0
  else
    cut_off = 1
    l = seq_array[0].size
    rate = seq_array.size * error_rate
    count_mut = variant_for_poisson(seq_array)
    max_count = count_mut.keys.max
    poisson_hash = poisson_distribution(rate, max_count)

    poisson_hash.each do |k,v|
      cal = l * v
      obs = count_mut[k] ? count_mut[k] : 0
      if obs >= fold_cutoff * cal
        cut_off = k
        break
      end
    end
    return cut_off
  end
end

#calculate binomial 95% confidence interval by R. refer to R function binom.test
#input number x and n, return low and high confidence intervals.
def r_binom_CI(x= 0, n= 0,temp_r_dir = File.dirname($0))
  temp_r_result = temp_r_dir + "/temp_r_result"
  print `Rscript -e 'binom.test(#{x},#{n})$conf.int[1];binom.test(#{x},#{n})$conf.int[2]' >#{temp_r_result}`
  lines = File.readlines(temp_r_result)
  low = lines[0].chomp[4..-1].to_f
  high = lines[1].chomp[4..-1].to_f
  File.unlink(temp_r_result)
  return [low.round(5), high.round(5)]
end

#input sequence hash, and Poisson cutoff for minority variants.
#HIV-1 PR region SDRM based on HIVDB.stanford.edu
#only for MPID-DR MiSeq sequences, PR codon 1-99
#return [substitution rate with 95% CI, halpotype abundance with 95% CI, amino acid sequence report spreadsheet]
def sdrm_pr_bulk(sequences, cutoff = 0, temp_r_dir = File.dirname($0))
  region = "PR"
  rf_label = 0
  start_codon_number = 1
  n_seq = sequences.size
  mut = {}
  mut_com = []
  aa = {}
  point_mutation_list = []
  sequences.each do |name,seq|
    s = Sequence.new(name,seq)
    s.get_aa_array(rf_label)
    aa_seq = s.aa_array
    aa[name] = aa_seq.join("")
    record = hiv_protease(aa_seq)
    mut_com << record
    record.each do |position,mutation|
      if mut[position]
        mut[position][1] << mutation[1]
      else
        mut[position] = [mutation[0],[]]
        mut[position][1] << mutation[1]
      end
    end
  end
  mut.each do |position,mutation|
    wt = mutation[0]
    mut_list = mutation[1]
    count_mut_list = count(mut_list)
    count_mut_list.each do |m,number|
      ci = r_binom_CI(number, n_seq, temp_r_dir)
      label = number < cutoff ? "*" : ""
      point_mutation_list << [region, n_seq, position, wt, m, number, (number/n_seq.to_f).round(5), ci[0], ci[1], label]
    end
  end
  point_mutation_list.sort_by! {|record| record[2]}

  link = count(mut_com)
  link2 = {}
  link.each do |k,v|
    pattern = []
    if k.size == 0
      pattern = ['WT']
    else
      k.each do |p,m|
        pattern << (m[0] + p.to_s + m[1])
      end
    end
    link2[pattern.join("+")] = v
  end
  linkage_list = []
  link2.sort_by{|_key,value|value}.reverse.to_h.each do |k,v|
    ci = r_binom_CI(v, n_seq, temp_r_dir)
    label = v < cutoff ? "*" : ""
    linkage_list << [region, n_seq, k, v, (v/n_seq.to_f).round(5), ci[0], ci[1], label]
  end

  report_list = []

  div_aa = {}
  aa_start = start_codon_number

  aa_size = aa.values[0].size - 1

  (0..aa_size).to_a.each do |p|
    aas = []
    aa.values.each do |r1|
      aas << r1[p]
    end
    count_aas = count(aas)
    div_aa[aa_start] = count_aas.sort_by{|k,v|v}.reverse.to_h
    aa_start += 1
  end

  div_aa.each do |k,v|
    record = [region, k, n_seq]
    $amino_acid_list.each do |amino_acid|
      aa_count = v[amino_acid]
      record << (aa_count.to_f/n_seq*100).round(4)
    end
    report_list << record
  end

  return [point_mutation_list, linkage_list, report_list]
end


#input sequence hash, and Poisson cutoff for minority variants.
#HIV-1 RT region SDRM based on HIVDB.stanford.edu
#only for MPID-DR MiSeq sequences
#RT codon 34-122, 152-236 two regions are linked.
#return [substitution rate with 95% CI, halpotype abundance with 95% CI, amino acid sequence report spreadsheet]
def sdrm_rt_bulk(sequences, cutoff = 0, temp_r_dir = File.dirname($0))
  region = "RT"
  rf_label = 1
  start_codon_number = 34
  gap = "AGACTTCAGGAAGTATACTGCATTTACCATACCTAGTATAAACAATGAGACACCAGGGATTAGATATCAGTACAATGTGCTTCCAC"

  n_seq = sequences.size
  mut_nrti = {}
  mut_nnrti = {}
  mut_com = []
  r1_aa = {}
  r2_aa = {}
  point_mutation_list = []
  sequences.each do |name,seq|
    r1 = seq[0,267]
    r2 = seq[267..-1]
    seq = r1 + gap + r2
    s = Sequence.new(name,seq)
    s.get_aa_array(rf_label)
    aa_seq = s.aa_array

    r1_aa[name] = aa_seq[0,89].join("")
    r2_aa[name] = aa_seq[85..-1].join("")
    nrti = sdrm_nrti(aa_seq,start_codon_number)
    nnrti = sdrm_nnrti(aa_seq,start_codon_number)
    mut_com << (nrti.merge(nnrti))

    nrti.each do |position,mutation|
      if mut_nrti[position]
        mut_nrti[position][1] << mutation[1]
      else
        mut_nrti[position] = [mutation[0],[]]
        mut_nrti[position][1] << mutation[1]
      end
    end
    nnrti.each do |position,mutation|
      if mut_nnrti[position]
        mut_nnrti[position][1] << mutation[1]
      else
        mut_nnrti[position] = [mutation[0],[]]
        mut_nnrti[position][1] << mutation[1]
      end
    end
  end

  mut_nrti.each do |position,mutation|
    wt = mutation[0]
    mut_list = mutation[1]
    count_mut_list = count(mut_list)
    count_mut_list.each do |m,number|
      ci = r_binom_CI(number, n_seq, temp_r_dir)
      label = number < cutoff ? "*" : ""
      point_mutation_list << ["NRTI", n_seq, position, wt, m, number, (number/n_seq.to_f).round(5), ci[0], ci[1], label]
    end
  end

  mut_nnrti.each do |position,mutation|
    wt = mutation[0]
    mut_list = mutation[1]
    count_mut_list = count(mut_list)
    count_mut_list.each do |m,number|
      ci = r_binom_CI(number, n_seq, temp_r_dir)
      label = number < cutoff ? "*" : ""
      point_mutation_list << ["NNRTI", n_seq, position, wt, m, number, (number/n_seq.to_f).round(5), ci[0], ci[1], label]
    end
  end
  point_mutation_list.sort_by! {|record| record[2]}

  link = count(mut_com)
  link2 = {}
  link.each do |k,v|
    pattern = []
    if k.size == 0
      pattern = ['WT']
    else
      k.each do |p,m|
        pattern << (m[0] + p.to_s + m[1])
      end
    end
    link2[pattern.join("+")] = v
  end
  linkage_list = []
  link2.sort_by{|_key,value|value}.reverse.to_h.each do |k,v|
    ci = r_binom_CI(v, n_seq, temp_r_dir)
    label = v < cutoff ? "*" : ""
    linkage_list << [region, n_seq, k, v, (v/n_seq.to_f).round(5), ci[0], ci[1], label]
  end

  report_list = []

  div_aa = {}
  r1_aa_start = 34
  r2_aa_start = 152

  r1_aa_size = r1_aa.values[0].size - 1
  r2_aa_size = r2_aa.values[0].size - 1

  (0..r1_aa_size).to_a.each do |p|
    aas = []
    r1_aa.values.each do |r1|
      aas << r1[p]
    end
    count_aas = count(aas)
    div_aa[r1_aa_start] = count_aas.sort_by{|k,v|v}.reverse.to_h
    r1_aa_start += 1
  end

  (0..r2_aa_size).to_a.each do |p|
    aas = []
    r2_aa.values.each do |r1|
      aas << r1[p]
    end
    count_aas = count(aas)
    div_aa[r2_aa_start] = count_aas.sort_by{|k,v|v}.reverse.to_h
    r2_aa_start += 1
  end

  div_aa.each do |k,v|
    record = [region, k, n_seq]
    $amino_acid_list.each do |amino_acid|
      aa_count = v[amino_acid]
      record << (aa_count.to_f/n_seq*100).round(4)
    end
    report_list << record
  end

  return [point_mutation_list, linkage_list, report_list]
end

#input sequence hash, and Poisson cutoff for minority variants.
#HIV-1 IN region SDRM based on HIVDB.stanford.edu
#only for MPID-DR MiSeq sequences
#IN codon 53-174
#return [substitution rate with 95% CI, halpotype abundance with 95% CI, amino acid sequence report spreadsheet]
def sdrm_in_bulk(sequences, cutoff = 0, temp_r_dir = File.dirname($0))
  region = "IN"
  rf_label = 2
  start_codon_number = 53
  n_seq = sequences.size
  mut = {}
  mut_com = []
  aa = {}
  point_mutation_list = []
  sequences.each do |name,seq|
    s = Sequence.new(name,seq)
    s.get_aa_array(rf_label)
    aa_seq = s.aa_array
    aa[name] = aa_seq.join("")
    record = sdrm_int(aa_seq, start_codon_number)
    mut_com << record
    record.each do |position,mutation|
      if mut[position]
        mut[position][1] << mutation[1]
      else
        mut[position] = [mutation[0],[]]
        mut[position][1] << mutation[1]
      end
    end
  end
  mut.each do |position,mutation|
    wt = mutation[0]
    mut_list = mutation[1]
    count_mut_list = count(mut_list)
    count_mut_list.each do |m,number|
      ci = r_binom_CI(number, n_seq, temp_r_dir)
      label = number < cutoff ? "*" : ""
      point_mutation_list << [region, n_seq, position, wt, m, number, (number/n_seq.to_f).round(5), ci[0], ci[1], label]
    end
  end
  point_mutation_list.sort_by! {|record| record[2]}

  link = count(mut_com)
  link2 = {}
  link.each do |k,v|
    pattern = []
    if k.size == 0
      pattern = ['WT']
    else
      k.each do |p,m|
        pattern << (m[0] + p.to_s + m[1])
      end
    end
    link2[pattern.join("+")] = v
  end
  linkage_list = []
  link2.sort_by{|_key,value|value}.reverse.to_h.each do |k,v|
    ci = r_binom_CI(v, n_seq, temp_r_dir)
    label = v < cutoff ? "*" : ""
    linkage_list << [region, n_seq, k, v, (v/n_seq.to_f).round(5), ci[0], ci[1], label]
  end

  report_list = []

  div_aa = {}
  aa_start = start_codon_number

  aa_size = aa.values[0].size - 1

  (0..aa_size).to_a.each do |p|
    aas = []
    aa.values.each do |r1|
      aas << r1[p]
    end
    count_aas = count(aas)
    div_aa[aa_start] = count_aas.sort_by{|k,v|v}.reverse.to_h
    aa_start += 1
  end

  div_aa.each do |k,v|
    record = [region, k, n_seq]
    $amino_acid_list.each do |amino_acid|
      aa_count = v[amino_acid]
      record << (aa_count.to_f/n_seq*100).round(4)
    end
    report_list << record
  end

  return [point_mutation_list, linkage_list, report_list]
end


#calculate Shannon's entropy, Euler's number as the base of logarithm
#https://en.wikipedia.org/wiki/Entropy_(information_theory)

def shannons_entropy(sequences)
  entropy_hash = {}
  seq_l = sequences[0].size
  seq_size = sequences.size
  (0..(seq_l - 1)).each do |position|
    element = []
    sequences.each do |seq|
      element << seq[position]
    end
    entropy = 0
    count(element).each do |k,v|
      p = v/seq_size.to_f
      entropy += (-p * Math.log(p))
    end
    entropy_hash[(position + 1)] = entropy
  end
  return entropy_hash
end

#TN93 distance function. Input: sequence array, output hash: diff => counts
def TN93(sequence_array = [])
  diff = []
  seq_hash = count(sequence_array)
  seq_hash.values.each do |v|
    comb = v * (v - 1) / 2
    comb.times {diff << 0}
  end

  seq_hash.keys.combination(2).to_a.each do |pair|
    s1 = pair[0]
    s2 = pair[1]
    diff_temp = compare_two_seq(s1,s2)
    comb = seq_hash[s1] * seq_hash[s2]
    comb.times {diff << diff_temp}
  end

  count_diff = count(diff)
  out_hash = Hash.new(0)
  Hash[count_diff.sort_by{|k,v|k}].each do |k,v|
    out_hash[k] = v
  end
  return out_hash
end


# gap strip from a sequence alignment

def gap_strip(sequence_alignment)
  new_seq_hash = {}
  seq_size = sequence_alignment.values[0].size
  seq_matrix = {}
  (0..(seq_size - 1)).each do |p|
    seq_matrix[p] = []
    sequence_alignment.values.each do |s|
      seq_matrix[p] << s[p]
    end
  end

  seq_matrix.delete_if do |p, list|
    list.include?("-")
  end

  sequence_alignment.each do |n,s|
    new_s = ""
    seq_matrix.keys.each {|p| new_s += s[p]}
    new_seq_hash[n] = new_s
  end
  return new_seq_hash
end

# gap strip from a sequence alignment, only strip the gaps at the ends of the alignment

def gap_strip_ends(sequence_alignment)
  new_seq_hash = {}
  seq_size = sequence_alignment.values[0].size
  seq_matrix = {}
  (0..(seq_size - 1)).each do |p|
    seq_matrix[p] = []
    sequence_alignment.values.each do |s|
      seq_matrix[p] << s[p]
    end
  end
  n1 = 0
  n2 = 0
  seq_matrix.each do |p, list|
    if list.include?("-")
      n1 += 1
    else
      break
    end
  end

  seq_matrix.keys.reverse.each do |p|
    list = seq_matrix[p]
    if list.include?("-")
      n2 += 1
    else
      break
    end
  end

  sequence_alignment.each do |n,s|
    new_s = s[n1..(- n2 - 1)]
    new_seq_hash[n] = new_s
  end
  return new_seq_hash
end

# return taxa names from a NWK tree, input is a nwk tree string

def taxa_name_nwk(nwk_line)
  taxa = []
  nwk_line.tr("()","").split(",").each {|n| taxa << n.split(":")[0].tr("\'","")}
  return taxa
end


# input paired-end sequence hash format seq_name => [r1_seq, r2_seq]
# overlap is pre-determined
def join_pid_seq_overlap(seq_pair_hash, diff = 0.0, overlap)
  joined_seq_hash = {}
  seq_pair_hash.each do |seq_name, seq_pair|
    r1_seq = seq_pair[0]
    r2_seq = seq_pair[1]
    if overlap.zero?
      joined_seq_hash[seq_name] = r1_seq + r2_seq
    elsif compare_two_seq(r1_seq[-overlap..-1], r2_seq[0,overlap]) <= (overlap * diff)
      joined_seq_hash[seq_name] = r1_seq + r2_seq[overlap..-1]
    else
      next
    end
  end
  return joined_seq_hash
end


# overlap is not predetermined
# model 1: overlap is determined based on consensus, all sequence pairs are supposed to have the same overlap size
# model 2: overlap is determined for each sequence pair, sequence pairs can have different size of overlap
def join_pid_seq_no_overlap(seq_pair_hash, diff = 0.0, model = 1)
  begin
    if model == 1
      overlap = determine_overlap_pid_pair(seq_pair_hash, diff)
      return join_pid_seq_overlap(seq_pair_hash, diff, overlap)
    elsif model == 2
      joined_seq_hash = {}
      seq_pair_hash.each do |seq_name, seq_pair|
        overlap_list = []
        overlap_matrix(seq_pair[0], seq_pair[1]).each do |overlap1, diff_nt|
          cut_off_base = overlap1 * diff
          overlap_list << overlap1 if diff_nt <= cut_off_base
        end
        if overlap_list.empty?
          joined_seq_hash[seq_name] = seq_pair[0] + seq_pair[1]
        else
          overlap = overlap_list.max
          joined_seq_hash[seq_name] = seq_pair[0] + seq_pair[1][overlap..-1]
        end
      end
      return joined_seq_hash
    else
      raise ArgumentError.new("Error::Wrong Overlap Model Argument. Given \'#{model}\', expected '1' or '2'.")
    end
  rescue ArgumentError => e
    puts e
  end
end


def determine_overlap_pid_pair(seq_pair_hash, diff = 0.0)
  overlaps = []
  seq_pair_hash.each do |_seq_name, seq_pair|
    overlap_list = []
    matrix = overlap_matrix(seq_pair[0], seq_pair[1])
    matrix.each do |overlap, diff_nt|
      cut_off_base = overlap * diff
      overlap_list << overlap if diff_nt <= cut_off_base
    end
    if overlap_list.empty?
      overlaps << 0
    else
      overlaps << overlap_list.max
    end
  end
  count_overlaps = count(overlaps)
  max_value = count_overlaps.values.max
  max_overlap_list = []
  count_overlaps.each {|overlap, counts| max_overlap_list << overlap if counts == max_value}
  max_overlap_list.max
end


def overlap_matrix(sequence1, sequence2)
  min_overlap = 6
  max_overlap = [sequence1.size, sequence2.size].max
  matrix_hash = {}
  (min_overlap..max_overlap).each do |overlap|
    matrix_hash[overlap] = compare_two_seq(sequence1[-overlap..-1], sequence2[0, overlap])
  end
  return matrix_hash
end

# input a directory with r1 and r2 sequences, return a hash :seq_name => [r1_seq, r2_seq]
# r1 and r2 file names should contain "r1" and "r2" respectively
# the sequence taxa should only differ by last 3 characters to distinguish r1 and r2 sequence.
def pair_fasta_to_hash(indir)
  files = Dir[indir + "/*"]
  r1_file = ""
  r2_file = ""
  files.each do |f|
    if File.basename(f) =~ /r1/
      r1_file = f
    elsif File.basename(f) =~ /r2/
      r2_file = f
    end
  end

  seq1 = fasta_to_hash(r1_file)
  seq2 = fasta_to_hash(r2_file)

  new_seq1 = seq1.each_with_object({}) {|(k, v), h| h[k[0..-4]] = v}
  new_seq2 = seq2.each_with_object({}) {|(k, v), h| h[k[0..-4]] = v}

  seq_pair_hash = {}

  new_seq1.each do |seq_name,seq|
    seq_pair_hash[seq_name] = [seq, new_seq2[seq_name]]
  end
  return seq_pair_hash
end

# batch quality check of HIV sequences based on :sequence_locator
# input a sequence hash, start nt position(s) in an array, end nt position(s) in a array
# and allow the sequence to contain indels
# return a hash of filtered sequences

def qc_hiv_seq_check(seq_hash, start_nt, end_nt, indel=true)
  seq_hash_unique = seq_hash.values.uniq
  seq_hash_unique_pass = []
  seq_hash_unique.each do |seq|
    loc = sequence_locator(seq)
    if start_nt.include?(loc[0]) && end_nt.include?(loc[1])
      if indel
        seq_hash_unique_pass << seq
      elsif loc[3] == nil
        seq_hash_unique_pass << seq
      end
    end
  end
  seq_pass = {}
  seq_hash_unique_pass.each do |seq|
    seq_hash.each do |seq_name, orginal_seq|
      if orginal_seq == seq
        seq_pass[seq_name] =  seq
        seq_hash.delete(seq_name)
      end
    end
  end
  return seq_pass
end
