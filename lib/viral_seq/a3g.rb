# viral_seq/a3g
# APOBEC3g/f hypermutation function including
# ViralSeq::a3g_hypermut_seq_hash
# ViralSeq::apobec3gf

# APOBEC3g/f G to A hypermutation
# APOBEC3G/F pattern: GRD -> ARD
# control pattern: G[YN|RC] -> A[YN|RC]
# use the sample consensus to determine potential a3g sites

# Two criteria to identify hypermutation
# 1. Fisher's exact test on the frequencies of G to A mutation at A3G positons vs. non-A3G positions
# 2. Poisson distribution of G to A mutations at A3G positions, outliers sequences
# note:  criteria 2 only applies on a sequence file containing more than 20 sequences
#        b/c Poisson model does not do well on small sample size.

# ViralSeq.a3g_hypermut_seq_hash(sequence_hash)
# sequence_hash is a Hash object for sequences. {:name => :sequence, ...}
# return array [hypermutation_hash, statistic_info]
# hypermutation_hash is a Hash object for sequences
# statistic_info is a hash object of [sequence_name, stats],
# in which stats String object in csv format (separated by ',') containing
#   sequence tag
#   G to A mutation numbers at potential a3g positions
#   total potential a3g G positions
#   G to A mutation numbers at non a3g positions
#   total non a3g G positions
#   a3g G to A mutation rate / non-a3g G to A mutation rate
#   Fishers Exact P-value
#
# =USAGE
#   # example 1
#   sequences = ViralSeq.fasta_to_hash('spec/sample_files/sample_a3g_sequence1.fasta')
#   hypermut = ViralSeq.a3g_hypermut_seq_hash(sequences)
#   hypermut[0].keys
#   => [">Seq7", ">Seq14"]
#   stats = hypermut[1]
#   stats.values
#   => [">Seq7,23,68,1,54,18.26,4.308329383112348e-06", ">Seq14,45,68,9,54,3.97,5.2143571971582974e-08"]
#
#   # example 2
#   sequences = ViralSeq.fasta_to_hash('spec/sample_files/sample_a3g_sequence2.fasta')
#   hypermut = ViralSeq.a3g_hypermut_seq_hash(sequences)
#   stats = hypermut[1]
#   stats = values
#   => [">CTAACACTCA_134_a3g-sample2,4,35,0,51,Infinity,0.02465676660128911", ">ATAGTGCCCA_60_a3g-sample2,4,35,1,51,5.83,0.1534487353839561"]
#   # notice sequence ">ATAGTGCCCA_60_a3g-sample2" has a p value at 0.15, greater than 0.05, but it is still called as hypermutation sequence b/c it's Poisson outlier sequence.


# ViralSeq.apobec3gf(sequence)
# APOBEC3G/F pattern: GRD -> ARD
# control pattern: G[YN|RC] -> A[YN|RC]
# input a sequence String object
# return all two arrays of position numbers of
#   a3g G positions (a3g)
#   non-a3g G positions (control)


module ViralSeq
  def ViralSeq.a3g_hypermut_seq_hash(seq_hash)
    # mut_hash number of apobec3g/f mutations per sequence
    mut_hash = {}
    hm_hash = {}
    out_hash = {}

    # total G->A mutations at apobec3g/f positions.
    total = 0

    # make consensus sequence for the input sequence hash
    ref = ViralSeq.consensus(seq_hash.values)

    # obtain apobec3g positions and control positions
    apobec = ViralSeq.apobec3gf(ref)
    mut = apobec[0]
    control = apobec[1]

    seq_hash.each do |k,v|
      a = 0 # muts
      b = 0 # potential mut sites
      c = 0 # control muts
      d = 0 # potenrial controls
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

    if seq_hash.size > 20
      rate = total.to_f/(seq_hash.size)

      count_mut = ViralSeq.count(mut_hash.values)
      maxi_count = count_mut.values.max

      poisson_hash = ViralSeq.poisson_distribution(rate,maxi_count)

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
    end

    hm_seq_hash = {}
    hm_hash.keys.each do |k|
      hm_seq_hash[k] = seq_hash[k]
    end
    return [hm_seq_hash,hm_hash]
  end

  # APOBEC3G/F mutation position identification
  # APOBEC3G/F pattern: GRD -> ARD
  # control pattern: G[YN|RC] -> A[YN|RC]

  def self.apobec3gf(seq = "")
    seq.tr!("-", "")
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

end
