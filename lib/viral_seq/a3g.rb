# APOBEC3g/f G to A hypermutation
# APOBEC3G/F pattern: GRD -> ARD
# control pattern: G[YN|RC] -> A[YN|RC]
#
# Two ways to identify hypermutation
# 1. Fisher's exact test on the frequencies of G to A mutation at A3G positons vs. non-A3G positions
# 2. Poisson distribution of G to A mutations at A3G positions, outliers sequences
#
# Input sequences hash
# Output [hypermutation sequence hash, lines for .csv file ("Sequence,Muts,Out of,Controls,Out of,Rate Ratio,Fishers Exact P-value")

module ViralSeq
  def ViralSeq.a3g_hypermut_seq_hash(seq_hash)
    # mut_hash number of apobec3g/f mutations per sequence
    mut_hash = {}
    hm_hash = {}
    out_hash = {}

    # total G->A mutations at apobec3g/f positions.
    total = 0

    #m ake consensus sequence for the input sequence hash
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

    rate = total.to_f/(seq_hash.size)

    count_mut = Viralseq.count(mut_hash.values)
    maxi_count = count_mut.values.max

    poisson_hash = Viralseq.poisson_distribution(rate,maxi_count)

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

  # APOBEC3G/F mutation position identification
  # APOBEC3G/F pattern: GRD -> ARD
  # control pattern: G[YN|RC] -> A[YN|RC]

  def self.apobec3gf(seq = "")
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
