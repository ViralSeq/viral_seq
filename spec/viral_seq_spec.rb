RSpec.describe ViralSeq do
  it "has a version number" do
    expect(ViralSeq::VERSION).not_to be nil
  end

  it "has reverse complement function" do
    expect("ACTG".rc).to eq "CAGT"
  end

  it "sequence_locator function PR example works" do
    seq = "CCTCAGATCACTCTTTGGCAACGACCCCTAGTTACAATAAGGGTAGGGGGGCAACTAAAGGAAGCCCTATTAGATACAGGAGCAGATGATACAGTATTAGAAGAAATAAATTTGCCAGGAAGATGGAAACCAAAAATGATAGGGGGAATTGGAGGTTTTATCAAAGTAAGACAATATGATCAGATACCCATAGAAATTTGTGGACATGAAGCTATAGGTACAGTATTAGTGGGACCTACACCTGTCAACATAATTGGGAGAAATCTGTTGACTCAGATTGGTTGCACTCTAAATTTT"
    s = ViralSeq::Sequence.new('my_sequence', seq)
    loc = s.locator
    expect(loc[0]).to eq 2253
    expect(loc[1]).to eq 2549
  end

  it "locator function V1V3 example works" do
    seq = "AAATTAACCCCACTCTGTGTTACCTTAAATTGCACTAACGCGGCCAACAGGACCAATAATGTGACCACTGAGACCAATGTGACCACTGAGACCAGAATTTACCCAGACATGATAGGTGAAATAAAAAATTGCTCTTTCAATACCTCCACAAACCTAGTAGGTAAGGATCAGAAAAATTATGCACTGTTTCGCAGCCTTGATATAGTACCAATAGAAGATAATAAGAGTAGTAATAGTAGTAATTTTACCAGCTATATGCTGACAAGCACAGTGCAATGTACACATGGAATTAGGCCAGTAGTGTCCACTCAACTGCTGTTAAATGGTAGTCTAGCAGAAGAAGACATAGTAATTAGGTCTGAGAACATCACAAATAATGTTAAAAACATAATAGTGCACCTGAATGAATCTGTAGAGATTAATTGTACGAGACCAGGCAACAATACAAGAAAAAGTATAACTATAGGACCAGGGAGAGCATTTTATGCCACAGGAGATATAATAGGAGATATAAGAAAA"
    s = ViralSeq::Sequence.new('my_sequence', seq)
    loc = s.locator
    expect(loc[1]).to eq 7208
    expect(loc[0]).to eq 6585
    expect(loc[2]).to eq 58.8
  end

  it "has sequence_clip function" do
    seq = "CCTCAGATCACTCTTTGGCAACGACCCCTAGTTACAATAAGGGTAGGGGGGCAACTAAAGGAAGCCCTATTAGATACAGGAGCAGATGATACAGTATTAGAAGAAATAAATTTGCCAGGAAGATGGAAACCAAAAATGATAGGGGGAATTGGAGGTTTTATCAAAGTAAGACAATATGATCAGATACCCATAGAAATTTGTGGACATGAAGCTATAGGTACAGTATTAGTGGGACCTACACCTGTCAACATAATTGGGAGAAATCTGTTGACTCAGATTGGTTGCACTCTAAATTTT"
    s = ViralSeq::Sequence.new('my_seq', seq)
    clip = s.sequence_clip(2333, 2433, :HXB2).dna
    expect(clip).to eq 'AGCAGATGATACAGTATTAGAAGAAATAAATTTGCCAGGAAGATGGAAACCAAAAATGATAGGGGGAATTGGAGGTTTTATCAAAGTAAGACAATATGATC'
  end

  it "can do mutiple sequence alignment with MUSCLE" do
    seq = %w{>CCCGCGGTTA_159_PD001-1_V1V3_combined AAATTAACCCCACTCTGTGTTACCTTAAATTGCACTAATGCGGCCAACAGAACCAATGTGACCACTGAGACCAGAATTTACCCAGACATGATAGGTGAAATAAAAAACTGCTCTTTCAATACCTCCACAGGCCTAGTAGGTAAGGATCAGAAAAATTATGCACTGTTTCGCAGCTTTGATGTAGTACCAATAGAACATAATAATAGTAGTAATTTTACCAGCTATATGCTGACAAGTTGTAACACCTCAGTCATTAAACAGGCCTGCACAGTGCAATGTACACATGGAATTAGGCCAGTAGTGTCCACTCAACTGCTGTTAAATGGTAGTCTAGCAGAAGAAGACATAGTAATTAGGTCTGAGAACATCACAAATAATGTTAAAAACATAATAGTGCACCTGAATGAATCTGTAGAGATTAAATGTATGAGACCAGGCAACAATACAAGAAAAAGTATAACTATAGGACCAGGGAGAGCATTTTATGCCACAGGAGATATAATAGGAGATATAAGAAAA >TGATCTGTGC_195_PD001-1_V1V3_combined AAATTAACCCCACTCTGTGTTACCTTAAATTGCACTAATGCGGCCAACAGAACCAATGTGACCACTGAGACCAGAATTTACCCAGACATGATAGGTGAAATAAAAAACTGCTCTTTCAATACCTCCACAGGCCTAGTAGGTAAGGATCAGAAAAATTATGCACTGTTTCGCAGCTTTGATGTAGTACCAATAGAACATAATAATAGTAGTAATTTTACCAGCTATATGCTGACAAGTTGTAACACCTCAGTCATTAAACAGGCCTGCACAGTGCAATGTACACATGGAATTAGGCCAGTAGTGTCCACTCAACTGCTGTTAAATGGTAGTCTAGCAGAAGAAGACATAGTAATTAGGTCTGAGAACATCACAAATAATGTTAAAAACATAATAGTGCACCTGAATGAATCTGTAGAGATTAAATGTATGAGACCAGGCAACAATACAAGAAAAAGTATAACTATAGGACCAGGGAGAGCATTTTATGCCACAGGAGATATAATAGGAGATATAAGAAAA >TAGAGGACTT_35_PD001-8_V1V3_combined AAGCTAACTCCACTCTGTGTTACCTTAAATTGCACTGACTATGTGGGGAATAATACTAAGAATGCCACTAAGAGTAAGGAAGAAATAGAAATGAAAAACTGCTCTTTCAATGTCACTGAAGTCATAAGGGATAAGGTGCAGAAAGAATATGCACTGTTTTATAAACTTGATATAGTACCAATAGATGAAGGTGGTCTTAACAAGACTGTTAATAATACCACATATAGGTTGATAAGTTGTAACACCTCAGTTATTAGACAGGCCTGCACAGTACAATGTACACATGGAATTAGGCCAGTAGTGTCAACACAATTGCTATTAAATGGTACTCTAGCAAAAGATAAGGTAGTAATTAGATCTGAAAATTTCACAGACAATGCAAAAACTATAATAGTACAGCTGAACGAATCTGTAGAAATTCACTGTGCAAGACCCAACAACAATACAAGAAAAAGTATACATATAGCACCAGGAAGAGCATTTTATGCAACAGGACAAATAATAGGAGATATAAGAAAA >AGTGAAAGTC_91_PD001-8_V1V3_combined AAGCTAACTCCACTCTGTGTTACTTTAAATTGCACTGACTATGTGGGGAATAATACTAAGAATGCCACTAAGAGTAAGGAAGAAATAGAAATGAAAAACTGCTCTTTCAATGTCACTGAAGTCATAAGAGATAAGGTGCAGAAAGAATATGCACTGTTTTATAAACTTGATATAGTGCCAATAGAGGAAGAGGGTCTTAACAAGACTGTTAATAATACCACATATAGGTTGATAAGTTGTAACACCTCAGTCATTAGACAGGCATGCACAGTACAATGTACACATGGAATTAGGCCAGTAGTGTCAACACAACTGCTATTAAATGGTACTCTAGCAAAAGATAAGGTAGTAATTAGATCTGAAAATTTCACAGACAATGCAAAAACTATAATAGTACAGCTGAACGAGTCTGTAGAAATTCACTGTGCAAGACCCAACAACAATACAAGAAAAAGTATACATATAGCACCAGGGAGAGCATTTTATGCAACAGGACAAATAATAGGAGATATAAGAAAA}
    seq_hash = ViralSeq::SeqHash.new(Hash[*seq])
    aln_seq = seq_hash.align
    expect(aln_seq.dna_hash.size).to eq 4
  end

  it "has a function to make consensus sequences" do
    seq_array = %w{

      ATTTTTTTTT
      AATTTTTTTT
      AAATTTTTTT
      AAAATTTTTT
      AAAAATTTTT
      AAAAAATTTT
      AAAAAAATTT
      AAAAAAAATT
      AAAAAAAAAT
      AAAAAAAAAA

     }
     my_seqhash = ViralSeq::SeqHash.array(seq_array)
     expect(my_seqhash.consensus).to eq 'AAAAAWTTTT'
     expect(my_seqhash.consensus(0.6)).to eq 'AAAAANTTTT'
     expect(my_seqhash.consensus(0.7)).to eq 'AAAANNNTTT'
     expect(my_seqhash.consensus(0.8)).to eq 'AAANNNNNTT'
     expect(my_seqhash.consensus(0.9)).to eq 'AANNNNNNNT'
     expect(my_seqhash.consensus(0.4)).to eq 'AAAAWWWTTT'
  end

  it "has a function to calculate Primer ID consensus cut-off" do
    expect(ViralSeq::Math.calculate_pid_cut_off(1000, 0.015)).to eq 17
    expect(ViralSeq::Math.calculate_pid_cut_off(10)).to eq 2
  end

  it "can filter HIV sequences by locations" do
    sequence_hash = ViralSeq::SeqHash.fa('spec/sample_files/sample_seq.fasta')
    filtered_sequence_hash = sequence_hash.hiv_seq_qc(4384, 4751, false, :HXB2, 'muscle')
    expect(filtered_sequence_hash.dna_hash.size).to eq 4
  end

  it "can read paired fasta sequences as a paired sequence hash" do
    expect(ViralSeq::SeqHashPair.fa('spec/sample_files/sample_paired_seq').dna_hash.size).to eq 29
  end

  it "can generate Primer ID pool by given Primer ID length" do
    expect(ViralSeq::PID.generate_pool(10).size).to eq 1048576
  end

  it "can filter sequences with residual offspring Primer IDs" do
    seqhash = ViralSeq::SeqHash.fa('spec/sample_files/sample_pid_filter.fasta')
    expect(seqhash.filter_similar_pid.dna_hash.size).to eq 3
  end

  it "can parse the nucleotide sequences as a String object and return a Regexp object for possible matches" do
    expect("ATRWCG".nt_parser.to_s).to eq "(?-mix:AT[A|G][A|T]CG)"
  end

  it "can calculate Poisson probability of k number of events given λ" do
    expect(ViralSeq::Math::PoissonDist.new(0.005).poisson_hash[2]).to eq 1.243765598990853e-05
  end

  it "can compare two sequences as String object and return number of differences" do
    seq1 = 'AAGGCGTAGGAC'
    seq2 = 'AAGCTTAGGACG'
    aligned_seqs = ViralSeq::Muscle.align(seq1,seq2)
    expect(seq1.compare_with(seq2)).to eq 8
    expect(aligned_seqs[0].compare_with(aligned_seqs[1])).to eq 4
  end

  it "has a gap strip function for a sequence alignment" do
      sequence_hash = {'>seq1' => 'AACCGGTT',
                       '>seq2' => 'A-CCGGTT',
                       '>seq3' => 'AAC-GGTT',
                       '>seq4' => 'AACCG-TT',
                       '>seq5' => 'AACCGGT-'}
      my_seqhash = ViralSeq::SeqHash.new(sequence_hash)
      expected_hash = {">seq1"=>"ACGT", ">seq2"=>"ACGT", ">seq3"=>"ACGT", ">seq4"=>"ACGT", ">seq5"=>"ACGT"}
      expect(my_seqhash.gap_strip.dna_hash).to eq expected_hash
  end

  it "has a gap strip function for a sequence alignment only at both ends" do
      sequence_hash = {'>seq1' => 'AACCGGTT',
                       '>seq2' => 'A-CCGGTT',
                       '>seq3' => 'AAC-GGTT',
                       '>seq4' => 'AACCG-TT',
                       '>seq5' => 'AACCGGT-'}
      my_seqhash = ViralSeq::SeqHash.new(sequence_hash)
      expected_hash = {">seq1"=>"AACCGGT", ">seq2"=>"A-CCGGT", ">seq3"=>"AAC-GGT", ">seq4"=>"AACCG-T", ">seq5"=>"AACCGGT"}
      expect(my_seqhash.gap_strip_ends.dna_hash).to eq expected_hash
  end

  paired_seqs = {">pair1"=>["GGGGGGGGGGAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA",
                            "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAATTTTTTTTTT"],
                ">pair2"=>["GGGGGGGGGGAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA",
                           "AAAAAAAAAAGAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAATTTTTTTTTT"],
                ">pair3"=>["GGGGGGGGGGAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA",
                           "AAAAAAAAAAGGAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAATTTTTTTTTT"]
                }
  p1 = ViralSeq::SeqHashPair.new(paired_seqs)

  paired_seq2 = {">pair4" => ["AAAGGGGGGG", "GGGGGGGTT"],
                 ">pair5" => ["AAAAAAGGGG", "GGGGTTTTT"],
                 ">pair6" => ["AAACAAGGGG", "GGGGTTTTT"]
                 }
  p2 = ViralSeq::SeqHashPair.new(paired_seq2)

  it "has a function to join paired-end reads with KNOWN overlap" do
    expect(p1.join1(100).dna_hash.keys).to eq [">pair1"]
    expect(p1.join1(100, 0.01).dna_hash.keys).to eq [">pair1", ">pair2"]
    expect(p1.join1(100, 0.02).dna_hash.keys).to eq [">pair1", ">pair2", ">pair3"]
  end

  it "has a function to join paired-end reads with UNKNOWN overlap" do
    expected_hash1 = {">pair4" => "AAAGGGGGGGGGGTT", ">pair5" => "AAAAAAGGGGTTTTT", ">pair6" => "AAACAAGGGGTTTTT"}
    expected_hash2 = {">pair4"=>"AAAGGGGGGGTT", ">pair5"=>"AAAAAAGGGGTTTTT", ">pair6"=>"AAACAAGGGGTTTTT"}
    expect(p2.join2.dna_hash).to eq expected_hash1
    expect(p2.join2(model: :indiv).dna_hash).to eq expected_hash2
  end

  it "has a function to find APOBEC3g/f hypermutation sequences" do
    sequences1 = ViralSeq::SeqHash.fa('spec/sample_files/sample_a3g_sequence1.fasta')
    sequences2 = ViralSeq::SeqHash.fa('spec/sample_files/sample_a3g_sequence2.fasta')
    hypermut1 = sequences1.a3g
    hypermut2 = sequences2.a3g
    expect(hypermut1[0].dna_hash.keys).to eq [">Seq7", ">Seq14"]
    expect(hypermut2[0].dna_hash.keys).to eq [">CTAACACTCA_134_a3g-sample2", ">ATAGTGCCCA_60_a3g-sample2"]
  end

  it "has a function to identify HCV NS5A drug resistance mutations given a sequence" do
    seq = 'GGCTAAGGGCCAAGCTCATGCCACAATTGCCCGGGATCCCTTTTGTGTCCTGCCAACGCGGATATAGGGGGGTCTGGAAGGGAGATGGCATTATGCACACTCGCTGCCACTGCGGAGCTGAGATCACTGGACATGTCAAGAACGGGACGATGAGGATCGCCGGTCCTAAGACCTGCAGAAACATGTGGAGTGGGACCTTCCCCATCAACGCCTGCACCACGGGCCCCTGTACCCCCCTTCCCGCGCCGAACTATACGTTCGCGTTGTGGAGGGTGTCTGCGGAGGAATACGTGGAAATAAGGCGGGTGGGAGACTTCCACTACGTAACGGGC'
    s = ViralSeq::Sequence.new('seq', seq)
    s.translate(2)
    expected_resistance = {30=>["L", "Q"], 44=>["K", "R"], 78=>["R", "K"], 83=>["T", "M"], 93=>["Y", "C"], 107=>["K", "T"]}
    expect(s.sdrm(:hcv_ns5a,23)).to eq expected_resistance
  end

  it "has a function to calculate Shannon's entropy from a sequence alignment" do
    sequence_file = 'spec/sample_files/sample_sequence_alignment_for_entropy.fasta'
    sequence_hash = ViralSeq::SeqHash.fa(sequence_file)
    entropy_hash = sequence_hash.shannons_entropy
    expect(entropy_hash[3].zero?).to be true
    expect(entropy_hash[14].round(3)).to be 0.639
    expect(entropy_hash[46].round(3)).to be 0.325
  end

  it "has a function to calculate nucleotide pairwise diversity (π)" do
    sequences = %w{ AAGGCCTT ATGGCCTT AAGGCGTT AAGGCCTT AACGCCTT AAGGCCAT }
    my_seqhash = ViralSeq::SeqHash.array(sequences)
    expect(my_seqhash.pi).to be 0.16667
  end

  it "has a function to tabulate pairwise comparison (TN93)" do
    sequences = %w{ AAGGCCTT ATGGCCTT AAGGCGTT AAGGCCTT AACGCCTT AAGGCCAT }
    my_seqhash = ViralSeq::SeqHash.array(sequences)
    expect(my_seqhash.tn93[2]).to be 6
  end

  it "has a function for Poisson minority cut-off" do
    sequence_file = 'spec/sample_files/sample_sequence_for_poisson.fasta'
    sequences = ViralSeq::SeqHash.fa(sequence_file)
    expect(sequences.pm).to be 2
  end

  it "has a function to make unique sequence Hash from a sequence Hash" do
    sequences = {'>seq1' => 'AAAA','>seq2' => 'AAAA', '>seq3' => 'AAAA',
                 '>seq4' => 'CCCC', '>seq5' => 'CCCC',
                 '>seq6' => 'TTTT' }
    uniq_sequence = {">sequence_1_3"=>"AAAA", ">sequence_2_2"=>"CCCC", ">sequence_3_1"=>"TTTT"}
    seq_hash = ViralSeq::SeqHash.new(sequences)
    expect(seq_hash.uniq_dna_hash.dna_hash).to eq uniq_sequence
  end

  it "has a function to detect SDRMs from a sequence Hash" do
    pr_sequence = ViralSeq::SeqHash.fa('spec/sample_files/sample_dr_sequences/pr.fasta')
    pr_p_cut_off = pr_sequence.pm
    pr_sdrm = pr_sequence.sdrm_hiv_pr(pr_p_cut_off)
    expect(pr_sdrm[0][0][5]).to eq 247

    rt_sequence = ViralSeq::SeqHash.fa('spec/sample_files/sample_dr_sequences/rt.fasta')
    rt_p_cut_off = rt_sequence.pm
    rt_sdrm = rt_sequence.sdrm_hiv_rt(rt_p_cut_off)
    expect(rt_sdrm[0][1][5]).to eq 52

    in_sequence = ViralSeq::SeqHash.fa('spec/sample_files/sample_dr_sequences/in.fasta')
    in_p_cut_off = in_sequence.pm
    in_sdrm = in_sequence.sdrm_hiv_in(in_p_cut_off)
    expect(in_sdrm[1][0][3]).to eq 452
  end

  it "can do sequence locator on a SeqHash object" do
    my_seqhash = ViralSeq::SeqHash.fa('spec/sample_files/sample_seq.fasta')
    loc = my_seqhash.loc
    expect(loc[0][4]).to eq 4384
    expect(loc[4][6]).to eq 80.3
  end
end
