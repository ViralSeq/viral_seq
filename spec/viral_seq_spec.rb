RSpec.describe ViralSeq do
  it "has a version number" do
    expect(ViralSeq::VERSION).not_to be nil
  end

  it "has reverse complement function" do
    expect("ACTG".rc).to eq "CAGT"
  end

  it "sequence_locator function PR example works" do
    seq = "CCTCAGATCACTCTTTGGCAACGACCCCTAGTTACAATAAGGGTAGGGGGGCAACTAAAGGAAGCCCTATTAGATACAGGAGCAGATGATACAGTATTAGAAGAAATAAATTTGCCAGGAAGATGGAAACCAAAAATGATAGGGGGAATTGGAGGTTTTATCAAAGTAAGACAATATGATCAGATACCCATAGAAATTTGTGGACATGAAGCTATAGGTACAGTATTAGTGGGACCTACACCTGTCAACATAATTGGGAGAAATCTGTTGACTCAGATTGGTTGCACTCTAAATTTT"
    loc = ViralSeq.sequence_locator(seq, :HXB2, 'muscle')
    expect(loc[0]).to eq 2253
    expect(loc[1]).to eq 2549
  end

  it "locator function V1V3 example works" do
    seq = "AAATTAACCCCACTCTGTGTTACCTTAAATTGCACTAACGCGGCCAACAGGACCAATAATGTGACCACTGAGACCAATGTGACCACTGAGACCAGAATTTACCCAGACATGATAGGTGAAATAAAAAATTGCTCTTTCAATACCTCCACAAACCTAGTAGGTAAGGATCAGAAAAATTATGCACTGTTTCGCAGCCTTGATATAGTACCAATAGAAGATAATAAGAGTAGTAATAGTAGTAATTTTACCAGCTATATGCTGACAAGCACAGTGCAATGTACACATGGAATTAGGCCAGTAGTGTCCACTCAACTGCTGTTAAATGGTAGTCTAGCAGAAGAAGACATAGTAATTAGGTCTGAGAACATCACAAATAATGTTAAAAACATAATAGTGCACCTGAATGAATCTGTAGAGATTAATTGTACGAGACCAGGCAACAATACAAGAAAAAGTATAACTATAGGACCAGGGAGAGCATTTTATGCCACAGGAGATATAATAGGAGATATAAGAAAA"
    loc = ViralSeq.sequence_locator(seq, :HXB2, 'muscle')
    expect(loc[1]).to eq 7208
    expect(loc[0]).to eq 6585
    expect(loc[2]).to eq 58.8
  end

  it "has sequence_clip function" do
    seq = "CCTCAGATCACTCTTTGGCAACGACCCCTAGTTACAATAAGGGTAGGGGGGCAACTAAAGGAAGCCCTATTAGATACAGGAGCAGATGATACAGTATTAGAAGAAATAAATTTGCCAGGAAGATGGAAACCAAAAATGATAGGGGGAATTGGAGGTTTTATCAAAGTAAGACAATATGATCAGATACCCATAGAAATTTGTGGACATGAAGCTATAGGTACAGTATTAGTGGGACCTACACCTGTCAACATAATTGGGAGAAATCTGTTGACTCAGATTGGTTGCACTCTAAATTTT"
    clip = ViralSeq.sequence_clip(seq, 2333, 2433, :HXB2, 'muscle')
    expect(clip).to eq 'AGCAGATGATACAGTATTAGAAGAAATAAATTTGCCAGGAAGATGGAAACCAAAAATGATAGGGGGAATTGGAGGTTTTATCAAAGTAAGACAATATGATC'
  end

  it "can do mutiple sequence alignment with MUSCLE" do
    seq = %w{>CCCGCGGTTA_159_PD001-1_V1V3_combined AAATTAACCCCACTCTGTGTTACCTTAAATTGCACTAATGCGGCCAACAGAACCAATGTGACCACTGAGACCAGAATTTACCCAGACATGATAGGTGAAATAAAAAACTGCTCTTTCAATACCTCCACAGGCCTAGTAGGTAAGGATCAGAAAAATTATGCACTGTTTCGCAGCTTTGATGTAGTACCAATAGAACATAATAATAGTAGTAATTTTACCAGCTATATGCTGACAAGTTGTAACACCTCAGTCATTAAACAGGCCTGCACAGTGCAATGTACACATGGAATTAGGCCAGTAGTGTCCACTCAACTGCTGTTAAATGGTAGTCTAGCAGAAGAAGACATAGTAATTAGGTCTGAGAACATCACAAATAATGTTAAAAACATAATAGTGCACCTGAATGAATCTGTAGAGATTAAATGTATGAGACCAGGCAACAATACAAGAAAAAGTATAACTATAGGACCAGGGAGAGCATTTTATGCCACAGGAGATATAATAGGAGATATAAGAAAA >TGATCTGTGC_195_PD001-1_V1V3_combined AAATTAACCCCACTCTGTGTTACCTTAAATTGCACTAATGCGGCCAACAGAACCAATGTGACCACTGAGACCAGAATTTACCCAGACATGATAGGTGAAATAAAAAACTGCTCTTTCAATACCTCCACAGGCCTAGTAGGTAAGGATCAGAAAAATTATGCACTGTTTCGCAGCTTTGATGTAGTACCAATAGAACATAATAATAGTAGTAATTTTACCAGCTATATGCTGACAAGTTGTAACACCTCAGTCATTAAACAGGCCTGCACAGTGCAATGTACACATGGAATTAGGCCAGTAGTGTCCACTCAACTGCTGTTAAATGGTAGTCTAGCAGAAGAAGACATAGTAATTAGGTCTGAGAACATCACAAATAATGTTAAAAACATAATAGTGCACCTGAATGAATCTGTAGAGATTAAATGTATGAGACCAGGCAACAATACAAGAAAAAGTATAACTATAGGACCAGGGAGAGCATTTTATGCCACAGGAGATATAATAGGAGATATAAGAAAA >TAGAGGACTT_35_PD001-8_V1V3_combined AAGCTAACTCCACTCTGTGTTACCTTAAATTGCACTGACTATGTGGGGAATAATACTAAGAATGCCACTAAGAGTAAGGAAGAAATAGAAATGAAAAACTGCTCTTTCAATGTCACTGAAGTCATAAGGGATAAGGTGCAGAAAGAATATGCACTGTTTTATAAACTTGATATAGTACCAATAGATGAAGGTGGTCTTAACAAGACTGTTAATAATACCACATATAGGTTGATAAGTTGTAACACCTCAGTTATTAGACAGGCCTGCACAGTACAATGTACACATGGAATTAGGCCAGTAGTGTCAACACAATTGCTATTAAATGGTACTCTAGCAAAAGATAAGGTAGTAATTAGATCTGAAAATTTCACAGACAATGCAAAAACTATAATAGTACAGCTGAACGAATCTGTAGAAATTCACTGTGCAAGACCCAACAACAATACAAGAAAAAGTATACATATAGCACCAGGAAGAGCATTTTATGCAACAGGACAAATAATAGGAGATATAAGAAAA >AGTGAAAGTC_91_PD001-8_V1V3_combined AAGCTAACTCCACTCTGTGTTACTTTAAATTGCACTGACTATGTGGGGAATAATACTAAGAATGCCACTAAGAGTAAGGAAGAAATAGAAATGAAAAACTGCTCTTTCAATGTCACTGAAGTCATAAGAGATAAGGTGCAGAAAGAATATGCACTGTTTTATAAACTTGATATAGTGCCAATAGAGGAAGAGGGTCTTAACAAGACTGTTAATAATACCACATATAGGTTGATAAGTTGTAACACCTCAGTCATTAGACAGGCATGCACAGTACAATGTACACATGGAATTAGGCCAGTAGTGTCAACACAACTGCTATTAAATGGTACTCTAGCAAAAGATAAGGTAGTAATTAGATCTGAAAATTTCACAGACAATGCAAAAACTATAATAGTACAGCTGAACGAGTCTGTAGAAATTCACTGTGCAAGACCCAACAACAATACAAGAAAAAGTATACATATAGCACCAGGGAGAGCATTTTATGCAACAGGACAAATAATAGGAGATATAAGAAAA}
    seq_hash = Hash[*seq]
    aln_seq = ViralSeq.muscle_align_multi(seq_hash)
    expect(aln_seq.size).to eq 4
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
     expect(ViralSeq.consensus(seq_array)).to eq 'AAAAAWTTTT'
     expect(ViralSeq.consensus(seq_array, 0.6)).to eq 'AAAAANTTTT'
     expect(ViralSeq.consensus(seq_array, 0.7)).to eq 'AAAANNNTTT'
     expect(ViralSeq.consensus(seq_array, 0.8)).to eq 'AAANNNNNTT'
     expect(ViralSeq.consensus(seq_array, 0.9)).to eq 'AANNNNNNNT'
     expect(ViralSeq.consensus(seq_array, 0.4)).to eq 'AAAAWWWTTT'
  end

  it "has a function to calculate Primer ID consensus cut-off" do
    expect(ViralSeq.calculate_pid_cut_off(1000, 0.015)).to eq 17
    expect(ViralSeq.calculate_pid_cut_off(10)).to eq 2
  end

  it "can filter HIV sequences by locations" do
    sequence_hash = ViralSeq.fasta_to_hash('spec/sample_files/sample_seq.fasta') # load the .fasta file as a sequence hash
    filtered_sequence_hash = ViralSeq.qc_hiv_seq_check(sequence_hash, 4384, 4751, false, :HXB2, 'muscle')
    expect(filtered_sequence_hash.size).to eq 4
  end

  it "can read paired fasta sequences as a paired sequence hash" do
    expect(ViralSeq.pair_fasta_to_hash('spec/sample_files/sample_paired_seq').size).to eq 29
  end

  it "can generate Primer ID pool by given Primer ID length" do
    expect(ViralSeq.generate_primer_id_pool(10).size).to eq 1048576
  end

  it "can filter sequences with residual offspring Primer IDs" do
    expect(ViralSeq.filter_similar_pid('spec/sample_files/sample_pid_filter.fasta').size).to eq 3
  end

  it "can parse the nucleotide sequences as a String object and return a Regexp object for possible matches" do
    expect("ATRWCG".nt_parser.to_s).to eq "(?-mix:AT[A|G][A|T]CG)"
  end

  it "can calculate Poisson probability of k number of events given λ" do
    expect(ViralSeq.poisson_distribution(0.005)[2]).to eq 1.243765598990853e-05
  end

  it "can compare two sequences as String object and return number of differences" do
    seq1 = 'AAGGCGTAGGAC'
    seq2 = 'AAGCTTAGGACG'
    aligned_seqs = ViralSeq.muscle_align(seq1,seq2)
    expect(ViralSeq.compare_two_seq(seq1, seq2)).to eq 8
    expect(ViralSeq.compare_two_seq(aligned_seqs[0], aligned_seqs[1])).to eq 4
  end

  it "has a gap strip function for a sequence alignment" do
      sequence_hash = {'>seq1' => 'AACCGGTT',
                       '>seq2' => 'A-CCGGTT',
                       '>seq3' => 'AAC-GGTT',
                       '>seq4' => 'AACCG-TT',
                       '>seq5' => 'AACCGGT-'}
      expected_hash = {">seq1"=>"ACGT", ">seq2"=>"ACGT", ">seq3"=>"ACGT", ">seq4"=>"ACGT", ">seq5"=>"ACGT"}
      expect(ViralSeq.gap_strip(sequence_hash)).to eq expected_hash
  end

  it "has a gap strip function for a sequence alignment only at both ends" do
      sequence_hash = {'>seq1' => 'AACCGGTT',
                       '>seq2' => 'A-CCGGTT',
                       '>seq3' => 'AAC-GGTT',
                       '>seq4' => 'AACCG-TT',
                       '>seq5' => 'AACCGGT-'}
      expected_hash = {">seq1"=>"AACCGGT", ">seq2"=>"A-CCGGT", ">seq3"=>"AAC-GGT", ">seq4"=>"AACCG-T", ">seq5"=>"AACCGGT"}
      expect(ViralSeq.gap_strip_ends(sequence_hash)).to eq expected_hash
  end


  paired_seqs = {">pair1"=>["GGGGGGGGGGAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA",
                            "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAATTTTTTTTTT"],
                ">pair2"=>["GGGGGGGGGGAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA",
                           "AAAAAAAAAAGAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAATTTTTTTTTT"],
                ">pair3"=>["GGGGGGGGGGAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA",
                           "AAAAAAAAAAGGAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAATTTTTTTTTT"]
                }

  paired_seq2 = {">pair4" => ["AAAGGGGGGG", "GGGGGGGTT"],
                 ">pair5" => ["AAAAAAGGGG", "GGGGTTTTT"],
                 ">pair6" => ["AAACAAGGGG", "GGGGTTTTT"]
                 }

  it "has a function to join paired-end reads with KNOWN overlap" do
    expect(ViralSeq.paired_join1(paired_seqs, 100, 0.00).keys).to eq [">pair1"]
    expect(ViralSeq.paired_join1(paired_seqs, 100, 0.01).keys).to eq [">pair1", ">pair2"]
    expect(ViralSeq.paired_join1(paired_seqs, 100, 0.02).keys).to eq [">pair1", ">pair2", ">pair3"]
  end

  it "has a function to join paired-end reads with UNKNOWN overlap" do
    expected_hash1 = {">pair4" => "AAAGGGGGGGGGGTT", ">pair5" => "AAAAAAGGGGTTTTT", ">pair6" => "AAACAAGGGGTTTTT"}
    expected_hash2 = {">pair4"=>"AAAGGGGGGGTT", ">pair5"=>"AAAAAAGGGGTTTTT", ">pair6"=>"AAACAAGGGGTTTTT"}
    expect(ViralSeq.paired_join2(paired_seq2, 1)).to eq expected_hash1
    expect(ViralSeq.paired_join2(paired_seq2, 2)).to eq expected_hash2
  end

  it "has a function to find APOBEC3g/f hypermutation sequences" do
    sequences1 = ViralSeq.fasta_to_hash('spec/sample_files/sample_a3g_sequence1.fasta')
    sequences2 = ViralSeq.fasta_to_hash('spec/sample_files/sample_a3g_sequence2.fasta')
    hypermut1 = ViralSeq.a3g_hypermut_seq_hash(sequences1)
    hypermut2 = ViralSeq.a3g_hypermut_seq_hash(sequences2)
    expect(hypermut1[0].keys).to eq [">Seq7", ">Seq14"]
    expect(hypermut2[0].keys).to eq [">CTAACACTCA_134_a3g-sample2", ">ATAGTGCCCA_60_a3g-sample2"]
  end

  it "has a function to identify HCV NS5A drug resistance mutations given a sequence" do
    seq = 'GGCTAAGGGCCAAGCTCATGCCACAATTGCCCGGGATCCCTTTTGTGTCCTGCCAACGCGGATATAGGGGGGTCTGGAAGGGAGATGGCATTATGCACACTCGCTGCCACTGCGGAGCTGAGATCACTGGACATGTCAAGAACGGGACGATGAGGATCGCCGGTCCTAAGACCTGCAGAAACATGTGGAGTGGGACCTTCCCCATCAACGCCTGCACCACGGGCCCCTGTACCCCCCTTCCCGCGCCGAACTATACGTTCGCGTTGTGGAGGGTGTCTGCGGAGGAATACGTGGAAATAAGGCGGGTGGGAGACTTCCACTACGTAACGGGC'
    s = ViralSeq::Sequence.new('seq', seq)
    s.get_aa_array(2)
    aa_seq = s.aa_array
    expected_resistance = {30=>["L", "Q"], 44=>["K", "R"], 78=>["R", "K"], 83=>["T", "M"], 93=>["Y", "C"], 107=>["K", "T"]}
    expect(ViralSeq.hcv_ns5a(aa_seq, 23)).to eq expected_resistance
  end
end
