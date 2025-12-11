
module ViralSeq
  class SeqHash

    # functions to identify SDRMs from a ViralSeq::SeqHash object at HIV PR region.
    #   works for MPID-DR protocol (dx.doi.org/10.17504/protocols.io.useewbe)
    #   PR codon 1-99
    #   RT codon 34-122 (HXB2 2649-2914) and 152-236(3001-3257)
    #   IN codon 53-174 (HXB2 4384-4751)
    # @param cutoff [Integer] cut-off for minimal abundance of a mutation to be called as valid mutation,
    #   can be obtained using ViralSeq::SeqHash#poisson_minority_cutoff function
    # @param fdr [Hash] hash of events => (false detecton rate)
    #   can be obtained using ViralSeq::SeqHash#fdr
    #
    # @return [Array] three elements `[point_mutation_list, linkage_list, report_list]`
    #
    #   # point_mutation_list: two demensional array for the following information,
    #     # [region,tcs_number,position,wildtype,mutation,count,%,CI_low,CI_high,fdr,label]
    #   # linkage_list: two demensional array for the following information,
    #     # [region,tcs_number,linkage,count,%,CI_low,CI_high,label]
    #   # report_list: two demensional array for the following information,
    #     # [position,codon,tcs_number,A,C,D,E,F,G,H,I,K,L,M,N,P,Q,R,S,T,V,W,Y,*]
    # @example identify SDRMs from a FASTA sequence file of HIV PR sequences obtained after MPID-DR sequencing
    #   my_seqhash = ViralSeq::SeqHash.fa('spec/sample_files/sample_dr_sequences/pr.fasta')
    #   p_cut_off = my_seqhash.pm
    #   fdr_hash = my_seqhash.fdr
    #   pr_sdrm = my_seqhash.sdrm_hiv_pr(p_cut_off, fdr_hash)
    #   puts "region,tcs_number,position,wildtype,mutation,count,%,CI_low,CI_high,fdr,label"; pr_sdrm[0].each {|n| puts n.join(',')}
    #   => region,tcs_number,position,wildtype,mutation,count,%,CI_low,CI_high,fdr,label
    #   => PR,396,30,D,N,247,0.62374,0.57398,0.67163,0,
    #   => PR,396,50,I,V,1,0.00253,6.0e-05,0.01399,0.18905,*
    #   => PR,396,88,N,D,246,0.62121,0.57141,0.66919,0,
    #
    #   puts "region,tcs_number,linkage,count,%,CI_low,CI_high,label"; pr_sdrm[1].each {|n| puts n.join(',')}
    #   => region,tcs_number,linkage,count,%,CI_low,CI_high,label
    #   => PR,396,D30N+N88D,245,0.61869,0.56884,0.66674,
    #   => PR,396,WT,149,0.37626,0.32837,0.42602,
    #   => PR,396,D30N,1,0.00253,6.0e-05,0.01399,*
    #   => PR,396,D30N+I50V+N88D,1,0.00253,6.0e-05,0.01399,*
    #
    #   puts "position,codon,tcs_number," + ViralSeq::AMINO_ACID_LIST.join(","); pr_sdrm[2].each {|n|puts n.join(",")}
    #   => position,codon,tcs_number,A,C,D,E,F,G,H,I,K,L,M,N,P,Q,R,S,T,V,W,Y,*
    #   => PR,1,396,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,100.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0
    #   => PR,2,396,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,100.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0
    #   => PR,3,396,0.0,0.0,0.0,0.0,0.0,0.0,0.0,100.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0
    #   => PR,4,396,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,100.0,0.0,0.0,0.0,0.0
    #   => PR,5,396,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,100.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0
    #   => PR,6,396,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,100.0,0.0,0.0
    #   => PR,7,396,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,100.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0
    #   => PR,8,396,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,100.0,0.0,0.0,0.0,0.0,0.0,0.0
    #   => PR,9,396,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,100.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0
    #   => PR,10,396,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,100.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0
    #   => PR,11,396,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,100.0,0.0,0.0,0.0
    #   => PR,12,396,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,37.8788,62.1212,0.0,0.0,0.0,0.0
    #   => PR,13,396,0.0,0.0,0.0,0.0,0.0,0.0,0.0,38.1313,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,61.8687,0.0,0.0,0.0
    #   => PR,14,396,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,100.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0
    #   => PR,15,396,0.0,0.0,0.0,0.0,0.0,0.0,0.0,62.3737,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,37.6263,0.0,0.0,0.0
    #   => PR,16,396,0.0,0.0,0.0,0.0,0.0,100.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0
    #   => PR,17,396,0.0,0.0,0.0,0.0,0.0,100.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0
    #   => PR,18,396,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,99.4949,0.5051,0.0,0.0,0.0,0.0,0.0,0.0
    #   => PR,19,396,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,100.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0
    #   => PR,20,396,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,100.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0
    #   => PR,21,396,0.0,0.0,0.0,100.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0
    #   => PR,22,396,100.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0
    #   => PR,23,396,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,100.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0
    #   => PR,24,396,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,100.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0
    #   => PR,25,396,0.0,0.0,100.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0
    #   => PR,26,396,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,100.0,0.0,0.0,0.0,0.0
    #   => PR,27,396,0.0,0.0,0.0,0.0,0.0,100.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0
    #   => PR,28,396,100.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0
    #   => PR,29,396,0.0,0.0,100.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0
    #   => PR,30,396,0.0,0.0,37.6263,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,62.3737,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0
    #   => PR,31,396,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,100.0,0.0,0.0,0.0,0.0
    #   => PR,32,396,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,100.0,0.0,0.0,0.0
    #   => PR,33,396,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,99.7475,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.2525,0.0,0.0,0.0
    #   => PR,34,396,0.0,0.0,0.0,100.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0
    #   => PR,35,396,0.0,0.0,62.1212,37.6263,0.0,0.0,0.0,0.0,0.2525,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0
    #   => PR,36,396,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,99.7475,0.0,0.0,0.0,0.0,0.0,0.0,0.2525,0.0,0.0,0.0
    #   => PR,37,396,0.0,0.0,37.8788,61.8687,0.0,0.0,0.0,0.0,0.2525,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0
    #   => PR,38,396,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,100.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0
    #   => PR,39,396,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,99.4949,0.0,0.0,0.5051,0.0,0.0,0.0,0.0,0.0
    #   => PR,40,396,0.0,0.0,0.0,0.0,0.0,100.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0
    #   => PR,41,396,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,37.8788,0.0,0.0,0.0,0.0,0.0,62.1212,0.0,0.0,0.0,0.0,0.0,0.0
    #   => PR,42,396,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,100.0,0.0,0.0
    #   => PR,43,396,0.0,0.0,0.0,0.2525,0.0,0.0,0.0,0.0,99.7475,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0
    #   => PR,44,396,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,100.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0
    #   => PR,45,396,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,100.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0
    #   => PR,46,396,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,100.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0
    #   => PR,47,396,0.0,0.0,0.0,0.0,0.0,0.0,0.0,100.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0
    #   => PR,48,396,0.0,0.0,0.0,0.0,0.0,100.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0
    #   => PR,49,396,0.0,0.0,0.0,0.0,0.0,100.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0
    #   => PR,50,396,0.0,0.0,0.0,0.0,0.0,0.0,0.0,99.7475,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.2525,0.0,0.0,0.0
    #   => PR,51,396,0.0,0.0,0.0,0.0,0.0,100.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0
    #   => PR,52,396,0.0,0.0,0.0,0.0,0.0,100.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0
    #   => PR,53,396,0.0,0.0,0.0,0.0,100.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0
    #   => PR,54,396,0.0,0.0,0.0,0.0,0.0,0.0,0.0,100.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0
    #   => PR,55,396,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,100.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0
    #   => PR,56,396,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,100.0,0.0,0.0,0.0
    #   => PR,57,396,0.0,0.0,0.0,0.0,0.0,0.2525,0.0,0.0,0.2525,0.0,0.0,0.0,0.0,0.0,99.4949,0.0,0.0,0.0,0.0,0.0,0.0
    #   => PR,58,396,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,100.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0
    #   => PR,59,396,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,100.0,0.0
    #   => PR,60,396,0.0,0.0,100.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0
    #   => PR,61,396,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,100.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0
    #   => PR,62,396,0.0,0.0,0.0,0.0,0.0,0.0,0.0,100.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0
    #   => PR,63,396,0.0,0.0,0.0,0.0,0.0,0.0,0.2525,0.0,0.0,37.8788,0.0,0.0,61.8687,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0
    #   => PR,64,396,0.0,0.0,0.0,0.0,0.0,0.0,0.0,62.1212,0.0,37.8788,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0
    #   => PR,65,396,0.0,0.0,0.0,100.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0
    #   => PR,66,396,0.0,0.0,0.0,0.0,0.0,0.0,0.0,100.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0
    #   => PR,67,396,0.0,100.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0
    #   => PR,68,396,0.0,0.0,0.0,0.0,0.0,100.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0
    #   => PR,69,396,0.0,0.0,0.0,0.0,0.0,0.0,100.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0
    #   => PR,70,396,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,100.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0
    #   => PR,71,396,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,62.1212,37.8788,0.0,0.0,0.0
    #   => PR,72,396,0.0,0.0,0.0,37.8788,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,62.1212,0.0,0.0,0.0,0.0
    #   => PR,73,396,0.0,0.0,0.0,0.0,0.0,100.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0
    #   => PR,74,396,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,100.0,0.0,0.0,0.0,0.0
    #   => PR,75,396,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,100.0,0.0,0.0,0.0
    #   => PR,76,396,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,100.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0
    #   => PR,77,396,0.0,0.0,0.0,0.0,0.0,0.0,0.0,99.7475,0.0,0.0,0.2525,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0
    #   => PR,78,396,0.0,0.0,0.0,0.0,0.0,100.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0
    #   => PR,79,396,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,100.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0
    #   => PR,80,396,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,100.0,0.0,0.0,0.0,0.0
    #   => PR,81,396,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,100.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0
    #   => PR,82,396,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,100.0,0.0,0.0,0.0
    #   => PR,83,396,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,99.4949,0.0,0.0,0.0,0.5051,0.0,0.0,0.0,0.0,0.0
    #   => PR,84,396,0.0,0.0,0.0,0.0,0.0,0.0,0.0,100.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0
    #   => PR,85,396,0.0,0.0,0.0,0.0,0.0,0.0,0.0,100.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0
    #   => PR,86,396,0.0,0.0,0.0,0.5051,0.0,99.4949,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0
    #   => PR,87,396,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,100.0,0.0,0.0,0.0,0.0,0.0,0.0
    #   => PR,88,396,0.0,0.0,62.1212,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,37.8788,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0
    #   => PR,89,396,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,100.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0
    #   => PR,90,396,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,100.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0
    #   => PR,91,396,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,100.0,0.0,0.0,0.0,0.0
    #   => PR,92,396,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,100.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0
    #   => PR,93,396,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,100.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0
    #   => PR,94,396,0.0,0.0,0.0,0.0,0.0,100.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0
    #   => PR,95,396,0.0,100.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0
    #   => PR,96,396,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,100.0,0.0,0.0,0.0,0.0
    #   => PR,97,396,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,100.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0
    #   => PR,98,396,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,99.7475,0.0,0.0,0.0,0.2525,0.0,0.0,0.0,0.0,0.0
    #   => PR,99,396,0.0,0.0,0.0,0.0,100.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0

    def sdrm_hiv_pr(cutoff = 0, fdr_hash = Hash.new(0))
      sequences = self.dna_hash
      region = "PR"
      rf_label = 0
      start_codon_number = 1
      n_seq = sequences.size
      mut = {}
      mut_com = []
      aa = {}
      point_mutation_list = []
      sequences.each do |name,seq|
        s = ViralSeq::Sequence.new(name,seq)
        s.translate(rf_label)
        aa[name] = s.aa_string
        record = s.sdrm(:PI)
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
        count_mut_list = mut_list.count_freq
        count_mut_list.each do |m,number|
          ci = ViralSeq::Math::BinomCI.new(number, n_seq)
          fdr = fdr_hash[number].round(5)
          label = number < cutoff ? "*" : ""
          point_mutation_list << [region, n_seq, position, wt, m, number, ci.mean.round(5), ci.lower.round(5), ci.upper.round(5), fdr, label]
        end
      end
      point_mutation_list.sort_by! {|record| record[2]}

      link = mut_com.count_freq
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
        ci = ViralSeq::Math::BinomCI.new(v, n_seq)
        label = v < cutoff ? "*" : ""
        linkage_list << [region, n_seq, k, v, ci.mean.round(5), ci.lower.round(5), ci.upper.round(5), label]
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
        count_aas = aas.count_freq
        div_aa[aa_start] = count_aas.sort_by{|_k,v|v}.reverse.to_h
        aa_start += 1
      end

      div_aa.each do |k,v|
        record = [region, k, n_seq]
        ViralSeq::AMINO_ACID_LIST.each do |amino_acid|
          aa_count = v[amino_acid]
          record << (aa_count.to_f/n_seq*100).round(4)
        end
        report_list << record
      end

      return [point_mutation_list, linkage_list, report_list]
    end


    # functions to identify SDRMs from a ViralSeq::SeqHash object at HIV RT region.
    #   works for MPID-DR protocol (dx.doi.org/10.17504/protocols.io.useewbe)
    #   RT codon 34-122, 152-236, two regions are linked
    # @param (see #sdrm_hiv_pr)
    # @return (see #sdrm_hiv_pr)

    def sdrm_hiv_rt(cutoff = 0, fdr_hash = Hash.new(0))
      sequences = self.dna_hash
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
        s = ViralSeq::Sequence.new(name,seq)
        s.translate(rf_label)

        r1_aa[name] = s.aa_string[0,89]
        r2_aa[name] = s.aa_string[-85..-1]
        nrti = s.sdrm(:nrti, start_codon_number)
        nnrti = s.sdrm(:nnrti, start_codon_number)
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
        count_mut_list = mut_list.count_freq
        count_mut_list.each do |m,number|
          ci = ViralSeq::Math::BinomCI.new(number, n_seq)
          fdr = fdr_hash[number].round(5)
          label = number < cutoff ? "*" : ""
          point_mutation_list << ["NRTI", n_seq, position, wt, m, number, ci.mean.round(5), ci.lower.round(5), ci.upper.round(5), fdr, label]
        end
      end

      mut_nnrti.each do |position,mutation|
        wt = mutation[0]
        mut_list = mutation[1]
        count_mut_list = mut_list.count_freq
        count_mut_list.each do |m,number|
          ci = ViralSeq::Math::BinomCI.new(number, n_seq)
          fdr = fdr_hash[number].round(5)
          label = number < cutoff ? "*" : ""
          point_mutation_list << ["NNRTI", n_seq, position, wt, m, number, ci.mean.round(5), ci.lower.round(5), ci.upper.round(5), fdr, label]
        end
      end

      point_mutation_list.sort_by! {|record| record[2]}

      link = mut_com.count_freq
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
        ci = ViralSeq::Math::BinomCI.new(v, n_seq)
        label = v < cutoff ? "*" : ""
        linkage_list << [region, n_seq, k, v, ci.mean.round(5), ci.lower.round(5), ci.upper.round(5), label]
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
        count_aas = aas.count_freq
        div_aa[r1_aa_start] = count_aas.sort_by{|_k,v|v}.reverse.to_h
        r1_aa_start += 1
      end

      (0..r2_aa_size).to_a.each do |p|
        aas = []
        r2_aa.values.each do |r1|
          aas << r1[p]
        end
        count_aas = aas.count_freq
        div_aa[r2_aa_start] = count_aas.sort_by{|_k,v|v}.reverse.to_h
        r2_aa_start += 1
      end

      div_aa.each do |k,v|
        record = [region, k, n_seq]
        ViralSeq::AMINO_ACID_LIST.each do |amino_acid|
          aa_count = v[amino_acid]
          record << (aa_count.to_f/n_seq*100).round(4)
        end
        report_list << record
      end

      return [point_mutation_list, linkage_list, report_list]
    end

    # functions to identify SDRMs from a ViralSeq::SeqHash object at HIV IN region.
    #   works for MPID-DR protocol (dx.doi.org/10.17504/protocols.io.useewbe)
    #   IN codon 53-174
    # @param (see #sdrm_hiv_pr)
    # @return (see #sdrm_hiv_pr)

    def sdrm_hiv_in(cutoff = 0, fdr_hash = Hash.new(0))
      sequences = self.dna_hash
      region = "IN"
      rf_label = 2
      start_codon_number = 53
      n_seq = sequences.size
      mut = {}
      mut_com = []
      aa = {}
      point_mutation_list = []
      sequences.each do |name,seq|
        s = ViralSeq::Sequence.new(name,seq)
        s.translate(rf_label)
        aa[name] = s.aa_string
        record = s.sdrm(:INSTI, start_codon_number)
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
        count_mut_list = mut_list.count_freq
        count_mut_list.each do |m,number|
          ci = ViralSeq::Math::BinomCI.new(number, n_seq)
          fdr = fdr_hash[number].round(5)
          label = number < cutoff ? "*" : ""
          point_mutation_list << [region, n_seq, position, wt, m, number, ci.mean.round(5), ci.lower.round(5), ci.upper.round(5), fdr, label]
        end
      end
      point_mutation_list.sort_by! {|record| record[2]}

      link = mut_com.count_freq
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
        ci = ViralSeq::Math::BinomCI.new(v, n_seq)
        label = v < cutoff ? "*" : ""
        linkage_list << [region, n_seq, k, v, ci.mean.round(5), ci.lower.round(5), ci.upper.round(5), label]
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
        count_aas = aas.count_freq
        div_aa[aa_start] = count_aas.sort_by{|_k,v|v}.reverse.to_h
        aa_start += 1
      end

      div_aa.each do |k,v|
        record = [region, k, n_seq]
        ViralSeq::AMINO_ACID_LIST.each do |amino_acid|
          aa_count = v[amino_acid]
          record << (aa_count.to_f/n_seq*100).round(4)
        end
        report_list << record
      end

      return [point_mutation_list, linkage_list, report_list]
    end


    # wrapper function for #a3g_hypermut and #stop_codon with ViralSeq::DrmRegionConfig as a param.

    def filter_for_drm(region_config)
      seq_coord = region_config.seq_coord
      reading_frame_number = region_config.get_reading_frame_number

      if !seq_coord["gap"]

        a3g_check = self.a3g
        a3g_seqs = a3g_check[:a3g_seq]
        a3g_filtered_seqs = a3g_check[:filtered_seq]

        stop_codon_check = a3g_filtered_seqs.stop_codon(reading_frame_number[0])
        stop_codon_seqs = stop_codon_check[:with_stop_codon]
        filtered_seqs = stop_codon_check[:without_stop_codon]

        return {
          filtered_seq: filtered_seqs,
          a3g_seq: a3g_seqs,
          stop_codon_seq: stop_codon_seqs
        }

      else

        r1_length, r2_length = region_config.r1_r2_length.values

        r1_seqs = {}
        r2_seqs = {}

        self.dna_hash.each do |k,v|
            r1_seqs[k] = v[0,r1_length]
            r2_seqs[k] = v[-r2_length..-1] # to ensure the length from the end. Sometimes the platform will return sequence with one extra base.
        end

        r1_sh = ViralSeq::SeqHash.new(r1_seqs)
        r2_sh = ViralSeq::SeqHash.new(r2_seqs)

        a3g_seqs_r1 = r1_sh.a3g[:a3g_seq]
        a3g_seqs_r2 = r2_sh.a3g[:a3g_seq]

        stop_codon_r1 = r1_sh.stop_codon(reading_frame_number[0])[:with_stop_codon]
        stop_codon_r2 = r2_sh.stop_codon(reading_frame_number[1])[:with_stop_codon]

        a3g_seq_keys = (a3g_seqs_r1.dna_hash.keys | a3g_seqs_r2.dna_hash.keys)
        a3g_seqs = ViralSeq::SeqHash.new(self.dna_hash.select {|k, _v| a3g_seq_keys.include? k})

        stop_codon_keys = (stop_codon_r1.dna_hash.keys | stop_codon_r2.dna_hash.keys)
        stop_codon_seqs = ViralSeq::SeqHash.new(self.dna_hash.select {|k, _v| stop_codon_keys.include? k})

        reject_keys = (a3g_seq_keys | stop_codon_keys)

        filtered_seqs = ViralSeq::SeqHash.new(self.dna_hash.reject { |k, _v| reject_keys.include? k })

        return {
          filtered_seq: filtered_seqs,
          a3g_seq: a3g_seqs,
          stop_codon_seq: stop_codon_seqs
        }

      end

    end # end of #filter_for_drm


    # insert the partial genome into the whole gene for HIV resistance analysis


    def complete_with_ref(region_config)
      complete_seqs = {}
      seq_coord = region_config.seq_coord

      ref = ViralSeq::RefSeq.get(region_config.ref_info["ref_type"].to_sym)
      a = region_config.ref_info["ref_coord"][0]
      b = region_config.ref_info["ref_coord"][1]
      c = seq_coord["minimum"]
      d = seq_coord["maximum"]

      if seq_coord["gap"]
        e = seq_coord["gap"]["minimum"]
        f = seq_coord["gap"]["maximum"]

        self.dna_hash.each do |k,v|
          complete_seqs[k] = ref[(a-1)..(c-2)] + v[0,(e-c)] + ref[(e-1)..(f-1)] + v[(e-c)..-1] + ref[d..(b-1)]
        end
      else
        self.dna_hash.each do |k,v|
          complete_seqs[k] = ref[(a-1)..(c-2)] + v + ref[d..(b-1)]
        end
      end

      return ViralSeq::SeqHash.new(complete_seqs)
    end #end of #complete_with_ref


    # function to interpret HIV drms with ViralSeq::DrmRegionConfig as a param.

    def drm(region_config)
      region = region_config.region
      fdr_hash = self.fdr # must run fdr before the completion of the sequences

      complete_gene = self.complete_with_ref(region_config)
      sequences = complete_gene.dna_hash

      n_seq = sequences.size
      aa = {}
      mut = {}
      mut_com = []
      point_mutation_list = []

      drm_list = region_config.drm_list

      sequences.each do |name, seq|
        s = ViralSeq::Sequence.new(name, seq)
        s.translate
        aa[name] = s.aa_string

        records_per_seq = {}

        drm_list.each do |drm_class, list|

          mut[drm_class] = {}  if !mut[drm_class]

          record = s.check_drm(list)
          records_per_seq = records_per_seq.merge(record)

          record.each do |position, mutation|
            if !mut[drm_class][position]
              mut[drm_class][position] = [mutation[0],[]]
            end
            mut[drm_class][position][1] << mutation[1]
          end
        end

        mut_com << records_per_seq.sort.to_h
      end

      mut.each do |drm_class, mutations|
        mutations.each do |position, mutation|
          wt = mutation[0]
          mut_list = mutation[1]
          count_mut_list = mut_list.count_freq
          count_mut_list.each do |m,number|
            ci = ViralSeq::Math::BinomCI.new(number, n_seq)
            fdr = fdr_hash[number].round(5)
            label = fdr >= 0.05 ? "*" : ""
            point_mutation_list << [drm_class, n_seq, position, wt, m, number, ci.mean.round(5), ci.lower.round(5), ci.upper.round(5), fdr, label]
          end
        end
      end

      point_mutation_list.sort_by! {|record| record[2]}

      link = mut_com.count_freq
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
        ci = ViralSeq::Math::BinomCI.new(v, n_seq)
        label = ""
        linkage_list << [region, n_seq, k, v, ci.mean.round(5), ci.lower.round(5), ci.upper.round(5), label]
      end

      report_list = []

      div_aa = {}
      aa_start = 1

      aa_size = aa.values[0].size - 1

      (0..aa_size).to_a.each do |p|
        aas = []
        aa.values.each do |r1|
          aas << r1[p]
        end
        count_aas = aas.count_freq
        div_aa[aa_start] = count_aas.sort_by{|_k,v|v}.reverse.to_h
        aa_start += 1
      end

      div_aa.each do |k,v|
        record = [region, k, n_seq]
        ViralSeq::AMINO_ACID_LIST.each do |amino_acid|
          aa_count = v[amino_acid]
          record << (aa_count.to_f/n_seq*100).round(4)
        end
        report_list << record
      end

      return [point_mutation_list, linkage_list, report_list]

    end

  end # end of ViralSeq::SeqHash

end # end of ViralSeq
