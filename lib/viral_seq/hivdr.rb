
module ViralSeq
  class SDRM

    # functions to identify SDRMs from a ViralSeq::SeqHash object at HIV PR region.
    #   works for MPID-DR protocol (dx.doi.org/10.17504/protocols.io.useewbe)
    #   PR codon 1-99
    #   RT codon 34-122 (HXB2 2650-2914) and 152-236(3001-3257)
    #   IN codon 53-174 (HXB2 4384-4751)
    # @param cutoff [Integer] cut-off for minimal abundance of a mutation to be called as valid mutation,
    #   can be obtained using ViralSeq::SeqHash#poisson_minority_cutoff function
    # @return [Array] three elements `[point_mutation_list, linkage_list, report_list]`
    #
    #   # point_mutation_list: two demensional array for the following information,
    #     # [region,tcs_number,position,wildtype,mutation,count,%,CI_low,CI_high,label]
    #   # linkage_list: two demensional array for the following information,
    #     # [region,tcs_number,linkage,count,%,CI_low,CI_high,label]
    #   # report_list: two demensional array for the following information,
    #     # [position,codon,tcs_number,A,C,D,E,F,G,H,I,K,L,M,N,P,Q,R,S,T,V,W,Y,*]
    # @example identify SDRMs from a FASTA sequence file of HIV PR sequences obtained after MPID-DR sequencing
    #   my_seqhash = ViralSeq::SeqHash.fa('spec/sample_files/sample_dr_sequences/pr.fasta')
    #   p_cut_off = my_seqhash.pm
    #   pr_sdrm = my_seqhash.sdrm_hiv_pr(p_cut_off)
    #   puts "region,tcs_number,position,wildtype,mutation,count,%,CI_low,CI_high,label"; pr_sdrm[0].each {|n| puts n.join(',')}
    #   => region,tcs_number,position,wildtype,mutation,count,%,CI_low,CI_high,label
    #   => PR,396,30,D,N,247,0.62374,0.57398,0.67163,
    #   => PR,396,50,I,V,1,0.00253,6.0e-05,0.01399,*
    #   => PR,396,88,N,D,246,0.62121,0.57141,0.66919,
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

    def sdrm_hiv_pr(cutoff = 0)
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
        record = s.sdrm(:hiv_pr)
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
          label = number < cutoff ? "*" : ""
          point_mutation_list << [region, n_seq, position, wt, m, number, ci.mean.round(5), ci.lower.round(5), ci.upper.round(5), label]
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

    def sdrm_hiv_rt(cutoff = 0)
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
          label = number < cutoff ? "*" : ""
          point_mutation_list << ["NRTI", n_seq, position, wt, m, number, ci.mean.round(5), ci.lower.round(5), ci.upper.round(5), label]
        end
      end

      mut_nnrti.each do |position,mutation|
        wt = mutation[0]
        mut_list = mutation[1]
        count_mut_list = mut_list.count_freq
        count_mut_list.each do |m,number|
          ci = ViralSeq::Math::BinomCI.new(number, n_seq)
          label = number < cutoff ? "*" : ""
          point_mutation_list << ["NNRTI", n_seq, position, wt, m, number, ci.mean.round(5), ci.lower.round(5), ci.upper.round(5), label]
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

    def sdrm_hiv_in(cutoff = 0)
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
        record = s.sdrm(:hiv_in, start_codon_number)
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
          label = number < cutoff ? "*" : ""
          point_mutation_list << [region, n_seq, position, wt, m, number, ci.mean.round(5), ci.lower.round(5), ci.upper.round(5), label]
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

  end # end of ViralSeq::SeqHash
end # end of ViralSeq
