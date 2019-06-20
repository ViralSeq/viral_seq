# viral_seq/sdrm_core.rb
# core functions for HIV SDRM analysis using MPID-DR protocol.
# More details for HIV Surveillance Drug Resistance Mutation (SDRM) can be found at
# https://hivdb.stanford.edu/pages/surveillance.html

# Including methods as:
#   ViralSeq::sdrm_nrti
#   ViralSeq::sdrm_nnrti
#   ViralSeq::hiv_protease
#   ViralSeq::sdrm_int
#   ViralSeq::sdrm_pr_bulk
#   ViralSeq::sdrm_rt_bulk
#   ViralSeq::sdrm_in_bulk

# ViralSeq.sdrm_nrti(aa_arry, start_aa)
# ViralSeq.sdrm_nnrti(aa_arry, start_aa)
# ViralSeq.hiv_protease(aa_arry, start_aa)
# ViralSeq.sdrm_int(aa_arry, start_aa)
#   # funtions to identify SDRMs from a given sequence in an Array object
#   # function names indicate which HIV drug resistance mutations it can identify
#   # input an Array object for amino acid sequence ['A', 'M', 'L', ...]
#   # start_aa is an Integer to indicate codon number of the 1st amino acid sequence in the input aa_array
#   # return a Hash object for SDRMs identified. {:posiiton =>[:wildtype_codon, :mutation_codon]}

# ViralSeq.sdrm_pr_bulk(sequence_hash, minority_cut_off)
# ViralSeq.sdrm_rt_bulk(sequence_hash, minority_cut_off)
# ViralSeq.sdrm_in_bulk(sequence_hash, minority_cut_off)
#   # functions to identify SDRMs from a sequence hash object.
#   # name of the functions indicate which region it works on
#   # works for MPID-DR protocol (dx.doi.org/10.17504/protocols.io.useewbe)
#   # PR codon 1-99
#   # RT codon 34-122, 152-236, two regions are linked
#   # IN codon 53-174
#   # sequence_hash is a Hash object of sequences {:name => :sequence, ...}
#   # sequences usually need to be QCed (remove sequences with stop codon and a3g hypermutations) first
#   # minority_cut_off is the Integer cut-off for minimal abundance of a mutation to be called as valid mutation
#   # minority_cut_off can be obtained using ViralSeq::poisson_minority_cutoff function
#   # return [point_mutation_list, linkage_list, report_list]
# =USAGE
#   # example (example files from ID:VS053118-0566)
#   sequence = ViralSeq.fasta_to_hash('spec/sample_files/sample_dr_sequences/pr.fasta')
#   p_cut_off = ViralSeq.poisson_minority_cutoff(sequences)
#   pr_sdrm = ViralSeq.sdrm_pr_bulk(sequence, p_cut_off)
#   puts "region,tcs_number,position,wildtype,mutation,count,%,CI_low,CI_high,label"
#   pr_sdrm[0].each {|n| puts n.join(',')}
#   => region,tcs_number,position,wildtype,mutation,count,%,CI_low,CI_high,label
#   => PR,396,30,D,N,247,0.62374,0.57398,0.67163,
#   => PR,396,50,I,V,1,0.00253,6.0e-05,0.01399,*
#   => PR,396,88,N,D,246,0.62121,0.57141,0.66919,
#
#   puts "region,tcs_number,linkage,count,%,CI_low,CI_high,label"
#   pr_sdrm[1].each {|n| puts n.join(',')}
#   => region,tcs_number,linkage,count,%,CI_low,CI_high,label
#   => PR,396,D30N+N88D,245,0.61869,0.56884,0.66674,
#   => PR,396,WT,149,0.37626,0.32837,0.42602,
#   => PR,396,D30N,1,0.00253,6.0e-05,0.01399,*
#   => PR,396,D30N+I50V+N88D,1,0.00253,6.0e-05,0.01399,*
#
#   puts "position,codon,tcs_number," + ViralSeq::AMINO_ACID_LIST.join(",")
#   pr_sdrm[2].each {|n|puts n.join(",")}
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



module ViralSeq

  # drug resistant mutation summary. input: amino acid array and starting codon, output, hash of summary
  def self.sdrm_nrti(aa_array,start_aa=1)
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

  def self.sdrm_nnrti(aa_array,start_aa=1)
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

  #HIV protease surveillance mutations

  def self.hiv_protease(aa_array,start_aa=1)
    out_hash = {}
    sdrm = {}
    sdrm[23] = ['L',['I']]
    sdrm[24] = ['L',['I']]
    sdrm[30] = ['D',['N']]
    sdrm[32] = ['V',['I']]
    sdrm[46] = ['M',['I','L','V']] # M46V not on the SDRM list but we still include it.
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

  def self.sdrm_int(aa_array,start_aa=1)
    out_hash = {}
    sdrm = {}
    sdrm[66] = ['T',['A','I','K']]
    sdrm[74] = ['L',['M']]
    sdrm[92] = ['E',['Q']]
    sdrm[95] = ['Q',['K']]
    sdrm[97] = ['T',['A']]
    sdrm[121] = ['F',['Y']]
    sdrm[140] = ['G',['A','S','C']]
    sdrm[143] = ["Y",["C","H","R"]]
    sdrm[147] = ['S',['G']]
    sdrm[148] = ['Q',['H','K','R']]
    sdrm[155] = ['N',['S','H']]
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

  # input sequence hash, and Poisson cutoff for minority variants.
  # HIV-1 PR region SDRM based on HIVDB.stanford.edu
  # only for MPID-DR MiSeq sequences, PR codon 1-99
  # return [substitution rate with 95% CI, halpotype abundance with 95% CI, amino acid sequence report spreadsheet]
  def self.sdrm_pr_bulk(sequences, cutoff = 0)
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
      s.get_aa_array(rf_label)
      aa_seq = s.aa_array
      aa[name] = aa_seq.join("")
      record = ViralSeq.hiv_protease(aa_seq)
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
      count_mut_list = ViralSeq.count(mut_list)
      count_mut_list.each do |m,number|
        ci = ViralSeq.r_binom_CI(number, n_seq)
        label = number < cutoff ? "*" : ""
        point_mutation_list << [region, n_seq, position, wt, m, number, (number/n_seq.to_f).round(5), ci[0], ci[1], label]
      end
    end
    point_mutation_list.sort_by! {|record| record[2]}

    link = ViralSeq.count(mut_com)
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
      ci = ViralSeq.r_binom_CI(v, n_seq)
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
      count_aas = ViralSeq.count(aas)
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


  #input sequence hash, and Poisson cutoff for minority variants.
  #HIV-1 RT region SDRM based on HIVDB.stanford.edu
  #only for MPID-DR MiSeq sequences
  #RT codon 34-122, 152-236 two regions are linked.
  #return [substitution rate with 95% CI, halpotype abundance with 95% CI, amino acid sequence report spreadsheet]
  def self.sdrm_rt_bulk(sequences, cutoff = 0)
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
      s.get_aa_array(rf_label)
      aa_seq = s.aa_array

      r1_aa[name] = aa_seq[0,89].join("")
      r2_aa[name] = aa_seq[-85..-1].join("")
      nrti = ViralSeq.sdrm_nrti(aa_seq,start_codon_number)
      nnrti = ViralSeq.sdrm_nnrti(aa_seq,start_codon_number)
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
      count_mut_list = ViralSeq.count(mut_list)
      count_mut_list.each do |m,number|
        ci = ViralSeq.r_binom_CI(number, n_seq)
        label = number < cutoff ? "*" : ""
        point_mutation_list << ["NRTI", n_seq, position, wt, m, number, (number/n_seq.to_f).round(5), ci[0], ci[1], label]
      end
    end

    mut_nnrti.each do |position,mutation|
      wt = mutation[0]
      mut_list = mutation[1]
      count_mut_list = ViralSeq.count(mut_list)
      count_mut_list.each do |m,number|
        ci = ViralSeq.r_binom_CI(number, n_seq)
        label = number < cutoff ? "*" : ""
        point_mutation_list << ["NNRTI", n_seq, position, wt, m, number, (number/n_seq.to_f).round(5), ci[0], ci[1], label]
      end
    end
    point_mutation_list.sort_by! {|record| record[2]}

    link = ViralSeq.count(mut_com)
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
      ci = ViralSeq.r_binom_CI(v, n_seq)
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
      count_aas = ViralSeq.count(aas)
      div_aa[r1_aa_start] = count_aas.sort_by{|_k,v|v}.reverse.to_h
      r1_aa_start += 1
    end

    (0..r2_aa_size).to_a.each do |p|
      aas = []
      r2_aa.values.each do |r1|
        aas << r1[p]
      end
      count_aas = ViralSeq.count(aas)
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

  #input sequence hash, and Poisson cutoff for minority variants.
  #HIV-1 IN region SDRM based on HIVDB.stanford.edu
  #only for MPID-DR MiSeq sequences
  #IN codon 53-174
  #return [substitution rate with 95% CI, halpotype abundance with 95% CI, amino acid sequence report spreadsheet]
  def self.sdrm_in_bulk(sequences, cutoff = 0)
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
      s.get_aa_array(rf_label)
      aa_seq = s.aa_array
      aa[name] = aa_seq.join("")
      record = ViralSeq.sdrm_int(aa_seq, start_codon_number)
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
      count_mut_list = ViralSeq.count(mut_list)
      count_mut_list.each do |m,number|
        ci = ViralSeq.r_binom_CI(number, n_seq)
        label = number < cutoff ? "*" : ""
        point_mutation_list << [region, n_seq, position, wt, m, number, (number/n_seq.to_f).round(5), ci[0], ci[1], label]
      end
    end
    point_mutation_list.sort_by! {|record| record[2]}

    link = ViralSeq.count(mut_com)
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
      ci = ViralSeq.r_binom_CI(v, n_seq)
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
      count_aas = ViralSeq.count(aas)
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

end
