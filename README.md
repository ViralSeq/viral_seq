# ViralSeq

A Ruby Gem containing bioinformatics tools for processing viral NGS data.

Specifically for Primer-ID sequencing and HIV drug resistance analysis.

## Installation

    $ gem install viral_seq

## Usage

#### Load all ViralSeq classes by requiring 'viral_seq.rb'

    #!/usr/bin/env ruby
    require 'viral_seq'
    
#### Use executable `locator` to get the coordinates of the sequences on HIV/SIV reference genome from a FASTA file through a terminal

    $ locator -i sequence.fasta -o sequence.fasta.csv

## Some Examples

#### Load nucleotide sequences from a FASTA format sequence file

    my_seqhash = ViralSeq::SeqHash.fa('my_seq_file.fasta')

#### Make an alignment (using MUSCLE)

    aligned_seqhash = my_seqhash.align

#### Filter nucleotide sequences with the reference coordinates (HIV Protease)

    qc_seqhash = aligned_seqhash.hiv_seq_qc(2253, 2549, false, :HXB2)

#### Further filter out sequences with Apobec3g/f hypermutations

    qc_seqhash = qc_seqhash.a3g

#### Calculate nucleotide diveristy π

    qc_seqhash.pi

#### Calculate cut-off for minority variants based on Poisson model

    cut_off = qc_seqhash.pm

#### Examine for drug resistance mutations for HIV PR region

    qc_seqhash.sdrm_hiv_pr(cut_off)

## Updates

Version 1.0.6-07232019:

    1. Several methods added to ViralSeq::SeqHash, including
        ViralSeq::SeqHash#size
        ViralSeq::SeqHash#+
        ViralSeq::SeqHash#write_nt_fa
        ViralSeq::SeqHash#mutation
    2. Update documentations and rspec samples. 

Version 1.0.5-07112019:

    1. Update ViralSeq::SeqHash#sequence_locator.
       Program will try to determine the direction (`+` or `-` of the query sequence)
    2. update executable `locator` to have a column of `direction` in output .csv file

Version 1.0.4-07102019:

    1. Use home directory (Dir.home) instead of the directory of the script file for temp MUSCLE file.
    2. Fix bugs in bin `locator`

Version 1.0.3-07102019:

    1. Bug fix.

Version 1.0.2-07102019:

    1. Fixed a gem loading issue.

Version 1.0.1-07102019:

    1. Add keyword argument :model to ViralSeq::SeqHashPair#join2.
    2. Add method ViralSeq::SeqHash#sequence_locator (also: #loc), a function to locate sequences on HIV/SIV reference genomes, as HIV Sequence Locator from LANL.
    3. Add executable 'locator'. An HIV/SIV sequence locator tool similar to LANL Sequence Locator.
    4. update documentations

Version 1.0.0-07092019:

    1. Rewrote the whole ViralSeq gem, grouping methods into modules and classes under main Module::ViralSeq

## Development

Bug reports and pull requests are welcome on GitHub at https://github.com/ViralSeq/viral_seq. This project is intended to be a safe, welcoming space for collaboration, and contributors are expected to adhere to the [Contributor Covenant](http://contributor-covenant.org) code of conduct.

## License

The gem is available as open source under the terms of the [MIT License](https://opensource.org/licenses/MIT).

## Code of Conduct

Everyone interacting in the viral_seq project’s codebases, issue trackers, chat rooms and mailing lists is expected to follow the [code of conduct](https://github.com/ViralSeq/viral_seq/blob/master/CODE_OF_CONDUCT.md).
