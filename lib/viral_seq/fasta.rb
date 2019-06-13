# fasta.rb
# methods for converting sequence formats, including
#   ViralSeq::fasta_to_hash
#   ViralSeq::fastq_to_fasta
#   ViralSeq::fastq_to_hash
#   ViralSeq::fasta_hash_to_rsphylip
#   ViralSeq::pair_fasta_to_hash

# =USAGE
#   sequence_fasta_hash = ViralSeq.fasta_to_hash(input_fasta_file)
#   # input a sequence file in fasta format, read as a sequence hash
#   # {:sequence_name1 => sequence1, ...}

#   sequence_fasta_hash = ViralSeq.fastq_to_fasta(input_fastq_file)
#   # input a sequence file in fastq format, read as a sequence hash
#   # discard sequence quality score

#   sequence_fastq_hash = ViralSeq.fasta_to_hash(input_fastq_file)
#   # input a sequence file in fastq format, read as a sequence hash
#   # keep sequence quality score
#   # {:sequence_name1 => [sequence1, quality1], ...}

#   phylip_hash = ViralSeq.fasta_hash_to_rsphylip(sequence_fasta_hash)
#   # convert a aligned fasta sequence hash into relaxed sequencial phylip format

#   paired_sequence_hash = ViralSeq.pair_fasta_to_hash(directory_of_paired_fasta)
#   # input a directory containing paired sequence files in the fasta format
#   # ├───lib1
#         │     lib1_r1.txt
#         │     lib1_r2.txt
#   # paired sequence files need to have "r1" and "r2" in their file names
#   # the sequence taxa should only differ by last 3 characters to distinguish r1 and r2 sequence.
#   # return a paired sequence hash :seq_name => [r1_seq, r2_seq]

module ViralSeq

  def self.fasta_to_hash(infile)
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


  # fastq file to fasta, discard quality, return a sequence hash

  def self.fastq_to_fasta(fastq_file)
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
      Hash[*sequence_a]
  end

  # fastq file to hash, including quality. {:seq_name => [seq,quality]}

  def self.fastq_to_hash(fastq_file)
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

  # fasta sequence hash to relaxed sequencial phylip format

  def self.fasta_hash_to_rsphylip(seqs)
    outline = "\s" + seqs.size.to_s + "\s" + seqs.values[0].size.to_s + "\n"
    names = seqs.keys
    max_name_l = (names.max.size - 1)
    max_name_l > 10 ? name_block_l = max_name_l : name_block_l = 10
    seqs.each do |k,v|
      outline += k[1..-1] + "\s" * (name_block_l - k.size + 2) + v.scan(/.{1,10}/).join("\s") + "\n"
    end
    return outline
  end

  # input a directory with r1 and r2 sequences, return a hash :seq_name => [r1_seq, r2_seq]
  # r1 and r2 file names should contain "r1" and "r2" respectively
  # the sequence taxa should only differ by last 3 characters to distinguish r1 and r2 sequence.
  def self.pair_fasta_to_hash(indir)
    files = Dir[indir + "/*"]
    r1_file = ""
    r2_file = ""
    files.each do |f|
      if File.basename(f) =~ /r1/i
        r1_file = f
      elsif File.basename(f) =~ /r2/i
        r2_file = f
      end
    end

    seq1 = ViralSeq.fasta_to_hash(r1_file)
    seq2 = ViralSeq.fasta_to_hash(r2_file)

    new_seq1 = seq1.each_with_object({}) {|(k, v), h| h[k[0..-4]] = v}
    new_seq2 = seq2.each_with_object({}) {|(k, v), h| h[k[0..-4]] = v}

    seq_pair_hash = {}

    new_seq1.each do |seq_name,seq|
      seq_pair_hash[seq_name] = [seq, new_seq2[seq_name]]
    end
    return seq_pair_hash
  end
end
