RSpec.describe ViralSeq do
  it "has a version number" do
    expect(ViralSeq::VERSION).not_to be nil
  end

  it "reverse complement" do
    expect("ACTG".rc).to eq "CAGT"
  end

  it "amino acid list" do
    expect(ViralSeq::AMINO_ACID_LIST[2]).to eq "D"
  end
end
