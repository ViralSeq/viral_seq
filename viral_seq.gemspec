
lib = File.expand_path("../lib", __FILE__)
$LOAD_PATH.unshift(lib) unless $LOAD_PATH.include?(lib)
require "viral_seq/version"

Gem::Specification.new do |spec|
  spec.name          = "viral_seq"
  spec.version       = ViralSeq::VERSION
  spec.authors       = ["Shuntai Zhou", "Michael Clark"]
  spec.email         = ["shuntai.zhou@gmail.com", "clarkmu@gmail.com"]

  spec.summary       = %q{A Ruby Gem containing bioinformatics tools for processing viral NGS data.}
  spec.description   = %q{A Ruby Gem with bioinformatics tools for processing viral NGS data.
                          Specifically for Primer-ID sequencing and HIV drug resistance analysis.}
  spec.homepage      = "https://github.com/ViralSeq/viral_seq"
  spec.license       = "MIT"

  # Specify which files should be added to the gem when it is released.
  # The `git ls-files -z` loads the files in the RubyGem that have been added into git.
  spec.files         = Dir.chdir(File.expand_path('..', __FILE__)) do
    `git ls-files -z`.split("\x0").reject { |f| f.match(%r{^(test|spec|features)/}) }
  end
  spec.bindir        = "bin"
  spec.executables   = spec.files.grep(%r{^bin/}) { |f| File.basename(f) }
  spec.require_paths = ["lib"]
  spec.post_install_message = "Thanks for installing!"

  spec.add_development_dependency "bundler", "~> 2.5"
  spec.add_development_dependency "rake", "~> 13.0"
  spec.add_development_dependency "rspec", "~> 3.0"

  # This gem will work with Ruby version 2.5.0 or greater...
  spec.required_ruby_version = '>= 2.5'
  # This gem will work with RubyGem version 1.3.6 or greater...
  spec.required_rubygems_version = '>= 1.3.6'

  # muscle_bio gem required
  spec.add_runtime_dependency "muscle_bio", "= 0.4"

  # virust-locator-ruby required
  spec.add_runtime_dependency "virust-locator-ruby", "~> 0.3"

  # colorize gem required
  spec.add_runtime_dependency "colorize", "~> 0.1"

  spec.add_runtime_dependency "prawn", "~> 2.3", '>= 2.3.0'

  spec.add_runtime_dependency "prawn-table", "~> 0.2", '>= 0.2.0'

  spec.add_runtime_dependency "combine_pdf", "~> 1.0", '>= 1.0.0'

  spec.requirements << 'R required for some functions'

  spec.add_dependency "shellwords", "~> 0.2"
end
