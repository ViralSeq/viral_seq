# Copyright (c) 2020 Shuntai Zhou (shuntai.zhou@gmail.com)
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.

module ViralSeq; end

# load all classes
require_relative "viral_seq/constant"
require_relative "viral_seq/enumerable"
require_relative "viral_seq/hash"
require_relative "viral_seq/hivdr"
require_relative "viral_seq/math"
require_relative "viral_seq/muscle"
require_relative "viral_seq/pid"
require_relative "viral_seq/ref_seq"
require_relative "viral_seq/rubystats"
require_relative "viral_seq/seq_hash"
require_relative "viral_seq/seq_hash_pair"
require_relative "viral_seq/sequence"
require_relative "viral_seq/string"
require_relative "viral_seq/version"

require "muscle_bio"
