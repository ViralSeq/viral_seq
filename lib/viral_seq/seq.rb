module ViralSeq
  #reverse complement
  class String
      def rc
          self.reverse.tr("ACTG","TGAC")
      end
  end
end
