module ViralSeq

  # array for all amino acid one letter abbreviations
  AMINO_ACID_LIST = ["A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y", "*"]

  # R script to check if the required libraries are installed, if not, install the missing libraries
  
  R_SCRIPT_CHECK_PACKAGES = 'dir.create(path = Sys.getenv("R_LIBS_USER"), showWarnings = FALSE, recursive = TRUE);' +
                            'packages <- c("ggplot2", "phangorn", "ape", "scales", "ggforce", "cowplot", "magrittr", "gridExtra");' +
                            'install.packages(setdiff(packages, rownames(installed.packages())), lib = Sys.getenv("R_LIBS_USER"), repos = "https://cran.rstudio.com/")'


  # R script for tcs_sdrm script

  R_SCRIPT = 'setwd("PATH_TO_FASTA")
              library(phangorn)
              library(ape)
              library(ggplot2)
              library(scales)
              library(ggforce)
              library(cowplot)
              library(magrittr)
              library(gridExtra)
              pdf("OUTPUT_PDF", onefile=T, width=11, height=8.5)
              fileNames <- list.files()
              for (fileName in fileNames) {
              dna <- read.dna(fileName, format="fasta")
              class(dna)
              D<- dist.dna(dna, model="raw")
              pi <- mean(D)
              dist20 <- quantile(D, prob=c(0.20))
              alldist <- data.frame(File=fileName, pi, dist20)
              write.table(alldist,"OUTPUT_CSV",append=TRUE, sep = ",", row.names = FALSE, col.names=FALSE)
              D2 <- dist.dna(dna, model="TN93")*100
              def.par <- par(no.readonly = TRUE)
              par(mfrow=c(1,2))
              hist<-hist(D, main=fileName, xlab="% Pairwise Distance", ylab="Frequency", col="gray")
              abline(v=dist20, col="royalblue",lwd=2)
              abline(v=pi, col="red", lwd=2)
              legend(x="topright", c("dist20", "pi"), col = c("royalblue", "red"), lwd = c(2,2), cex=0.5)
              njtree<-NJ(D2)
              njtreeplot <- plot(njtree, show.tip.label=F, "unrooted", main=fileName)
              add.scale.bar(cex=0.7, font=2, col="red")
              }
              dev.off()'


end
