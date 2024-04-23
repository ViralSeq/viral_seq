setwd("PATH_TO_FASTA")
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
dev.off()
