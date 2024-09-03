library(ggplot2)

setwd("/Users/corinna/Documents/Work/Beinart_Lab/Mollusk_symbioses_bottleneck/A_boucheti/host")

ABEh <- read.table("A_boucheti-ABE-Allio.final.summary",header=T)
TMh <- read.table("A_boucheti-TM-Allio.final.summary",header=T)

pdf("StairwayPlot.FINAL.reduced.pdf",height = 6.25,width = 10)

options(scipen=999)

plot(TMh$year, log(TMh$Ne_median), type="n", xlab="Time (years)", ylab="log(Ne)", xlim=c(0,1000), ylim=c(5,15))

lines(TMh$year,log(TMh$Ne_median),type="s",col="darkgoldenrod1",lwd = 5)
lines(TMh$year,log(TMh$Ne_2.5.),type="s",col="darkgoldenrod1",lty=3)
lines(TMh$year,log(TMh$Ne_97.5.),type="s",col="darkgoldenrod1",lty=3)

lines(ABEh$year,log(ABEh$Ne_median),type="s",col="royalblue3",lwd = 5)
lines(ABEh$year,log(ABEh$Ne_2.5.),type="s",col="royalblue3",lty=3)
lines(ABEh$year,log(ABEh$Ne_97.5.),type="s",col="royalblue3",lty=3)


legend("topright",legend = c("Tu'i Malila","ABE","95% CI"),col=c("darkgoldenrod1","royalblue3","black"),lty=c(1,1,3),lwd=c(5,5,1),cex=0.8)

dev.off()
