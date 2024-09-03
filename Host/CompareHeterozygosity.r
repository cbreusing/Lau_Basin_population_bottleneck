library(ggplot2)
library(cowplot)
library(patchwork)
library(respR)

setwd("/Users/corinna/Documents/Work/Beinart_Lab/Mollusk_symbioses_bottleneck/A_boucheti/host")

data <- read.table("Het-After.txt", header=T)

Hexp <- subset(data, data$Type=='Hexp')
HexpS <- subsample(Hexp, length.out=100, plot=FALSE)
Hobs <- subset(data, data$Type=='Hobs')
HobsS <- subsample(Hobs, length.out=100, plot=FALSE)

df <- rbind(HexpS, HobsS)

sink("Het_Statistics-After.txt")
shapiro.test(df$Het)
pairwise.t.test(df$Het, df$Type, p.adjust.method = "fdr", pool.sd = FALSE)
pairwise.wilcox.test(df$Het, df$Type, p.adjust.method = "fdr", pool.sd = FALSE)
sink()

data <- read.table("Het-Before.txt", header=T)

Hexp <- subset(data, data$Type=='Hexp')
HexpS <- subsample(Hexp, length.out=100, plot=FALSE)
Hobs <- subset(data, data$Type=='Hobs')
HobsS <- subsample(Hobs, length.out=100, plot=FALSE)

df <- rbind(HexpS, HobsS)

sink("Het_Statistics-Before.txt")
shapiro.test(df$Het)
pairwise.t.test(df$Het, df$Type, p.adjust.method = "fdr", pool.sd = FALSE)
pairwise.wilcox.test(df$Het, df$Type, p.adjust.method = "fdr", pool.sd = FALSE)
sink()
