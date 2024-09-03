library(ggplot2)
library(cowplot)
library(patchwork)
library(respR)

setwd("/Users/corinna/Documents/Work/Beinart_Lab/Mollusk_symbioses_bottleneck/A_boucheti/host")

data <- read.table("Thetas_per_site-BW.txt", header=T)
data$Theta <- exp(data$Theta)

p <- ggplot(data, aes(x=Vent, y=Theta, fill=Vent)) + ylim(0, 0.6) + geom_boxplot(outlier.shape=NA) + theme_classic() + scale_fill_manual(values=c("maroon4", "cornflowerblue")) + labs(x="Time", y="Theta per site") + scale_x_discrete(limits=c("Before", "After")) + theme(legend.position = "none", text = element_text(size = 15))

pdf("Theta_per_site-BW.pdf",  width=5, height=5)
p
dev.off()

Before <- subset(data, data$Vent=='Before')
Befores <- subsample(Before, length.out=100, plot=FALSE)
After <- subset(data, data$Vent=='After')
Afters <- subsample(After, length.out=100, plot=FALSE)

df <- rbind(Befores, Afters)

sink("Theta_Site_Statistics-BW.txt")
shapiro.test(df$Theta)
pairwise.t.test(df$Theta, df$Vent, p.adjust.method = "fdr", pool.sd = FALSE)
pairwise.wilcox.test(df$Theta, df$Vent, p.adjust.method = "fdr", pool.sd = FALSE)
sink()

data <- read.table("TajimasD-BW.txt", header=T)

p <- ggplot(data, aes(x=Vent, y=Tajima, fill=Vent)) + ylim(-2.5, 4) + geom_violin(aes(fill = Vent), size = 0.5, bw = 0.2) + geom_boxplot(fill = "white", size = 0.5, width = 0.2, outlier.shape=NA) + theme_classic() + scale_fill_manual(values=c("maroon4", "cornflowerblue")) + labs(x="Time", y="Tajima's D") + scale_x_discrete(limits=c("Before", "After")) + theme(legend.position = "none", text = element_text(size = 15))

pdf("TajimasD-BW.pdf",  width=5, height=5)
p
dev.off()

Before <- subset(data, data$Vent=='Before')
Befores <- subsample(Before, length.out=100, plot=FALSE)
After <- subset(data, data$Vent=='After')
Afters <- subsample(After, length.out=100, plot=FALSE)

df <- rbind(Befores, Afters)

sink("Tajima_Statistics-BW.txt")
shapiro.test(df$Tajima)
pairwise.t.test(df$Tajima, df$Vent, p.adjust.method = "fdr", pool.sd = FALSE)
pairwise.wilcox.test(df$Tajima, df$Vent, p.adjust.method = "fdr", pool.sd = FALSE)
sink()

data <- read.table("Het-BW.txt", header=T)

p <- ggplot(data, aes(x=Vent, y=Hobs, fill=Vent)) + ylim(0, 1) + geom_boxplot(outlier.shape=NA) + geom_point(data = aggregate(Hexp ~ Vent, data = data, median), aes(x = Vent, y = Hexp), color = "black", size = 3) + theme_classic() + scale_fill_manual(values=c("maroon4", "cornflowerblue")) + labs(x="Time", y="Heterozygosity") + scale_x_discrete(limits=c("Before", "After")) + theme(legend.position = "none", text = element_text(size = 15))

pdf("Het-BW.pdf",  width=5, height=5)
p
dev.off()

Before <- subset(data, data$Vent=='Before')
Befores <- subsample(Before, length.out=100, plot=FALSE)
After <- subset(data, data$Vent=='After')
Afters <- subsample(After, length.out=100, plot=FALSE)

df <- rbind(Befores, Afters)

sink("Het_Statistics-BW.txt")
shapiro.test(df$Hobs)
pairwise.t.test(df$Hobs, df$Vent, p.adjust.method = "fdr", pool.sd = FALSE)
pairwise.wilcox.test(df$Hobs, df$Vent, p.adjust.method = "fdr", pool.sd = FALSE)
sink()

