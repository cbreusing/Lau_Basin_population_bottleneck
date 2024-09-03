library(ggplot2)
library(cowplot)
library(patchwork)
library(respR)

setwd("/Users/corinna/Documents/Work/Beinart_Lab/Mollusk_symbioses_bottleneck/A_boucheti/symbiont/basin-wide_analyses")

data <- read.table("Inter-Pi-BW.txt", header=T)

p <- ggplot(data, aes(x=site, y=inter_pi, fill=site)) + geom_boxplot(outlier.shape=NA) + theme_classic() + scale_fill_manual(values=c("maroon4", "cornflowerblue")) + labs(x="Time", y="Inter-host nucleotide diversity") + scale_x_discrete(limits=c("Before", "After")) + theme(legend.position = "none", text = element_text(size = 15))

pdf("Boxplot_inter-pi-BW.pdf", width=5, height=5)
p
dev.off()

sink("InterPi_Statistics.txt")
shapiro.test(data$inter_pi)
pairwise.t.test(data$inter_pi, data$site, p.adjust.method = "fdr", pool.sd = FALSE)
pairwise.wilcox.test(data$inter_pi, data$site, p.adjust.method = "fdr", pool.sd = FALSE)
sink()

data <- read.table("Site-Pi-BW.txt", header=T)

p <- ggplot(data, aes(x=site, y=pi, fill=site)) + ylim(0, 0.6) + geom_boxplot(outlier.shape=NA) + theme_classic() + scale_fill_manual(values=c("maroon4", "cornflowerblue")) + labs(x="Time", y="Nucleotide diversity") + scale_x_discrete(limits=c("Before", "After")) + theme(legend.position = "none", text = element_text(size = 15))

pdf("Boxplot_site-pi-BW.pdf", width=5, height=5)
p
dev.off()

sink("SitePi_Statistics.txt")
shapiro.test(data$pi)
pairwise.t.test(data$pi, data$site, p.adjust.method = "fdr", pool.sd = FALSE)
pairwise.wilcox.test(data$pi, data$site, p.adjust.method = "fdr", pool.sd = FALSE)
sink()

data <- read.table("TajimaD-BW.txt", header=T)

p <- ggplot(data, aes(x=site, y=tajima, fill=site)) + ylim(-3, 4) + geom_violin(aes(fill = site), size = 0.5, bw = 0.2) + geom_boxplot(fill = "white", size = 0.5, width = 0.2, outlier.shape=NA) + theme_classic() + scale_fill_manual(values=c("maroon4", "cornflowerblue")) + labs(x="Time", y="Tajima's D") + scale_x_discrete(limits=c("Before", "After")) + theme(legend.position = "none", text = element_text(size = 15))

pdf("TajimaD-BW.pdf", width=5, height=5)
p
dev.off()

sink("Tajima_Statistics.txt")
shapiro.test(data$tajima)
pairwise.t.test(data$tajima, data$site, p.adjust.method = "fdr", pool.sd = FALSE)
pairwise.wilcox.test(data$tajima, data$site, p.adjust.method = "fdr", pool.sd = FALSE)
sink()
