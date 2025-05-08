library(vegan)
library(ggplot2)
library(ape)
library(tidyverse)
library(poppr)
library(pegas)
library(adegenet)
library(ade4)
library(ggrepel)

setwd("/Users/corinna/Documents/Work/Beinart_Lab/Mollusk_symbioses_bottleneck/A_boucheti/symbiont/basin-wide_analyses")

data <- read.table("Epsilon.Freebayes.FINAL.CT.FORMAT", header=T, row.names=1)
com <- t(apply(data, 2, function(i) i/sum(i, na.rm=TRUE)))
data2 <- t(read.table("Epsilon.Freebayes.FINAL.GT.FORMAT", header=T, row.names=1))
pop <- read.table("Epsilon.clst", header=T)

dist <- vegdist(com, method="bray", binary=F, na.rm=TRUE)
matrix <- df2genind(data2, ploidy = 1, ncode = 1)
dist2 <- vegdist(matrix, method="euclidean", binary=F, na.rm=TRUE)
dist2[is.na(dist2)] = 0
pcoa <- pcoa(dist, correction = "cailliez")
pcoa2 <- pcoa(dist2, correction = "cailliez")

nmds1 <- metaMDS(dist, trymax=1000, autotransform=FALSE, k=2)
nmds2 <- metaMDS(dist2, trymax=1000, autotransform=FALSE, k=2)

Axis1 <- pcoa$vectors[,1]
Axis2 <- pcoa$vectors[,2]
var1 <- round(pcoa$values[1,2], digits=4)*100
var2 <- round(pcoa$values[2,2], digits=4)*100
x_axis <- paste("Axis1"," ","[",var1,"%]",sep="")
y_axis <- paste("Axis2"," ","[",var2,"%]",sep="")
Axis1b <- pcoa2$vectors[,1]
Axis2b <- pcoa2$vectors[,2]
var1b <- round(pcoa2$values[1,2], digits=4)*100
var2b <- round(pcoa2$values[2,2], digits=4)*100
x_axisb <- paste("Axis1"," ","[",var1b,"%]",sep="")
y_axisb <- paste("Axis2"," ","[",var2b,"%]",sep="")
scores1 <- as.data.frame(scores(nmds1))
stress1 <- nmds1$stress
scores2 <- as.data.frame(scores(nmds2))
stress2 <- nmds2$stress

plot1 <- cbind(com, Axis1, Axis2, pop)
plot2 <- cbind(matrix, Axis1b, Axis2b, pop)
plot3 <- cbind(com, scores1, pop)
plot4 <- cbind(matrix, scores2, pop)

plot1$Site <- factor(plot1$Site, levels = c("Kilo Moana", "Tow Cam", "ABE", "Tu'i Malila"))
plot2$Site <- factor(plot2$Site, levels = c("Kilo Moana", "Tow Cam", "ABE", "Tu'i Malila"))
plot3$Site <- factor(plot3$Site, levels = c("Kilo Moana", "Tow Cam", "ABE", "Tu'i Malila"))
plot4$Site <- factor(plot4$Site, levels = c("Kilo Moana", "Tow Cam", "ABE", "Tu'i Malila"))

col <- c("After" = "maroon4", "Before" = "cornflowerblue")

pdf("Epsilon_PCoA.12.Eruption.pdf")
p <- ggplot(plot1, aes(Axis1, Axis2, color=Eruption, shape=Site)) + geom_point(position=position_jitter(width=0, height=0), size=6) + theme_classic() + theme(text = element_text(size = 15)) + scale_shape_manual(values = c(8, 18, 16, 17), name="Site") + scale_color_manual(values = col, name="Eruption") + labs(x = x_axis, y = y_axis)
p
dev.off()

pdf("Epsilon_PCoA.geno.12.Eruption.pdf")
p <- ggplot(plot2, aes(Axis1b, Axis2b, color=Eruption, shape=Site)) + geom_point(position=position_jitter(width=0, height=0), size=6) + theme_classic() + theme(text = element_text(size = 15)) + scale_shape_manual(values = c(8, 18, 16, 17), name="Site") + scale_color_manual(values = col, name="Eruption") + labs(x = x_axisb, y = y_axisb)
p
dev.off()

pdf("Epsilon_NMDS.12.Eruption.pdf")
p <- ggplot(plot3, aes(NMDS1, NMDS2, color=Eruption, shape=Site)) + geom_point(position=position_jitter(width=0, height=0), size=6) + theme_classic() + theme(text = element_text(size = 15)) + scale_shape_manual(values = c(8, 18, 16, 17), name="Site") + scale_color_manual(values = col, name="Eruption")
p
dev.off()

pdf("Epsilon_NMDS.geno.12.Eruption.pdf")
p <- ggplot(plot4, aes(NMDS1, NMDS2, color=Eruption, shape=Site)) + geom_point(position=position_jitter(width=0, height=0), size=6) + theme_classic() + theme(text = element_text(size = 15)) + scale_shape_manual(values = c(8, 18, 16, 17), name="Site") + scale_color_manual(values = col, name="Eruption")
p
dev.off()

df <- data.frame(pop)

sink("PERMANOVA.txt")
# Adonis test
adonis2(dist ~ Eruption, data = df, add = "cailliez", permutations = 999, by = NULL)
adonis2(dist2 ~ Eruption, data = df, add = "cailliez", permutations = 999, by = NULL)

# Dispersion test
beta1 <- betadisper(dist, df$Eruption)
permutest(beta1)
beta2 <- betadisper(dist2, df$Eruption)
permutest(beta2)
sink()