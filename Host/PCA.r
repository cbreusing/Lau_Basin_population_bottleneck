library(vegan)
library(ggplot2)
library(ape)
library(tidyverse)
library(poppr)
library(pegas)
library(adegenet)
library(ade4)
library(goeveg)

setwd("/Users/corinna/Documents/Work/Beinart_Lab/Mollusk_symbioses_bottleneck/A_boucheti/host")

#IBS: Identical by descent/state, average genotype distance between individuals, ignores allele frequencies, used in MDS
data <- as.matrix(read.table("A_boucheti.ibsMat"))
data[is.na(data)] = 0 
#COV: Covariance matrix, weighting scheme based on allele frequencies, used in PCA
data2 <- as.matrix(read.table("A_boucheti.covMat"))
data2[is.na(data2)] = 0
pop <- read.table("A_boucheti.clst", header=T)
pcoa <- pcoa(data, correction = "cailliez")
pca <- prcomp(data2, scale=TRUE, center=TRUE)

pdf("Screeplot.pdf")
screeplot_NMDS(data, distance = "euclidean", k = 6, trymax = 100, autotransform = FALSE)
dev.off()

nmds <- metaMDS(as.dist(data), trymax=1000, autotransform=FALSE, k=4)

Axis1 <- pcoa$vectors[,1]
Axis2 <- pcoa$vectors[,2]
var1 <- round(pcoa$values[1,2], digits=4)*100
var2 <- round(pcoa$values[2,2], digits=4)*100
x_axis <- paste("Axis1"," ","[",var1,"%]",sep="")
y_axis <- paste("Axis2"," ","[",var2,"%]",sep="")
PC1 <- pca$rotation[,1]
PC2 <- pca$rotation[,2]
eigs <- pca$sdev^2
pc1_var <- round(eigs[1]/sum(eigs), digits=4)*100
pc2_var <- round(eigs[2]/sum(eigs), digits=4)*100
pc1_axis <- paste("PC1"," ","[",pc1_var,"%]",sep="")
pc2_axis <- paste("PC2"," ","[",pc2_var,"%]",sep="")
scores <- as.data.frame(scores(nmds))
stress <- nmds$stress

plot1 <- cbind(data, Axis1, Axis2, pop)
plot2 <- cbind(data2, PC1, PC2, pop)
plot3 <- cbind(data, scores, pop)

plot1$Vent <- factor(plot1$Vent, levels = c("Kilo Moana", "Tow Cam", "ABE", "Tu'i Malila"))
plot2$Vent <- factor(plot2$Vent, levels = c("Kilo Moana", "Tow Cam", "ABE", "Tu'i Malila"))
plot3$Vent <- factor(plot3$Vent, levels = c("Kilo Moana", "Tow Cam", "ABE", "Tu'i Malila"))

col <- c("After" = "maroon4", "Before" = "cornflowerblue")

pdf("A_boucheti_PCoA.12.Eruption.pdf")
p <- ggplot(plot1, aes(Axis1, Axis2, color=Eruption, shape=Vent)) + geom_point(position=position_jitter(width=0, height=0), size=6) + theme_classic() + theme(text = element_text(size = 15)) + scale_shape_manual(values = c(8, 18, 16, 17), name="Vent") + scale_color_manual(values = col, name="Eruption") + labs(x = x_axis, y = y_axis)
p
dev.off()

pdf("A_boucheti_PCA.12.Eruption.pdf")
p <- ggplot(plot2, aes(PC1, PC2, color=Eruption, shape=Vent)) + geom_point(position=position_jitter(width=0, height=0), size=6) + theme_classic() + theme(text = element_text(size = 15)) + scale_shape_manual(values = c(8, 18, 16, 17), name="Vent") + scale_color_manual(values = col, name="Eruption") + labs(x = pc1_axis, y = pc2_axis)
p
dev.off()

pdf("A_boucheti_NMDS.12.Eruption.pdf")
p <- ggplot(plot3, aes(NMDS1, NMDS2, color=Eruption, shape=Vent)) + geom_point(position=position_jitter(width=0, height=0), size=6) + theme_classic() + theme(text = element_text(size = 15)) + scale_shape_manual(values = c(8, 18, 16, 17), name="Vent") + scale_color_manual(values = col, name="Eruption")
p
dev.off()

dist <- as.dist(data)
dist2 <- dist(data2, method="euclidean")
df <- data.frame(pop)

sink("PERMANOVA.txt")
# Adonis test
adonis2(dist ~ Eruption, data = df, add = "cailliez", permutations = 999, by = NULL)
adonis2(dist2 ~ Eruption + Vent, data = df, add = "cailliez", permutations = 999, by = NULL)

# Dispersion test
beta1 <- betadisper(dist, df$Eruption)
permutest(beta1)
beta2 <- betadisper(dist2, df$Eruption)
permutest(beta2)
sink()