library(ggplot2)
library(plyr)
library(dplyr)
library(devtools)
library(phyloseq)
library(ape)
library(vegan)
library(rbiom)
library(metagMisc)

setwd("/users/corinna/Documents/Work/Beinart_Lab/Mollusk_symbioses_bottleneck/Lau_16S_FL_amplicons")

# Import final ASV biom and mapping files
biomFile <- import_biom("Lau_16S_exact_zotus_tax.biom", parseFunction = parse_taxonomy_default)
mapFile <- read.table("Lau_16S_metadata.txt", header=TRUE, row.names=1)
sampleData = sample_data(mapFile)
biomMapFile <- merge_phyloseq(biomFile, sampleData)

# Import representative sequences and remove non-ASV information
repsetFile <- Biostrings::readDNAStringSet("Lau_16S_zotus.fa")
names(repsetFile) <- gsub("\\s.+$", "", names(repsetFile))

# Create full phyloseq object
phyloseq <- merge_phyloseq(biomMapFile, repsetFile)
colnames(tax_table(phyloseq)) = c("Domain", "Phylum", "Class", "Order", "Family", "Genus")

# Remove low abundance ASVs
marmic = filter_taxa(phyloseq, function (x) sum(x) > 0, TRUE)
marmic = filter_taxa(marmic, function (x) {sum(x > 0) > 1}, TRUE)
marmic = phyloseq_filter_prevalence(marmic, prev.trh = 0.05, abund.trh = 10, threshold_condition = "AND", abund.type = "total")

write.table(otu_table(marmic), file="Lau_16S_exact_zotus_clean.txt", sep="\t")