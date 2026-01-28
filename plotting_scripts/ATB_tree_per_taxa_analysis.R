library(phytools)
library(viridis)
library(svglite)
library(ggplot2)
library(MCMCglmm)
library(dplyr)
library(stringr)
library(ape)
library(coda)
library(tools)
library(tidyr)
library(paletteer)
library(ggsci)
library(ggpubr)

overall_indir <- ""
indir <- paste0(overall_indir, "ATB_tree/")
tree.file <- paste0(indir, "data/", "reps_ppmod_tree.treefile")
summary.file <- paste0(indir, "data/", "ppmod_summary.csv")
index.file <- paste0(indir, "data/", "sampled_alignments.tsv")
gtdb.file <- paste0(indir, "data/", "gtdbtk.bac120.summary.tsv")
outpref <- paste0(indir, "figures/")

summary.df <- read.csv(summary.file, row.names = 1)

values.NA <- summary.df[is.na(summary.df$order),]

index.df <- read.csv(index.file, sep="\t", header = FALSE)
colnames(index.df) <- c("taxa", "tip", "file_path")

meta <- merge(summary.df, index.df, by = "taxa")

# read in taxonomy
gtdb.df <- read.csv(gtdb.file, sep="\t", header = TRUE)
gtdb.df <- select(gtdb.df, c('user_genome', 'classification'))
gtdb.df <- gtdb.df |> separate_wider_delim(classification, delim = ";", names = c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species"))
gtdb.df$Domain <- gsub(".*__", "", gtdb.df$Domain)
gtdb.df$Phylum <- gsub(".*__", "", gtdb.df$Phylum)
gtdb.df$Class <- gsub(".*__", "", gtdb.df$Class)
gtdb.df$Order <- gsub(".*__", "", gtdb.df$Order)
gtdb.df$Family <- gsub(".*__", "", gtdb.df$Family)
gtdb.df$Genus <- gsub(".*__", "", gtdb.df$Genus)
gtdb.df$Species <- gsub(".*__", "", gtdb.df$Species)

gtdb.df <- as.data.frame(gtdb.df)

meta <- merge(meta, gtdb.df, by.x = "tip", by.y = "user_genome")

traits <- meta[, c("rate_genes1_median", "core_mu_median", "prop_genes2_median", "Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")]
colnames(traits) <- c("Basal gene turnover rate", "Core mutation rate", "Proportion fast genes", "Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")

traits.stacked <- traits %>%
  pivot_longer(
    cols = c("Basal gene turnover rate", "Core mutation rate", "Proportion fast genes"),
    names_to = "Trait",
    values_to = "Value"
  )

base_cols <- paletteer_d("ggsci::category20_d3")

tax <- meta$Phylum
names(tax) <- rownames(meta)
tax <- factor(tax)

tax_cols <- setNames(
  c(base_cols[1:length(levels(tax))]),
  levels(tax)
)

# generate summary statistics of traits
p <- ggplot(traits.stacked, aes(x=Phylum, y = Value, colour = Phylum)) + 
  theme_light() + 
  facet_grid(Trait ~ ., scales = "free_y") +
  geom_boxplot() +
  scale_y_log10() +
  xlab("Phylum") + ylab("Trait Value") + 
  scale_color_manual(values = tax_cols) +
  stat_compare_means(method = "kruskal.test") +
  theme(axis.text.x = element_text(size = 18, angle = 45, hjust = 1, vjust=1), axis.text.y = element_text(size = 18), axis.title=element_text(size=20,face="bold"), legend.position = "none", strip.text.y = element_text(size=16, face="bold"))
p
ggsave(paste(outpref, "per_taxa_trait_values.png", sep=""), plot = p, height = 12, width = 10)
