library(ggplot2)
library(stringr)
library(dplyr)
library(ggpubr)
library(ggsci)
library(readxl)
library(Rtsne)
library(umap)
library(tools)

df <- read.csv("/Users/shorsfield/Library/Mobile Documents/com~apple~CloudDocs/Work/Postdoc/Analysis/PopPUNK-mod/publication_figures/inter_species_comparisons/poppunk_js_distances.tsv", sep = "\t")
str <- str_split(df$File1, "_distances_")

df$File1 <- basename(file_path_sans_ext(df$File1))
df$File2 <- basename(file_path_sans_ext(df$File2))

df %>% mutate(
  str_remove(File1, "_distances_.*")
) -> tmp.df
df$File1 <- tmp.df$`str_remove(File1, "_distances_.*")`

df %>% mutate(
  str_remove(File2, "_distances_.*")
) -> tmp.df

df$File2 <- tmp.df$`str_remove(File2, "_distances_.*")`

#rename files
{
  df$File1[df$File1 == "ab"] <- "A. baumanii"
  df$File1[df$File1 == "bp"] <- "B. pertussis"
  df$File1[df$File1 == "cj"] <- "C. jejuni"
  df$File1[df$File1 == "GBS_full"] <- "S. agalactiae"
  df$File1[df$File1 == "GPSv4"] <- "S. pneumoniae"
  df$File1[df$File1 == "e_feacalis"] <- "E. faecalis"
  df$File1[df$File1 == "e_feacium"] <- "E. faecium"
  df$File1[df$File1 == "ec"] <- "E. coli"
  df$File1[df$File1 == "hflu"] <- "H. influenzae"
  df$File1[df$File1 == "hp"] <- "H. pylori"
  df$File1[df$File1 == "kp"] <- "K. pneumoniae"
  df$File1[df$File1 == "lm"] <- "L. monocytogenes"
  df$File1[df$File1 == "lp"] <- "L. pneumophila"
  df$File1[df$File1 == "ma"] <- "M. abscessus"
  df$File1[df$File1 == "mcat"] <- "M. catarrhalis"
  df$File1[df$File1 == "mtb"] <- "M. tuberculosis"
  df$File1[df$File1 == "ngon"] <- "N. gonorrhoeae"
  df$File1[df$File1 == "nm"] <- "N. meningitidis"
  df$File1[df$File1 == "pa"] <- "P. aeruginosa"
  df$File1[df$File1 == "sal"] <- "S. enterica"
  df$File1[df$File1 == "se"] <- "S. equisimilis"
  df$File1[df$File1 == "sm"] <- "S. maltophilia"
  df$File1[df$File1 == "staph_aureus"] <- "S. aureus"
  df$File1[df$File1 == "GAS"] <- "S. pyogenes"
  
  df$File2[df$File2 == "ab"] <- "A. baumanii"
  df$File2[df$File2 == "bp"] <- "B. pertussis"
  df$File2[df$File2 == "cj"] <- "C. jejuni"
  df$File2[df$File2 == "GBS_full"] <- "S. agalactiae"
  df$File2[df$File2 == "GPSv4"] <- "S. pneumoniae"
  df$File2[df$File2 == "e_feacalis"] <- "E. faecalis"
  df$File2[df$File2 == "e_feacium"] <- "E. faecium"
  df$File2[df$File2 == "ec"] <- "E. coli"
  df$File2[df$File2 == "hflu"] <- "H. influenzae"
  df$File2[df$File2 == "hp"] <- "H. pylori"
  df$File2[df$File2 == "kp"] <- "K. pneumoniae"
  df$File2[df$File2 == "lm"] <- "L. monocytogenes"
  df$File2[df$File2 == "lp"] <- "L. pneumophila"
  df$File2[df$File2 == "ma"] <- "M. abscessus"
  df$File2[df$File2 == "mcat"] <- "M. catarrhalis"
  df$File2[df$File2 == "mtb"] <- "M. tuberculosis"
  df$File2[df$File2 == "ngon"] <- "N. gonorrhoeae"
  df$File2[df$File2 == "nm"] <- "N. meningitidis"
  df$File2[df$File2 == "pa"] <- "P. aeruginosa"
  df$File2[df$File2 == "sal"] <- "S. enterica"
  df$File2[df$File2 == "se"] <- "S. equisimilis"
  df$File2[df$File2 == "sm"] <- "S. maltophilia"
  df$File2[df$File2 == "staph_aureus"] <- "S. aureus"
  df$File2[df$File2 == "GAS"] <- "S. pyogenes"
}

# TODO: generate embedding from pairwise distances
# Get all unique item names
items <- unique(c(df$File1, df$File2))
n <- length(items)

# Initialize distance matrix
dist_mat <- matrix(0, nrow = n, ncol = n, dimnames = list(items, items))

# Fill in the distances
for (i in 1:nrow(df)) {
  a <- df$File1[i]
  b <- df$File2[i]
  d <- df$Distance[i]
  dist_mat[a, b] <- d
  dist_mat[b, a] <- d  # make symmetric
}
dist_obj <- as.dist(dist_mat)

#compare with Bobay & Ochman 2018 https://pmc.ncbi.nlm.nih.gov/articles/PMC6186134/
{
  # get NE data
  ne_dataframe <- read_excel("Bobay_Ochman_2018_Ne_estimates.xlsx", sheet="thesis_uni")
  
  ne_dataframe$Species[ne_dataframe$Abbreviation == "ab"] <- "A. baumanii"
  ne_dataframe$Species[ne_dataframe$Abbreviation == "bp"] <- "B. pertussis"
  ne_dataframe$Species[ne_dataframe$Abbreviation == "cj"] <- "C. jejuni"
  ne_dataframe$Species[ne_dataframe$Abbreviation == "GBS_full"] <- "S. agalactiae"
  ne_dataframe$Species[ne_dataframe$Abbreviation == "GPSv4"] <- "S. pneumoniae"
  ne_dataframe$Species[ne_dataframe$Abbreviation == "e_feacalis"] <- "E. faecalis"
  ne_dataframe$Species[ne_dataframe$Abbreviation == "e_feacium"] <- "E. faecium"
  ne_dataframe$Species[ne_dataframe$Abbreviation == "ec"] <- "E. coli"
  ne_dataframe$Species[ne_dataframe$Abbreviation == "hflu"] <- "H. influenzae"
  ne_dataframe$Species[ne_dataframe$Abbreviation == "hp"] <- "H. pylori"
  ne_dataframe$Species[ne_dataframe$Abbreviation == "kp"] <- "K. pneumoniae"
  ne_dataframe$Species[ne_dataframe$Abbreviation == "lm"] <- "L. monocytogenes"
  ne_dataframe$Species[ne_dataframe$Abbreviation == "lp"] <- "L. pneumophila"
  ne_dataframe$Species[ne_dataframe$Abbreviation == "ma"] <- "M. abscessus"
  ne_dataframe$Species[ne_dataframe$Abbreviation == "mcat"] <- "M. catarrhalis"
  ne_dataframe$Species[ne_dataframe$Abbreviation == "mtb"] <- "M. tuberculosis"
  ne_dataframe$Species[ne_dataframe$Abbreviation == "ngon"] <- "N. gonorrhoeae"
  ne_dataframe$Species[ne_dataframe$Abbreviation == "nm"] <- "N. meningitidis"
  ne_dataframe$Species[ne_dataframe$Abbreviation == "pa"] <- "P. aeruginosa"
  ne_dataframe$Species[ne_dataframe$Abbreviation == "sal"] <- "S. enterica"
  ne_dataframe$Species[ne_dataframe$Abbreviation == "se"] <- "S. equisimilis"
  ne_dataframe$Species[ne_dataframe$Abbreviation == "sm"] <- "S. maltophilia"
  ne_dataframe$Species[ne_dataframe$Abbreviation == "staph_aureus"] <- "S. aureus"
  ne_dataframe$Species[ne_dataframe$Abbreviation == "GAS"] <- "S. pyogenes"
  ne_dataframe$Ne <- as.numeric(ne_dataframe$Ne)
  ne_dataframe$rm <- as.numeric(ne_dataframe$rm)
  ne_dataframe$hm <- as.numeric(ne_dataframe$hm)
  ne_dataframe$dNdS <- as.numeric(ne_dataframe$dNdS)
  ne_dataframe$dNdS_all <- as.numeric(ne_dataframe$dNdS_all)
  ne_dataframe$Phi <- as.numeric(ne_dataframe$Phi)
  ne_dataframe$Adjusted_pangenome <- as.numeric(ne_dataframe$Adjusted_pangenome)
}


# use linear transformation
embedding <- cmdscale(dist_obj, k = 2)
embedding_df <- as.data.frame(embedding)
embedding_df$Species <- rownames(embedding_df)

embedding_df <- merge(embedding_df, ne_dataframe, by="Species")

#use PCA
pc <- prcomp(dist_obj, center=TRUE, scale. = TRUE)
pc_coords <-as.data.frame(pc$x)
pc_coords$Species <- row.names(pc_coords)
embedding_df <- merge(embedding_df, pc_coords, by="Species")

# use tSNE
tsne_result <- Rtsne(dist_obj, is_distance = TRUE, dims = 2, perplexity=5)
tsne_coords <- as.data.frame(tsne_result$Y)
colnames(tsne_coords) <- c("tSNE1", "tSNE2")
tsne_coords$Species <- embedding_df$Species  # same order as dist_obj
embedding_df <- merge(embedding_df, tsne_coords, by="Species")

# use UMAP
mds_highdim <- as.data.frame(cmdscale(dist_obj, k = min(n - 1, 10)))  # up to n-1 dimensions

# Run UMAP
umap_result <- umap(mds_highdim)
umap_coords <- as.data.frame(umap_result$layout)
colnames(umap_coords) <- c("UMAP1", "UMAP2")
umap_coords$Species <- embedding_df$Species
embedding_df <- merge(embedding_df, umap_coords, by="Species")

# colour by lifestyle
{
  # colour by lifestyle, plot linear
  p <- ggplot(embedding_df, aes(x = V1, y = V2, label = Species, color = Life_Style)) +
    geom_point(size = 3) +
    geom_text(size=2,nudge_y = 0.02) +
    theme_light() +
    labs(x = "Dimension 1", y = "Dimension 2", color = "Life Style") +
    theme(axis.text.x = element_text(size = 18), axis.text.y = element_text(size = 18), axis.title=element_text(size=16,face="bold"), strip.text.x = element_text(size = 12), legend.title = element_text(size = 16), legend.text = element_text(size = 14))
  p
  ggsave(paste("./inter_species_comparisons/pairwise_embeddings_linear_label.svg", sep=""), plot = p, width=10, height = 6)
  ggsave(paste("./inter_species_comparisons/pairwise_embeddings_linear_label.png", sep=""), plot = p, width=10, height = 6)
  
  # colour by lifestyle, plot linear
  p <- ggplot(embedding_df, aes(x = PC1, y = PC2, label = Species, color = Life_Style)) +
    geom_point(size = 3) +
    geom_text(size=2,nudge_y = 0.2) +
    theme_light() +
    labs(x = "PC1", y = "PC2", color = "Life Style") +
    theme(axis.text.x = element_text(size = 18), axis.text.y = element_text(size = 18), axis.title=element_text(size=16,face="bold"), strip.text.x = element_text(size = 12), legend.title = element_text(size = 16), legend.text = element_text(size = 14))
  p
  ggsave(paste("./inter_species_comparisons/pairwise_embeddings_PCA_label.svg", sep=""), plot = p, width=10, height = 6)
  ggsave(paste("./inter_species_comparisons/pairwise_embeddings_PCA_label.png", sep=""), plot = p, width=10, height = 6)
  
  
  # colour by lifestyle, plot UMAP
  p <- ggplot(embedding_df, aes(x = UMAP1, y = UMAP1, color = Life_Style)) +
    geom_point(size = 3) +
    theme_light() +
    labs(x = "Dimension 1", y = "Dimension 2", color = "Life Style") +
    theme(axis.text.x = element_text(size = 18), axis.text.y = element_text(size = 18), axis.title=element_text(size=16,face="bold"), strip.text.x = element_text(size = 12), legend.title = element_text(size = 16), legend.text = element_text(size = 14))
  p
  ggsave(paste("./inter_species_comparisons/pairwise_embeddings_UMAP.svg", sep=""), plot = p, width=10, height = 6)
  ggsave(paste("./inter_species_comparisons/pairwise_embeddings_UMAP.png", sep=""), plot = p, width=10, height = 6)

  # colour by lifestyle, plot tsne
  p <- ggplot(embedding_df, aes(x = tSNE1, y = tSNE2, color = Life_Style)) +
    geom_point(size = 3) +
    theme_light() +
    labs(x = "Dimension 1", y = "Dimension 2", color = "Life Style") +
    theme(axis.text.x = element_text(size = 18), axis.text.y = element_text(size = 18), axis.title=element_text(size=16,face="bold"), strip.text.x = element_text(size = 12), legend.title = element_text(size = 16), legend.text = element_text(size = 14))
  p
  ggsave(paste("./inter_species_comparisons/pairwise_embeddings_tSNE.svg", sep=""), plot = p, width=10, height = 6)
  ggsave(paste("./inter_species_comparisons/pairwise_embeddings_tSNE.png", sep=""), plot = p, width=10, height = 6)
}

# colour by species
{
  # colour by lifestyle, plot linear
  p <- ggplot(embedding_df, aes(x = V1, y = V2, color = Species)) +
    geom_point(size = 3) +
    theme_light() +
    labs(x = "Dimension 1", y = "Dimension 2", color = " Species") +
    theme(axis.text.x = element_text(size = 18), axis.text.y = element_text(size = 18), axis.title=element_text(size=16,face="bold"), strip.text.x = element_text(size = 12), legend.title = element_text(size = 16), legend.text = element_text(size = 14))
  p
  ggsave(paste("./inter_species_comparisons/pairwise_embeddings_linear_species.svg", sep=""), plot = p, width=10, height = 6)
  ggsave(paste("./inter_species_comparisons/pairwise_embeddings_linear_species.png", sep=""), plot = p, width=10, height = 6)
  
  # colour by lifestyle, plot UMAP
  p <- ggplot(embedding_df, aes(x = UMAP1, y = UMAP1, color = Specie)) +
    geom_point(size = 3) +
    theme_light() +
    labs(x = "Dimension 1", y = "Dimension 2", color = "Specie") +
    theme(axis.text.x = element_text(size = 18), axis.text.y = element_text(size = 18), axis.title=element_text(size=16,face="bold"), strip.text.x = element_text(size = 12), legend.title = element_text(size = 16), legend.text = element_text(size = 14))
  p
  ggsave(paste("./inter_species_comparisons/pairwise_embeddings_UMAP.svg", sep=""), plot = p, width=10, height = 6)
  ggsave(paste("./inter_species_comparisons/pairwise_embeddings_UMAP.png", sep=""), plot = p, width=10, height = 6)

  # colour by lifestyle, plot tsne
  p <- ggplot(embedding_df, aes(x = tSNE1, y = tSNE2, color = Specie)) +
    geom_point(size = 3) +
    theme_light() +
    labs(x = "Dimension 1", y = "Dimension 2", color = "Specie") +
    theme(axis.text.x = element_text(size = 18), axis.text.y = element_text(size = 18), axis.title=element_text(size=16,face="bold"), strip.text.x = element_text(size = 12), legend.title = element_text(size = 16), legend.text = element_text(size = 14))
  p
  ggsave(paste("./inter_species_comparisons/pairwise_embeddings_tSNE.svg", sep=""), plot = p, width=10, height = 6)
  ggsave(paste("./inter_species_comparisons/pairwise_embeddings_tSNE.png", sep=""), plot = p, width=10, height = 6)
}

