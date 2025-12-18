library(tidyverse)
library(ggplot2)
library(ggsci)
library(ggpmisc)
library(ggrepel)

overall_indir <- ""
indir <- paste0(overall_indir, "pangenome_size_comparison/")
pangenome.file <- paste0(indir, "data/", "ppmod_real_input.txt")
fits.file <- paste0(overall_indir, "real_genome_fits/figures/exp_core_mu_prop_genes2_rate_genes1_merged.csv")
outpref <- paste0(indir, "figures/")

pangenome.df <- read.csv(pangenome.file, sep="\t")
fits.df <- read.csv(fits.file)

overall.df <- merge(fits.df, pangenome.df, by.x = "Species", by.y = "name")

overall.df$core_prop <- overall.df$core_genes / overall.df$pan_genes
overall.df$inter_prop <- overall.df$intermediate_genes / overall.df$pan_genes
overall.df$rare_prop <- overall.df$rare_genes / overall.df$pan_genes

subset_df_init <- subset(overall.df, Fitted_param_short == "prop_genes2")

subset_df_whole <- subset_df_init %>%
  pivot_longer(
    cols = -c(Species, species_short, Experiment, Fitted_param_short, Fitted_param_long, Abbreviation, Life_Style, combo, Species_ordered, dist_file, Median, UCI, LCI),           # all category columns
    names_to = "category",     # new column for category names
    values_to = "value"        # new column for values
  )

subset_df_whole$category_long <- NA
subset_df_whole$category_long[subset_df_whole$category == "core_prop"] <- "No. core genes /\n Total pangenome size"
subset_df_whole$category_long[subset_df_whole$category == "inter_prop"] <- "No. intermediate genes /\n Total pangenome size"
subset_df_whole$category_long[subset_df_whole$category == "rare_prop"] <- "No. rare genes /\n Total pangenome size"
subset_df_whole$category_long[subset_df_whole$category == "core_genes"] <- "No. core genes"
subset_df_whole$category_long[subset_df_whole$category == "intermediate_genes"] <- "No. intermediate genes"
subset_df_whole$category_long[subset_df_whole$category == "rare_genes"] <- "No. rare genes"
subset_df_whole$category_long[subset_df_whole$category == "pan_genes"] <- "Total pangenome size"

subset_df_whole$category_long[subset_df_whole$category == "rm"] <- "r/m"
subset_df_whole$category_long[subset_df_whole$category == "dNdS_all"] <- "dNdS (all genes)"
subset_df_whole$category_long[subset_df_whole$category == "Ne"] <- "Effective population size (Ne)"
subset_df_whole$category_long[subset_df_whole$category == "Phi"] <- "Genome fluidity"

#subset_df <- subset(subset_df_whole, category == "core_genes" | category == "intermediate_genes" | category == "rare_genes")
subset_df <- subset(subset_df_whole, category == "core_prop" | category == "inter_prop" | category == "rare_prop")

# plot
{
  # 1. Define your desired Life_Style order
  subset_df$Life_Style <- factor(
    subset_df$Life_Style,
    levels = c("Animal pathogen", "Commensal", "Free living")
  )
  
  # 2. Create the combined ordering key (Life_Style first, then species)
  subset_df <- subset_df %>%
    mutate(
      combo = paste(Life_Style, species_short, sep = "_")
    ) %>%
    arrange(Life_Style, species_short) %>%  # <-- this does the ordering logic
    mutate(
      Species_ordered = factor(combo, levels = unique(combo))
    )
  
  # 3. Apply that ordering to species_short (so you can plot species_short directly)
  species_levels <- subset_df %>%
    arrange(Life_Style, species_short) %>%
    distinct(species_short, .keep_all = TRUE) %>%
    pull(species_short)
  
  subset_df$species_short <- factor(subset_df$species_short, levels = species_levels)
  
  
  fitted_params <- unique(subset_df$Fitted_param_long)
  subset_df$Species <- factor(subset_df$Species)
  subset_df$value <- as.numeric(subset_df$value)
  
  p <- ggplot(data = subset_df, aes(x = species_short, y=value, colour = Life_Style)) + geom_point(size = 3) +
    scale_y_log10() +
    facet_wrap( ~ category, ncol=1, scales = "free_y", strip.position = "left") +
    scale_color_npg() +
    theme_light() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, face = "italic")) +
    theme(axis.text.x = element_text(size = 14), axis.text.y = element_text(size = 14), axis.title=element_text(size=16,face="bold"), strip.text.y = element_text(size = 16, colour = "black", face = "bold"), 
          legend.text = element_text(size = 12), legend.title = element_text(size = 14, face="bold"),
          strip.placement = "outside", strip.background = element_blank()) +
    xlab("Species") + ylab("") + labs(color = "Life style")
  p
  
  fname <- paste0(outpref, "pangenome_stats")
  ggsave(paste0(fname, ".png"), p, width = 9, height = 10)
  ggsave(paste0(fname, ".svg"), p, width = 9, height = 10)
}

subset_df <- subset(subset_df_whole, category == "core_prop" | category == "inter_prop" | category == "rare_prop" | category == "pan_genes")
subset_df$category_long <- factor(subset_df$category_long, levels = c("No. core genes /\n Total pangenome size", "No. intermediate genes /\n Total pangenome size", "No. rare genes /\n Total pangenome size", "Total pangenome size"))

{
  p <- ggplot(data = subset_df, aes(x = Median, y=value, group = category, colour = Life_Style, label = species_short)) + geom_point(size = 3) +
    geom_smooth(method='lm', formula= y~x, se=FALSE, colour = "black") +
    stat_poly_eq(
      formula = y ~ x,
      aes(label = paste(..rr.label.., ..p.value.label.., sep = "~~~")),
      parse = TRUE
    ) + 
    scale_y_log10() +
    facet_wrap( ~ category_long, ncol=1, scales = "free_y", strip.position = "left") +
    scale_color_npg() +
    geom_text_repel(size=2) +
    theme_light() +
    theme(axis.text.x = element_text(size = 14), axis.text.y = element_text(size = 14), axis.title=element_text(size=16,face="bold"), strip.text.y = element_text(size = 14, colour = "black", face = "bold"), 
          legend.text = element_text(size = 12), legend.title = element_text(size = 14, face="bold"),
          strip.placement = "outside", strip.background = element_blank()) +
    xlab("Proportion of fast genes (Median)") + ylab("") + labs(color = "Life style")
  p
  fname <- paste0(outpref, "median_v_pangenome_stats")
  ggsave(paste0(fname, ".png"), p, width = 9, height = 14)
  ggsave(paste0(fname, ".svg"), p, width = 9, height = 14)
}

subset_df <- subset(subset_df_whole, category == "rm" | category == "dNdS_all" | category == "Ne" | category == "Phi")

{
  p <- ggplot(data = subset_df, aes(x = Median, y=value, group = category, colour = Life_Style, label = species_short)) + geom_point(size = 3) +
    geom_smooth(method='lm', formula= y~x, se=FALSE, colour = "black") +
    stat_poly_eq(
      formula = y ~ x,
      aes(label = paste(..rr.label.., ..p.value.label.., sep = "~~~")),
      parse = TRUE
    ) + 
    #scale_y_log10() +
    facet_wrap( ~ category_long, ncol=1, scales = "free_y", strip.position = "left") +
    scale_color_npg() +
    geom_text_repel(size=2) +
    theme_light() +
    theme(axis.text.x = element_text(size = 14), axis.text.y = element_text(size = 14), axis.title=element_text(size=16,face="bold"), strip.text.y = element_text(size = 14, colour = "black", face = "bold"), 
          legend.text = element_text(size = 12), legend.title = element_text(size = 14, face="bold"),
          strip.placement = "outside", strip.background = element_blank()) +
    xlab("Proportion of fast genes (Median)") + ylab("") + labs(color = "Life style")
  p
  fname <- paste0(outpref, "median_v_population_stats")
  ggsave(paste0(fname, ".png"), p, width = 9, height = 14)
  ggsave(paste0(fname, ".svg"), p, width = 9, height = 14)
}

