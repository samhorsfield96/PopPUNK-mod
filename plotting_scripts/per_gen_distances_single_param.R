library(ggplot2)
library(dplyr)
library(ggsci)
library(GGally)
library(tools)
library(ggExtra)
library(purrr)
library(ggpubr)
library(scales)
library(cowplot)

search_string <- function(filename, string){
  search_rg <- paste0(".*", string, "_(-?[0-9.]+(?:[eE][-+]?[0-9]+)?)(?:_|\\.|$).*")
  label = if (grepl(string, filename)) {
    as.numeric(sub(search_rg, "\\1", filename))
  } else {
    NA
  }
  return(label)
}

# v placehold stands for value
parse_filename <- function(filename) {
  list(
    rate_genes2_v = search_string(filename, "rate_genes2"),
    core_genes_v = search_string(filename, "core_genes"),
    core_nuc_v = search_string(filename, "core_nuc"),
    pan_size_v = search_string(filename, "pan_size"),
    n_gen_v = search_string(filename, "n_gen"),
    avg_gene_freq_v = search_string(filename, "avg_gene_freq"),
    pop_size_v = search_string(filename, "pop_size"),
    prop_positive_v = search_string(filename, "prop_positive"),
    HR_rate_v = search_string(filename, "HR_rate"),
    HGT_rate_v = search_string(filename, "HGT_rate"),
    rate_genes1_v = search_string(filename, "rate_genes1"),
    prop_genes2_v = search_string(filename, "prop_genes2"),
    core_mu_v = search_string(filename, "core_mu"),
    competition_strength_v = search_string(filename, "competition_strength"),
    pos_lambda_v = search_string(filename, "pos_lambda"),
    neg_lambda_v = search_string(filename, "neg_lambda")
  )
}

parse_results <- function(df_paths, ori.col.names, col.names) 
{
  j <- 1
  for (j in 1:length(df_paths))
  {
    filename <- df_paths[j]
    df <- read.table(filename, sep = "\t", comment.char = "", header=FALSE)
    colnames(df) <- c("avg_core", "avg_acc", "std_core", "std_acc")
    
    # add metadata
    parsed_values <- parse_filename(filename)
    name <- "n_gen_v"
    for (name in names(parsed_values)) {
      matched <- match(sub("_v$", "", name), ori.col.names)
      matched.real <- ori.col.names[matched]
      df[[matched.real]] <- parsed_values[[name]]
    }
    
    # remove values that aren't present
    df <- df[,colSums(is.na(df))<nrow(df)]
    
    df$generation <- seq(from = 1, to = nrow(df))
    
    base.name <- file_path_sans_ext(basename(filename))
    
    # create df_all if first iteration
    if (j == 1)
    {
      df_all <- df
    } else {
      df_all <- rbind(df_all, df)
    }
  }
  df_all
}

plot_per_gen <- function(df_all, outpref){
  # subset to ignore results for ngen =10000
  df_all <- subset(df_all, n_gen != 10000)
  
  df_all$min_core <- df_all$avg_core - df_all$std_core
  df_all$max_core <- df_all$avg_core + df_all$std_core
  df_all$min_acc <- df_all$avg_acc - df_all$std_acc
  df_all$max_acc <- df_all$avg_acc + df_all$std_acc
  
  df_all$min_core[df_all$min_core<0] <- 0
  df_all$min_acc[df_all$min_acc<0] <- 0
  
  df_all$n_gen <- factor(df_all$n_gen)
  df_all$pop_size <- factor(df_all$pop_size)
  
  p.acc <- ggplot(df_all, aes(x = generation, y = avg_acc, colour = pop_size, group = pop_size)) + geom_line() +
    #geom_ribbon(aes(ymin = min_acc, ymax = max_acc, colour = n_gen), alpha = 0.2) +
    labs(x = "", y = "Avg. Accessory Distance") +
    facet_wrap( ~ n_gen, scales = "free_x", nrow = 1) + 
    theme_light() + scale_colour_brewer(palette = "YlOrRd") +
    theme(strip.text = element_text(color = "white", size = 12),
          #strip.text.y.right = element_text(angle = 0),
          axis.title=element_text(size=14,face="bold"),
          axis.text=element_text(size=10),
          legend.title=element_text(size=12, face="bold"), 
          legend.text=element_text(size=12))  +
    labs(color = "Population size")
  p.acc
  
  p.acc.std <- ggplot(df_all) + geom_line(aes(x=generation, y=std_acc, colour = pop_size, group = pop_size)) +
    labs(x = "", y = "StdDev. Accessory Distance") +
    facet_wrap( ~ n_gen, scales = "free_x", nrow = 1) + 
    theme_light() + scale_colour_brewer(palette = "YlOrRd") +
    theme(strip.text = element_text(color = "white", size = 12),
          #strip.text.y.right = element_text(angle = 0),
          axis.title=element_text(size=14,face="bold"),
          axis.text=element_text(size=10),
          legend.title=element_text(size=12, face="bold"), 
          legend.text=element_text(size=12)) +
    labs(color = "Population size")
  p.acc.std
  
  p.core <- ggplot(df_all, aes(x = generation, y = avg_core, colour = pop_size, group = pop_size)) + geom_line() +
    #geom_ribbon(aes(ymin = min_core, ymax = max_core), fill = "#0c589c", alpha = 0.2) +
    labs(x = "", y = "Avg. Core Distance") +
    facet_wrap( ~ n_gen, scales = "free_x", nrow = 1) + 
    theme_light() + scale_colour_brewer(palette = "YlOrRd") +
    theme(strip.text = element_text(color = "white", size = 12),
          #strip.text.y.right = element_text(angle = 0),
          axis.title=element_text(size=14,face="bold"),
          axis.text=element_text(size=10),
          legend.title=element_text(size=12, face="bold"), 
          legend.text=element_text(size=12)) +
    labs(color = "Population size")
  p.core
  
  p.core.std <- ggplot(df_all) + geom_line(aes(x=generation, y=std_core, colour = pop_size, group = pop_size)) +
    labs(x = "", y = "StdDev. Core Distance") +
    facet_wrap( ~ n_gen, scales = "free_x", nrow = 1) + 
    theme_light() + scale_colour_brewer(palette = "YlOrRd") +
    theme(strip.text = element_text(color = "white", size = 12),
          #strip.text.y.right = element_text(angle = 0),
          axis.title=element_text(size=14,face="bold"),
          axis.text=element_text(size=10),
          legend.title=element_text(size=12, face="bold"), 
          legend.text=element_text(size=12)) +
    labs(color = "Population size")
  p.core.std
  
  avg_plots <- ggarrange(
    p.core, p.acc,
    ncol = 1,
    align = "hv",
    #labels="AUTO",
    legend="right",
    common.legend = TRUE)
  
  
  avg_plots <- annotate_figure(avg_plots,
                                    bottom = text_grob("Generation", color = "black", face = "bold", size = 18))
  avg_plots
  
  std_plots <- ggarrange(
    p.core.std, p.acc.std,
    ncol = 1,
    align = "hv",
    #labels="AUTO",
    legend="right",
    common.legend = TRUE)
  
  
  std_plots <- annotate_figure(std_plots,
                               bottom = text_grob("Generation", color = "black", face = "bold", size = 18))
  std_plots
  
  all_plots <- ggarrange(
    p.core, p.acc, p.core.std, p.acc.std,
    ncol = 1,
    align = "hv",
    labels="AUTO",
    legend="right",
    common.legend = TRUE)
  
  all_plots <- annotate_figure(all_plots,
                               bottom = text_grob("Generation", color = "black", face = "bold", size = 18))
  all_plots
  
  
  ggsave(paste(outpref, "avg_distance_per_gen.png", sep=""), plot = avg_plots, width=16, height = 6.5)
  ggsave(paste(outpref, "stddev_distance_per_gen.png", sep=""), plot = std_plots, width=16, height = 6.5)
  ggsave(paste(outpref, "avg_distance_per_gen.svg", sep=""), plot = avg_plots, width=16, height = 6.5)
  ggsave(paste(outpref, "stddev_distance_per_gen.svg", sep=""), plot = std_plots, width=16, height = 6.5)
  
  ggsave(paste(outpref, "all_per_gen.png", sep=""), plot = all_plots, width=16, height = 14)
  ggsave(paste(outpref, "all_per_gen.svg", sep=""), plot = all_plots, width=16, height = 14)
}

indir <- "per_gen_distances_single_param/data/"
outpref <- "per_gen_distances_single_param/figures/"
df_paths <- Sys.glob(paste0(indir, "*_per_gen.tsv"))
ori.col.names <- c("rate_genes2", "core_genes", "core_nuc", 
                   "pan_size", "n_gen", "avg_gene_freq", "pop_size",
                   "rate_genes1", "core_mu", "prop_positive", 
                   "HR_rate", "HGT_rate", "prop_genes2", 
                   "competition_strength", "pos_lambda", "neg_lambda")
col.names <- c("Fast gene turnover rate", "Core genome size (genes)", "Core genome size (nucleotides)", 
               "Pangenome size (genes)", "No. generations", "Average proportion of pangenome per genome",  "Population size",
               "Basal gene turnover rate", "Core mutation rate", "Proportion +ve selected genes",
               "HR rate", "HGT rate", "Proportion fast genes",
               "Competition strength", "Gene +ve selection lambda", "Gene -ve selection lambda")

df_all <- parse_results(df_paths = df_paths, ori.col.names = ori.col.names, col.names = col.names)

plot_per_gen(df_all, outpref)


