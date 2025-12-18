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
  in_string <- grepl(string, filename)
  label = if (in_string) {
    as.numeric(sub(search_rg, "\\1", filename))
  } else {
    NA
  }
  
  # check if lable present but no numerical value
  if (in_string & is.na(label)){
    label <- 1
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

parse_results <- function(df_paths, ori.col.names, col.names, base.values) 
{
  statistic <- c("avg_core", "avg_acc", "std_core", "std_acc")
  j <- 1
  for (j in 1:length(df_paths))
  {
    filename <- df_paths[j]
    df <- read.table(filename, sep = "\t", comment.char = "", header=FALSE)
    colnames(df) <- statistic
    
    # add metadata
    parsed_values <- parse_filename(filename)
    name <- "n_gen_v"
    for (name in names(parsed_values)) {
      matched <- match(sub("_v$", "", name), ori.col.names)
      matched.real <- ori.col.names[matched]
      matched.val <- base.values[matched]
      if (is.na(parsed_values[[name]]))
      {
        real.val <- matched.val
      } else {
        real.val <- parsed_values[[name]]
      }
      df[[matched.real]] <- real.val
    }
    
    # Identify constant columns (the others)
    constant_cols <- setdiff(names(df), statistic)
    
    # Compute means and standard deviations
    means <- colMeans(df[statistic], na.rm = TRUE)
    sds   <- sapply(df[statistic], sd, na.rm = TRUE)
    
    # Combine results into a data frame
    summary_stats <- data.frame(
      variable = statistic,
      mean = means,
      sd = sds,
      row.names = NULL
    )
    
    # Take the constants from the first row
    constants <- df[1, constant_cols, drop = FALSE]
    
    # Combine everything — each varying column has mean/sd, constants repeated once if desired
    df <- cbind(statistic, t(as.data.frame(rbind(means, sds))), constants[rep(1, length(statistic)), ])
    row.names(df) <- NULL 
    
    # remove values that aren't present
    #df <- df[,colSums(is.na(df))<nrow(df)]
    
    #df$generation <- seq(from = 1, to = nrow(df))
    
    base.name <- file_path_sans_ext(basename(filename))
    df$file <- base.name
    
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

plot_per_gen <- function(df_all, outpref, base.values, ori.col.names){
  
  # subset to ignore results for ngen =10000
  df_all <- subset(df_all, n_gen != 10000)
  
  # Base values mapping
  base.map <- setNames(base.values, ori.col.names)
  
  # Columns to ignore
  ignore_cols <- c("n_gen", "pop_size")
  
  # Use only columns that are in both df_all and base.map (minus ignored)
  compare_cols <- intersect(setdiff(names(base.map), ignore_cols), names(df_all))
  
  # Generate parameter strings per row
  i <- 5
  param_strings <- sapply(seq_len(nrow(df_all)), function(i) {
    row <- df_all[i, ]
    
    col <- "core_mu"
    diffs <- compare_cols[sapply(compare_cols, function(col) {
      val <- row[[col]]
      base <- base.map[[col]]
      
      # Safe numeric comparison with tolerance
      if (is.numeric(val) && is.numeric(base)) {
        is.na(val) || abs(val - base) > 1e-8
      } else {
        as.character(val) != as.character(base)
      }
    })]
    
    if (length(diffs) == 0) {
      "baseline"
    } else {
      paste(paste0(diffs, ": ", row[diffs]), collapse = ", ")
    }
  })
  
  # Attach results to dataframe
  df_all$param_string <- param_strings
  
  df_all$n_gen <- factor(df_all$n_gen)
  df_all$pop_size <- factor(df_all$pop_size)
  df_all$param_string <- factor(df_all$param_string)
  
  df_all$min <- df_all$means - df_all$sds
  df_all$max <- df_all$means + df_all$sds
  df_all$min[df_all$min<0] <- 0
  
  subset.df <- subset(df_all, statistic == "avg_acc")
  p.acc <- ggplot(subset.df, aes(x = param_string, y = means, colour = pop_size, group = pop_size)) + geom_point(position=position_dodge(.9)) +
    geom_errorbar(aes(ymin=min, ymax=max), width=.2,
                  position=position_dodge(.9)) +
    labs(x = "", y = "Avg. Accessory Distance") +
    facet_wrap( ~ n_gen, nrow = 1) + 
    theme_light() + scale_colour_brewer(palette = "YlOrRd") +
    theme(strip.text = element_text(color = "white", size = 12),
          #strip.text.y.right = element_text(angle = 0),
          axis.title=element_text(size=14,face="bold"),
          axis.text=element_text(size=10),
          legend.title=element_text(size=12, face="bold"), 
          legend.text=element_text(size=12)) +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
    theme(axis.title.x = element_blank(),
          axis.text.x  = element_blank()) +
    labs(color = "Population size")
  p.acc
  
  subset.df <- subset(df_all, statistic == "avg_core")
  p.core <- ggplot(subset.df, aes(x = param_string, y = means, colour = pop_size, group = pop_size)) + geom_point(position=position_dodge(.9)) +
    geom_errorbar(aes(ymin=min, ymax=max), width=.2,
                  position=position_dodge(.9)) +
    labs(x = "", y = "Avg. Core Distance") +
    facet_wrap( ~ n_gen, nrow = 1) + 
    theme_light() + scale_colour_brewer(palette = "YlOrRd") +
    theme(strip.text = element_text(color = "white", size = 12),
          #strip.text.y.right = element_text(angle = 0),
          axis.title=element_text(size=14,face="bold"),
          axis.text=element_text(size=10),
          legend.title=element_text(size=12, face="bold"), 
          legend.text=element_text(size=12)) +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
    theme(axis.title.x = element_blank(),
          axis.text.x  = element_blank()) +
    labs(color = "Population size")
  p.core
  
  subset.df <- subset(df_all, statistic == "std_core")
  p.core.std <- ggplot(subset.df, aes(x = param_string, y = means, colour = pop_size, group = pop_size)) + geom_point(position=position_dodge(.9)) +
    geom_errorbar(aes(ymin=min, ymax=max), width=.2,
                  position=position_dodge(.9)) +
    labs(x = "", y = "StdDev. Core Distance") +
    facet_wrap( ~ n_gen, nrow = 1) + 
    theme_light() + scale_colour_brewer(palette = "YlOrRd") +
    theme(strip.text = element_text(color = "white", size = 12),
          #strip.text.y.right = element_text(angle = 0),
          axis.title=element_text(size=14,face="bold"),
          axis.text=element_text(size=10),
          legend.title=element_text(size=12, face="bold"), 
          legend.text=element_text(size=12)) +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
    theme(axis.title.x = element_blank(),
          axis.text.x  = element_blank()) +
    labs(color = "Population size")
  p.core.std
  
  subset.df <- subset(df_all, statistic == "std_acc")
  p.acc.std <- ggplot(subset.df, aes(x = param_string, y = means, colour = pop_size, group = pop_size)) + geom_point(position=position_dodge(.9)) +
    geom_errorbar(aes(ymin=min, ymax=max), width=.2,
                  position=position_dodge(.9)) +
    labs(x = "", y = "StdDev. Accessory Distance") +
    facet_wrap( ~ n_gen, nrow = 1) + 
    theme_light() + scale_colour_brewer(palette = "YlOrRd") +
    theme(strip.text = element_text(color = "white", size = 12),
          #strip.text.y.right = element_text(angle = 0),
          axis.title=element_text(size=14,face="bold"),
          axis.text=element_text(size=10),
          legend.title=element_text(size=12, face="bold"), 
          legend.text=element_text(size=12)) +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
    labs(color = "Population size")
  p.acc.std
  
  # avg_plots <- ggarrange(
  #   p.core, p.acc,
  #   ncol = 1,
  #   align = "hv",
  #   #labels="AUTO",
  #   legend="right",
  #   common.legend = TRUE)
  # 
  # 
  # avg_plots <- annotate_figure(avg_plots,
  #                                   bottom = text_grob("Parameter combination", color = "black", face = "bold", size = 18))
  # avg_plots
  # 
  # std_plots <- ggarrange(
  #   p.core.std, p.acc.std,
  #   ncol = 1,
  #   align = "hv",
  #   #labels="AUTO",
  #   legend="right",
  #   common.legend = TRUE)
  # 
  # 
  # std_plots <- annotate_figure(std_plots,
  #                              bottom = text_grob("Parameter combination", color = "black", face = "bold", size = 18))
  # std_plots

  all_plots <- ggarrange(
    p.core, p.acc, p.core.std, p.acc.std,
    ncol = 1,
    align = "v",
    labels="AUTO",
    legend="right",
    heights = c(1 , 1,  1, 1.5),
    common.legend = TRUE)

  all_plots <- annotate_figure(all_plots,
                               bottom = text_grob("Parameter combination", color = "black", face = "bold", size = 18))
  all_plots

  # all_plots <- ggplot(df_all, aes(x = param_string, y = means, colour = pop_size, group = pop_size)) + geom_point(position=position_dodge(.9)) +
  #   geom_errorbar(aes(ymin=min, ymax=max), width=.2,
  #                 position=position_dodge(.9)) +
  #   labs(x = "", y = "Avg. Core Distance") +
  #   facet_grid(statistic ~ n_gen) +
  #   theme_light() + scale_colour_brewer(palette = "YlOrRd") +
  #   theme(strip.text = element_text(color = "white", size = 12),
  #         #strip.text.y.right = element_text(angle = 0),
  #         axis.title=element_text(size=14,face="bold"),
  #         axis.text=element_text(size=10),
  #         legend.title=element_text(size=12, face="bold"),
  #         legend.text=element_text(size=12)) +
  #   theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
  #   labs(color = "Population size")
  # all_plots
  
  
  #ggsave(paste(outpref, "avg_distance_per_gen.png", sep=""), plot = avg_plots, width=16, height = 6.5)
  #ggsave(paste(outpref, "stddev_distance_per_gen.png", sep=""), plot = std_plots, width=16, height = 6.5)
  #ggsave(paste(outpref, "avg_distance_per_gen.svg", sep=""), plot = avg_plots, width=16, height = 6.5)
  #ggsave(paste(outpref, "stddev_distance_per_gen.svg", sep=""), plot = std_plots, width=16, height = 6.5)
  
  ggsave(paste(outpref, "all_per_gen.png", sep=""), plot = all_plots, width=16, height = 16)
  ggsave(paste(outpref, "all_per_gen.svg", sep=""), plot = all_plots, width=16, height = 16)
}

indir <- "per_gen_distances/data/"
outpref <- "per_gen_distances/figures/"
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
base.values <- c(10.0, 1000, 1000 * 1000, 
                 5000, 100, 0.5, 1000,
                 0.001, 0.00005, -0.1,
                 0, 0, 0,
                 0, 100.0, 100.0)

df_all <- parse_results(df_paths = df_paths, ori.col.names = ori.col.names, col.names = col.names, base.values = base.values)

plot_per_gen(df_all, outpref, base.values, ori.col.names)


