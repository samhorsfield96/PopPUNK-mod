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

parse_results <- function(df_paths, ori.col.names, col.names, base.values, outpref) 
{
  statistic <- c("Core", "Accessory")
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
    #constants <- df[1, constant_cols, drop = FALSE]
    
    # Combine everything — each varying column has mean/sd, constants repeated once if desired
    #df <- cbind(statistic, t(as.data.frame(rbind(means, sds))), constants[rep(1, length(statistic)), ])
    #row.names(df) <- NULL 
    
    # remove values that aren't present
    #df <- df[,colSums(is.na(df))<nrow(df)]
    
    #df$generation <- seq(from = 1, to = nrow(df))
    
    base.name <- file_path_sans_ext(basename(filename))
    df$file <- base.name
    
    p <- ggplot() + 
      geom_point(data = df, aes(x = Core, y = Accessory), colour="#4DBBD5FF", alpha=0.05) +
      theme_light() +
      #scale_colour_npg() +
      scale_colour_brewer(palette = "YlOrRd") +
      theme_light() + ylab("Accessory Distance") + xlab("Core Distance") +
      theme(axis.text.x = element_text(size = 14), axis.text.y = element_text(size = 14), axis.title=element_text(size=16,face="bold"), strip.text.x = element_text(size = 12), 
            legend.text = element_text(size = 12), legend.title = element_text(size = 14, face="bold"))
    
    p <- p + geom_density_2d(
      data = df,
      aes(x = Core, y = Accessory),
      bins = 20,
      color = "black",
      alpha = 0.5
    )
    p
    
    ggsave(paste0(outpref, base.name, ".png"), plot = p, width=10, height = 6)
    
    # # create df_all if first iteration
    # if (j == 1)
    # {
    #   df_all <- df
    # } else {
    #   df_all <- rbind(df_all, df)
    # }
  }
}

indir <- "per_gen_distances/data_distances/"
outpref <- "per_gen_distances/figures_distances/"
df_paths <- Sys.glob(paste0(indir, "*.tsv"))
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

parse_results(df_paths = df_paths, ori.col.names = ori.col.names, col.names = col.names, base.values = base.values, outpref = outpref)



