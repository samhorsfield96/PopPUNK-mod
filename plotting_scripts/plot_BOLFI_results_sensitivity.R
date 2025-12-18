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

get_summary <- function(df, ori.col.names, col.names, experiment_idx)
{
  idx_keep <- which(!grepl("_v$", colnames(df)))
  idx_absent <- which(grepl("_v$", colnames(df)))
  col.names.ori <- colnames(df)[idx_keep]
  col.names.absent <- colnames(df)[idx_absent]
  
  # strip "_v" and compare
  unmatched <- col.names.absent[!sub("_v$", "", col.names.absent) %in% col.names.ori]
  unmatched.real <- sub("_v$", "", unmatched)
  names_unmatched = match(unmatched.real, ori.col.names)
  unmatched.name <- col.names[names_unmatched]
  
  names_keep <- match(col.names.ori, ori.col.names)
  col.names.new <- col.names[names_keep]
  
  # get unmatched value
  unmatched.value <- max(df[[unmatched]])
  
  # Generate marginal plots for each pair
  idx <- 1
  
  # here, need to generate a row table of one row, with the unmatched paramter, its value, and all subsequent parameters changed with their true value, median and 95% credible
  # interval, as well as whether the true parameter lies within
  sum_df <- sapply(idx_keep, function(idx) {
    xvar <- idx_keep[idx]
    
    # get original column names
    x.name.ori <- col.names.ori[idx]
    
    # get new column names
    x.name <- col.names.new[xvar]
    
    # Corresponding real values
    x_real <- unique(df[[paste0(x.name.ori, "_v")]])
    
    x.median <- median(df[[x.name.ori]])
    
    x.5quantile <- quantile(df[[x.name.ori]], 0.025)
    x.95quantile <- quantile(df[[x.name.ori]], 0.975)
    
    within.range <- if (x_real >= x.5quantile && x_real <= x.95quantile) {
      "True"
    } else {
      "False"
    }
    
    range.size <- (x.95quantile - x.5quantile) / abs(x.median)
    
    c(experiment_idx, unmatched.name, unmatched.value, x.name, x_real, x.5quantile, x.median, x.95quantile, within.range, range.size)
  })
  
  sum_df <- as.data.frame(t(sum_df))
  
  sum_df
}

parse_results <- function(df_paths, ori.col.names, col.names, outpref)
{
  j <- 18
  for (j in 1:length(df_paths))
  {
    filename <- df_paths[j]
    df <- read.table(filename, sep = ",", comment.char = "", header=TRUE)
    
    # add metadata
    parsed_values <- parse_filename(filename)
    for (name in names(parsed_values)) {
      df[[name]] <- parsed_values[[name]]
    }
    
    # remove values that aren't present
    df <- df[,colSums(is.na(df))<nrow(df)]
    
    base.name <- file_path_sans_ext(basename(filename))
    
    df <- get_summary(df, ori.col.names, col.names, j)
    
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

#TODO: plot results here, showing bars for each condition where value is within or outside CIs
plot_sensitivity <- function(summary_df, outpref)
{
  # create 2D grid of parameters and values
  dependent.params <- unique(summary_df$Affected_param)
  independent.params <- unique(summary_df$Changed_param)
  
  counts <- summary_df %>%
    group_by(Affected_param, Changed_param) %>%
    summarise(TrueCount = sum(Affected_in_range == "True"),
              TotalCount = n(),
              .groups = "drop")

  counts$PercentTrue <- (counts$TrueCount / counts$TotalCount) * 100
  
  Affected_param_order <- c("Basal gene turnover rate", "Core mutation rate", 
                            "Proportion fast genes", "HGT rate", "HR rate",
                            "Proportion +ve selected genes", "Gene +ve selection lambda", "Gene -ve selection lambda",
                            "Competition strength")
  counts$Affected_param <- factor(counts$Affected_param, levels = rev(Affected_param_order))
  
  p <- ggplot(counts, aes(x = Changed_param, y = Affected_param, fill = PercentTrue)) +
    geom_tile(color = "white") +
    scale_fill_gradient(low = "white", high = "red") +
    theme_light() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
    labs(x = "Fixed Parameter",
         y = "Variable",
         fill = "% True variable\nvalue in CI")
  p
  ggsave(paste0(outpref, "/sensitivity_heatmap.png"), plot = p, width=6.5, height = 6.5)
  ggsave(paste0(outpref, "/sensitivity_heatmap.svg"), plot = p, width=6.5, height = 6.5)
  
  summary_df$Changed_value <- factor(as.numeric(summary_df$Changed_value))
  p <- ggplot(summary_df, aes(x= Changed_value, fill=Affected_in_range)) + geom_bar(position = "fill") +
    scale_y_continuous(labels = scales::percent) +
    facet_wrap(. ~ Changed_param, scales = "free_x") +
    theme_minimal() +
    labs(x = "Variable value", y = "Percentage", fill = "True value in\nrange") +
    theme(axis.text.x = element_text(size = 12), axis.text.y = element_text(size = 12), axis.title=element_text(size=16,face="bold"), strip.text.x = element_text(size = 12), legend.title = element_text(size = 16), legend.text = element_text(size = 14))
  p
  ggsave(paste0(outpref, "/sensitivity_plot.png"), plot = p, width=14, height = 6.5)
}

# all
{
  indir <- "BOLFI_gridsearch/data/sensitivity_analysis/"
  df_paths <- Sys.glob(paste0(indir, "*_mcmc_posterior_samples.csv"))
  # df_paths <- c(Sys.glob(paste0(indir, "*prop_genes2*_mcmc_posterior_samples.csv")),
  #               Sys.glob(paste0(indir, "*HGT_rate*_mcmc_posterior_samples.csv")),
  #               Sys.glob(paste0(indir, "*HR_rate*_mcmc_posterior_samples.csv")),
  #               Sys.glob(paste0(indir, "*prop_positive*_mcmc_posterior_samples.csv")),
  #               Sys.glob(paste0(indir, "*pos_lambda*_mcmc_posterior_samples.csv")),
  #               Sys.glob(paste0(indir, "*neg_lambda*_mcmc_posterior_samples.csv")),
  #               Sys.glob(paste0(indir, "*competition_strength*_mcmc_posterior_samples.csv")))
  outpref <- "BOLFI_gridsearch/figures/sensitivity_analysis/"
  # get columns which don't have value placeholders and name if required
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
  #max.values <- list(c(NA, NA), c(NA, NA))
  summary_df <- parse_results(df_paths, ori.col.names, col.names, outpref)
  colnames(summary_df) <- c("Experiment", "Changed_param", "Changed_value", "Affected_param", "Affected_value", "LCI", "Median", "UCI", "Affected_in_range", "CI_size")
  
  plot_sensitivity(summary_df, outpref)
  write.csv(summary_df, file=paste0(outpref, "sensitivity_analysis.csv"))
}




