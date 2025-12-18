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
    0
  }
  return(label)
}

# v placehold stands for value
parse_filename <- function(filename) {
  list(
    prop_positive_v = search_string(filename, "prop_positive"),
    HR_rate_v = search_string(filename, "HR_rate"),
    HGT_rate_v = search_string(filename, "HGT_rate"),
    rate_genes1_v = search_string(filename, "rate_genes1"),
    rate_genes2_v = search_string(filename, "rate_genes2"),
    prop_genes2_v = search_string(filename, "prop_genes2"),
    core_mu_v = search_string(filename, "core_mu"),
    competition_strength_v = search_string(filename, "competition_strength"),
    pos_lambda_v = search_string(filename, "pos_lambda"),
    neg_lambda_v = search_string(filename, "neg_lambda")
  )
}

scatter_with_contour <- function(data, mapping, ...) {
  ggplot(data = data, mapping = mapping) +
    geom_point(alpha = 0.5, colour = "#808180FF", size = 0.5) +
    geom_density_2d(alpha = 0.7, colour = "#1B1919FF", bins=10)  # adds contour lines
}

custom_diag_hist <- function(data, mapping, ...) {
  ggplot(data = data, mapping = mapping) +
    geom_histogram(
      fill = "#808180FF",   # fill color
      color = "black", # outline color
      bins = 20,          # number of bins
      alpha = 0.7
    ) +
    theme_light()
}

plot_posterior <- function(df, ori.col.names, col.names, outpref, max.values)
{
  idx_keep <- which(!grepl("_v$", colnames(df)))
  col.names.ori <- colnames(df)[idx_keep]
  
  names_keep <- match(col.names.ori, ori.col.names)
  col.names.new <- col.names[names_keep]
  max.values.new <- max.values[names_keep]
  
  colnames(df)[idx_keep] <- col.names.new
  
  # Generate marginal plots for each pair
  idx <- 1
  no_log <- c("prop_positive", "prop_genes2", "competition_strength")
  
  plots <- lapply(idx_keep, function(idx) {
    xvar <- idx_keep[idx]
    
    x.limits <- max.values.new[[idx]]
    
    # get original column names
    x.name.ori <- col.names.ori[idx]
    
    # get new column names
    x.name <- col.names.new[xvar]
    
    # Corresponding real values
    x_real <- unique(df[[paste0(x.name.ori, "_v")]])
    
    x.median <- median(df[[x.name]])
    
    x.5quantile <- quantile(df[[x.name]], 0.025)
    x.95quantile <- quantile(df[[x.name]], 0.975)
    
    v_line_color <- if (x_real >= x.5quantile && x_real <= x.95quantile) {
      "#00A087FF"
    } else {
      "#E64B35FF"
    }
    
    p <- ggplot(df, aes(x = .data[[x.name]])) +
      geom_density(aes(y = ..scaled..), alpha = 0.5, size=0.5, colour = "#808180FF") +
      
      geom_vline(xintercept = x.median, linetype = "solid", color = "#3C5488FF") +
      
      geom_vline(xintercept = x_real, linetype = "dashed", color = v_line_color) +
      
      # apply limits
      xlim(x.limits[1], x.limits[2]) +
      
      xlab(x.name) + 
      ylab("Density") +
      theme_minimal()
    p <- p + annotate(
      "rect",
      xmin = x.5quantile, xmax = x.95quantile,   # x range of box
      ymin = 0, ymax = Inf,   # y range of box
      alpha = 0.2,            # transparency
      fill = "#4DBBD5FF",           # color
      #color = "black",
      #linetype = "solid",
      linewidth = 0.5
    ) + ylab("")
    p <- p + theme(axis.text.x = element_text(size = 14), axis.text.y = element_text(size = 14), axis.title=element_text(size=16,face="bold"), strip.text.x = element_text(size = 12), 
                  legend.text = element_text(size = 12), legend.title = element_text(size = 14, face="bold"))
    # only log if required
    to_log <- x.name.ori %in% no_log
    if (!to_log) {
      p <- p + scale_x_log10(limits=x.limits, labels = trans_format("log10", math_format(10^.x)))
    }
    p
  })
  
  plots[[2]]
  
  # merge all graphs together
  
  # combine plots
  merged_plots <- ggarrange(
    plotlist=plots,
    nrow = 1,
    align = "hv")
  
  
  merged_plots
  # Save each plot with informative filenames
  fname <- paste0(outpref, "_merged.png")
  ggsave(fname, annotate_figure(merged_plots, left = text_grob("Density", color = "black", face = "bold", rot = 90, size = 16)), 
         width = 12, height = 6)
  
  merged_plots
}

parse_results <- function(df_paths, ori.col.names, col.names, outpref, max.values)
{
  j <- 1
  for (j in 1:length(df_paths))
  {
    filename <- df_paths[j]
    df <- read.table(filename, sep = ",", comment.char = "", header=TRUE)
    
    # add metadata
    parsed_values <- parse_filename(filename)
    for (name in names(parsed_values)) {
      df[[name]] <- parsed_values[[name]]
    }
    
    base.name <- file_path_sans_ext(basename(filename))
    
    plots <- plot_posterior(df, ori.col.names, col.names, paste(outpref, base.name, sep = "/"), max.values)
    
    # create df_all if first iteration
    if (j == 1)
    {
      #df_all <- df
      plots_all <- list()
      plots_all <- append(plots_all, list(plots))
    } else {
      #df_all <- rbind(df_all, df)
      plots_all <- append(plots_all, list(plots))
    }
  }
  
  plots_all
}

# all
{
  indir <- "BOLFI_gridsearch/data/param_step/"
  #df_paths <- Sys.glob(paste0(indir, "*_mcmc_posterior_samples.csv"))
  df_paths <- c(Sys.glob(paste0(indir, "*prop_genes2*_mcmc_posterior_samples.csv")),
                Sys.glob(paste0(indir, "*HGT_rate*_mcmc_posterior_samples.csv")),
                Sys.glob(paste0(indir, "*HR_rate*_mcmc_posterior_samples.csv")),
                Sys.glob(paste0(indir, "*prop_positive*_mcmc_posterior_samples.csv")),
                Sys.glob(paste0(indir, "*pos_lambda*_mcmc_posterior_samples.csv")),
                Sys.glob(paste0(indir, "*neg_lambda*_mcmc_posterior_samples.csv")),
                Sys.glob(paste0(indir, "*competition_strength*_mcmc_posterior_samples.csv")))
  outpref <- "BOLFI_gridsearch/figures/param_step/"
  # get columns which don't have value placeholders and name if required
  ori.col.names <- c("rate_genes1", "core_mu", "prop_positive", 
                     "HR_rate", "HGT_rate", "rate_genes2", "prop_genes2", 
                     "competition_strength", "pos_lambda", "neg_lambda")
  col.names <- c("Basal gene turnover rate", "Core mutation rate", "Proportion +ve selected genes",
                 "HR rate", "HGT rate", "Fast gene turnover rate", "Proportion fast genes",
                 "Competition strength", "Gene +ve selection lambda", "Gene -ve selection lambda")
  max.values <- list(c(1e-5, 1.0), c(1e-10, 1e-2), c(0.0, 1.0), 
                     c(1e-1, 1000), c(1e-7, 10.0), c(10.0, 10.0), c(0.0, 1.0), 
                     c(0.0, 1000), c(1e-2, 100), c(1e-2, 100))
  #max.values <- list(c(NA, NA), c(NA, NA))
  plot_list <- parse_results(df_paths, ori.col.names, col.names, outpref, max.values)
  
  plot_list[1]
  # combine plots

  
  combined_plots <- ggarrange(
    plotlist=plot_list,
    ncol = 1,
    align = "hv",
    labels="AUTO")
  
  combined_plots <- annotate_figure(combined_plots,
                                         left = text_grob("Density", color = "black", face = "bold", rot = 90, size = 20))
  
  combined_plots
  ggsave(paste(outpref, "/param_step_BOLFI_fit.png", sep=""), plot = combined_plots, width=16, height = 16)
  ggsave(paste(outpref, "/param_step_BOLFI_fit.svg", sep=""), plot = combined_plots, width=16, height = 16)
  
}




