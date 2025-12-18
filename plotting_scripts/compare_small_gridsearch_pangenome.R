library(ggplot2)
library(stringr)
library(dplyr)
library(readxl)
library(ggpubr)
library(ggsci)
library(RColorBrewer)

parse_filename <- function(filename) {
  list(
    prop_positive = as.numeric(sub(".*prop_positive_(-?[0-9.]+)(?:_|\\.|$).*", "\\1", filename)),
    HR_rate = as.numeric(sub(".*HR_rate_(-?[0-9.]+)(?:_|\\.|$).*", "\\1", filename)),
    HGT_rate = as.numeric(sub(".*HGT_rate_(-?[0-9.]+)(?:_|\\.|$).*", "\\1", filename)),
    rate_genes1 = as.numeric(sub(".*rate_genes1_(-?[0-9.]+)(?:_|\\.|$).*", "\\1", filename)),
    prop_genes2 = as.numeric(sub(".*prop_genes2_(-?[0-9.]+)(?:_|\\.|$).*", "\\1", filename)),
    competition = if (grepl("competition_", filename)) {
      as.numeric(sub(".*competition_(-?[0-9.]+)(?:_|\\.|$).*", "\\1", filename))
    } else {
      -0.1
    },
    pos_lambda = if (grepl("pos_lambda_", filename)) {
      as.numeric(sub(".*pos_lambda_(-?[0-9.]+)(?:_|\\.|$).*", "\\1", filename))
    } else {
      -0.1
    },
    neg_lambda = if (grepl("neg_lambda_", filename)) {
      as.numeric(sub(".*neg_lambda_(-?[0-9.]+)(?:_|\\.|$).*", "\\1", filename))
    } else {
      -0.1
    },
    baseline = ifelse(grepl("baseline", filename), 1, 0)
  )
}

plot_results <- function(df_paths, parameter, outpref, sample.num, plot_line, contours, facets, facet_bins, label, facet_align, bins)
{
  bin_breaks <- seq(0, 100, by = 100/bins)
  
  j <- 1
  for (j in 1:length(df_paths))
  {
    filename <- df_paths[j]
    df <- read.table(filename, sep = ",", comment.char = "")
    
    sum_cols <- colSums(df)
    perc_cols <- (sum_cols / nrow(df)) * 100
    perc_cols <- perc_cols[perc_cols > 0]
    
    # Use hist() to define breaks and counts
    h <- hist(perc_cols, breaks = bin_breaks, plot = FALSE)
    
    # Build a tidy data frame
    new.df <- data.frame(
      bin_start = head(h$breaks, -1),
      bin_end   = tail(h$breaks, -1),
      count     = h$counts
    )
    
    # Add a readable label column like "0–10"
    new.df$bin <- sprintf("%.1f–%.1f", new.df$bin_start, new.df$bin_end)
    
    # Order by bin_start just to be sure
    new.df <- new.df[order(new.df$bin_start), ]
    
    # add metadata
    parsed_values <- parse_filename(filename)
    for (name in names(parsed_values)) {
      new.df[[name]] <- parsed_values[[name]]
    }
    
    new.df$group <- filename
    
    if (j == 1)
    {
      df_all = new.df
    } else {
      df_all <- rbind(df_all, new.df)
    }
  }
  
  df_all$group <- as.factor(df_all[[parameter]])
  df_all$bin_mid <- (df_all$bin_start + df_all$bin_end) / 2
  
  df_split <- df_all %>% group_split(group)
  
  p <- ggplot() + 
    geom_line(data = df_all, aes(x = bin_mid, y = count, group = group, colour = group), linewidth=1.2) +
    theme_light() +
    #scale_colour_npg() +
    scale_colour_brewer(palette = "YlOrRd") + 
    labs(color = label, facet_align) +
    theme_light() + ylab("Count") + xlab("Gene Frequency (%)") +
    theme(axis.text.x = element_text(size = 14), axis.text.y = element_text(size = 14), axis.title=element_text(size=16,face="bold"), strip.text.x = element_text(size = 12), 
          legend.text = element_text(size = 12), legend.title = element_text(size = 14, face="bold"))
  p
  
  if (facets == TRUE)
  {
    if (facet_align == "v"){
      p <- p + facet_grid(rows = vars(group)) + 
        theme(
          strip.background = element_blank(),
          strip.text.x = element_blank()
        )
    } else if (facet_align == "h") {
      p <- p + facet_grid(cols = vars(group)) + 
        theme(
          strip.background = element_blank(),
          strip.text.x = element_blank()
        )
    }
  }
  
  ggsave(paste(outpref, ".svg", sep=""), plot = p, width=8, height = 4)
  ggsave(paste(outpref, ".png", sep=""), plot = p, width=8, height = 4)
  
  return(p)
}


files <- Sys.glob("pansim_runs/pangenome/*.csv")
indir <- "pansim_runs/pangenome/"
outdir <- "pansim_runs/figures_pangenome"
facet_align = "v"
bins <- 20

# prop_genes2
{
  parameter = "prop_genes2"
  label = "Proportion\nfast genes"
  plot_line = FALSE
  contours = FALSE
  facets = FALSE
  facet_bins = NA
  df_paths <- c(paste(indir, "prop_positive_-0.1_HR_rate_0.0_HGT_rate_0.0_rate_genes1_0.001_prop_genes2_0.0_baseline_pangenome.csv", sep = ""),
                paste(indir, "prop_positive_-0.1_HR_rate_0.0_HGT_rate_0.0_rate_genes1_0.001_prop_genes2_0.1_pangenome.csv", sep = ""),
                paste(indir, "prop_positive_-0.1_HR_rate_0.0_HGT_rate_0.0_rate_genes1_0.001_prop_genes2_0.5_pangenome.csv", sep = ""),
                paste(indir, "prop_positive_-0.1_HR_rate_0.0_HGT_rate_0.0_rate_genes1_0.001_prop_genes2_1.0_pangenome.csv", sep = ""))
  #paste(indir, "prop_positive_-0.1_HR_rate_0.0_HGT_rate_0.0_rate_genes1_0.001_prop_genes2_0.0_competition_1000.0_pangenome.csv", sep = ""),
  #paste(indir, "prop_positive_-0.1_HR_rate_0.0_HGT_rate_0.0_rate_genes1_0.001_prop_genes2_0.0_competition_100000.0_pangenome.csv", sep = ""))
  outpref <- paste(outdir, "comparisons_pangenome_", parameter, sep = "")
  p <- plot_results(df_paths, parameter, outpref, sample.num, plot_line, contours, facets, facet_bins, label, facet_align, bins)
  assign(parameter, p)
}

# rate_genes1
{
  parameter = "rate_genes1"
  label = "Basal gene\nturnover rate"
  plot_line = FALSE
  contours = FALSE
  facets = FALSE
  facet_bins = NA
  df_paths <- c(paste(indir, "prop_positive_-0.1_HR_rate_0.0_HGT_rate_0.0_rate_genes1_0.001_prop_genes2_0.0_baseline_pangenome.csv", sep = ""),
                paste(indir, "prop_positive_-0.1_HR_rate_0.0_HGT_rate_0.0_rate_genes1_0.01_prop_genes2_0.0_pangenome.csv", sep = ""),
                paste(indir, "prop_positive_-0.1_HR_rate_0.0_HGT_rate_0.0_rate_genes1_0.1_prop_genes2_0.0_pangenome.csv", sep = ""),
                paste(indir, "prop_positive_-0.1_HR_rate_0.0_HGT_rate_0.0_rate_genes1_1.0_prop_genes2_0.0_pangenome.csv", sep = ""))
  #paste(indir, "prop_positive_-0.1_HR_rate_0.0_HGT_rate_0.0_rate_genes1_0.001_prop_genes2_0.0_competition_1000.0_pangenome.csv", sep = ""),
  #paste(indir, "prop_positive_-0.1_HR_rate_0.0_HGT_rate_0.0_rate_genes1_0.001_prop_genes2_0.0_competition_100000.0_pangenome.csv", sep = ""))
  outpref <- paste(outdir, "comparisons_pangenome_", parameter, sep = "")
  p <- plot_results(df_paths, parameter, outpref, sample.num, plot_line, contours, facets, facet_bins, label, facet_align, bins)
  assign(parameter, p)
}

# HGT rate
{
  parameter = "HGT_rate"
  label = "HGT rate"
  plot_line = FALSE
  contours = FALSE
  facets = FALSE
  facet_bins = NA
  df_paths <- c(paste(indir, "prop_positive_-0.1_HR_rate_0.0_HGT_rate_0.0_rate_genes1_0.001_prop_genes2_0.0_baseline_pangenome.csv", sep = ""),
                paste(indir, "prop_positive_-0.1_HR_rate_0.0_HGT_rate_0.1_rate_genes1_0.001_prop_genes2_0.0_pangenome.csv", sep = ""),
                paste(indir, "prop_positive_-0.1_HR_rate_0.0_HGT_rate_1.0_rate_genes1_0.001_prop_genes2_0.0_pangenome.csv", sep = ""),
                paste(indir, "prop_positive_-0.1_HR_rate_0.0_HGT_rate_10.0_rate_genes1_0.001_prop_genes2_0.0_pangenome.csv", sep = ""))
  outpref <- paste(outdir, "comparisons_pangenome_", parameter, sep = "")
  p <- plot_results(df_paths, parameter, outpref, sample.num, plot_line, contours, facets, facet_bins, label, facet_align, bins)
  assign(parameter, p)
}

# prop_postive rate
{
  parameter = "prop_positive"
  label = "Proportion +ve\nselected genes"
  plot_line = FALSE
  contours = TRUE
  facets = TRUE
  facet_bins = 10
  df_paths <- c(paste(indir, "prop_positive_-0.1_HR_rate_0.0_HGT_rate_0.0_rate_genes1_0.001_prop_genes2_0.0_baseline_pangenome.csv", sep = ""),
                paste(indir, "prop_positive_0.0_HR_rate_0.0_HGT_rate_0.0_rate_genes1_0.001_prop_genes2_0.0_pangenome.csv", sep = ""),
                paste(indir, "prop_positive_0.1_HR_rate_0.0_HGT_rate_0.0_rate_genes1_0.001_prop_genes2_0.0_pangenome.csv", sep = ""),
                #paste(indir, "prop_positive_0.5_HR_rate_0.0_HGT_rate_0.0_rate_genes1_0.001_prop_genes2_0.0_pangenome.csv", sep = ""),
                paste(indir, "prop_positive_1.0_HR_rate_0.0_HGT_rate_0.0_rate_genes1_0.001_prop_genes2_0.0_pangenome.csv", sep = ""))
  outpref <- paste(outdir, "comparisons_pangenome_", parameter, sep = "")
  p <- plot_results(df_paths, parameter, outpref, sample.num, plot_line, contours, facets, facet_bins, label, facet_align, bins)
  assign(parameter, p)
}

# HR_rate
{
  parameter = "HR_rate"
  label = "HR rate"
  plot_line = FALSE
  contours = FALSE
  facets = FALSE
  facet_bins = NA
  df_paths <- c(paste(indir, "prop_positive_-0.1_HR_rate_0.0_HGT_rate_0.0_rate_genes1_0.001_prop_genes2_0.0_baseline_pangenome.csv", sep = ""),
                #paste(indir, "prop_positive_-0.1_HR_rate_0.1_HGT_rate_0.0_rate_genes1_0.001_prop_genes2_0.0_pangenome.csv", sep = ""),
                #paste(indir, "prop_positive_-0.1_HR_rate_1.0_HGT_rate_0.0_rate_genes1_0.001_prop_genes2_0.0_pangenome.csv", sep = ""),
                #paste(indir, "prop_positive_-0.1_HR_rate_10.0_HGT_rate_0.0_rate_genes1_0.001_prop_genes2_0.0_pangenome.csv", sep = ""),
                paste(indir, "prop_positive_-0.1_HR_rate_100.0_HGT_rate_0.0_rate_genes1_0.001_prop_genes2_0.0_pangenome.csv", sep = ""),
                paste(indir, "prop_positive_-0.1_HR_rate_1000.0_HGT_rate_0.0_rate_genes1_0.001_prop_genes2_0.0_pangenome.csv", sep = ""))
  outpref <- paste(outdir, "comparisons_pangenome_", parameter, sep = "")
  p <- plot_results(df_paths, parameter, outpref, sample.num, plot_line, contours, facets, facet_bins, label, facet_align, bins)
  assign(parameter, p)
}

# High/low positive
{
  parameter = "pos_lambda"
  label = "Gene +ve\nselection\nlambda"
  plot_line = FALSE
  contours = TRUE
  facets = TRUE
  facet_bins = 10
  df_paths <- c(paste(indir, "prop_positive_-0.1_HR_rate_0.0_HGT_rate_0.0_rate_genes1_0.001_prop_genes2_0.0_baseline_pangenome.csv", sep = ""),
                paste(indir, "prop_positive_0.5_HR_rate_0.0_HGT_rate_0.0_rate_genes1_0.001_prop_genes2_0.0_pos_lambda_0.1_neg_lambda_100.0_pangenome.csv", sep = ""),
                paste(indir, "prop_positive_0.5_HR_rate_0.0_HGT_rate_0.0_rate_genes1_0.001_prop_genes2_0.0_pos_lambda_1.0_neg_lambda_100.0_pangenome.csv", sep = ""),
                paste(indir, "prop_positive_0.5_HR_rate_0.0_HGT_rate_0.0_rate_genes1_0.001_prop_genes2_0.0_pos_lambda_10.0_neg_lambda_100.0_pangenome.csv", sep = ""))
  outpref <- paste(outdir, "comparisons_pangenome_", parameter, sep = "")
  p <- plot_results(df_paths, parameter, outpref, sample.num, plot_line, contours, facets, facet_bins, label, facet_align, bins)
  assign(parameter, p)
}

# High/low negative
{
  parameter = "neg_lambda"
  label = "Gene -ve\nselection\nlambda"
  plot_line = FALSE
  contours = TRUE
  facets = TRUE
  facet_bins = 10
  df_paths <- c(paste(indir, "prop_positive_-0.1_HR_rate_0.0_HGT_rate_0.0_rate_genes1_0.001_prop_genes2_0.0_baseline_pangenome.csv", sep = ""),
                paste(indir, "prop_positive_0.5_HR_rate_0.0_HGT_rate_0.0_rate_genes1_0.001_prop_genes2_0.0_pos_lambda_100.0_neg_lambda_0.1_pangenome.csv", sep = ""),
                paste(indir, "prop_positive_0.5_HR_rate_0.0_HGT_rate_0.0_rate_genes1_0.001_prop_genes2_0.0_pos_lambda_100.0_neg_lambda_1.0_pangenome.csv", sep = ""),
                paste(indir, "prop_positive_0.5_HR_rate_0.0_HGT_rate_0.0_rate_genes1_0.001_prop_genes2_0.0_pos_lambda_100.0_neg_lambda_10.0_pangenome.csv", sep = ""))
  outpref <- paste(outdir, "comparisons_pangenome_", parameter, sep = "")
  p <- plot_results(df_paths, parameter, outpref, sample.num, plot_line, contours, facets, facet_bins, label, facet_align, bins)
  assign(parameter, p)
}

# competition
{
  parameter = "competition"
  label = "Competition\nstrength"
  plot_line = FALSE
  contours = TRUE
  facets = TRUE
  facet_bins = 10
  facet_align = "h"
  df_paths <- c(paste(indir, "prop_positive_-0.1_HR_rate_0.0_HGT_rate_0.0_rate_genes1_0.001_prop_genes2_0.0_baseline_pangenome.csv", sep = ""),
                paste(indir, "prop_positive_-0.1_HR_rate_0.0_HGT_rate_0.0_rate_genes1_0.001_prop_genes2_0.0_competition_1.0_pangenome.csv", sep = ""),
                paste(indir, "prop_positive_-0.1_HR_rate_0.0_HGT_rate_0.0_rate_genes1_0.001_prop_genes2_0.0_competition_10.0_pangenome.csv", sep = ""))
                #paste(indir, "prop_positive_-0.1_HR_rate_0.0_HGT_rate_0.0_rate_genes1_0.001_prop_genes2_0.0_competition_100.0_pangenome.csv", sep = ""))
                #paste(indir, "prop_positive_-0.1_HR_rate_0.0_HGT_rate_0.0_rate_genes1_0.001_prop_genes2_0.0_competition_1000.0_pangenome.csv", sep = ""))
  outpref <- paste(outdir, "comparisons_pangenome_", parameter, sep = "")
  p <- plot_results(df_paths, parameter, outpref, sample.num, plot_line, contours, facets, facet_bins, label, facet_align, bins)
  assign(parameter, p)
}

# combine plots
non_selection_plots <- ggarrange(
  rate_genes1 + xlab("") + ylab(""),
  prop_genes2 + xlab("") + ylab(""),
  HGT_rate + xlab("") + ylab(""),
  HR_rate + xlab("") + ylab(""),
  align = "hv",
  labels = "AUTO")

non_selection_plots <- annotate_figure(non_selection_plots,
                bottom = text_grob("Gene Frequency (%)", color = "black", face = "bold", size = 18),
                left = text_grob("Count", color = "black", face = "bold", rot = 90, size = 18))

non_selection_plots
ggsave(paste(outdir, "non_selection_plots_pangenome.png", sep=""), plot = non_selection_plots, width=10, height = 6)
ggsave(paste(outdir, "non_selection_plots_pangenome.svg", sep=""), plot = non_selection_plots, width=10, height = 6)

# run plot_selection_landscape before running this

selection_plots_a <- ggarrange(
  prop_positive + xlab("") + ylab("") + theme(legend.position = "none"), 
  pos_lambda + xlab("") + ylab("") + theme(legend.position = "none"),
  neg_lambda + xlab("") + ylab("") + theme(legend.position = "none"),
  nrow = 3,
  ncol = 1,
  align = "hv",
  labels = "AUTO",
  font.label = list(size = 18))

selection_plots_a <- annotate_figure(selection_plots_a,
                                       bottom = text_grob("Gene Frequency (%)", color = "black", face = "bold", size = 18),
                                       left = text_grob("Count", color = "black", face = "bold", rot = 90, size = 18))

selection_plots_all <- ggarrange(
  selection_plots_a,
  selection_plots_b,
  nrow = 1,
  ncol = 2,
  align = "hv")

selection_plots_all
ggsave(paste(outdir, "selection_plots_pangenome.png", sep=""), plot = selection_plots_all, width=14, height = 18)
ggsave(paste(outdir, "selection_plots_pangenome.svg", sep=""), plot = selection_plots_all, width=14, height = 18)

