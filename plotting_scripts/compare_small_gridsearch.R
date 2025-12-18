library(ggplot2)
library(stringr)
library(dplyr)
library(readxl)
library(ggpubr)
library(ggsci)
library(reticulate)
library(RColorBrewer)

path_to_python="/Users/shorsfield/miniforge3/envs/biopython/bin/python"
use_python(path_to_python)
py_run_string("
import numpy as np
from scipy.optimize import curve_fit

# fit asymptotic curve using exponential decay
def negative_exponential(x, b0, b1, b2): # based on https://isem-cueb-ztian.github.io/Intro-Econometrics-2017/handouts/lecture_notes/lecture_10/lecture_10.pdf and https://www.statforbiology.com/articles/usefulequations/
    return b0 - (b0 - b1) * np.exp(-b2 * x)
  
def fit_curve(x, y):
  popt, pcov = curve_fit(negative_exponential, x, y, p0=[1.0, 0.0, 1.0], bounds=([0.0, 0.0, 0.0], [1.0, 1.0, np.inf]))
  b0_err, b1_err, b2_err = np.sqrt(np.diag(pcov))
  return popt.tolist(), [b0_err, b1_err, b2_err]
")

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

plot_results <- function(df_paths, parameter, outpref, sample.num, plot_line, contours, facets, facet_bins, label, facet_align)
{
  df_all <- data.frame(Core = c(), Accessory = c(), prop_positive = c(), HR_rate = c(), HGT_rate = c(), rate_genes1 = c(), prop_genes2 = c(), competition = c(), baseline = c())
  params_all <- data.frame(b0 = c(), b1 = c(), b2 = c(), b0_err = c(), b1_err = c(), b2_err = c(), prop_positive = c(), HR_rate = c(), HGT_rate = c(), rate_genes1 = c(), prop_genes2 = c(), competition = c(), baseline = c())
  
  for (j in 1:length(df_paths))
  {
    filename <- df_paths[j]
    df <- read.table(filename, sep = "\t", comment.char = "")
    colnames(df) <- c("Core", "Accessory")
    
    sample.num.temp <- sample.num
    
    if (nrow(df) < sample.num.temp)
    {
      sample.num.temp <- nrow(df)
    }
    
    df = sample_n(df, sample.num.temp)
    
    params <- py$fit_curve(df$Core, df$Accessory)
    b0 <- params[[1]][[1]]
    b1 <- params[[1]][[2]]
    b2 <- params[[1]][[3]]
    b0_err <- params[[2]][[1]]
    b1_err <- params[[2]][[2]]
    b2_err <- params[[2]][[3]]
    
    params_df <- data.frame(b0, b1, b2, b0_err, b1_err, b2_err)
    
    # add metadata
    parsed_values <- parse_filename(filename)
    for (name in names(parsed_values)) {
      df[[name]] <- parsed_values[[name]]
      params_df[[name]] <- parsed_values[[name]]
    }
    
    df_all <- rbind(df_all, df)
    params_all <- rbind(params_all, params_df)
  }
  
  df_all$group <- as.factor(df_all[[parameter]])
  params_all$group<- as.factor(params_all[[parameter]])
  x_values <- seq(0, max(df_all$Core), length.out = 200)
  
  lines_df <- params_all %>% 
    rowwise() %>% 
    do({
      data.frame(group = .$group,
                 x = x_values,
                 y = .$b0 - (.$b0 - .$b1) * exp(-.$b2 * x_values))
    }) %>%
    bind_rows()
  
  df_split <- df_all %>% group_split(group)
  
  p <- ggplot() + 
    geom_point(data = df_all, aes(x = Core, y = Accessory, group = group, colour = group)) +
    theme_light() +
    #scale_colour_npg() +
    scale_colour_brewer(palette = "YlOrRd") + 
    labs(color = label, facet_align) +
    theme_light() + ylab("Accessory Distance") + xlab("Core Distance") +
    theme(axis.text.x = element_text(size = 14), axis.text.y = element_text(size = 14), axis.title=element_text(size=16,face="bold"), strip.text.x = element_text(size = 12), 
          legend.text = element_text(size = 12), legend.title = element_text(size = 14, face="bold"))
  
  if (plot_line == TRUE) {
    p <- p + geom_line(data = lines_df, aes(x = x, y = y, colour = group, group = group), linewidth=1.5)
  }
  
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
  
  if (contours == TRUE)
  {
    #p <- p + geom_density_2d(data = df_all, aes(x = Core, y = Accessory, group = group), colour = "black", bins = facet_bins, alpha = 0.5)
    #p <- p + geom_density_2d_filled(data = df_all, aes(x = Core, y = Accessory, fill = group), fill = "black", bins = 500)
    # Add one density layer per group
    for (df_group in df_split) {
      p <- p + geom_density_2d(
        data = df_group,
        aes(x = Core, y = Accessory),
        bins = facet_bins,
        color = "black",
        alpha = 0.5
      )
    }
    max_x = max(df_all$Core)
    max_y = max(df_all$Accessory)
    p <- p + scale_x_continuous(limits = c(0, max_x * 1.1)) + scale_y_continuous(limits = c(0, max_y * 1.1))
  }
  
  ggsave(paste(outpref, ".svg", sep=""), plot = p, width=10, height = 6)
  ggsave(paste(outpref, ".png", sep=""), plot = p, width=10, height = 6)
  
  return(p)
}


files <- Sys.glob("pansim_runs/core_vs_acc/*.tsv")
sample.num <- 30000
indir <- "pansim_runs/core_vs_acc/"
outdir <- "pansim_runs/figures_distances"
facet_align = "v"

# prop_genes2
{
  parameter = "prop_genes2"
  label = "Proportion\nfast genes"
  plot_line = FALSE
  contours = FALSE
  facets = FALSE
  facet_bins = NA
  df_paths <- c(paste(indir, "prop_positive_-0.1_HR_rate_0.0_HGT_rate_0.0_rate_genes1_0.001_prop_genes2_0.0_baseline.tsv", sep = ""),
                paste(indir, "prop_positive_-0.1_HR_rate_0.0_HGT_rate_0.0_rate_genes1_0.001_prop_genes2_0.1.tsv", sep = ""),
                paste(indir, "prop_positive_-0.1_HR_rate_0.0_HGT_rate_0.0_rate_genes1_0.001_prop_genes2_0.5.tsv", sep = ""),
                paste(indir, "prop_positive_-0.1_HR_rate_0.0_HGT_rate_0.0_rate_genes1_0.001_prop_genes2_1.0.tsv", sep = ""))
  #paste(indir, "prop_positive_-0.1_HR_rate_0.0_HGT_rate_0.0_rate_genes1_0.001_prop_genes2_0.0_competition_1000.0.tsv", sep = ""),
  #paste(indir, "prop_positive_-0.1_HR_rate_0.0_HGT_rate_0.0_rate_genes1_0.001_prop_genes2_0.0_competition_100000.0.tsv", sep = ""))
  outpref <- paste(outdir, "comparisons_", parameter, sep = "")
  p <- plot_results(df_paths, parameter, outpref, sample.num, plot_line, contours, facets, facet_bins, label, facet_align)
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
  df_paths <- c(paste(indir, "prop_positive_-0.1_HR_rate_0.0_HGT_rate_0.0_rate_genes1_0.001_prop_genes2_0.0_baseline.tsv", sep = ""),
                paste(indir, "prop_positive_-0.1_HR_rate_0.0_HGT_rate_0.0_rate_genes1_0.01_prop_genes2_0.0.tsv", sep = ""),
                paste(indir, "prop_positive_-0.1_HR_rate_0.0_HGT_rate_0.0_rate_genes1_0.1_prop_genes2_0.0.tsv", sep = ""),
                paste(indir, "prop_positive_-0.1_HR_rate_0.0_HGT_rate_0.0_rate_genes1_1.0_prop_genes2_0.0.tsv", sep = ""))
  #paste(indir, "prop_positive_-0.1_HR_rate_0.0_HGT_rate_0.0_rate_genes1_0.001_prop_genes2_0.0_competition_1000.0.tsv", sep = ""),
  #paste(indir, "prop_positive_-0.1_HR_rate_0.0_HGT_rate_0.0_rate_genes1_0.001_prop_genes2_0.0_competition_100000.0.tsv", sep = ""))
  outpref <- paste(outdir, "comparisons_", parameter, sep = "")
  p <- plot_results(df_paths, parameter, outpref, sample.num, plot_line, contours, facets, facet_bins, label, facet_align)
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
  df_paths <- c(paste(indir, "prop_positive_-0.1_HR_rate_0.0_HGT_rate_0.0_rate_genes1_0.001_prop_genes2_0.0_baseline.tsv", sep = ""),
                paste(indir, "prop_positive_-0.1_HR_rate_0.0_HGT_rate_0.1_rate_genes1_0.001_prop_genes2_0.0.tsv", sep = ""),
                paste(indir, "prop_positive_-0.1_HR_rate_0.0_HGT_rate_1.0_rate_genes1_0.001_prop_genes2_0.0.tsv", sep = ""),
                paste(indir, "prop_positive_-0.1_HR_rate_0.0_HGT_rate_10.0_rate_genes1_0.001_prop_genes2_0.0.tsv", sep = ""))
  outpref <- paste(outdir, "comparisons_", parameter, sep = "")
  p <- plot_results(df_paths, parameter, outpref, sample.num, plot_line, contours, facets, facet_bins, label, facet_align)
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
  df_paths <- c(paste(indir, "prop_positive_-0.1_HR_rate_0.0_HGT_rate_0.0_rate_genes1_0.001_prop_genes2_0.0_baseline.tsv", sep = ""),
                paste(indir, "prop_positive_0.0_HR_rate_0.0_HGT_rate_0.0_rate_genes1_0.001_prop_genes2_0.0.tsv", sep = ""),
                paste(indir, "prop_positive_0.1_HR_rate_0.0_HGT_rate_0.0_rate_genes1_0.001_prop_genes2_0.0.tsv", sep = ""),
                #paste(indir, "prop_positive_0.5_HR_rate_0.0_HGT_rate_0.0_rate_genes1_0.001_prop_genes2_0.0.tsv", sep = ""),
                paste(indir, "prop_positive_1.0_HR_rate_0.0_HGT_rate_0.0_rate_genes1_0.001_prop_genes2_0.0.tsv", sep = ""))
  outpref <- paste(outdir, "comparisons_", parameter, sep = "")
  p <- plot_results(df_paths, parameter, outpref, sample.num, plot_line, contours, facets, facet_bins, label, facet_align)
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
  df_paths <- c(paste(indir, "prop_positive_-0.1_HR_rate_0.0_HGT_rate_0.0_rate_genes1_0.001_prop_genes2_0.0_baseline.tsv", sep = ""),
                #paste(indir, "prop_positive_-0.1_HR_rate_0.1_HGT_rate_0.0_rate_genes1_0.001_prop_genes2_0.0.tsv", sep = ""),
                #paste(indir, "prop_positive_-0.1_HR_rate_1.0_HGT_rate_0.0_rate_genes1_0.001_prop_genes2_0.0.tsv", sep = ""),
                #paste(indir, "prop_positive_-0.1_HR_rate_10.0_HGT_rate_0.0_rate_genes1_0.001_prop_genes2_0.0.tsv", sep = ""),
                paste(indir, "prop_positive_-0.1_HR_rate_100.0_HGT_rate_0.0_rate_genes1_0.001_prop_genes2_0.0.tsv", sep = ""),
                paste(indir, "prop_positive_-0.1_HR_rate_1000.0_HGT_rate_0.0_rate_genes1_0.001_prop_genes2_0.0.tsv", sep = ""))
  outpref <- paste(outdir, "comparisons_", parameter, sep = "")
  p <- plot_results(df_paths, parameter, outpref, sample.num, plot_line, contours, facets, facet_bins, label, facet_align)
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
  df_paths <- c(paste(indir, "prop_positive_-0.1_HR_rate_0.0_HGT_rate_0.0_rate_genes1_0.001_prop_genes2_0.0_baseline.tsv", sep = ""),
                paste(indir, "prop_positive_0.5_HR_rate_0.0_HGT_rate_0.0_rate_genes1_0.001_prop_genes2_0.0_pos_lambda_0.1_neg_lambda_100.0.tsv", sep = ""),
                paste(indir, "prop_positive_0.5_HR_rate_0.0_HGT_rate_0.0_rate_genes1_0.001_prop_genes2_0.0_pos_lambda_1.0_neg_lambda_100.0.tsv", sep = ""),
                paste(indir, "prop_positive_0.5_HR_rate_0.0_HGT_rate_0.0_rate_genes1_0.001_prop_genes2_0.0_pos_lambda_10.0_neg_lambda_100.0.tsv", sep = ""))
  outpref <- paste(outdir, "comparisons_", parameter, sep = "")
  p <- plot_results(df_paths, parameter, outpref, sample.num, plot_line, contours, facets, facet_bins, label, facet_align)
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
  df_paths <- c(paste(indir, "prop_positive_-0.1_HR_rate_0.0_HGT_rate_0.0_rate_genes1_0.001_prop_genes2_0.0_baseline.tsv", sep = ""),
                paste(indir, "prop_positive_0.5_HR_rate_0.0_HGT_rate_0.0_rate_genes1_0.001_prop_genes2_0.0_pos_lambda_100.0_neg_lambda_0.1.tsv", sep = ""),
                paste(indir, "prop_positive_0.5_HR_rate_0.0_HGT_rate_0.0_rate_genes1_0.001_prop_genes2_0.0_pos_lambda_100.0_neg_lambda_1.0.tsv", sep = ""),
                paste(indir, "prop_positive_0.5_HR_rate_0.0_HGT_rate_0.0_rate_genes1_0.001_prop_genes2_0.0_pos_lambda_100.0_neg_lambda_10.0.tsv", sep = ""))
  outpref <- paste(outdir, "comparisons_", parameter, sep = "")
  p <- plot_results(df_paths, parameter, outpref, sample.num, plot_line, contours, facets, facet_bins, label, facet_align)
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
  df_paths <- c(paste(indir, "prop_positive_-0.1_HR_rate_0.0_HGT_rate_0.0_rate_genes1_0.001_prop_genes2_0.0_baseline.tsv", sep = ""),
                paste(indir, "prop_positive_-0.1_HR_rate_0.0_HGT_rate_0.0_rate_genes1_0.001_prop_genes2_0.0_competition_1.0.tsv", sep = ""),
                paste(indir, "prop_positive_-0.1_HR_rate_0.0_HGT_rate_0.0_rate_genes1_0.001_prop_genes2_0.0_competition_10.0.tsv", sep = ""))
                #paste(indir, "prop_positive_-0.1_HR_rate_0.0_HGT_rate_0.0_rate_genes1_0.001_prop_genes2_0.0_competition_100.0.tsv", sep = ""))
                #paste(indir, "prop_positive_-0.1_HR_rate_0.0_HGT_rate_0.0_rate_genes1_0.001_prop_genes2_0.0_competition_1000.0.tsv", sep = ""))
  outpref <- paste(outdir, "comparisons_", parameter, sep = "")
  p <- plot_results(df_paths, parameter, outpref, sample.num, plot_line, contours, facets, facet_bins, label, facet_align)
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
                bottom = text_grob("Core distance", color = "black", face = "bold", size = 18),
                left = text_grob("Accessory distance", color = "black", face = "bold", rot = 90, size = 18))

non_selection_plots
ggsave(paste(outdir, "non_selection_plots.png", sep=""), plot = non_selection_plots, width=10, height = 6)
ggsave(paste(outdir, "non_selection_plots.svg", sep=""), plot = non_selection_plots, width=10, height = 6)

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
                                       bottom = text_grob("Core distance", color = "black", face = "bold", size = 18),
                                       left = text_grob("Accessory distance", color = "black", face = "bold", rot = 90, size = 18))

selection_plots_all <- ggarrange(
  selection_plots_a,
  selection_plots_b,
  nrow = 1,
  ncol = 2,
  align = "hv")

selection_plots_all
ggsave(paste(outdir, "selection_plots.png", sep=""), plot = selection_plots_all, width=14, height = 18)
ggsave(paste(outdir, "selection_plots.svg", sep=""), plot = selection_plots_all, width=14, height = 18)

