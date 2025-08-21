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

files <- Sys.glob("/Users/shorsfield/Library/Mobile Documents/com~apple~CloudDocs/Work/Postdoc/Analysis/PopPUNK-mod/publication_figures/pansim_runs/selection/*.tsv")
indir <- "/Users/shorsfield/Library/Mobile Documents/com~apple~CloudDocs/Work/Postdoc/Analysis/PopPUNK-mod/publication_figures/pansim_runs/selection/"
outdir <- "/Users/shorsfield/Library/Mobile Documents/com~apple~CloudDocs/Work/Postdoc/Analysis/PopPUNK-mod/publication_figures/pansim_runs/"

plot_results_selection <- function(df_paths, parameter, outpref, label)
{
  df_all <- data.frame(Selection_coefficient = c(), prop_positive = c(), HR_rate = c(), HGT_rate = c(), rate_genes1 = c(), prop_genes2 = c(), competition = c(), baseline = c())
  params_all <- data.frame(prop_positive = c(), HR_rate = c(), HGT_rate = c(), rate_genes1 = c(), prop_genes2 = c(), competition = c(), baseline = c())

  j <- 1
  for (j in 1:length(df_paths))
  {
    filename <- df_paths[j]
    df <- read.table(filename, sep = "\t", comment.char = "")
    colnames(df) <- c("Selection_coefficient")
    
    params_df <- data.frame(remove = 1)
    
    # add metadata
    parsed_values <- parse_filename(filename)
    for (name in names(parsed_values)) {
      df[[name]] <- parsed_values[[name]]
      params_df[[name]] <- parsed_values[[name]]
    }
    
    params_df <- params_df[,!names(params_df) %in% c("remove")]
    df_all <- rbind(df_all, df)
    params_all <- rbind(params_all, params_df)
  }
  
  # Separate text and numeric elements
  df_all[[parameter]][df_all[[parameter]] == -0.1] <- "No selection"
  unique_elements <- unique(df_all[[parameter]])
  is_num <- suppressWarnings(!is.na(as.numeric(unique_elements)))
  text_levels <- sort(unique_elements[!is_num])
  numeric_levels <- sort(as.numeric(unique_elements[is_num]))
  
  # Combine text and numeric levels (convert numeric levels back to character)
  all_levels <- c(text_levels, as.character(numeric_levels))
  
  # Create factor with specified level order
  df_all$group <- factor(df_all[[parameter]], levels = all_levels)
  
  params_all$group<- as.factor(params_all[[parameter]])
  
  p <- ggplot(data = df_all, aes(x = Selection_coefficient, colour = group, group = group, fill = group, after_stat(ndensity))) + 
    geom_histogram(colour = "black", bins = 20, position = "stack") +
    #geom_density() +
    #geom_freqpoly(bins=50) +
    geom_vline(xintercept = 0, linewidth=1.0, colour = "black", linetype = 3) +
    theme_light() +
    #scale_fill_npg() +
    scale_fill_brewer(palette = "YlOrRd") + 
    #scale_colour_npg() +
    labs(fill=label) +
    facet_grid(rows = vars(group)) +
    theme_light() + ylab("Density") + xlab("Gene selection coefficient (s)") +
    theme(axis.text.x = element_text(size = 14), axis.text.y = element_text(size = 14), axis.title=element_text(size=16,face="bold"), strip.text.x = element_text(size = 12), legend.text = element_text(size = 12), legend.title = element_text(size = 16, face="bold"), strip.background = element_blank(), strip.text.y = element_blank()) +
    #coord_cartesian(xlim = c(-1, 5))
    xlim(c(-1, 4))
  p
  
  ggsave(paste(outpref, ".svg", sep=""), plot = p, width=10, height = 6)
  ggsave(paste(outpref, ".png", sep=""), plot = p, width=10, height = 6)
  return(p)
}

# prop_postive rate
{
  parameter = "prop_positive"
  label = "Proportion +ve\nselected genes"
  plot_line = FALSE
  contours = TRUE
  facets = TRUE
  facet_bins = 10
  df_paths <- c(paste(indir, "prop_positive_-0.1_HR_rate_0.0_HGT_rate_0.0_rate_genes1_0.001_prop_genes2_0.0_baseline_selection.tsv", sep = ""),
                paste(indir, "prop_positive_0.0_HR_rate_0.0_HGT_rate_0.0_rate_genes1_0.001_prop_genes2_0.0_selection.tsv", sep = ""),
                paste(indir, "prop_positive_0.1_HR_rate_0.0_HGT_rate_0.0_rate_genes1_0.001_prop_genes2_0.0_selection.tsv", sep = ""),
                paste(indir, "prop_positive_1.0_HR_rate_0.0_HGT_rate_0.0_rate_genes1_0.001_prop_genes2_0.0_selection.tsv", sep = ""))
  outpref <- paste(outdir, "selection_", parameter, sep = "")
  p <- plot_results_selection(df_paths, parameter, outpref, label)
  assign(paste(parameter, "selection", sep="_"), p)
}

# High/low positive
{
  parameter = "pos_lambda"
  label = "Gene +ve\nselection\nlambda"
  plot_line = FALSE
  contours = TRUE
  facets = TRUE
  facet_bins = 10
  df_paths <- c(paste(indir, "prop_positive_-0.1_HR_rate_0.0_HGT_rate_0.0_rate_genes1_0.001_prop_genes2_0.0_baseline_selection.tsv", sep = ""),
                paste(indir, "prop_positive_0.5_HR_rate_0.0_HGT_rate_0.0_rate_genes1_0.001_prop_genes2_0.0_pos_lambda_0.1_neg_lambda_100.0_selection.tsv", sep = ""),
                paste(indir, "prop_positive_0.5_HR_rate_0.0_HGT_rate_0.0_rate_genes1_0.001_prop_genes2_0.0_pos_lambda_1.0_neg_lambda_100.0_selection.tsv", sep = ""),
                paste(indir, "prop_positive_0.5_HR_rate_0.0_HGT_rate_0.0_rate_genes1_0.001_prop_genes2_0.0_pos_lambda_10.0_neg_lambda_100.0_selection.tsv", sep = ""))
  outpref <- paste(outdir, "selection_", parameter, sep = "")
  p <- plot_results_selection(df_paths, parameter, outpref, label)
  assign(paste(parameter, "selection", sep="_"), p)
}

# High/low negative
{
  parameter = "neg_lambda"
  label = "Gene -ve\nselection\nlambda"
  plot_line = FALSE
  contours = TRUE
  facets = TRUE
  facet_bins = 10
  df_paths <- c(paste(indir, "prop_positive_-0.1_HR_rate_0.0_HGT_rate_0.0_rate_genes1_0.001_prop_genes2_0.0_baseline_selection.tsv", sep = ""),
                paste(indir, "prop_positive_0.5_HR_rate_0.0_HGT_rate_0.0_rate_genes1_0.001_prop_genes2_0.0_pos_lambda_100.0_neg_lambda_0.1_selection.tsv", sep = ""),
                paste(indir, "prop_positive_0.5_HR_rate_0.0_HGT_rate_0.0_rate_genes1_0.001_prop_genes2_0.0_pos_lambda_100.0_neg_lambda_1.0_selection.tsv", sep = ""),
                paste(indir, "prop_positive_0.5_HR_rate_0.0_HGT_rate_0.0_rate_genes1_0.001_prop_genes2_0.0_pos_lambda_100.0_neg_lambda_10.0_selection.tsv", sep = ""))
  outpref <- paste(outdir, "selection_", parameter, sep = "")
  p <- plot_results_selection(df_paths, parameter, outpref, label)
  assign(paste(parameter, "selection", sep="_"), p)
}

selection_plots_b <- ggarrange(
  prop_positive_selection + xlab("") + ylab(""), 
  pos_lambda_selection + xlab("") + ylab(""),
  neg_lambda_selection + xlab("") + ylab(""),
  nrow = 3,
  ncol = 1,
  align = "hv")

selection_plots_b <- annotate_figure(selection_plots_b,
                                   bottom = text_grob("Gene selection coefficient (s)", color = "black", face = "bold", size = 18),
                                   left = text_grob("Density", color = "black", face = "bold", rot = 90, size = 18))
selection_plots_b
