library(ggplot2)
library(dplyr)
library(ggsci)
library(GGally)
library(tools)
library(ggExtra)
library(purrr)
library(ggpubr)
library(scales)
library(reporter)

add_headers <- function(plots, nrow, ncol, row_headers=NULL, col_headers=NULL) {
  stopifnot(length(plots) == nrow * ncol)
  
  # --- Column headers (must align with row headers setup) ---
  if (!is.null(col_headers)) {
    col_grobs <- lapply(col_headers, function(lbl) {
      ggplot() + 
        annotate("text", x=0.5, y=0.5, label=lbl, 
                 size=5, fontface="bold", vjust=0.9, hjust=0.1) +
        theme_void() +
        theme(plot.margin = margin(0,0,0,0))
    })
    
    # add an empty cell on the right for row labels
    header_row <- ggarrange(
      plotlist = c(col_grobs, list(NULL)), 
      ncol = ncol + 1,
      widths = c(rep(1, ncol), 0.2)
    )
  } else {
    header_row <- NULL
  }
  
  # --- Rows with row headers on the right ---
  row_list <- lapply(seq_len(nrow), function(i) {
    row_plots <- plots[((i-1)*ncol+1):(i*ncol)]
    
    if (!is.null(row_headers)) {
      row_lab <- ggplot() + 
        annotate("text", x=0.5, y=0.5, label=row_headers[i], 
                 angle=90, size=5, fontface="bold", vjust=0.75, hjust=0.1) +
        theme_void() +
        theme(plot.margin = margin(0,0,0,0))
      ggarrange(
        plotlist = c(row_plots, list(row_lab)), 
        ncol = ncol + 1,
        widths = c(rep(1, ncol), 0.2)
      )
    } else {
      ggarrange(plotlist=row_plots, ncol=ncol)
    }
  })
  
  # --- Stack headers + rows ---
  if (!is.null(header_row)) {
    ggarrange(
      header_row, 
      ggarrange(plotlist=row_list, ncol=1),
      ncol=1,
      heights=c(0.15, 1)  # adjust header spacing
    )
  } else {
    ggarrange(plotlist=row_list, ncol=1)
  }
}


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
    competition_v = search_string(filename, "competition"),
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

plot_posterior <- function(df, col.names, outpref, max.values)
{
  idx_keep <- which(!grepl("_v$", colnames(df)))
  col.names.ori <- colnames(df)[idx_keep]
  colnames(df)[idx_keep] <- col.names
  
  idx_pairs <- combn(idx_keep, 2, simplify = FALSE)
  
  # p <- ggpairs(
  #   df[, idx_keep],
  #   diag = list(continuous = custom_diag_hist),   # histograms
  #   lower = list(continuous = scatter_with_contour),   # scatterplots
  #   upper = list(continuous = "blank")     # blank upper triangle
  # )
  # p + theme_light()
  # p
  #col.names.pairs <- combn(col.names, 2, simplify = FALSE)
  #cols.name.keep.pairs <- combn(cols.name.keep, 2, simplify = FALSE)
  
  # p <- ggscatmat(
  #   df,
  #   columns = cols_keep,
  #   alpha = 0.01,          # transparency of points
  # )
  # 
  # p <- p + 
  #   theme_light(base_size = 12) + 
  #   theme(
  #     axis.title.x = element_text(face = "bold"),
  #     axis.title.y = element_text(face = "bold")
  #   ) + xlab("X-values") + ylab("Y-values")
  # p
  
  # Generate marginal plots for each pair
  idx <- c(1,2)
  plots <- lapply(idx_pairs, function(idx) {
    xvar <- idx_keep[idx[1]]
    yvar <- idx_keep[idx[2]]
    
    x.limits <- max.values[[idx[1]]]
    y.limits <- max.values[[idx[2]]]
    
    # get original column names
    x.name.ori <- col.names.ori[[xvar]]
    y.name.ori <- col.names.ori[[yvar]]
    
    # get new column names
    x.name <- col.names[[xvar]]
    y.name <- col.names[[yvar]]
    
    # Corresponding real values
    x_real <- unique(df[[paste0(x.name.ori, "_v")]])
    y_real <- unique(df[[paste0(y.name.ori, "_v")]])
    
    x.median <- median(df[[x.name]])
    y.median <- median(df[[y.name]])
    
    x.5quantile <- quantile(df[[x.name]], 0.025)
    y.5quantile <- quantile(df[[y.name]], 0.025)
    x.95quantile <- quantile(df[[x.name]], 0.975)
    y.95quantile <- quantile(df[[y.name]], 0.975)
    
    v_line_color <- if (x_real >= x.5quantile && x_real <= x.95quantile) {
      "#00A087FF"
    } else {
      "#E64B35FF"
    }
    
    h_line_color <- if (y_real >= y.5quantile && y_real <= y.95quantile) {
      "#00A087FF"
    } else {
      "#E64B35FF"
    }
    
    
    p <- ggplot(df, aes(x = .data[[x.name]], y = .data[[y.name]])) +
      geom_point(alpha = 0.5, size=0.5, colour = "#808180FF") +
      geom_density_2d(color = "#1B1919FF", alpha = 0.6) +
      
      geom_vline(xintercept = x.median, linetype = "solid", color = "#3C5488FF") +
      geom_hline(yintercept = y.median, linetype = "solid", color = "#3C5488FF") +
      
      geom_vline(xintercept = x_real, linetype = "dashed", color = v_line_color) +
      geom_hline(yintercept = y_real, linetype = "dashed", color = h_line_color) +
      
      # apply limits
      #xlim(x.limits[1], x.limits[2]) + ylim(y.limits[1], y.limits[2]) +
      scale_y_log10(limits=y.limits, labels = trans_format("log10", math_format(10^.x))) + 
      scale_x_log10(limits=x.limits, labels = trans_format("log10", math_format(10^.x))) +
      
      xlab(x.name) + 
      ylab(y.name) +
      theme_minimal()
    p <- p + annotate(
      "rect",
      xmin = x.5quantile, xmax = x.95quantile,   # x range of box
      ymin = y.5quantile, ymax = y.95quantile,   # y range of box
      alpha = 0.1,            # transparency
      fill = "#4DBBD5FF",           # color
      color = "black",
      linetype = "solid",
      linewidth = 0.5
    ) + xlab("") + ylab("")
    
    p
    ggMarginal(p, type = "density", fill="#808180FF")
  })
  
  plots[[1]]
  
  # Save each plot with informative filenames
  mapply(function(plot, idx) {
    xvar <- idx_keep[idx[1]]
    yvar <- idx_keep[idx[2]]
    
    # get original column names
    x.name.ori <- col.names.ori[[xvar]]
    y.name.ori <- col.names.ori[[yvar]]
    
    fname <- paste0(outpref, "_", x.name.ori, "_vs_", y.name.ori, ".png")
    ggsave(fname, plot, width = 5, height = 5)
  }, plots, idx_pairs)
  
  plots
}

parse_results <- function(df_paths, col.names, outpref, max.values)
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
    
    plots <- plot_posterior(df, col.names, paste(outpref, base.name, sep = "/"), max.values)
    
    # create df_all if first iteration
    if (j == 1)
    {
      df_all <- df
      plots_all <- list()
      plots_all <- append(plots_all, plots)
    } else {
      df_all <- rbind(df_all, df)
      plots_all <- append(plots_all, plots)
    }
  }
  return_list = list(df_all, plots_all)
  
  return(return_list)
}

# baseline
{
  indir <- "BOLFI_gridsearch/data/baseline/"
  #df_paths <- Sys.glob(paste0(indir, "*_mcmc_posterior_samples.csv"))
  df_paths <- c(Sys.glob(paste0(indir, "*core_mu_1e-3*rate_genes1_1e-5*_mcmc_posterior_samples.csv")),
                Sys.glob(paste0(indir, "*core_mu_1e-3*rate_genes1_1e-4*_mcmc_posterior_samples.csv")),
                Sys.glob(paste0(indir, "*core_mu_1e-3*rate_genes1_1e-3*_mcmc_posterior_samples.csv")),
                Sys.glob(paste0(indir, "*core_mu_1e-3*rate_genes1_1e-2*_mcmc_posterior_samples.csv")),
                Sys.glob(paste0(indir, "*core_mu_1e-3*rate_genes1_1e-1*_mcmc_posterior_samples.csv")),
                Sys.glob(paste0(indir, "*core_mu_1e-4*rate_genes1_1e-5*_mcmc_posterior_samples.csv")),
                Sys.glob(paste0(indir, "*core_mu_1e-4*rate_genes1_1e-4*_mcmc_posterior_samples.csv")),
                Sys.glob(paste0(indir, "*core_mu_1e-4*rate_genes1_1e-3*_mcmc_posterior_samples.csv")),
                Sys.glob(paste0(indir, "*core_mu_1e-4*rate_genes1_1e-2*_mcmc_posterior_samples.csv")),
                Sys.glob(paste0(indir, "*core_mu_1e-4*rate_genes1_1e-1*_mcmc_posterior_samples.csv")),
                Sys.glob(paste0(indir, "*core_mu_1e-5*rate_genes1_1e-5*_mcmc_posterior_samples.csv")),
                Sys.glob(paste0(indir, "*core_mu_1e-5*rate_genes1_1e-4*_mcmc_posterior_samples.csv")),
                Sys.glob(paste0(indir, "*core_mu_1e-5*rate_genes1_1e-3*_mcmc_posterior_samples.csv")),
                Sys.glob(paste0(indir, "*core_mu_1e-5*rate_genes1_1e-2*_mcmc_posterior_samples.csv")),
                Sys.glob(paste0(indir, "*core_mu_1e-5*rate_genes1_1e-1*_mcmc_posterior_samples.csv")),
                Sys.glob(paste0(indir, "*core_mu_1e-6*rate_genes1_1e-5*_mcmc_posterior_samples.csv")),
                Sys.glob(paste0(indir, "*core_mu_1e-6*rate_genes1_1e-4*_mcmc_posterior_samples.csv")),
                Sys.glob(paste0(indir, "*core_mu_1e-6*rate_genes1_1e-3*_mcmc_posterior_samples.csv")),
                Sys.glob(paste0(indir, "*core_mu_1e-6*rate_genes1_1e-2*_mcmc_posterior_samples.csv")),
                Sys.glob(paste0(indir, "*core_mu_1e-6*rate_genes1_1e-1*_mcmc_posterior_samples.csv")),
                Sys.glob(paste0(indir, "*core_mu_1e-7*rate_genes1_1e-5*_mcmc_posterior_samples.csv")),
                Sys.glob(paste0(indir, "*core_mu_1e-7*rate_genes1_1e-4*_mcmc_posterior_samples.csv")),
                Sys.glob(paste0(indir, "*core_mu_1e-7*rate_genes1_1e-3*_mcmc_posterior_samples.csv")),
                Sys.glob(paste0(indir, "*core_mu_1e-7*rate_genes1_1e-2*_mcmc_posterior_samples.csv")),
                Sys.glob(paste0(indir, "*core_mu_1e-7*rate_genes1_1e-1*_mcmc_posterior_samples.csv")))
  
  outpref <- "BOLFI_gridsearch/figures/baseline"
  # get columns which don't have value placeholders and name if required
  col.names <- c("Basal gene turnover rate", "Core mutation rate")
  max.values <- list(c(1e-5, 1.0), c(1e-10, 1e-2))
  #max.values <- list(c(NA, NA), c(NA, NA))
  return_list <- parse_results(df_paths, col.names, outpref, max.values)
  
  row_headers = c(paste0('10',supsc('-3')), paste0('10',supsc('-4')), paste0('10',supsc('-5')), paste0('10',supsc('-6')), paste0('10',supsc('-7')))
  col_headers = c(paste0('10',supsc('-5')), paste0('10',supsc('-4')), paste0('10',supsc('-3')), paste0('10',supsc('-2')), paste0('10',supsc('-1')))
  
  # -------------------------------
  # Example usage with your setup
  # -------------------------------
  baseline_plots <- add_headers(
    plots = return_list[[2]], 
    nrow = length(row_headers),
    ncol = length(col_headers),
    row_headers = row_headers,
    col_headers = col_headers
  )
  #grid_with_headers
  #ggsave(paste(outpref, "/baseline_BOLFI_fit_grid.png", sep=""), plot = baseline_plots, width=8, height = 6.5)
  
  # combine plots
  # baseline_plots <- ggarrange(
  #   plotlist=return_list[[2]],
  #   ncol = 5,
  #   nrow = 5,
  #   align = "hv")
  
  baseline_plots <- annotate_figure(baseline_plots,
                                         bottom = text_grob("Basal gene turnover rate", color = "black", face = "bold", size = 18),
                                         left = text_grob("Core mutation rate", color = "black", face = "bold", rot = 90, size = 18))
  
  #baseline_plots
  ggsave(paste(outpref, "/baseline_BOLFI_fit.png", sep=""), plot = baseline_plots, width=8, height = 6.5)
  #ggsave(paste(outpref, "/baseline_BOLFI_fit.svg", sep=""), plot = baseline_plots, width=8, height = 6.5)
  
}




