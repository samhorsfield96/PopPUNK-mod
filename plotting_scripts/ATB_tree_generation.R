library(phytools)
library(viridis)
library(svglite)
library(ggplot2)
library(MCMCglmm)
library(dplyr)
library(stringr)
library(ape)
library(coda)
library(tools)
library(tidyr)
library(paletteer)

overall_indir <- ""
indir <- paste0(overall_indir, "ATB_tree/")
tree.file <- paste0(indir, "data/", "reps_ppmod_tree.treefile")
summary.file <- paste0(indir, "data/", "ppmod_summary.csv")
index.file <- paste0(indir, "data/", "sampled_alignments.tsv")
gtdb.file <- paste0(indir, "data/", "gtdbtk.bac120.summary.tsv")
outpref <- paste0(indir, "figures/")

summary.df <- read.csv(summary.file, row.names = 1)

values.NA <- summary.df[is.na(summary.df$order),]

index.df <- read.csv(index.file, sep="\t", header = FALSE)
colnames(index.df) <- c("taxa", "tip", "file_path")

meta <- merge(summary.df, index.df, by = "taxa")

# read in taxonomy
gtdb.df <- read.csv(gtdb.file, sep="\t", header = TRUE)
gtdb.df <- select(gtdb.df, c('user_genome', 'classification'))
gtdb.df <- gtdb.df |> separate_wider_delim(classification, delim = ";", names = c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species"))
gtdb.df$Domain <- gsub(".*__", "", gtdb.df$Domain)
gtdb.df$Phylum <- gsub(".*__", "", gtdb.df$Phylum)
gtdb.df$Class <- gsub(".*__", "", gtdb.df$Class)
gtdb.df$Order <- gsub(".*__", "", gtdb.df$Order)
gtdb.df$Family <- gsub(".*__", "", gtdb.df$Family)
gtdb.df$Genus <- gsub(".*__", "", gtdb.df$Genus)
gtdb.df$Species <- gsub(".*__", "", gtdb.df$Species)

gtdb.df <- as.data.frame(gtdb.df)

meta <- merge(meta, gtdb.df, by.x = "tip", by.y = "user_genome")

# generate CI range
meta$rate_genes1_CI_range <- meta$rate_genes1_upper_CI - meta$rate_genes1_lower_CI
meta$core_mu_CI_range <- meta$core_mu_upper_CI - meta$core_mu_lower_CI
meta$prop_genes2_CI_range <- meta$prop_genes2_upper_CI - meta$prop_genes2_lower_CI

# log rate_genes1_median
# meta$rate_genes1_median <- log(meta$rate_genes1_median)
# meta$rate_genes1_lower_CI <- log(meta$rate_genes1_lower_CI)
# meta$rate_genes1_upper_CI <- log(meta$rate_genes1_upper_CI)

rownames(meta) <- meta$tip

# Choose interval to analyse
#traits <- meta[, c("rate_genes1_median", "core_mu_median", "prop_genes2_median")]
#traits <- meta[, c("rate_genes1_lower_CI", "core_mu_lower_CI", "prop_genes2_lower_CI")]
#traits <- meta[, c("rate_genes1_upper_CI", "core_mu_upper_CI", "prop_genes2_upper_CI")]
traits <- meta[, c("rate_genes1_CI_range", "core_mu_CI_range", "prop_genes2_CI_range")]

colnames(traits) <- c("Basal gene turnover rate", "Core mutation rate", "Proportion fast genes")

# -------------------------------
# Load tree
# -------------------------------
tree <- read.tree(tree.file)

# -------------------------------
# Taxonomy vector (named!)
# -------------------------------
tax <- meta$Phylum
names(tax) <- rownames(meta)
tax <- factor(tax)

# -------------------------------
# Drop tips with no taxonomy
# -------------------------------
tips_with_data <- rownames(traits)   # or meta$tip
tree <- drop.tip(
  tree,
  setdiff(tree$tip.label, tips_with_data)
)

# -------------------------------
# Root tree
# -------------------------------
tree <- midpoint.root(tree)

# -------------------------------
# Reorder taxonomy to tree tips
# -------------------------------
tax <- tax[tree$tip.label]
tax <- droplevels(tax)
traits <- traits[tree$tip.label, , drop = FALSE]

# -------------------------------
# HARD SAFETY CHECKS
# -------------------------------
stopifnot(
  length(tax) == length(tree$tip.label),
  all(names(tax) == tree$tip.label),
  !any(is.na(tax)),
  nrow(traits) == length(tree$tip.label),
  all(rownames(traits) == tree$tip.label),
  !any(is.na(traits))
)

# calculate pagel's lambda
lambda_results <- lapply(seq_len(ncol(traits)), function(i) {
  phylosig(
    tree,
    traits[, i],
    method = "lambda",
    test = TRUE
  )
})

names(lambda_results) <- colnames(traits)

# individual trees, internal branches annotated
{
  png(
    paste0(outpref, "phylo_internal_branches_indiv_trees_CI_range.png"),
    units = "in",
    width = 7.25, height = 8, res = 150
  )

  # svglite(
  #   filename = paste0(outpref, "phylo_internal_branches_indiv_trees_median.svg"),
  #   width = 7.25,    # width in inches
  #   height = 8,  # height in inches (scale up as needed)
  # )

  
  # Layout: 3 rows for plots + 1 column for legend on the right
  # The layout matrix defines 2 columns: left for plots (panels 1,2,3 stacked), right for legend
  layout_matrix <- matrix(c(1,4,
                            2,4,
                            3,4), ncol = 2, byrow = TRUE)
  layout(layout_matrix, widths = c(5, 2))
  
  # Common margins for plots (no right margin, since legend is separate)
  par(mar = c(1, 1, 1, 1), oma = c(1, 1, 1, 1))
  
  panel_labels <- c("A", "B", "C")
  
  i <- 2
  for (i in seq_len(ncol(traits))) {
    
    trait <- colnames(traits)[i]
    sig   <- lambda_results[[trait]]
    
    ## ---- extract + name trait vector ----
    vals <- traits[, i]
    names(vals) <- rownames(traits)
    
    ## ---- drop tips without data ----
    keep <- intersect(tree$tip.label, names(vals))
    vals <- vals[keep]
    tree_i <- drop.tip(tree, setdiff(tree$tip.label, keep))
    
    ## ---- numeric + finite ----
    vals <- as.numeric(vals)
    names(vals) <- keep
    
    good <- is.finite(vals)
    vals <- vals[good]
    tree_i <- drop.tip(tree_i, setdiff(tree_i$tip.label, names(vals)))
    
    ## ---- skip invariant traits ----
    if (length(unique(vals)) < 2) next
    
    ## ---- reorder EXACTLY to tree ----
    vals <- vals[tree_i$tip.label]
    
    ## ---- scale per trait ----
    vals <- (vals - min(vals)) / (max(vals) - min(vals))
    
    ## ---- build contMap ----
    cm <- contMap(tree_i, vals, plot = FALSE)
    #colours <- c("#33608C", "#F5B355","#B81840")
    #colours <- c("#2E5A87","white", "#B81840")
    colours <- c("#72BCDC","#FD8E3F","white")
    cm <- setMap(cm, colors = colours)
    
    cm$border <- adjustcolor("black", alpha.f = 0.1)  # make outlines transparent
  
    
    ## ---- plot THIS panel ----
    plot(
      cm,
      type  = "fan",
      ftype = "off",
      lwd   = 1.1,
      legend = FALSE,
      outline = FALSE
    )
    
    ## ---- panel label (A/B/C) ----
    mtext(
      panel_labels[i],
      side = 3,
      adj  = 0,
      line = -1.2,
      cex  = 1.4,
      font = 2
    )
    
    ## ---- get plotting coordinates ----
    pp <- get("last_plot.phylo", envir = .PlotPhyloEnv)

    ntip <- length(tree_i$tip.label)

    x <- pp$xx[1:ntip]
    y <- pp$yy[1:ntip]

    angles <- atan2(y, x)
    tip_r  <- sqrt(x^2 + y^2)
    
    # --- fixed base radius for all bars (start radius) ---
    base_radius <- max(tip_r) + 0.1  # slightly outside tips

    ## --- ring geometry ---
    ring_inner <- max(tip_r) * 1.05
    ring_width <- 0.08
    dtheta     <- pi / ntip

    #base_cols <- RColorBrewer::brewer.pal(8, "Set3")
    base_cols <- paletteer_d("ggsci::category20_d3")
    #base_cols <- paletteer_d("colorBlindness::SteppedSequential5Steps")
    #base_cols <- paletteer_d("colorBlindness::Blue2DarkOrange12Steps")
    #base_cols <- paletteer_d("ggsci::nrc_npg")

    # tax_cols <- setNames(
    #   colorRampPalette(base_cols)(length(levels(tax))),
    #   levels(tax)
    # )
    
    tax_cols <- setNames(
      c(base_cols[1:length(levels(tax))]),
      levels(tax)
    )

    ## --- draw blocks ---
    for (i in seq_len(ntip)) {
      theta <- angles[i]
      col   <- tax_cols[tax[i]]

      polygon(
        x = c(
          ring_inner * cos(theta - dtheta),
          (ring_inner + ring_width) * cos(theta - dtheta),
          (ring_inner + ring_width) * cos(theta + dtheta),
          ring_inner * cos(theta + dtheta)
        ),
        y = c(
          ring_inner * sin(theta - dtheta),
          (ring_inner + ring_width) * sin(theta - dtheta),
          (ring_inner + ring_width) * sin(theta + dtheta),
          ring_inner * sin(theta + dtheta)
        ),
        col = col,
        border = NA
      )
    }
    
    # --- draw connecting lines from tip to bars ---
    for (j in 1:ntip) {
      segments(
        x0 = x[j], y0 = y[j],
        x1 = base_radius * cos(angles[j]),
        y1 = base_radius * sin(angles[j]),
        col = "grey40",
        lwd = 0.1
      )
    }
    
    ## ---- caption to the right ----
    x_cap <- max(pp$xx) * 1.15
    y_cap <- mean(range(pp$yy))
    
    text(
      x = x_cap,
      y = y_cap,
      labels = sprintf(
        "%s\nPagel's \u03BB = %.2f\np = %.3g",
        trait, sig$lambda, sig$P
      ),
      adj = c(0, 0.5),
      cex = 0.9
    )
  }
  
  # Now plot the legend in the reserved right panel
  par(mar = c(5, 2, 4, 2))  # bigger margins for legend panel
  plot.new()
  legend(
    "center",
    legend = names(tax_cols),
    fill = tax_cols,
    border = NA,
    bty = "n",
    cex = 0.8,
    y.intersp = 1.2,
    x.intersp = 1,
    title = "Phyla"
  )

  dev.off()
}