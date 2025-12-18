library(phytools)
library(viridis)

overall_indir <- ""
indir <- paste0(overall_indir, "ATB_tree/")
tree.file <- paste0(indir, "data/", "reps_ppmod_tree.treefile")
summary.file <- paste0(indir, "data/", "ppmod_summary.csv")
index.file <- paste0(indir, "data/", "sampled_alignments.tsv")
outpref <- paste0(indir, "figures/")

summary.df <- read.csv(summary.file)
index.df <- read.csv(index.file, sep="\t", header = FALSE)
colnames(index.df) <- c("taxa", "tip", "file_path")

meta <- merge(summary.df, index.df, by = "taxa")
rownames(meta) <- meta$tip

traits <- meta[, c("rate_genes1_median", "core_mu_median", "prop_genes2_median")]
colnames(traits) <- c("Basal gene turnover rate", "Core mutation rate", "Proportion fast genes")

tree <- read.tree(tree.file)

# remove tips with no data
tips_with_data <- rownames(traits)   # or meta$tip
tree <- drop.tip(
  tree,
  setdiff(tree$tip.label, tips_with_data)
)
traits <- traits[tree$tip.label, , drop = FALSE]

# --- midpoint root the tree ---
tree <- midpoint.root(tree)

png(paste0(outpref, "phylo_radial_bars.png"), width = 800, height = 800, res = 150)

# 1. Plot the tree invisibly to get coordinates
plotTree(tree, type = "fan", ftype = "off", lwd = 1, plot = FALSE)  # plot=FALSE ignored by phytools, but let's try

# Instead, plot normally but capture last_plot
plotTree(tree, type = "fan", ftype = "off", lwd = 1)
pp <- get("last_plot.phylo", envir = .PlotPhyloEnv)

# 2. Calculate expanded limits manually
expand_factor <- 1.4
xlim <- range(pp$xx) * expand_factor
ylim <- range(pp$yy) * expand_factor

# 3. Replot with expanded limits
plotTree(
  tree,
  type = "fan",
  ftype = "off",
  lwd = 1,
  xlim = xlim,
  ylim = ylim
)

# --- extract tip angles and positions ---
pp <- get("last_plot.phylo", envir = .PlotPhyloEnv)
ntip <- Ntip(tree)

x <- pp$xx[1:ntip]
y <- pp$yy[1:ntip]

angles <- atan2(y, x)
tip_radius <- sqrt(x^2 + y^2)

# --- fixed base radius for all bars (start radius) ---
base_radius <- max(tip_radius) + 0.1  # slightly outside tips

# --- ring parameters ---
ring_gap   <- 0.02
ring_width <- 0.15
dtheta     <- pi / ntip

# --- Viridis palettes for rings ---
palettes <- c("A", "B", "C")

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


log <- FALSE
ring_colours <- c("#E64B35FF", "#4DBBD5FF", "#00A087FF")
trait_names <- colnames(traits)  # assuming colnames are trait names

for (i in 1:3) {
  
  vals <- traits[, i]
  vals <- vals / max(vals)  # scale 0 to 1, relative to max only
  
  # Add pseudocount to avoid log(0)
  pseudocount <- 1e-6
  if (log == TRUE)
  {
    vals <- log10(vals + pseudocount)
    
    # Scale 0 to 1 based on max logged value
    vals_scaled <- (vals_log - min(vals_log)) / (max(vals_log) - min(vals_log))
  } else
  {
    vals_scaled <- vals
  }
  
  cols <- rep(ring_colours[i], ntip)  # or any colors you like
  
  for (j in 1:ntip) {
    
    r0 <- base_radius + (i - 1) * (ring_width + ring_gap)
    r1 <- r0 + vals_scaled[j] * ring_width
    theta <- angles[j]
    
    polygon(
      x = c(
        r0 * cos(theta - dtheta),
        r1 * cos(theta - dtheta),
        r1 * cos(theta + dtheta),
        r0 * cos(theta + dtheta)
      ),
      y = c(
        r0 * sin(theta - dtheta),
        r1 * sin(theta - dtheta),
        r1 * sin(theta + dtheta),
        r0 * sin(theta + dtheta)
      ),
      col = cols[j],
      border = NA
    )
  }
}

# Add legend outside the plot (adjust x, y as needed)
legend(
  x = "topright",
  legend = trait_names,
  fill = ring_colours,
  border = NA,
  bty = "n",
  cex = 0.5,
  title = ""
)


dev.off()

# pagels lambda trait
lambda_test_results <- lapply(1:ncol(traits), function(i) {
  phylosig(tree, traits[, i], method = "lambda", test = TRUE)
})

names(lambda_test_results) <- colnames(traits)

# Access estimated lambda and p-value for each trait like:
for (trait in names(lambda_test_results)) {
  cat(trait, ": lambda =", lambda_test_results[[trait]]$lambda,
      ", p-value =", lambda_test_results[[trait]]$P, "\n")
}
