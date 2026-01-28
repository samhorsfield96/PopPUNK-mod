library(phytools)
library(viridis)
library(svglite)
library(ggplot2)
library(MCMCglmm)
library(dplyr)
library(stringr)
library(ape)
library(coda)
library(taxize)
library(tools)
library(tidyr)


# from https://github.com/AnnaEDewar/pangenome_lifestyle/blob/main/Code_S1.R
summary_mcmc_glmm <- function(mcmc_model) {
  summary <- summary(mcmc_model)
  summary_subset <- as.data.frame(summary$solutions)
  
  p_values <- summary_subset$pMCMC
  stars <- ifelse(p_values < 0.001, "***",
                  ifelse(p_values < 0.01, "**",
                         ifelse(p_values < 0.05, "*",
                                ifelse(p_values < 0.1, ".","ns"))))
  summary_subset$signif. <- stars
  
  # Extract numeric columns
  num_cols <- sapply(summary_subset, is.numeric)
  numeric_table <- summary_subset[, num_cols]
  
  # Format numeric columns to 4 significant figures
  formatted_table <- format(numeric_table, scientific = FALSE, digits = 4)
  
  # Replace original numeric columns with formatted columns
  summary_subset[, num_cols] <- formatted_table
  
  return(summary_subset)
}

overall_indir <- "/Users/shorsfield/Library/Mobile Documents/com~apple~CloudDocs/Work/Postdoc/Analysis/PopPUNK-mod/publication_figures/"
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
rownames(meta) <- meta$tip
traits <- meta[, c("rate_genes1_median", "core_mu_median", "prop_genes2_median")]
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

# add traits
i <- 1
for (i in seq_len(ncol(traits))) {
  
  trait <- colnames(traits)[i]
  vals <- traits[, i]
  
  mean <- mean(vals)
  sd <- sd(vals)
  
  # As a unitless value
  cv <- sd / mean
  
  Q3 <- quantile(vals, 0.75)
  Q2 <- quantile(vals, 0.5)
  Q1 <- quantile(vals, 0.25)
  
  IQR <- (Q3 - Q1)
  
  QCV = (Q3 - Q1) / (Q3 + Q1)
  
  # cat(trait, "\n",
  #     "Mean :", mean, "\n",
  #     "SD:", sd, "\n",
  #     "CV", cv, "\n",
  #     "Median :", Q2, "\n",
  #     "IQR :", IQR, "\n",
  #     "QCV :", QCV, "\n"
  # )
  
  temp.df <- data.frame(Trait = trait, Mean = mean, SD = sd, CV = cv, Median = Q2, IQR = IQR, QCV = round(QCV, digits=3), QCV_text = paste0("QCD: ", round(QCV, digits=3)))
  
  if (i == 1) {
    summary.stats.df <- temp.df
  } else {
    summary.stats.df <- rbind(summary.stats.df, temp.df)
  }
}

traits.stacked <- stack(traits, select = c("Basal gene turnover rate", "Core mutation rate", "Proportion fast genes"))
colnames(traits.stacked) <- c("Value", "Trait")

# generate summary statistics of traits
p <- ggplot(traits.stacked, aes(x=Trait, y = Value)) + 
  theme_light() + 
  geom_text(data = summary.stats.df, size = 5, aes(x = Trait, y = 10, label = QCV_text)) +
  facet_wrap(Trait ~ ., scales = "free_x") +
  geom_violin() + geom_jitter(alpha=0.2) +
  scale_y_log10() +
  xlab("Trait") + ylab("Trait Value") + 
  stat_summary(fun=median, colour="red", geom="text", show.legend = FALSE, vjust=-1, size = 5, aes( label=format(signif(after_stat(y), digits=3), scientific = TRUE))) + 
  stat_summary(fun.y=median, colour="red", geom="point", shape=18, size=2, show.legend= FALSE) +
  theme(axis.text.x = element_text(size = 18), axis.text.y = element_text(size = 18), axis.title=element_text(size=20,face="bold"), legend.position = "none", strip.background = element_blank(), strip.text.x = element_blank())
p
ggsave(paste(outpref, "trait_values.png", sep=""), plot = p, height = 8, width = 10)

# MCMCglmm analysis
# parse generalism analysis
traits.file <- paste0(indir, "data/", "metadata_with_traits.tsv")
traits.df <- read.csv(traits.file, sep="\t", header = TRUE)
to.keep <- c("taxonID.GTDB", "total.genomes", "majority.fraction", "species..GTDB.", "strain..GTDB.", "generalism_score", "habitat_count", "pangenome_openness", "intra_species_nucleotide_diversity", "generalist")

traits.df <- traits.df[, to.keep]
traits.df <- traits.df[!is.na(traits.df$generalism_score),]

traits.df.final <- traits.df %>% 
  group_by(species..GTDB.) %>% 
  summarize(generalism_score=mean(generalism_score), habitat_count=mean(habitat_count), pangenome_openness=mean(pangenome_openness), intra_species_nucleotide_diversity=mean(intra_species_nucleotide_diversity))

traits.df.final <- traits.df.final %>%
  mutate(label = str_to_lower(str_replace_all(species..GTDB., "\\s+", "_")))

meta$tip <- rownames(meta)

# match across datasets
merged.df <- merge(meta, traits.df.final, by.x = "Species", by.y = "species..GTDB.")
rownames(merged.df) <- merged.df$tip

# MCMCglmm
tree <- read.tree(tree.file)

# remove tips with no data
tips_with_data <- merged.df$tip   # or meta$tip
tree <- drop.tip(
  tree,
  setdiff(tree$tip.label, tips_with_data)
)

# --- midpoint root the tree ---
tree <- midpoint.root(tree)

# -------------------------------
# Reorder taxonomy to tree tips
# -------------------------------
merged.df <- merged.df[tree$tip.label, , drop = FALSE]

# -------------------------------
# HARD SAFETY CHECKS
# -------------------------------
stopifnot(
  nrow(merged.df) == length(tree$tip.label),
  all(rownames(merged.df) == tree$tip.label)
)

Ainv <- inverseA(tree, nodes = "TIPS", scale = FALSE)$Ainv

prior <- list(
  G = list(
    G1 = list(V = 1, nu = 1)
  ),
  R = list(
    R1 = list(V = 1, nu = 1)
  )
)

# with tree
{
  var <- "prop_genes2_median"
  mcmc_model <- MCMCglmm(generalism_score ~ prop_genes2_median, random=~tip, 
                                   data=merged.df, prior = prior,
                                   nitt=50000, ginverse = list(tip = Ainv),verbose = FALSE)
  
  mcmc_summary <- summary_mcmc_glmm(mcmc_model)
  
  effectiveSize(mcmc_model$Sol)
  effectiveSize(mcmc_model$VCV)
  
  autocorr.diag(mcmc_model$Sol)
  autocorr.diag(mcmc_model$VCV)
  
  summary(mcmc_model$Sol)
  summary(mcmc_model$VCV)
  
  # Calculate fixed effects predicted values variance
  beta_fixed <- colMeans(mcmc_model$Sol)
  X <- model.matrix(~ prop_genes2_median, data = merged.df)
  mFixed <- as.numeric(X %*% beta_fixed)
  mVarF <- var(mFixed)
  
  # Extract posterior mean variance components
  VCV_mean <- colMeans(mcmc_model$VCV)
  
  # Random effect variance (phylogeny)
  random_var <- VCV_mean["tip"]
  
  # Residual variance
  resid_var <- VCV_mean["units"]
  
  # Calculate R² values
  fixed_r2 <- mVarF / (mVarF + random_var + resid_var)
  random_r2 <- random_var / (mVarF + random_var + resid_var)
  conditional_r2 <- (mVarF + random_var) / (mVarF + random_var + resid_var)
  
  # Create a summary table
  R2_table <- data.frame(
    Component = c("Fixed effect", "Random effect (phylogeny)", "Conditional (total)"),
    R2 = round(c(fixed_r2, random_r2, conditional_r2), 4)
  )
  
  mcmc_summary$param <- var
  R2_table$param <- var
  
  mcmc_summary_all <- mcmc_summary
  R2_table_all <- R2_table
  
  # model checks
  png(filename=paste0(outpref, var, "_MCMC_variance_components.png"))
  plot(mcmc_model$VCV)    # variance components
  dev.off()
  
  
  png(filename=paste0(outpref, var, "_MCMC_fixed_components.png"))
  plot(mcmc_model$Sol)    # fixed effects
  dev.off()
}


{
  var <- "rate_genes1_median"
  mcmc_model <- MCMCglmm(generalism_score ~ rate_genes1_median, random=~tip, 
                                   data=merged.df, prior = prior,
                                   nitt=50000, ginverse = list(tip = Ainv),verbose = FALSE)
  
  
  mcmc_summary <- summary_mcmc_glmm(mcmc_model)
  
  effectiveSize(mcmc_model$Sol)
  effectiveSize(mcmc_model$VCV)
  
  autocorr.diag(mcmc_model$Sol)
  autocorr.diag(mcmc_model$VCV)
  
  summary_mcmc_model <- summary_mcmc_glmm(mcmc_model)
  
  summary(mcmc_model$Sol)
  summary(mcmc_model$VCV)
  
  # Calculate fixed effects predicted values variance
  beta_fixed <- colMeans(mcmc_model$Sol)
  X <- model.matrix(~ rate_genes1_median, data = merged.df)
  mFixed <- as.numeric(X %*% beta_fixed)
  mVarF <- var(mFixed)
  
  # Extract posterior mean variance components
  VCV_mean <- colMeans(mcmc_model$VCV)
  
  # Random effect variance (phylogeny)
  random_var <- VCV_mean["tip"]
  
  # Residual variance
  resid_var <- VCV_mean["units"]
  
  # Calculate R² values
  fixed_r2 <- mVarF / (mVarF + random_var + resid_var)
  random_r2 <- random_var / (mVarF + random_var + resid_var)
  conditional_r2 <- (mVarF + random_var) / (mVarF + random_var + resid_var)
  
  # Create a summary table
  R2_table <- data.frame(
    Component = c("Fixed effect", "Random effect (phylogeny)", "Conditional (total)"),
    R2 = round(c(fixed_r2, random_r2, conditional_r2), 4)
  )
  
  mcmc_summary$param <- var
  R2_table$param <- var
  
  mcmc_summary_all <- rbind(mcmc_summary_all, mcmc_summary)
  R2_table_all <- rbind(R2_table_all, R2_table)
  
  # model checks
  png(filename=paste0(outpref, var, "_MCMC_variance_components.png"))
  plot(mcmc_model$VCV)    # variance components
  dev.off()
  
  
  png(filename=paste0(outpref, var, "_MCMC_fixed_components.png"))
  plot(mcmc_model$Sol)    # fixed effects
  dev.off()
}

# core_mu_median
{
  var <- "core_mu_median"
  mcmc_model <- MCMCglmm(generalism_score ~ core_mu_median, random=~tip, 
                         data=merged.df, prior = prior,
                         nitt=50000, ginverse = list(tip = Ainv),verbose = FALSE)
  
  
  mcmc_summary <- summary_mcmc_glmm(mcmc_model)
  
  effectiveSize(mcmc_model$Sol)
  effectiveSize(mcmc_model$VCV)
  
  autocorr.diag(mcmc_model$Sol)
  autocorr.diag(mcmc_model$VCV)
  
  summary_mcmc_model <- summary_mcmc_glmm(mcmc_model)
  
  summary(mcmc_model$Sol)
  summary(mcmc_model$VCV)
  
  # Calculate fixed effects predicted values variance
  beta_fixed <- colMeans(mcmc_model$Sol)
  X <- model.matrix(~ core_mu_median, data = merged.df)
  mFixed <- as.numeric(X %*% beta_fixed)
  mVarF <- var(mFixed)
  
  # Extract posterior mean variance components
  VCV_mean <- colMeans(mcmc_model$VCV)
  
  # Random effect variance (phylogeny)
  random_var <- VCV_mean["tip"]
  
  # Residual variance
  resid_var <- VCV_mean["units"]
  
  # Calculate R² values
  fixed_r2 <- mVarF / (mVarF + random_var + resid_var)
  random_r2 <- random_var / (mVarF + random_var + resid_var)
  conditional_r2 <- (mVarF + random_var) / (mVarF + random_var + resid_var)
  
  # Create a summary table
  R2_table <- data.frame(
    Component = c("Fixed effect", "Random effect (phylogeny)", "Conditional (total)"),
    R2 = round(c(fixed_r2, random_r2, conditional_r2), 4)
  )
  
  mcmc_summary$param <- var
  R2_table$param <- var
  
  mcmc_summary_all <- rbind(mcmc_summary_all, mcmc_summary)
  R2_table_all <- rbind(R2_table_all, R2_table)
  
  # model checks
  png(filename=paste0(outpref, var, "_MCMC_variance_components.png"))
  plot(mcmc_model$VCV)    # variance components
  dev.off()
  
  
  png(filename=paste0(outpref, var, "_MCMC_fixed_components.png"))
  plot(mcmc_model$Sol)    # fixed effects
  dev.off()
}

# write outputs
write.csv(R2_table_all, file=paste0(outpref, "R2_table.csv"))
write.csv(mcmc_summary_all, file=paste0(outpref, "mcmc_summary.csv"))
merged.df.all.x <- merge(meta, traits.df.final, by.x = "taxa", by.y = "label", all.x = TRUE)

merged.df.all.x$taxa <- merged.df.all.x$taxa %>%
  str_replace_all("_", " ") %>%        # replace underscores with spaces
  str_to_sentence()                    # capitalize first letter only

to.remove <- c("habitat_count", "pangenome_openness", "intra_species_nucleotide_diversity", "intra_species_nucleotide_diversity", "species..GTDB.", "file_path")

merged.df.all.x <- merged.df.all.x[ , -which(names(merged.df.all.x) %in% to.remove)]
colnames(merged.df.all.x) <- c("Species", "Basal gene turnover rate median", "Basal gene turnover rate LCI", "Basal gene turnover rate UCI", 
                               "Core mutation rate median", "Core mutation rate LCI", "Core mutation rate UCI", 
                               "Proportion fast genes median", "Proportion fast genes LCI", "Proportion fast genes UCI",
                               "SampleID", "Generalism score")

write.csv(merged.df.all.x, file=paste0(outpref, "408_species_data.csv"), row.names = FALSE)

# match across datasets
merged.df <- merge(meta, traits.df.final, by.x = "taxa", by.y = "label")

# plot tree just for parameters tested
rownames(merged.df) <- merged.df$tip

traits <- merged.df[, c("rate_genes1_median", "core_mu_median", "prop_genes2_median", "generalism_score")]
colnames(traits) <- c("Basal gene turnover rate", "Core mutation rate", "Proportion fast genes", "Generalism score")

# redraw tree, with all datapoints
tree <- read.tree(tree.file)

# match across datasets
merged.df.all.x <- merge(meta, traits.df.final, by.x = "taxa", by.y = "label", all.x = TRUE)

# plot tree just for parameters tested
rownames(merged.df.all.x) <- merged.df.all.x$tip

traits <- merged.df.all.x[, c("rate_genes1_median", "core_mu_median", "prop_genes2_median", "generalism_score")]
colnames(traits) <- c("Basal gene turnover rate", "Core mutation rate", "Proportion fast genes", "Generalism score")

# remove tips with no data
tips_with_data <- rownames(traits)   # or meta$tip
tree <- drop.tip(
  tree,
  setdiff(tree$tip.label, tips_with_data)
)

tree <- midpoint.root(tree)

# -------------------------------
# Reorder taxonomy to tree tips
# -------------------------------
tax <- tax[tree$tip.label]
tax <- droplevels(tax)
traits <- traits[tree$tip.label, , drop = FALSE]

generalism_vals <- traits[, 4]
na_idx <- which(!is.na(generalism_vals))

traits[is.na(traits)] <- 0

{
  png(paste0(outpref, "phylo_radial_bars_generalism_all.png"), units = "in", width = 9, height = 5.5, res = 150)
  
  layout_matrix <- matrix(c(1,2), ncol = 2, byrow = TRUE)
  layout(layout_matrix, widths = c(5, 4))
  par(mar = c(1, 1, 1, 1), oma = c(1, 1, 1, 1))
  
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
  
  ## ---- get plotting coordinates ----
  pp <- get("last_plot.phylo", envir = .PlotPhyloEnv)
  
  ntip <- length(tree$tip.label)
  
  x <- pp$xx[1:ntip]
  y <- pp$yy[1:ntip]
  
  angles <- atan2(y, x)
  tip_r  <- sqrt(x^2 + y^2)
  
  # --- fixed base radius for all bars (start radius) ---
  base_radius <- max(tip_r) + 0.1  # slightly outside tips
  
  ## --- ring geometry ---
  ring_inner <- max(tip_r) * 1.05
  ring_width <- 0.08
  ring_gap   <- 0.02
  dtheta     <- pi / ntip
  
  base_cols <- paletteer_d("ggsci::category20_d3")
  
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
  
  # update base_radius
  base_radius <- base_radius + 1 * (ring_width + ring_gap) + 0.1
  
  # --- ring parameters ---
  ring_width <- 0.15
  dtheta     <- pi / ntip
  
  # --- Viridis palettes for rings ---
  palettes <- c("A", "B", "C", "D")
  
  log <- FALSE
  ring_colours <- c("#E64B35FF", "#4DBBD5FF", "#00A087FF", "#3C5488FF")
  trait_names <- colnames(traits)  # assuming colnames are trait names
  
  for (i in 1:4) {
    
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

  # add asterisks if missing
  if (length(na_idx) > 0) {

    # radius just outside the outermost ring
    outer_ring_radius <- base_radius +
      (4 - 1) * (ring_width + ring_gap) +
      ring_width + 0.03

    text(
      x = outer_ring_radius * cos(angles[na_idx]),
      y = outer_ring_radius * sin(angles[na_idx]),
      labels = "*",
      cex = 0.6,
      col = "black"
    )
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
    title = "Taxonomy"
  )
  
  
  dev.off()
}
