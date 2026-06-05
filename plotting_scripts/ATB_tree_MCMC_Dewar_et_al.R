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

generate_tree <- function(tree, dataframe){
  # remove tips with no data
  tips_with_data <- dataframe$tip   # or meta$tip
  new.tree <- drop.tip(
    tree,
    setdiff(tree$tip.label, tips_with_data)
  )
  
  # --- midpoint root the tree ---
  new.tree <- midpoint.root(new.tree)
  
  # -------------------------------
  # Reorder taxonomy to tree tips
  # -------------------------------
  new.dataframe <- dataframe[new.tree$tip.label, , drop = FALSE]
  
  # -------------------------------
  # HARD SAFETY CHECKS
  # -------------------------------
  stopifnot(
    nrow(new.dataframe) == length(new.tree$tip.label),
    all(rownames(new.dataframe) == new.tree$tip.label)
  )
  
  list(new.tree, new.dataframe)
}

indir <- "/Users/samhorsfield/Software/PopPUNK-mod/Dewar_et_al_analysis/"

tree.file <- paste0(indir, "reps_ppmod_tree.treefile")
summary.file <- paste0(indir, "ppmod_summary.csv")
index.file <- paste0(indir, "sampled_alignments.tsv")
gtdb.file <- paste0(indir, "gtdbtk.bac120.summary.tsv")
lifestyle.file <- paste0(indir, "pangenome_lifestyles.csv")
outpref <- paste0(indir, "/Dewar_analysis")

summary.df <- read.csv(summary.file, row.names = 1)

values.NA <- summary.df[is.na(summary.df$order),]

index.df <- read.csv(index.file, sep="\t", header = FALSE)
colnames(index.df) <- c("taxa", "tip", "file_path")

meta <- merge(summary.df, index.df, by = "taxa")

tree <- read.tree(tree.file)

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

# add Dewar et al. lifestyle analysis
pangenome_lifestyles <- read.csv(lifestyle.file)
## Order lifestyle traits in order of 'least' to 'most' variable 
pangenome_lifestyles$Host_or_free <- factor(pangenome_lifestyles$Host_or_free, levels=c("Host","Both","Free","Unknown"))
pangenome_lifestyles$Obligate_facultative <- factor(pangenome_lifestyles$Obligate_facultative, levels=c("Obligate","Facultative","Unknown",""))
pangenome_lifestyles$Intra_or_extracellular <- factor(pangenome_lifestyles$Intra_or_extracellular, levels=c("Intracellular","Both","Extracellular", "Unknown", ""))
pangenome_lifestyles$Effect_on_host <- factor(pangenome_lifestyles$Effect_on_host, levels=c("Pathogen","Both","Mutualist","Unknown", ""))
pangenome_lifestyles$Motility <- factor(pangenome_lifestyles$Motility, levels=c("Non-motile","Both","Motile"))
pangenome_lifestyles$Category_primary_env <- factor(pangenome_lifestyles$Category_primary_env, levels=c("Host","Free","Unknown"))

## Make data set a 'data.frame' (to avoid any errors from data format)
pangenome_lifestyles$pangenome_fluidity <- pangenome_lifestyles$genome_fluidity

pangenome_lifestyles <- as.data.frame(pangenome_lifestyles)

# Normalize species names: replace spaces with underscores, strip trailing GTDB letter suffix (e.g. " maltophilia_A" -> "_maltophilia")
meta$Species_norm <- gsub("_[A-Z]+$", "", gsub(" ", "_", meta$Species))
merged.df <- merge(meta, pangenome_lifestyles, by.x = "Species_norm", by.y = "Species")
rownames(merged.df) <- merged.df$tip

# Host_or_free vs. fitted params
pangenome_lifestyles_host_free_no_unknown <- merged.df %>%
  filter(Host_or_free != "Unknown")

pangenome_lifestyles_host_free_no_unknown <- as.data.frame(pangenome_lifestyles_host_free_no_unknown)

tree.df.list <- generate_tree(tree, pangenome_lifestyles_host_free_no_unknown)
new.tree <- tree.df.list[[1]]
new.df <- tree.df.list[[2]]

Ainv <- inverseA(new.tree, nodes = "TIPS", scale = FALSE)$Ainv

prior <- list(
  G = list(
    G1 = list(V = 1, nu = 1)
  ),
  R = list(
    R1 = list(V = 1, nu = 1)
  )
)

prior <- list(R=list(V = 1, nu = 0.002), G=list(G1=list(V=1,nu=0.002))) 

# --------------------------------------------------------------------------
# Generic MCMCglmm runner with phylogenetic correction
# --------------------------------------------------------------------------
run_mcmc_model <- function(y_var, x_vars, data, tree, prior, nitt = 50000) {
  # Drop rows where any categorical x variable is "Unknown" or ""
  filter_mask <- rep(TRUE, nrow(data))
  for (xv in x_vars) {
    if (xv %in% colnames(data) && is.factor(data[[xv]])) {
      filter_mask <- filter_mask & !(as.character(data[[xv]]) %in% c("Unknown", ""))
    }
  }
  sub_data <- droplevels(data[filter_mask, , drop = FALSE])

  # Prune tree to match filtered data and align row order
  tree_list <- generate_tree(tree, sub_data)
  sub_tree  <- tree_list[[1]]
  sub_df    <- tree_list[[2]]

  Ainv_sub <- inverseA(sub_tree, nodes = "TIPS", scale = FALSE)$Ainv

  model_formula <- as.formula(
    paste(y_var, "~", paste(x_vars, collapse = " + "))
  )

  model <- MCMCglmm(model_formula,
                    random   = ~tip,
                    data     = sub_df,
                    nitt     = nitt,
                    prior    = prior,
                    ginverse = list(tip = Ainv_sub),
                    verbose  = FALSE)

  # R2 (Nakagawa & Schielzeth style)
  mFixed    <- as.vector(model$X %*% apply(model$Sol, 2, mean))
  mVarF     <- var(mFixed)
  vcv_means <- apply(model$VCV, 2, mean)
  rand_var  <- sum(vcv_means[names(vcv_means) != "units"])

  fixed_r2       <- mVarF / (mVarF + sum(vcv_means))
  random_r2      <- rand_var / (sum(vcv_means) + mVarF)
  conditional_r2 <- (rand_var + mVarF) / (sum(vcv_means) + mVarF)

  R2_table <- matrix(c(fixed_r2, random_r2, conditional_r2), ncol = 1)
  R2_table <- format(R2_table, digits = 4)
  rownames(R2_table) <- c("Fixed effect", "Random effect", "Total model")
  colnames(R2_table) <- c("R-squared value")

  list(model = model, summary = summary_mcmc_glmm(model), R2 = R2_table)
}

# --------------------------------------------------------------------------
# Y and X variable definitions
# --------------------------------------------------------------------------
y_vars <- c(
  "rate_genes1_median",  "rate_genes1_lower_CI",  "rate_genes1_upper_CI",
  "core_mu_median",      "core_mu_lower_CI",      "core_mu_upper_CI",
  "prop_genes2_median",  "prop_genes2_lower_CI",  "prop_genes2_upper_CI"
)

x_var_sets <- list(
  Host_or_free           = c("Host_or_free"),
  Category_primary_env   = c("Category_primary_env"),
  Obligate_facultative   = c("Obligate_facultative"),
  Intra_or_extracellular = c("Intra_or_extracellular"),
  Effect_on_host         = c("Effect_on_host"),
  Motility               = c("Motility"),
  combined               = c("Obligate_facultative", "Intra_or_extracellular",
                              "Effect_on_host", "Motility")
)

# --------------------------------------------------------------------------
# Run all y ~ x combinations
# --------------------------------------------------------------------------
all_results <- list()
for (y in y_vars) {
  all_results[[y]] <- list()
  for (x_name in names(x_var_sets)) {
    message(sprintf("Fitting: %s ~ %s", y,
                    paste(x_var_sets[[x_name]], collapse = " + ")))
    all_results[[y]][[x_name]] <- run_mcmc_model(
      y_var  = y,
      x_vars = x_var_sets[[x_name]],
      data   = merged.df,
      tree   = tree,
      prior  = prior
    )
  }
}

# --------------------------------------------------------------------------
# Helper functions to extract results
# --------------------------------------------------------------------------
get_all_summaries <- function(results) {
  rows <- list()
  for (y in names(results)) {
    for (x in names(results[[y]])) {
      summ          <- results[[y]][[x]]$summary
      summ$y_var    <- y
      summ$x_model  <- x
      summ$term     <- rownames(summ)
      rownames(summ) <- NULL
      rows[[length(rows) + 1]] <- summ
    }
  }
  do.call(rbind, rows)
}

get_all_r2 <- function(results) {
  rows <- list()
  for (y in names(results)) {
    for (x in names(results[[y]])) {
      r2_mat <- results[[y]][[x]]$R2
      rows[[length(rows) + 1]] <- data.frame(
        y_var   = y,
        x_model = x,
        metric  = rownames(r2_mat),
        value   = as.numeric(r2_mat[, 1]),
        stringsAsFactors = FALSE
      )
    }
  }
  do.call(rbind, rows)
}

all_summaries <- get_all_summaries(all_results)
all_r2        <- get_all_r2(all_results)
