suppressPackageStartupMessages(library(optparse))

option_list <- list(
  make_option(c("-i", "--infile"),    type = "character", default = NULL,
              help = "Path to odds ratios TSV file [required]"),
  make_option(c("-o", "--outpref"),   type = "character", default = NULL,
              help = "Output prefix [default: derived from --infile]"),
  make_option(c("--min-freq"),        type = "double",    default = 0.0,
              help = "Minimum gene frequency [default: %default]"),
  make_option(c("--max-freq"),        type = "double",    default = 1.0,
              help = "Maximum gene frequency [default: %default]"),
  make_option(c("-p", "--p-cutoff"),  type = "double",    default = 0.05,
              help = "BH-adjusted p-value cutoff [default: %default]"),
  make_option(c("--product"),          type = "character", default = NULL,
              help = "Comma-separated regex pattern(s) to match against the Product field for enrichment analysis [optional]"),
  make_option(c("--eggnog"),           type = "character", default = NULL,
              help = "Path to eggNOG-mapper annotations file (.emapper.annotations); merged on Gene_Name = #query [optional]")
)

opt <- parse_args(OptionParser(option_list = option_list))

if (is.null(opt$infile)) {
  stop("--infile is required. Run with --help for usage.", call. = FALSE)
}

infile   <- opt$infile
min.freq <- opt[["min-freq"]]
max.freq <- opt[["max-freq"]]
p.cutoff         <- opt[["p-cutoff"]]
product.patterns <- if (!is.null(opt$product)) trimws(strsplit(opt$product, ",")[[1]]) else NULL
eggnog.file      <- opt$eggnog
outpref <- if (!is.null(opt$outpref)) opt$outpref else
              tools::file_path_sans_ext(infile)

outfile           <- paste0(outpref, "_go_enrichment.csv")
product.outfile   <- paste0(outpref, "_product_enrichment.csv")
go.rank.outfile   <- paste0(outpref, "_go_rank_test.csv")
prod.rank.outfile <- paste0(outpref, "_product_rank_test.csv")
cog.outfile       <- paste0(outpref, "_cog_enrichment.csv")
cog.rank.outfile  <- paste0(outpref, "_cog_rank_test.csv")

df <- read.csv(infile, sep = "\t")

if (!is.null(eggnog.file)) {
  egg <- read.csv(eggnog.file, sep = "\t", skip = 5,
                  header = FALSE, stringsAsFactors = FALSE)
  # Re-read just the header line (first non-comment line)
  hdr_line <- grep("^#query", readLines(eggnog.file, n = 50), value = TRUE)[1]
  egg_cols  <- strsplit(sub("^#", "", hdr_line), "\t")[[1]]
  colnames(egg) <- egg_cols
  egg <- egg[egg$query != "", ]
  # Keep only join key and annotation columns
  egg_keep <- intersect(c("query", "GOs", "Description", "COG_category",
                          "Preferred_name", "EC", "KEGG_ko",
                          "KEGG_Pathway", "PFAMs"), colnames(egg))
  egg <- egg[, egg_keep, drop = FALSE]
  colnames(egg)[colnames(egg) == "query"] <- "Gene_Name"

  df <- merge(df, egg, by = "Gene_Name", all.x = TRUE, suffixes = c("", ".eggnog"))

  # Backfill GO and Product from eggnog where the original values are absent
  if ("GOs" %in% colnames(df)) {
    missing_go <- is.na(df$GO) | df$GO == "" | df$GO == "-"
    df$GO[missing_go] <- df$GOs[missing_go]
    df$GOs <- NULL
  }
  if ("Description" %in% colnames(df)) {
    missing_prod <- is.na(df$Product) | df$Product == "" | df$Product == "-"
    df$Product[missing_prod] <- df$Description[missing_prod]
    df$Description <- NULL
  }
}
# subset.df <- subset(df, Gene_Frequency > min.freq & Gene_Frequency <= max.freq & Significance != "NS")
# 
# df.lo <- subset(subset.df[order(subset.df$Odds_Ratio, decreasing = TRUE),], Significance == "Lo")
# df.hi <- subset(subset.df[order(subset.df$Odds_Ratio, decreasing = FALSE),], Significance == "Hi")

assign_p_stars <- function(p_values)
{
  stars <- ifelse(p_values < 0.001, "***",
                  ifelse(p_values < 0.01, "**",
                         ifelse(p_values < 0.05, "*",
                                ifelse(p_values < 0.1, ".","ns"))))
  stars
}

# GO term enrichment analysis
# Tests whether GO terms are enriched in a given Significance group (Hi/Lo/NS)
# vs all other genes using Fisher's exact test with BH FDR correction.
# Args:
#   df          : data frame from the odds ratios TSV (full, unfiltered)
#   group       : one of "Hi", "Lo", or "NS"
#   min.freq    : minimum gene frequency threshold (default 0.0)
#   max.freq    : maximum gene frequency threshold (default 0.95)
#   p.adj.cutoff: adjusted p-value significance cutoff (default 0.05)
# Returns a data frame of enriched GO terms sorted by adjusted p-value.
go_enrichment <- function(df, group, min.freq = 0.0, max.freq = 0.95, p.adj.cutoff = 0.05) {
  df <- subset(df, Gene_Frequency > min.freq & Gene_Frequency <= max.freq)

  go_list <- strsplit(as.character(df$GO), ",")
  go_list <- lapply(go_list, function(x) trimws(x[nchar(trimws(x)) > 0]))

  in_group <- df$Significance == group

  all_go <- unique(unlist(go_list))
  all_go <- all_go[nchar(all_go) > 0]

  if (length(all_go) == 0) {
    message("No GO terms found in data.")
    return(NULL)
  }

  results <- lapply(all_go, function(go_term) {
    has_go <- vapply(go_list, function(x) go_term %in% x, logical(1))

    # 2x2 contingency table:
    #            in_group  other
    # has_go        a        b
    # no_go         c        d
    a <- sum(has_go &  in_group) + 0.5
    b <- sum(has_go & !in_group) + 0.5
    c <- sum(!has_go &  in_group) + 0.5
    d <- sum(!has_go & !in_group) + 0.5

    ft <- fisher.test(matrix(c(a, c, b, d), nrow = 2))

    data.frame(
      GO          = go_term,
      In_Group    = a,
      Other       = b,
      Odds_Ratio  = as.numeric(ft$estimate),
      P_Value     = ft$p.value,
      stringsAsFactors = FALSE
    )
  })

  results_df <- do.call(rbind, results)
  results_df$P_Adj <- p.adjust(results_df$P_Value, method = "BH")
  results_df$Group <- group
  results_df <- results_df[order(results_df$P_Adj, results_df$P_Value), ]
  rownames(results_df) <- NULL
  results_df$P_Stars <- assign_p_stars(results_df$P_Adj)
  results_df
}

# Product text enrichment analysis
# Tests whether genes whose Product field matches a regex pattern are enriched
# in each Significance group (Hi/Lo/NS) vs all other genes, using Fisher's
# exact test with BH FDR correction across all pattern x group combinations.
# Args:
#   df          : data frame from the odds ratios TSV (full, unfiltered)
#   patterns    : character vector of regex patterns to match against Product
#   min.freq    : minimum gene frequency threshold
#   max.freq    : maximum gene frequency threshold
#   p.adj.cutoff: adjusted p-value significance cutoff
# Returns a data frame of significant pattern x group combinations.
product_enrichment <- function(df, patterns, min.freq = 0.0, max.freq = 1.0, p.adj.cutoff = 0.05) {
  df <- subset(df, Gene_Frequency > min.freq & Gene_Frequency <= max.freq)

  groups <- sort(unique(df$Significance))

  results <- do.call(rbind, lapply(patterns, function(pat) {
    matches <- grepl(pat, df$Product, ignore.case = TRUE, perl = TRUE)

    do.call(rbind, lapply(groups, function(group) {
      in_group <- df$Significance == group

      a <- sum( matches &  in_group) + 0.5
      b <- sum( matches & !in_group) + 0.5
      c <- sum(!matches &  in_group) + 0.5
      d <- sum(!matches & !in_group) + 0.5

      ft <- fisher.test(matrix(c(a, c, b, d), nrow = 2))

      data.frame(
        Pattern    = pat,
        Group      = group,
        In_Group   = a,
        Other      = b,
        Odds_Ratio = as.numeric(ft$estimate),
        P_Value    = ft$p.value,
        stringsAsFactors = FALSE
      )
    }))
  }))

  results$P_Adj <- p.adjust(results$P_Value, method = "BH")
  results <- results[order(results$P_Adj, results$P_Value), ]
  rownames(results) <- NULL
  results$P_Stars <- assign_p_stars(results$P_Adj)
  results
}

# GO term rank test
# Wilcoxon rank-sum test: are Odds_Ratio values higher/lower for genes with a
# given GO term vs genes without? More sensitive than Fisher when signal is
# spread across the full odds-ratio continuum rather than just Hi/Lo/NS labels.
go_rank_test <- function(df, min.freq = 0.0, max.freq = 1.0, p.adj.cutoff = 0.05) {
  df <- subset(df, Gene_Frequency > min.freq & Gene_Frequency <= max.freq)

  go_list <- strsplit(as.character(df$GO), ",")
  go_list <- lapply(go_list, function(x) trimws(x[nchar(trimws(x)) > 0]))

  all_go <- unique(unlist(go_list))
  all_go <- all_go[nchar(all_go) > 0]

  if (length(all_go) == 0) {
    message("No GO terms found in data.")
    return(NULL)
  }

  results <- lapply(all_go, function(go_term) {
    has_go   <- vapply(go_list, function(x) go_term %in% x, logical(1))
    vals_in  <- df$Odds_Ratio[has_go]
    vals_out <- df$Odds_Ratio[!has_go]
    if (length(vals_in) < 2 || length(vals_out) < 2) return(NULL)
    wt <- wilcox.test(vals_in, vals_out, exact = FALSE)
    kt <- suppressWarnings(ks.test(vals_in, vals_out))
    data.frame(
      GO             = go_term,
      N_With         = length(vals_in),
      N_Without      = length(vals_out),
      Median_With    = median(vals_in),
      Median_Without = median(vals_out),
      W              = as.numeric(wt$statistic),
      W_P_Value      = wt$p.value,
      KS_D           = as.numeric(kt$statistic),
      KS_P_Value     = kt$p.value,
      stringsAsFactors = FALSE
    )
  })

  results <- Filter(Negate(is.null), results)
  if (length(results) == 0) return(NULL)
  results_df <- do.call(rbind, results)
  results_df$W_P_Adj  <- p.adjust(results_df$W_P_Value,  method = "BH")
  results_df$KS_P_Adj <- p.adjust(results_df$KS_P_Value, method = "BH")
  results_df <- results_df[order(results_df$W_P_Adj, results_df$W_P_Value), ]
  rownames(results_df) <- NULL
  results_df$W_P_Stars  <- assign_p_stars(results_df$W_P_Adj)
  results_df$KS_P_Stars <- assign_p_stars(results_df$KS_P_Adj)
  results_df
}

# Product term rank test
# Wilcoxon rank-sum test: are Odds_Ratio values higher/lower for genes matching
# each product pattern vs non-matching genes?
product_rank_test <- function(df, patterns, min.freq = 0.0, max.freq = 1.0, p.adj.cutoff = 0.05) {
  df <- subset(df, Gene_Frequency > min.freq & Gene_Frequency <= max.freq)

  results <- lapply(patterns, function(pat) {
    matches  <- grepl(pat, df$Product, ignore.case = TRUE, perl = TRUE)
    vals_in  <- df$Odds_Ratio[matches]
    vals_out <- df$Odds_Ratio[!matches]
    if (length(vals_in) < 2 || length(vals_out) < 2) return(NULL)
    wt <- wilcox.test(vals_in, vals_out, exact = FALSE)
    kt <- suppressWarnings(ks.test(vals_in, vals_out))
    data.frame(
      Pattern        = pat,
      N_With         = length(vals_in),
      N_Without      = length(vals_out),
      Median_With    = median(vals_in),
      Median_Without = median(vals_out),
      W              = as.numeric(wt$statistic),
      W_P_Value      = wt$p.value,
      KS_D           = as.numeric(kt$statistic),
      KS_P_Value     = kt$p.value,
      stringsAsFactors = FALSE
    )
  })

  results <- Filter(Negate(is.null), results)
  if (length(results) == 0) return(NULL)
  results_df <- do.call(rbind, results)
  results_df$W_P_Adj  <- p.adjust(results_df$W_P_Value,  method = "BH")
  results_df$KS_P_Adj <- p.adjust(results_df$KS_P_Value, method = "BH")
  results_df <- results_df[order(results_df$W_P_Adj, results_df$W_P_Value), ]
  rownames(results_df) <- NULL
  results_df$W_P_Stars  <- assign_p_stars(results_df$W_P_Adj)
  results_df$KS_P_Stars <- assign_p_stars(results_df$KS_P_Adj)
  results_df
}

# COG category enrichment analysis
# Each gene may carry multiple COG letters (e.g. "EGP"); each letter is tested
# independently against each Significance group using Fisher's exact test with
# BH FDR correction across all letter x group combinations.
cog_enrichment <- function(df, min.freq = 0.0, max.freq = 1.0, p.adj.cutoff = 0.05) {
  if (!"COG_category" %in% colnames(df)) {
    message("COG_category column not found; skipping COG enrichment.")
    return(NULL)
  }
  df <- subset(df, Gene_Frequency > min.freq & Gene_Frequency <= max.freq)

  cog_list <- lapply(as.character(df$COG_category), function(x) {
    x <- trimws(x)
    if (is.na(x) || x == "" || x == "-") return(character(0))
    strsplit(x, "")[[1]]
  })

  groups   <- sort(unique(df$Significance))
  all_cats <- sort(unique(unlist(cog_list)))
  all_cats <- all_cats[nchar(all_cats) > 0]

  if (length(all_cats) == 0) {
    message("No COG categories found in data.")
    return(NULL)
  }

  results <- do.call(rbind, lapply(all_cats, function(cat) {
    has_cat <- vapply(cog_list, function(x) cat %in% x, logical(1))
    do.call(rbind, lapply(groups, function(group) {
      in_group <- df$Significance == group
      a <- sum( has_cat &  in_group) + 0.5
      b <- sum( has_cat & !in_group) + 0.5
      c <- sum(!has_cat &  in_group) + 0.5
      d <- sum(!has_cat & !in_group) + 0.5
      ft <- fisher.test(matrix(c(a, c, b, d), nrow = 2))
      data.frame(
        COG_category = cat,
        Group        = group,
        In_Group     = a,
        Other        = b,
        Odds_Ratio   = as.numeric(ft$estimate),
        P_Value      = ft$p.value,
        stringsAsFactors = FALSE
      )
    }))
  }))

  results$P_Adj <- p.adjust(results$P_Value, method = "BH")
  results <- results[order(results$P_Adj, results$P_Value), ]
  rownames(results) <- NULL
  results$P_Stars <- assign_p_stars(results$P_Adj)
  results
}

# COG category rank test
# Wilcoxon rank-sum test: are Odds_Ratio values shifted for genes carrying each
# COG category letter vs genes without it?
cog_rank_test <- function(df, min.freq = 0.0, max.freq = 1.0, p.adj.cutoff = 0.05) {
  if (!"COG_category" %in% colnames(df)) {
    message("COG_category column not found; skipping COG rank test.")
    return(NULL)
  }
  df <- subset(df, Gene_Frequency > min.freq & Gene_Frequency <= max.freq)

  cog_list <- lapply(as.character(df$COG_category), function(x) {
    x <- trimws(x)
    if (is.na(x) || x == "" || x == "-") return(character(0))
    strsplit(x, "")[[1]]
  })

  all_cats <- sort(unique(unlist(cog_list)))
  all_cats <- all_cats[nchar(all_cats) > 0]

  if (length(all_cats) == 0) {
    message("No COG categories found in data.")
    return(NULL)
  }

  results <- lapply(all_cats, function(cat) {
    has_cat  <- vapply(cog_list, function(x) cat %in% x, logical(1))
    vals_in  <- df$Odds_Ratio[has_cat]
    vals_out <- df$Odds_Ratio[!has_cat]
    if (length(vals_in) < 2 || length(vals_out) < 2) return(NULL)
    wt <- wilcox.test(vals_in, vals_out, exact = FALSE)
    kt <- suppressWarnings(ks.test(vals_in, vals_out))
    data.frame(
      COG_category   = cat,
      N_With         = length(vals_in),
      N_Without      = length(vals_out),
      Median_With    = median(vals_in),
      Median_Without = median(vals_out),
      W              = as.numeric(wt$statistic),
      W_P_Value      = wt$p.value,
      KS_D           = as.numeric(kt$statistic),
      KS_P_Value     = kt$p.value,
      stringsAsFactors = FALSE
    )
  })

  results <- Filter(Negate(is.null), results)
  if (length(results) == 0) return(NULL)
  results_df <- do.call(rbind, results)
  results_df$W_P_Adj  <- p.adjust(results_df$W_P_Value,  method = "BH")
  results_df$KS_P_Adj <- p.adjust(results_df$KS_P_Value, method = "BH")
  results_df <- results_df[order(results_df$W_P_Adj, results_df$W_P_Value), ]
  rownames(results_df) <- NULL
  results_df$W_P_Stars  <- assign_p_stars(results_df$W_P_Adj)
  results_df$KS_P_Stars <- assign_p_stars(results_df$KS_P_Adj)
  results_df <- results_df[order(results_df$COG_category),]
  results_df
}

# GO term significance only if nothing else specified
if (is.null(product.patterns) & is.null(eggnog.file)) {
  total.go.df <- rbind(go_lo, go_ns, go_hi)
  go_hi <- go_enrichment(df, group = "Hi",  min.freq = min.freq, max.freq = max.freq, p.adj.cutoff = p.cutoff)
  go_lo <- go_enrichment(df, group = "Lo",  min.freq = min.freq, max.freq = max.freq, p.adj.cutoff = p.cutoff)
  go_ns <- go_enrichment(df, group = "NS",  min.freq = min.freq, max.freq = max.freq, p.adj.cutoff = p.cutoff)
  # Rank tests on continuous Odds_Ratio values
  go.rank.df <- go_rank_test(df, min.freq = min.freq, max.freq = max.freq, p.adj.cutoff = p.cutoff)
  if (!is.null(go.rank.df) && nrow(go.rank.df) > 0) {
    write.csv(go.rank.df, go.rank.outfile, row.names = FALSE)
    message("Written: ", go.rank.outfile)
  }
  
  write.csv(total.go.df, outfile, row.names = FALSE)
  message("Written: ", outfile)
}

if (!is.null(product.patterns)) {
  prod.df <- product_enrichment(df, product.patterns, min.freq = min.freq, max.freq = max.freq, p.adj.cutoff = p.cutoff)
  write.csv(prod.df, product.outfile, row.names = FALSE)
  message("Written: ", product.outfile)
}

if (!is.null(product.patterns)) {
  prod.rank.df <- product_rank_test(df, product.patterns, min.freq = min.freq, max.freq = max.freq, p.adj.cutoff = p.cutoff)
  if (!is.null(prod.rank.df) && nrow(prod.rank.df) > 0) {
    write.csv(prod.rank.df, prod.rank.outfile, row.names = FALSE)
    message("Written: ", prod.rank.outfile)
  }
}

if (!is.null(eggnog.file)) {
  cog.df <- cog_enrichment(df, min.freq = min.freq, max.freq = max.freq, p.adj.cutoff = p.cutoff)
  if (!is.null(cog.df) && nrow(cog.df) > 0) {
    write.csv(cog.df, cog.outfile, row.names = FALSE)
    message("Written: ", cog.outfile)
  }

  cog.rank.df <- cog_rank_test(df, min.freq = min.freq, max.freq = max.freq, p.adj.cutoff = p.cutoff)
  if (!is.null(cog.rank.df) && nrow(cog.rank.df) > 0) {
    write.csv(cog.rank.df, cog.rank.outfile, row.names = FALSE)
    message("Written: ", cog.rank.outfile)
  }
}
