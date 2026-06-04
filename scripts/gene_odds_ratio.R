suppressPackageStartupMessages(library(optparse))

option_list <- list(
  make_option(c("-i", "--infile"),    type = "character", default = NULL,
              help = "Path to odds ratios TSV file [required]"),
  make_option(c("-o", "--outfile"),   type = "character", default = NULL,
              help = "Output CSV path [default: <infile_basename>_go_enrichment.csv]"),
  make_option(c("--min-freq"),        type = "double",    default = 0.0,
              help = "Minimum gene frequency [default: %default]"),
  make_option(c("--max-freq"),        type = "double",    default = 1.0,
              help = "Maximum gene frequency [default: %default]"),
  make_option(c("-p", "--p-cutoff"),  type = "double",    default = 0.05,
              help = "BH-adjusted p-value cutoff [default: %default]"),
  make_option(c("--product"),          type = "character", default = NULL,
              help = "Comma-separated regex pattern(s) to match against the Product field for enrichment analysis [optional]")
)

opt <- parse_args(OptionParser(option_list = option_list))

if (is.null(opt$infile)) {
  stop("--infile is required. Run with --help for usage.", call. = FALSE)
}

infile   <- opt$infile
min.freq <- opt[["min-freq"]]
max.freq <- opt[["max-freq"]]
p.cutoff        <- opt[["p-cutoff"]]
product.patterns <- if (!is.null(opt$product)) trimws(strsplit(opt$product, ",")[[1]]) else NULL
outfile  <- if (!is.null(opt$outfile)) opt$outfile else
              sub("(\\.tsv|\\.csv)?$", "_go_enrichment.csv",
                  tools::file_path_sans_ext(infile), ignore.case = TRUE)
product.outfile <- sub("_go_enrichment\\.csv$", "_product_enrichment.csv", outfile)

df <- read.csv(infile, sep = "\t")
# subset.df <- subset(df, Gene_Frequency > min.freq & Gene_Frequency <= max.freq & Significance != "NS")
# 
# df.lo <- subset(subset.df[order(subset.df$Odds_Ratio, decreasing = TRUE),], Significance == "Lo")
# df.hi <- subset(subset.df[order(subset.df$Odds_Ratio, decreasing = FALSE),], Significance == "Hi")

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

  subset(results_df, P_Adj <= p.adj.cutoff)
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

  subset(results, P_Adj <= p.adj.cutoff)
}

# Example usage: run enrichment for each Significance group
go_hi <- go_enrichment(df, group = "Hi",  min.freq = min.freq, max.freq = max.freq, p.adj.cutoff = p.cutoff)
go_lo <- go_enrichment(df, group = "Lo",  min.freq = min.freq, max.freq = max.freq, p.adj.cutoff = p.cutoff)
go_ns <- go_enrichment(df, group = "NS",  min.freq = min.freq, max.freq = max.freq, p.adj.cutoff = p.cutoff)

total.go.df <- rbind(go_lo, go_ns, go_hi)

write.csv(total.go.df, outfile, row.names = FALSE)
message("Written: ", outfile)

if (!is.null(product.patterns)) {
  prod.df <- product_enrichment(df, product.patterns, min.freq = min.freq, max.freq = max.freq, p.adj.cutoff = p.cutoff)
  write.csv(prod.df, product.outfile, row.names = FALSE)
  message("Written: ", product.outfile)
}
