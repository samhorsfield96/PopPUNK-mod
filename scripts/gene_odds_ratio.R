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
              help = "BH-adjusted p-value cutoff [default: %default]")
)

opt <- parse_args(OptionParser(option_list = option_list))

if (is.null(opt$infile)) {
  stop("--infile is required. Run with --help for usage.", call. = FALSE)
}

infile   <- opt$infile
min.freq <- opt[["min-freq"]]
max.freq <- opt[["max-freq"]]
p.cutoff <- opt[["p-cutoff"]]
outfile  <- if (!is.null(opt$outfile)) opt$outfile else
              sub("(\\.tsv|\\.csv)?$", "_go_enrichment.csv",
                  tools::file_path_sans_ext(infile), ignore.case = TRUE)

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
    a <- sum(has_go &  in_group)
    b <- sum(has_go & !in_group)
    c <- sum(!has_go &  in_group)
    d <- sum(!has_go & !in_group)

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

# Example usage: run enrichment for each Significance group
go_hi <- go_enrichment(df, group = "Hi",  min.freq = min.freq, max.freq = max.freq, p.adj.cutoff = p.cutoff)
go_lo <- go_enrichment(df, group = "Lo",  min.freq = min.freq, max.freq = max.freq, p.adj.cutoff = p.cutoff)
go_ns <- go_enrichment(df, group = "NS",  min.freq = min.freq, max.freq = max.freq, p.adj.cutoff = p.cutoff)

total.go.df <- rbind(go_lo, go_ns, go_hi)

write.csv(total.go.df, outfile, row.names = FALSE)
message("Written: ", outfile)
