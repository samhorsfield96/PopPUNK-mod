library(ggplot2)
library(dplyr)
library(stringr)
library(tidyr)
library(ggrepel)

# ── COG category functional descriptions ──────────────────────────────────────
cog_descriptions <- c(
  "A" = "A: RNA processing\n& modification",
  "B" = "B: Chromatin structure\n& dynamics",
  "C" = "C: Energy production\n& conversion",
  "D" = "D: Cell cycle control\n& division",
  "E" = "E: Amino acid transport\n& metabolism",
  "F" = "F: Nucleotide transport\n& metabolism",
  "G" = "G: Carbohydrate transport\n& metabolism",
  "H" = "H: Coenzyme transport\n& metabolism",
  "I" = "I: Lipid transport\n& metabolism",
  "J" = "J: Translation &\nribosomal biogenesis",
  "K" = "K: Transcription",
  "L" = "L: Replication,\nrecombination & repair",
  "M" = "M: Cell wall/membrane\nbiogenesis",
  "N" = "N: Cell motility",
  "O" = "O: Post-translational\nmodification & chaperones",
  "P" = "P: Inorganic ion transport\n& metabolism",
  "Q" = "Q: Secondary metabolites\nbiosynthesis & transport",
  "R" = "R: General function\nprediction only",
  "S" = "S: Function unknown",
  "T" = "T: Signal transduction\nmechanisms",
  "U" = "U: Intracellular trafficking\n& secretion",
  "V" = "V: Defense mechanisms",
  "W" = "W: Extracellular structures",
  "X" = "X: Mobilome: prophages\n& transposons"
)

assign_p_stars <- function(p_values)
{
  stars <- ifelse(p_values < 0.001, "***",
                  ifelse(p_values < 0.01, "**",
                         ifelse(p_values < 0.05, "*",
                                ifelse(p_values < 0.1, ".","ns"))))
  stars
}

# ── Input directory (first CLI arg, or current dir) ───────────────────────────
args <- commandArgs(trailingOnly = TRUE)
input_dir <- if (length(args) > 0) args[1] else "."

# ── Read all *_cog_rank_test.csv files ────────────────────────────────────────
files <- Sys.glob(file.path(input_dir, "*_cog_rank_test.csv"))
if (length(files) == 0) {
  stop("No *_cog_rank_test.csv files found in: ", input_dir)
}

all_data <- bind_rows(lapply(files, function(f) {
  species <- sub("_cog_rank_test\\.csv$", "", basename(f))
  df <- read.csv(f)
  df$Species <- species
  df
}))

# Prettify species names (underscores → spaces)
all_data$Species <- gsub("_", " ", all_data$Species)

# ── Global BH adjustment across all species and both tests ────────────────────
all_data$W_P_Adj_global  <- p.adjust(all_data$W_P_Value,  method = "BH")
all_data$KS_P_Adj_global <- p.adjust(all_data$KS_P_Value, method = "BH")
all_data$W_P_Stars_global <- assign_p_stars(all_data$W_P_Adj_global)
all_data$KS_P_Stars_global <- assign_p_stars(all_data$KS_P_Adj_global)

# ── Attach COG labels and set factor order (alphabetical by letter) ───────────
all_data$COG_label <- cog_descriptions[all_data$COG_category]
# Fallback: categories not in lookup table keep their raw letter
missing <- is.na(all_data$COG_label)
all_data$COG_label[missing] <- all_data$COG_category[missing]

present_cats  <- intersect(names(cog_descriptions), unique(all_data$COG_category))
label_levels  <- cog_descriptions[present_cats]
all_data$COG_label <- factor(all_data$COG_label, levels = label_levels)

# ── Reshape to long format (one row per test per species per category) ─────────
plot_data <- all_data %>%
  select(COG_category, COG_label, Species, Median_With, Median_Without,
         W_P_Adj_global, KS_P_Adj_global) %>%
  pivot_longer(
    cols      = c(W_P_Adj_global, KS_P_Adj_global),
    names_to  = "Test",
    values_to = "P_Adj"
  ) %>%
  mutate(
    Test        = recode(Test,
                         "W_P_Adj_global"  = "Wilcoxon rank-sum",
                         "KS_P_Adj_global" = "Kolmogorov-Smirnov"),
    # log2 fold-change: scale-free, symmetric, 0 = no effect
    effect_size = log2(Median_With / Median_Without),
    neg_log10_p = -log10(P_Adj),
    sig_colour  = ifelse(P_Adj < 0.05, Species, "ns")
  )

# ── Plot helper ───────────────────────────────────────────────────────────────
sig_line  <- -log10(0.05)
n_species <- length(unique(plot_data$Species))
n_cog     <- length(levels(plot_data$COG_label))

# Use a palette that scales with the number of species
species_names <- unique(plot_data$Species)
raw_colours <- if (n_species <= 8) {
  RColorBrewer::brewer.pal(max(3, n_species), "Set1")[seq_len(n_species)]
} else {
  scales::hue_pal()(n_species)
}
species_colours <- setNames(raw_colours, species_names)

make_plot <- function(data, test_label, x_label) {
  ggplot(data, aes(x = effect_size, y = neg_log10_p,
                   colour = sig_colour, label = COG_category)) +
    geom_point(size = 2.5) +
    geom_text_repel(data  = filter(data, P_Adj < 0.05),
                    size  = 3,
                    show.legend = FALSE,
                    max.overlaps = Inf) +
    geom_hline(yintercept = sig_line,
               linetype   = "dashed",
               colour     = "black",
               linewidth  = 0.6) +
    geom_vline(xintercept = 0,
               linetype   = "solid",
               colour     = "grey40",
               linewidth  = 0.4) +
    annotate("text",
             x      = Inf,
             y      = sig_line,
             hjust  = 1,
             vjust  = -0.4,
             label  = "p = 0.05",
             size   = 3,
             colour = "black") +
    facet_grid(Species ~ ., scales = "free") +
    scale_colour_manual(
      values = c(species_colours, "ns" = "grey70"),
      breaks = names(species_colours),
      labels = names(species_colours)
    ) +
    labs(
      x       = x_label,
      y       = expression(-log[10](adjusted~italic(p)~value)),
      colour  = "Species",
    ) +
    theme_bw(base_size = 12) +
    theme(
      legend.text      = element_text(face = "italic"),
      legend.position  = "none",
      panel.grid.minor = element_blank(),
      strip.background = element_rect(fill = "grey90"),
      strip.text.y     = element_text(face = "bold.italic"),
      plot.caption     = element_text(hjust = 0, size = 8, colour = "grey40")
    )
}

p_wilcox <- make_plot(filter(plot_data, Test == "Wilcoxon rank-sum"),
                      "Wilcoxon rank-sum test",
                      expression(log[2](Median[with] / Median[without])))
p_ks     <- make_plot(filter(plot_data, Test == "Kolmogorov-Smirnov"),
                      "Kolmogorov-Smirnov test",
                      expression(log[2](Median[with] / Median[without])))
# ── Save ───────────────────────────────────────────────────────────────────────
plot_height <- 4 + n_species * 1.5   # scale height with number of species

for (pair in list(
  list(p = p_wilcox, stem = "cog_rank_test_wilcoxon"),
  list(p = p_ks,     stem = "cog_rank_test_ks")
)) {
  ggsave(file.path(input_dir, paste0(pair$stem, ".pdf")),
         plot = pair$p, width = 8, height = plot_height, device = "pdf")
  ggsave(file.path(input_dir, paste0(pair$stem, ".png")),
         plot = pair$p, width = 8, height = plot_height, dpi = 300)
  message("Saved: ", pair$stem)
}
