#install.packages("ggmagnify", repos = c("https://hughjonesd.r-universe.dev", 
                                        #"https://cloud.r-project.org"))
library(ggplot2)
library(stringr)
library(dplyr)
library(readxl)
library(ggpubr)
library(reticulate)
library(ggmagnify)

path_to_python="/Users/shorsfield/miniforge3/envs/biopython/bin/python"
use_python(path_to_python)
annotate_graphs = TRUE
trend_line = TRUE
plot_on_same_axis = TRUE
generate_plots = TRUE

py_run_string("
import numpy as np
from scipy.optimize import curve_fit

# fit asymptotic curve using exponential decay
def negative_exponential(x, b0, b1, b2): # based on https://isem-cueb-ztian.github.io/Intro-Econometrics-2017/handouts/lecture_notes/lecture_10/lecture_10.pdf and https://www.statforbiology.com/articles/usefulequations/
    return b0 - (b0 - b1) * np.exp(-b2 * x)
  
def fit_curve(x, y):
  popt, pcov = curve_fit(negative_exponential, x, y, p0=[1.0, 0.0, 1.0], bounds=([0.0, 0.0, 0.0], [1.0, 1.0, np.inf]))
  b0_err, b1_err, b2_err = np.sqrt(np.diag(pcov))
  return popt.tolist(), [b0_err, b1_err, b2_err]
")

#files <- Sys.glob("/Users/shorsfield/Library/Mobile Documents/com~apple~CloudDocs/Work/PhD/experiments/popPUNK-models/sample_distances/sample_wo_replacement/samples/*.txt")
files <- Sys.glob("/Users/shorsfield/Library/Mobile Documents/com~apple~CloudDocs/Work/PhD/experiments/popPUNK-models/sample_distances/sample_wo_replacement/samples/*.txt")
sample.replace <- FALSE
sample.num <- 30000

all_pairs <- data.frame(Species = c(), Core = c(), Accessory = c())

for (j in 1:length(files))
{
  filename <- files[j]
  df <- read.table(files[j], sep = "\t", comment.char = "")
  
  str <- str_split(filename, "_distances_")[[1]]
  species.name <- str_split(str[1], "/")[[1]]
  species.name <- species.name[length(species.name)]
  colnames(df) <- c("Species", "Sample", "Core", "Accessory")
  
  if (sample.replace == TRUE) {
    sample.name <- str_split(str[2], ".txt")[[1]]
    sample.name <- as.numeric(sample.name[1])
    df["Sample"] <- sample.name
  } else {
    drop <- c("Sample")
    df = df[,!(names(df) %in% drop)]
  }
  
  df["Species"] <- species.name
  
  sample.num.temp <- sample.num
  
  if (nrow(df) < sample.num.temp)
  {
    sample.num.temp <- nrow(df)
  }
  
  df = sample_n(df, sample.num.temp)
  
  all_pairs <- rbind(all_pairs, df)
}


# rename bacteria
all_pairs$Species[all_pairs$Species == "ab"] <- "A. baumanii"
all_pairs$Species[all_pairs$Species == "bp"] <- "B. pertussis"
all_pairs$Species[all_pairs$Species == "cj"] <- "C. jejuni"
all_pairs$Species[all_pairs$Species == "GBS_full"] <- "S. agalactiae"
all_pairs$Species[all_pairs$Species == "GPSv4"] <- "S. pneumoniae"
all_pairs$Species[all_pairs$Species == "e_feacalis"] <- "E. faecalis"
all_pairs$Species[all_pairs$Species == "e_feacium"] <- "E. faecium"
all_pairs$Species[all_pairs$Species == "ec"] <- "E. coli"
all_pairs$Species[all_pairs$Species == "hflu"] <- "H. influenzae"
all_pairs$Species[all_pairs$Species == "hp"] <- "H. pylori"
all_pairs$Species[all_pairs$Species == "kp"] <- "K. pneumoniae"
all_pairs$Species[all_pairs$Species == "lm"] <- "L. monocytogenes"
all_pairs$Species[all_pairs$Species == "ma"] <- "M. abscessus"
all_pairs$Species[all_pairs$Species == "mcat"] <- "M. catarrhalis"
all_pairs$Species[all_pairs$Species == "mtb"] <- "M. tuberculosis"
all_pairs$Species[all_pairs$Species == "ngon"] <- "N. gonorrhoeae"
all_pairs$Species[all_pairs$Species == "nm"] <- "N. meningitidis"
all_pairs$Species[all_pairs$Species == "pa"] <- "P. aeruginosa"
all_pairs$Species[all_pairs$Species == "sal"] <- "S. enterica"
all_pairs$Species[all_pairs$Species == "se"] <- "S. equisimilis"
all_pairs$Species[all_pairs$Species == "sm"] <- "S. maltophilia"
all_pairs$Species[all_pairs$Species == "staph_aureus"] <- "S. aureus"
all_pairs$Species[all_pairs$Species == "GAS"] <- "S. pyogenes"

all_species <- unique(all_pairs$Species)
max_core <- max(all_pairs$Core)
max_acc <- max(all_pairs$Accessory)
min_core <- min(all_pairs$Core)
min_acc <- min(all_pairs$Accessory)

all_coeffs <- data.frame(Species = c(), b0 = c(), b1 = c(), b2 = c(), b0_err = c(), b1_err = c(), b2_err = c())

for (i in 1:length(all_species))
{
  Species <- all_species[i]
  sample_points <- all_pairs[all_pairs$Species == Species,]
  
  params <- py$fit_curve(sample_points$Core, sample_points$Accessory)
  b0 <- params[[1]][[1]]
  b1 <- params[[1]][[2]]
  b2 <- params[[1]][[3]]
  b0_err <- params[[2]][[1]]
  b1_err <- params[[2]][[2]]
  b2_err <- params[[2]][[3]]
  
  drop <- c("Species","Sample")
  sample_points = sample_points[,!(names(sample_points) %in% drop)]

  if (generate_plots == TRUE)
  {
    p <- ggplot(sample_points, aes(x = Core, y = Accessory)) + geom_point(colour = "#2980B9", alpha=0.5) +
      geom_density_2d(aes(color = ..level..), bins=20) +
      theme_light() +
      scale_color_viridis_c() + theme_light() + ylab("Accessory Distance") + xlab("Core Distance") +
      theme(axis.text.x = element_text(size = 18), axis.text.y = element_text(size = 18), axis.title=element_text(size=16,face="bold"), strip.text.x = element_text(size = 12), legend.position="none")
    
    anno_pos_x <- min(min_core)
    anno_pos_y <- max(max_acc)
    
    if (trend_line == TRUE)
    {
      p <- p + stat_function(fun = function(x) b0 - (b0 - b1) * exp(-b2 * x),
                             color = "red", linewidth = 1)
      if (annotate_graphs == TRUE)
      {
        p <- p + annotate("text", x = anno_pos_x, y = anno_pos_y, hjust = 0, vjust = 1,
                          #label = sprintf("i = %.2f, j = %.2f, k = %.2f\ni err = %.2f, j err = %.2f, k err = %.2f", b0, b1, b2, b0_err, b1_err, b2_err),
                          label = sprintf("i = %.2f, j = %.2f, k = %.2f", b0, b1, b2),
                          size = 5)
      }
    }
    
    if (plot_on_same_axis == TRUE) {
      p <- p + coord_cartesian(xlim = c(0, max_core), ylim = c(0, max_acc))#scale_x_continuous(limits = c()) + scale_y_continuous(limits = c(0, max_acc)) # xlim(0, max_core) + ylim()
      anno_pos_x <- min_core
      anno_pos_y <- max_acc
      
      # magnify small plots
      if (max(sample_points$Core) <= (max_core / 6) | max(sample_points$Accessory) <= (max_acc / 6)) {
        #from <- c(xmin = min(sample_points$Core), xmax = max(sample_points$Core), ymin = min(sample_points$Accessory), ymax = max(sample_points$Accessory))
        to <- c(xmin = max_core * 0.35, xmax =  max_core * 0.95, ymin = max_acc * 0.25, ymax = max_acc * 0.85)
        
        p <- p + geom_magnify(aes(from = Core <= max(Core) & Accessory <= max(Accessory)), 
                              to = to, axes = "xy")
      }
    }
    
    
    #ggsave(paste("./species_distances/", Species, " contour_trend_line.svg", sep=""), plot = p, width=10, height = 6)
    ggsave(paste("./species_distances/", Species, " contour_trend_line.png", sep=""), plot = p, width=10, height = 6)
  }

  
  all_coeffs <- rbind(data.frame(Species, b0, b1, b2, b0_err, b1_err, b2_err), all_coeffs)
}

# plot specific species
{
  species <- c("M. tuberculosis", "S. pneumoniae", "E. coli", "L. monocytogenes")
  subset_all_pairs = subset(all_pairs, Species %in% species)
  subset_all_pairs$Species <- factor(subset_all_pairs$Species, levels = species)
  
  max_core <- max(subset_all_pairs$Core)
  max_acc <- max(subset_all_pairs$Accessory)
  min_core <- min(subset_all_pairs$Core)
  min_acc <- min(subset_all_pairs$Accessory)
  
  plot_list = list()
  
  for (spec in species)
  {
    sample_points <- all_pairs[all_pairs$Species == spec,]
    
    params <- py$fit_curve(sample_points$Core, sample_points$Accessory)
    b0 <- params[[1]][[1]]
    b1 <- params[[1]][[2]]
    b2 <- params[[1]][[3]]
    b0_err <- params[[2]][[1]]
    b1_err <- params[[2]][[2]]
    b2_err <- params[[2]][[3]]
    
    p <- ggplot(sample_points, aes(x = Core, y = Accessory)) + geom_point(colour = "#2980B9", alpha=0.5) +
      geom_density_2d(aes(color = ..level..), bins=20) +
      theme_light() +
      scale_color_viridis_c() + theme_light() + ylab("Accessory Distance") + xlab("Core Distance") +
      theme(axis.text.x = element_text(size = 16), axis.text.y = element_text(size = 16), axis.title=element_text(size=16,face="bold"), strip.text.x = element_text(size = 12), legend.position="none")
    
    anno_pos_x <- min(sample_points$Core)
    anno_pos_y <- max(sample_points$Accessory)
    
    if (trend_line == TRUE)
    {
      p <- p + stat_function(fun = function(x) b0 - (b0 - b1) * exp(-b2 * x),
                             color = "red", linewidth = 1)
      if (annotate_graphs == TRUE)
      {
        p <- p + annotate("text", x = anno_pos_x, y = anno_pos_y, hjust = 0, vjust = 1,
                          label = sprintf("b0 = %.2f, b1 = %.2f, b2 = %.2f\nb0_err = %.2f, b1_err = %.2f, b2_err = %.2f", b0, b1, b2, b0_err, b1_err, b2_err),
                          size = 5)
      }
    }
    
    if (plot_on_same_axis == TRUE) {
      p <- p + coord_cartesian(xlim = c(0, max_core), ylim = c(0, max_acc))#scale_x_continuous(limits = c()) + scale_y_continuous(limits = c(0, max_acc)) # xlim(0, max_core) + ylim()
      anno_pos_x <- min_core
      anno_pos_y <- max_acc
      
      # magnify small plots
      if (max(sample_points$Core) <= (max_core / 8) | max(sample_points$Accessory) <= (max_acc / 10)) {
        #from <- c(xmin = min(sample_points$Core), xmax = max(sample_points$Core), ymin = min(sample_points$Accessory), ymax = max(sample_points$Accessory))
        to <- c(xmin = max_core * 0.35, xmax =  max_core * 0.95, ymin = max_acc * 0.25, ymax = max_acc * 0.85)
        
        p <- p + geom_magnify(aes(from = Core <= max(Core) & Accessory <= max(Accessory)), 
                              to = to, axes = "xy")
      }
    }
    
    ggsave(paste("./species_distances/", spec, " contour_trend_line.svg", sep=""), plot = p, width=10, height = 6)
    ggsave(paste("./species_distances/", spec, " contour_trend_line.png", sep=""), plot = p, width=10, height = 6)
    p
    plot_list <- append(plot_list, list(p))
  }
  
  final_plot <- ggarrange(plotlist = plot_list, labels = c("A", "B", "C", "D"))
  final_plot
  
  annotate_figure(final_plot,
                  bottom = text_grob("Core distance", color = "black", face = "bold", size = 18),
                  left = text_grob("Accessory distance", color = "black", rot = 90, face = "bold", size = 18))
}


# compare with Bobay & Ochman 2018 https://pmc.ncbi.nlm.nih.gov/articles/PMC6186134/
{
  # get NE data
  ne_dataframe <- read_excel("Bobay_Ochman_2018_Ne_estimates.xlsx", sheet="thesis_uni")
  
  ne_dataframe$Species[ne_dataframe$Abbreviation == "ab"] <- "A. baumanii"
  ne_dataframe$Species[ne_dataframe$Abbreviation == "bp"] <- "B. pertussis"
  ne_dataframe$Species[ne_dataframe$Abbreviation == "cj"] <- "C. jejuni"
  ne_dataframe$Species[ne_dataframe$Abbreviation == "GBS_full"] <- "S. agalactiae"
  ne_dataframe$Species[ne_dataframe$Abbreviation == "GPSv4"] <- "S. pneumoniae"
  ne_dataframe$Species[ne_dataframe$Abbreviation == "e_feacalis"] <- "E. faecalis"
  ne_dataframe$Species[ne_dataframe$Abbreviation == "e_feacium"] <- "E. faecium"
  ne_dataframe$Species[ne_dataframe$Abbreviation == "ec"] <- "E. coli"
  ne_dataframe$Species[ne_dataframe$Abbreviation == "hflu"] <- "H. influenzae"
  ne_dataframe$Species[ne_dataframe$Abbreviation == "hp"] <- "H. pylori"
  ne_dataframe$Species[ne_dataframe$Abbreviation == "kp"] <- "K. pneumoniae"
  ne_dataframe$Species[ne_dataframe$Abbreviation == "lm"] <- "L. monocytogenes"
  ne_dataframe$Species[ne_dataframe$Abbreviation == "ma"] <- "M. abscessus"
  ne_dataframe$Species[ne_dataframe$Abbreviation == "mcat"] <- "M. catarrhalis"
  ne_dataframe$Species[ne_dataframe$Abbreviation == "mtb"] <- "M. tuberculosis"
  ne_dataframe$Species[ne_dataframe$Abbreviation == "ngon"] <- "N. gonorrhoeae"
  ne_dataframe$Species[ne_dataframe$Abbreviation == "nm"] <- "N. meningitidis"
  ne_dataframe$Species[ne_dataframe$Abbreviation == "pa"] <- "P. aeruginosa"
  ne_dataframe$Species[ne_dataframe$Abbreviation == "sal"] <- "S. enterica"
  ne_dataframe$Species[ne_dataframe$Abbreviation == "se"] <- "S. equisimilis"
  ne_dataframe$Species[ne_dataframe$Abbreviation == "sm"] <- "S. maltophilia"
  ne_dataframe$Species[ne_dataframe$Abbreviation == "staph_aureus"] <- "S. aureus"
  ne_dataframe$Species[ne_dataframe$Abbreviation == "GAS"] <- "S. pyogenes"
  ne_dataframe$Ne <- as.numeric(ne_dataframe$Ne)
  ne_dataframe$rm <- as.numeric(ne_dataframe$rm)
  ne_dataframe$hm <- as.numeric(ne_dataframe$hm)
  ne_dataframe$dNdS <- as.numeric(ne_dataframe$dNdS)
  ne_dataframe$dNdS_all <- as.numeric(ne_dataframe$dNdS_all)
  ne_dataframe$Phi <- as.numeric(ne_dataframe$Phi)
  ne_dataframe$Adjusted_pangenome <- as.numeric(ne_dataframe$Adjusted_pangenome)
  
  all_coeffs <- merge(all_coeffs, ne_dataframe, by="Species")
  
  my_comparisons <- list( c("animal pathogen", "commensal"), c("animal pathogen", "free living"), c("commensal", "free living") )
  
  spearmans.rho.bobay <- data.frame(Coeff = c(), Comparator = c(), Rho = c(), P.value = c())
  
  # boxplot comparisons with lifestyle
  {
    pbox_b0 <- ggplot(all_coeffs, aes(x=Life_Style, y = b0, color=Life_Style)) + theme_light() + geom_boxplot(outlier.shape = NA) + geom_jitter(size=2, alpha=0.4) + xlab("Lifestyle") + ylab("b0") + theme(axis.text.x = element_text(size = 18), axis.text.y = element_text(size = 18), axis.title=element_text(size=20,face="bold"), strip.text.x = element_text(size = 12), legend.position = "none") + stat_compare_means(comparisons = my_comparisons, method = "wilcox.test", aes(label = ..p.adj..), p.adjust.method = "bonferroni")
    pbox_b0
    pbox_b1 <- ggplot(all_coeffs, aes(x=Life_Style, y = b1, color=Life_Style)) + theme_light() + geom_boxplot(outlier.shape = NA) + geom_jitter(size=2, alpha=0.4) + xlab("Lifestyle") + ylab("b1") + theme(axis.text.x = element_text(size = 18), axis.text.y = element_text(size = 18), axis.title=element_text(size=20,face="bold"), strip.text.x = element_text(size = 12), legend.position = "none") + stat_compare_means(comparisons = my_comparisons, method = "wilcox.test", aes(label = ..p.adj..), p.adjust.method = "bonferroni")
    pbox_b1
    pbox_b2 <- ggplot(all_coeffs, aes(x=Life_Style, y = b2, color=Life_Style)) + theme_light() + geom_boxplot(outlier.shape = NA) + geom_jitter(size=2, alpha=0.4) + xlab("Lifestyle") + ylab("b2") + theme(axis.text.x = element_text(size = 18), axis.text.y = element_text(size = 18), axis.title=element_text(size=20,face="bold"), strip.text.x = element_text(size = 12), legend.position = "none") + stat_compare_means(comparisons = my_comparisons, method = "wilcox.test", aes(label = ..p.adj..), p.adjust.method = "bonferroni")
    pbox_b2
  }
  
  corr.function <- function(x, y, spearmans.rho, Coeff, Comparator)
  {
    corr <- cor.test(x=x, y=y, method = 'spearman')
    spearmans.rho <- rbind(spearmans.rho, data.frame(Coeff = Coeff, Comparator = Comparator, Rho = corr$statistic, P.value = corr$p.value))
    spearmans.rho
  }
  
  # NE comparisons
  {
    pNe_b0 <- ggplot(all_coeffs, aes(x=Ne, y = b0, color=Life_Style)) + theme_light() + geom_point(size=4) + xlab("Effective population size (Ne)") + ylab("b0") + theme(axis.text.x = element_text(size = 18), axis.text.y = element_text(size = 18), axis.title=element_text(size=20,face="bold"), strip.text.x = element_text(size = 12), legend.title=element_text(size=20,face="bold"), legend.text=element_text(size=16)) + guides(colour=guide_legend(title="Lifestyle"))# + scale_x_continuous(labels = label_scientific(digits = 2))
    pNe_b0
    pNe_b1 <- ggplot(all_coeffs, aes(x=Ne, y = b1, color=Life_Style)) + theme_light() + geom_point(size=4) + xlab("Effective population size (Ne)") + ylab("b1") + theme(axis.text.x = element_text(size = 18), axis.text.y = element_text(size = 18), axis.title=element_text(size=20,face="bold"), strip.text.x = element_text(size = 12), legend.title=element_text(size=20,face="bold"), legend.text=element_text(size=16)) + guides(colour=guide_legend(title="Lifestyle"))# + scale_x_continuous(labels = label_scientific(digits = 2))
    pNe_b1
    pNe_b2 <- ggplot(all_coeffs, aes(x=Ne, y = b2, color=Life_Style)) + theme_light() + geom_point(size=4) + xlab("Effective population size (Ne)") + ylab("b2") + theme(axis.text.x = element_text(size = 18), axis.text.y = element_text(size = 18), axis.title=element_text(size=20,face="bold"), strip.text.x = element_text(size = 12), legend.title=element_text(size=20,face="bold"), legend.text=element_text(size=16)) + guides(colour=guide_legend(title="Lifestyle"))# + scale_x_continuous(labels = label_scientific(digits = 2))
    pNe_b2
    spearmans.rho.bobay <- corr.function(all_coeffs$b0, all_coeffs$Ne, spearmans.rho.bobay, "b0", "Ne")
    spearmans.rho.bobay <- corr.function(all_coeffs$b1, all_coeffs$Ne, spearmans.rho.bobay, "b1", "Ne")
    spearmans.rho.bobay <- corr.function(all_coeffs$b2, all_coeffs$Ne, spearmans.rho.bobay, "b2", "Ne")
  }
  
  # genome fluidity comparisons
  {
    pPhi_b0 <- ggplot(all_coeffs, aes(x=Phi, y = b0, color=Life_Style)) + theme_light() + geom_point(size=4) + xlab("Genome fluidity") + ylab("b0") + theme(axis.text.x = element_text(size = 18), axis.text.y = element_text(size = 18), axis.title=element_text(size=20,face="bold"), strip.text.x = element_text(size = 12), legend.title=element_text(size=20,face="bold"), legend.text=element_text(size=16)) + guides(colour=guide_legend(title="Lifestyle"))# + scale_x_continuous(labels = label_scientific(digits = 2))
    pPhi_b0
    pPhi_b1 <- ggplot(all_coeffs, aes(x=Phi, y = b1, color=Life_Style)) + theme_light() + geom_point(size=4) + xlab("Genome fluidity") + ylab("b1") + theme(axis.text.x = element_text(size = 18), axis.text.y = element_text(size = 18), axis.title=element_text(size=20,face="bold"), strip.text.x = element_text(size = 12), legend.title=element_text(size=20,face="bold"), legend.text=element_text(size=16)) + guides(colour=guide_legend(title="Lifestyle"))# + scale_x_continuous(labels = label_scientific(digits = 2))
    pPhi_b1
    pPhi_b2 <- ggplot(all_coeffs, aes(x=Phi, y = b2, color=Life_Style)) + theme_light() + geom_point(size=4) + xlab("Genome fluidity") + ylab("b2") + theme(axis.text.x = element_text(size = 18), axis.text.y = element_text(size = 18), axis.title=element_text(size=20,face="bold"), strip.text.x = element_text(size = 12), legend.title=element_text(size=20,face="bold"), legend.text=element_text(size=16)) + guides(colour=guide_legend(title="Lifestyle"))# + scale_x_continuous(labels = label_scientific(digits = 2))
    pPhi_b2
    spearmans.rho.bobay <- corr.function(all_coeffs$b0, all_coeffs$Phi, spearmans.rho.bobay, "b0", "Phi")
    spearmans.rho.bobay <- corr.function(all_coeffs$b1, all_coeffs$Phi, spearmans.rho.bobay, "b1", "Phi")
    spearmans.rho.bobay <- corr.function(all_coeffs$b2, all_coeffs$Phi, spearmans.rho.bobay, "b2", "Phi")
  }
  
  # recombination comparisons
  {
    prm_b0 <- ggplot(all_coeffs, aes(x=rm, y = b0, color=Life_Style)) + theme_light() + geom_point(size=4) + xlab("R/M") + ylab("b0") + theme(axis.text.x = element_text(size = 18), axis.text.y = element_text(size = 18), axis.title=element_text(size=20,face="bold"), strip.text.x = element_text(size = 12), legend.title=element_text(size=20,face="bold"), legend.text=element_text(size=16)) + guides(colour=guide_legend(title="Lifestyle"))# + scale_x_continuous(labels = label_scientific(digits = 2))
    prm_b0
    prm_b1 <- ggplot(all_coeffs, aes(x=rm, y = b1, color=Life_Style)) + theme_light() + geom_point(size=4) + xlab("R/M") + ylab("b1") + theme(axis.text.x = element_text(size = 18), axis.text.y = element_text(size = 18), axis.title=element_text(size=20,face="bold"), strip.text.x = element_text(size = 12), legend.title=element_text(size=20,face="bold"), legend.text=element_text(size=16)) + guides(colour=guide_legend(title="Lifestyle"))# + scale_x_continuous(labels = label_scientific(digits = 2))
    prm_b1
    prm_b2 <- ggplot(all_coeffs, aes(x=rm, y = b2, color=Life_Style)) + theme_light() + geom_point(size=4) + xlab("R/M") + ylab("b2") + theme(axis.text.x = element_text(size = 18), axis.text.y = element_text(size = 18), axis.title=element_text(size=20,face="bold"), strip.text.x = element_text(size = 12), legend.title=element_text(size=20,face="bold"), legend.text=element_text(size=16)) + guides(colour=guide_legend(title="Lifestyle"))# + scale_x_continuous(labels = label_scientific(digits = 2))
    prm_b2
    spearmans.rho.bobay <- corr.function(all_coeffs$b0, all_coeffs$rm, spearmans.rho.bobay, "b0", "rm")
    spearmans.rho.bobay <- corr.function(all_coeffs$b1, all_coeffs$rm, spearmans.rho.bobay, "b1", "rm")
    spearmans.rho.bobay <- corr.function(all_coeffs$b2, all_coeffs$rm, spearmans.rho.bobay, "b2", "rm")
  }
  
  # hm comparisons
  {
    phm_b0 <- ggplot(all_coeffs, aes(x=hm, y = b0, color=Life_Style)) + theme_light() + geom_point(size=4) + xlab("homoplasic / non-homoplasic alleles") + ylab("b0") + theme(axis.text.x = element_text(size = 18), axis.text.y = element_text(size = 18), axis.title=element_text(size=20,face="bold"), strip.text.x = element_text(size = 12), legend.title=element_text(size=20,face="bold"), legend.text=element_text(size=16)) + guides(colour=guide_legend(title="Lifestyle"))# + scale_x_continuous(labels = label_scientific(digits = 2))
    phm_b0
    phm_b1 <- ggplot(all_coeffs, aes(x=hm, y = b1, color=Life_Style)) + theme_light() + geom_point(size=4) + xlab("homoplasic / non-homoplasic alleles") + ylab("b1") + theme(axis.text.x = element_text(size = 18), axis.text.y = element_text(size = 18), axis.title=element_text(size=20,face="bold"), strip.text.x = element_text(size = 12), legend.title=element_text(size=20,face="bold"), legend.text=element_text(size=16)) + guides(colour=guide_legend(title="Lifestyle"))# + scale_x_continuous(labels = label_scientific(digits = 2))
    phm_b1
    phm_b2 <- ggplot(all_coeffs, aes(x=hm, y = b2, color=Life_Style)) + theme_light() + geom_point(size=4) + xlab("homoplasic / non-homoplasic alleles") + ylab("b2") + theme(axis.text.x = element_text(size = 18), axis.text.y = element_text(size = 18), axis.title=element_text(size=20,face="bold"), strip.text.x = element_text(size = 12), legend.title=element_text(size=20,face="bold"), legend.text=element_text(size=16)) + guides(colour=guide_legend(title="Lifestyle"))# + scale_x_continuous(labels = label_scientific(digits = 2))
    phm_b2
    spearmans.rho.bobay <- corr.function(all_coeffs$b0, all_coeffs$hm, spearmans.rho.bobay, "b0", "hm")
    spearmans.rho.bobay <- corr.function(all_coeffs$b1, all_coeffs$hm, spearmans.rho.bobay, "b1", "hm")
    spearmans.rho.bobay <- corr.function(all_coeffs$b2, all_coeffs$hm, spearmans.rho.bobay, "b2", "hm")
  }
  
  # dnds comparisons
  {
    pdNdS_all_b0 <- ggplot(all_coeffs, aes(x=dNdS_all, y = b0, color=Life_Style)) + theme_light() + geom_point(size=4) + xlab("dN/dS (all genes)") + ylab("b0") + theme(axis.text.x = element_text(size = 18), axis.text.y = element_text(size = 18), axis.title=element_text(size=20,face="bold"), strip.text.x = element_text(size = 12), legend.title=element_text(size=20,face="bold"), legend.text=element_text(size=16)) + guides(colour=guide_legend(title="Lifestyle"))# + scale_x_continuous(labels = label_scientific(digits = 2))
    pdNdS_all_b0
    pdNdS_all_b1 <- ggplot(all_coeffs, aes(x=dNdS_all, y = b1, color=Life_Style)) + theme_light() + geom_point(size=4) + xlab("dN/dS (all genes)") + ylab("b1") + theme(axis.text.x = element_text(size = 18), axis.text.y = element_text(size = 18), axis.title=element_text(size=20,face="bold"), strip.text.x = element_text(size = 12), legend.title=element_text(size=20,face="bold"), legend.text=element_text(size=16)) + guides(colour=guide_legend(title="Lifestyle"))# + scale_x_continuous(labels = label_scientific(digits = 2))
    pdNdS_all_b1
    pdNdS_all_b2 <- ggplot(all_coeffs, aes(x=dNdS_all, y = b2, color=Life_Style)) + theme_light() + geom_point(size=4) + xlab("dN/dS (all genes)") + ylab("b2") + theme(axis.text.x = element_text(size = 18), axis.text.y = element_text(size = 18), axis.title=element_text(size=20,face="bold"), strip.text.x = element_text(size = 12), legend.title=element_text(size=20,face="bold"), legend.text=element_text(size=16)) + guides(colour=guide_legend(title="Lifestyle"))# + scale_x_continuous(labels = label_scientific(digits = 2))
    pdNdS_all_b2
    spearmans.rho.bobay <- corr.function(all_coeffs$b0, all_coeffs$dNdS_all, spearmans.rho.bobay, "b0", "dNdS_all")
    spearmans.rho.bobay <- corr.function(all_coeffs$b1, all_coeffs$dNdS_all, spearmans.rho.bobay, "b1", "dNdS_all")
    spearmans.rho.bobay <- corr.function(all_coeffs$b2, all_coeffs$dNdS_all, spearmans.rho.bobay, "b2", "dNdS_all")
  }
  
  # pangenome size comparisons
  {
    ppan_b0 <- ggplot(all_coeffs, aes(x=Adjusted_pangenome, y = b0, color=Life_Style)) + theme_light() + geom_point(size=4) + xlab("Pangenome size") + ylab("b0") + theme(axis.text.x = element_text(size = 18), axis.text.y = element_text(size = 18), axis.title=element_text(size=20,face="bold"), strip.text.x = element_text(size = 12), legend.title=element_text(size=20,face="bold"), legend.text=element_text(size=16)) + guides(colour=guide_legend(title="Lifestyle"))# + scale_x_continuous(labels = label_scientific(digits = 2))
    ppan_b0
    ppan_b1 <- ggplot(all_coeffs, aes(x=Adjusted_pangenome, y = b1, color=Life_Style)) + theme_light() + geom_point(size=4) + xlab("Pangenome size") + ylab("b1") + theme(axis.text.x = element_text(size = 18), axis.text.y = element_text(size = 18), axis.title=element_text(size=20,face="bold"), strip.text.x = element_text(size = 12), legend.title=element_text(size=20,face="bold"), legend.text=element_text(size=16)) + guides(colour=guide_legend(title="Lifestyle"))# + scale_x_continuous(labels = label_scientific(digits = 2))
    ppan_b1
    ppan_b2 <- ggplot(all_coeffs, aes(x=Adjusted_pangenome, y = b2, color=Life_Style)) + theme_light() + geom_point(size=4) + xlab("Pangenome size") + ylab("b2") + theme(axis.text.x = element_text(size = 18), axis.text.y = element_text(size = 18), axis.title=element_text(size=20,face="bold"), strip.text.x = element_text(size = 12), legend.title=element_text(size=20,face="bold"), legend.text=element_text(size=16)) + guides(colour=guide_legend(title="Lifestyle"))# + scale_x_continuous(labels = label_scientific(digits = 2))
    ppan_b2
    spearmans.rho.bobay <- corr.function(all_coeffs$b0, all_coeffs$Adjusted_pangenome, spearmans.rho.bobay, "b0", "Adjusted_pangenome")
    spearmans.rho.bobay <- corr.function(all_coeffs$b1, all_coeffs$Adjusted_pangenome, spearmans.rho.bobay, "b1", "Adjusted_pangenome")
    spearmans.rho.bobay <- corr.function(all_coeffs$b2, all_coeffs$Adjusted_pangenome, spearmans.rho.bobay, "b2", "Adjusted_pangenome")
  }
  
  spearmans.rho.bobay$P.value.adj <- p.adjust(spearmans.rho.bobay$P.value, method="BH")
}

write.csv(spearmans.rho.bobay, "spearmans_rho_bobay.csv", row.names=FALSE, quote=FALSE) 



