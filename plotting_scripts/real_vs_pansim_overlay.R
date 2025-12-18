library(tidyverse)
library(ggplot2)
library(ggsci)
library(ggpmisc)
library(ggrepel)
library(tools)
library(ggpubr)

rename_species <- function(df){
  df$species_short[df$taxa == "Acinetobacter_baumannii"] <- "A. baumanii"
  df$species_short[df$taxa == "Bordetalla_pertussis"] <- "B. pertussis"
  df$species_short[df$taxa == "Campylobacter_jejuni"] <- "C. jejuni"
  df$species_short[df$taxa == "Streptococcus_agalactiae"] <- "S. agalactiae"
  df$species_short[df$taxa == "Streptococcus_pneumoniae"] <- "S. pneumoniae"
  df$species_short[df$taxa == "Enterococcus_faecalis"] <- "E. faecalis"
  df$species_short[df$taxa == "Enterococcus_faecium"] <- "E. faecium"
  df$species_short[df$taxa == "Escherichia_coli"] <- "E. coli"
  df$species_short[df$taxa == "Haemophilus_influenzae"] <- "H. influenzae"
  df$species_short[df$taxa == "Helicobacter_pylori"] <- "H. pylori"
  df$species_short[df$taxa == "Klebsiella_pneumoniae"] <- "K. pneumoniae"
  df$species_short[df$taxa == "Listeria_monocytogenes"] <- "L. monocytogenes"
  df$species_short[df$taxa == "Mycobacterium_abscessus"] <- "M. abscessus"
  df$species_short[df$taxa == "Morexalla_catarrhalis"] <- "M. catarrhalis"
  df$species_short[df$taxa == "Mycobacterium_tuberculosis"] <- "M. tuberculosis"
  df$species_short[df$taxa == "Neisseria_gonorrhoeae"] <- "N. gonorrhoeae"
  df$species_short[df$taxa == "Neisseria_meningitidis"] <- "N. meningitidis"
  df$species_short[df$taxa == "Pseudomonas_aeruginosa"] <- "P. aeruginosa"
  df$species_short[df$taxa == "Salmonella_enterica"] <- "S. enterica"
  df$species_short[df$taxa == "Streptococcus_dysgalactiae_subspecies_equisimilis"] <- "S. equisimilis"
  df$species_short[df$taxa == "Stenotrophomonas_maltophilia"] <- "S. maltophilia"
  df$species_short[df$taxa == "Staphylococcus_aureus"] <- "S. aureus"
  df$species_short[df$taxa == "Streptococcus_pyogenes"] <- "S. pyogenes"
  
  df
}

sample.df <- function(df, sample.num) {
  sample.num.temp <- sample.num
  
  if (nrow(df) < sample.num.temp)
  {
    sample.num.temp <- nrow(df)
  }
  
  df = sample_n(df, sample.num.temp)
  df
}


overall_indir <- ""
indir <- paste0(overall_indir, "real_vs_pansim_overlay/")
outpref <- paste0(indir, "figures/")
sample.num <- 30000

# pansim running
{
  pangenome.file <- paste0(indir, "data/", "ppmod_real_input.txt")
  summary.file <- paste0(indir, "data/", "ppmod_summary_N24.csv")
  
  pangenome.df <- read.csv(pangenome.file, sep="\t")
  summary.df <- read.csv(summary.file)
  
  overall.df <- merge(summary.df, pangenome.df, by.x = "taxa", by.y = "name")
  
  write_tsv(overall.df, paste0(outpref, "pansim_run_params.tsv"))
}

# overlay plots
{
  real_genomes <- Sys.glob(paste0(indir, "data/real/*_real.txt"))
  sim_genomes <- Sys.glob(paste0(indir, "data/pansim/*_pansim.tsv"))
  
  real_genomes.df <- data.frame(file_real = real_genomes)
  real_genomes.df$taxa <- file_path_sans_ext(basename(real_genomes.df$file_real))
  real_genomes.df$taxa <- sub("_[^_]+$", "", real_genomes.df$taxa)
  
  sim_genomes.df <- data.frame(file_sim = sim_genomes)
  sim_genomes.df$taxa <- file_path_sans_ext(basename(sim_genomes.df$file_sim))
  sim_genomes.df$taxa <- sub("_[^_]+$", "", sim_genomes.df$taxa)
  
  merged.df <- merge(real_genomes.df, sim_genomes.df, by = "taxa")
  
  
  merged.df <- rename_species(merged.df)
  
  idx_list <- seq(1, nrow(merged.df))
  #idx_list <- seq(1, 4)
  idx <- 1
  plots <- lapply(idx_list, function(idx) {
    row <- merged.df[idx,]
    
    taxa <- row$taxa
    species_short <- row$species_short
  
    sim.dist <- row$file_sim
    real.dist <- row$file_real
    
    sim.df <- read.csv(sim.dist, sep = "\t", header = FALSE)
    real.df <- read.csv(real.dist, sep = "\t", header = FALSE)
    
    sim.df <- sample.df(sim.df, sample.num)
    real.df <- sample.df(real.df, sample.num)
    
    colnames <- c("Core", "Accessory")
    colnames(sim.df) <- colnames
    colnames(real.df) <- colnames
    
    anno_pos_x <- min(sim.df$Core, real.df$Core)
    max_pos_x <- max(sim.df$Core, real.df$Core)
    anno_pos_y <- max(sim.df$Accessory, real.df$Accessory)
    anno_range <- (max_pos_x - anno_pos_x) * 0.2
    
    anno_pos_x <- anno_pos_x + anno_range
    
    p <- ggplot(real.df, aes(x = Core, y = Accessory)) + geom_point(colour= "#2980B9", alpha=0.5) +
      geom_density_2d(data =real.df, colour= "#3A5FCD", bins=10) +
      geom_point(data=sim.df, aes(x = Core, y = Accessory), colour= "#FF0000", alpha=0.5) +
      geom_density_2d(data =sim.df,  color = "#8B0000", bins=10) +
      theme_light() + ylab("Accessory Distance") + xlab("Core Distance") +
      theme(axis.text.x = element_text(size = 18), axis.text.y = element_text(size = 18), axis.title=element_blank(), strip.text.x = element_text(size = 12), legend.position="none") +
      annotate("text", x=anno_pos_x, y=anno_pos_y, label= species_short, fontface = 'italic')
    p
  })
  
  plots[[2]]
  
  combined_plots <- ggarrange(
    plotlist=plots,
    #ncol = 3,
    #align = "v",
    labels=NULL)
  
  combined_plots <- annotate_figure(combined_plots,
                                    left = text_grob("Core Distance", color = "black", face = "bold", rot = 90, size = 20),
                                    bottom = text_grob("Accessory Distance", color = "black", face = "bold", size = 20))
  #combined_plots
  
  ggsave(paste(outpref, "/real_vs_pansim_overlay.png", sep=""), plot = combined_plots, width=20, height = 20)
  #ggsave(paste(outpref, "/param_step_BOLFI_fit.svg", sep=""), plot = combined_plots, width=16, height = 16)
  
}