library(ggplot2)
library(stringr)
library(dplyr)
library(ggpubr)
library(ggsci)

real.files <- Sys.glob("sample_distances/*.txt")
sim.files <- Sys.glob("pansim_distances/chosen_dists/*.tsv")
sample.replace <- FALSE
sample.num <- 30000
plot_on_same_axis = FALSE

all_pairs <- data.frame(Species = c(), Core = c(), Accessory = c())

for (j in 1:length(real.files))
{
  filename <- real.files[j]
  df <- read.table(real.files[j], sep = "\t", comment.char = "")
  
  str <- str_split(filename, "_distances_")[[1]]
  species.name <- str_split(str[1], "/")[[1]]
  species.name <- species.name[length(species.name)]
  colnames(df) <- c("Species", "Sample", "Core", "Accessory")
  
  df["Species"] <- species.name
  drop <- c("Sample")
  df = df[,!(names(df) %in% drop)]
  
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

sim_pairs <- data.frame(Species = c(), Core = c(), Accessory = c())

#j <- 1 
for (j in 1:length(sim.files))
{
  filename <- sim.files[j]
  df <- read.table(sim.files[j], sep = "\t", comment.char = "")
  
  str <- str_split(filename, "_task_")[[1]]
  species.name <- str_split(str[2], "_")[[1]]
  species.name <- str_split(species.name[[2]], ".tsv")[[1]][[1]]
  
  df["Species"] <- species.name
  colnames(df) <- c("Core", "Accessory", "Species")
  
  sample.num.temp <- sample.num
  
  if (nrow(df) < sample.num.temp)
  {
    sample.num.temp <- nrow(df)
  }
  
  df = sample_n(df, sample.num.temp)
  
  sim_pairs <- rbind(sim_pairs, df)
}

sim_pairs$Species[sim_pairs$Species == "pneumo"] <- "S. pneumoniae"
sim_pairs$Species[sim_pairs$Species == "ecoli"] <- "E. coli"
sim_pairs$Species[sim_pairs$Species == "lm"] <- "L. monocytogenes"
sim_pairs$Species[sim_pairs$Species == "MtB"] <- "M. tuberculosis"

all_pairs$Type <- "Real"
sim_pairs$Type <- "Sim"

merged.df <- rbind(all_pairs, sim_pairs)

all_species <- unique(merged.df$Species)
max_core <- max(merged.df$Core)
max_acc <- max(merged.df$Accessory)
min_core <- min(merged.df$Core)
min_acc <- min(merged.df$Accessory)

#i <- 1
for (i in 1:length(all_species))
{
  Species <- all_species[i]
  sample_points <- merged.df[merged.df$Species == Species,]
  sample_points_real <- sample_points[sample_points$Type == "Real",]
  sample_points_sim <- sample_points[sample_points$Type == "Sim",]
  
  max_x = max(sample_points$Core)
  max_y = max(sample_points$Accessory)
  
  p <- ggplot(sample_points_real, aes(x = Core, y = Accessory)) + geom_point(colour= "#2980B9", alpha=0.5) +
    geom_density_2d(data =sample_points_real, colour= "#3A5FCD", bins=10) +
    geom_point(data=sample_points_sim, aes(x = Core, y = Accessory), colour= "#FF0000", alpha=0.5) +
    geom_density_2d(data =sample_points_sim,  color = "#8B0000", bins=10) +
    theme_light() + ylab("Accessory Distance") + xlab("Core Distance") +
    theme(axis.text.x = element_text(size = 18), axis.text.y = element_text(size = 18), axis.title=element_text(size=16,face="bold"), strip.text.x = element_text(size = 12), legend.position="none") +
    scale_x_continuous(limits = c(0, max_x * 1.1)) + scale_y_continuous(limits = c(0, max_y * 1.1))
  p
  anno_pos_x <- min(sample_points$Core)
  anno_pos_y <- max(sample_points$Accessory)
  
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

  ggsave(paste("./pansim_distances/", Species, " pansim_contour.svg", sep=""), plot = p, width=10, height = 6)
  ggsave(paste("./pansim_distances/", Species, " pansim_contour.png", sep=""), plot = p, width=10, height = 6)
}

