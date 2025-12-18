library(ggplot2)
library(dplyr)
library(ggsci)
library(GGally)
library(tools)
library(ggExtra)
library(purrr)
library(ggpubr)
library(scales)
library(cowplot)
library(readxl)

species_list <- c(
  "Acinetobacter_baumannii",
  "Bordetalla_pertussis",
  "Campylobacter_jejuni",
  "Streptococcus_agalactiae",
  "Streptococcus_pneumoniae",
  "Enterococcus_faecalis",
  "Enterococcus_faecium",
  "Escherichia_coli",
  "Haemophilus_influenzae",
  "Helicobacter_pylori",
  "Klebsiella_pneumoniae",
  "Listeria_monocytogenes",
  "Mycobacterium_abscessus",
  "Morexalla_catarrhalis",
  "Mycobacterium_tuberculosis",
  "Neisseria_gonorrhoeae",
  "Neisseria_meningitidis",
  "Pseudomonas_aeruginosa",
  "Salmonella_enterica",
  "Streptococcus_dysgalactiae_subspecies_equisimilis",
  "Stenotrophomonas_maltophilia",
  "Staphylococcus_aureus",
  "Streptococcus_pyogenes"
)

#compare with Bobay & Ochman 2018 https://pmc.ncbi.nlm.nih.gov/articles/PMC6186134/
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
}

ori.col.names <- c("rate_genes2", "core_genes", "core_nuc", 
                   "pan_size", "n_gen", "avg_gene_freq", "pop_size",
                   "rate_genes1", "core_mu", "prop_positive", 
                   "HR_rate", "HGT_rate", "prop_genes2", 
                   "competition_strength", "pos_lambda", "neg_lambda", "all_params")
col.names <- c("Fast gene turnover rate", "Core genome size (genes)", "Core genome size (nucleotides)", 
               "Pangenome size (genes)", "No. generations", "Average proportion of pangenome per genome",  "Population size",
               "Basal gene turnover rate", "Core mutation rate", "Proportion +ve selected genes",
               "HR rate", "HGT rate", "Proportion fast genes",
               "Competition strength", "Gene +ve selection lambda", "Gene -ve selection lambda", "All parameters")

rename_species <- function(df){
  df$species_short[df$Species == "Acinetobacter_baumannii"] <- "A. baumanii"
  df$species_short[df$Species == "Bordetalla_pertussis"] <- "B. pertussis"
  df$species_short[df$Species == "Campylobacter_jejuni"] <- "C. jejuni"
  df$species_short[df$Species == "Streptococcus_agalactiae"] <- "S. agalactiae"
  df$species_short[df$Species == "Streptococcus_pneumoniae"] <- "S. pneumoniae"
  df$species_short[df$Species == "Enterococcus_faecalis"] <- "E. faecalis"
  df$species_short[df$Species == "Enterococcus_faecium"] <- "E. faecium"
  df$species_short[df$Species == "Escherichia_coli"] <- "E. coli"
  df$species_short[df$Species == "Haemophilus_influenzae"] <- "H. influenzae"
  df$species_short[df$Species == "Helicobacter_pylori"] <- "H. pylori"
  df$species_short[df$Species == "Klebsiella_pneumoniae"] <- "K. pneumoniae"
  df$species_short[df$Species == "Listeria_monocytogenes"] <- "L. monocytogenes"
  df$species_short[df$Species == "Mycobacterium_abscessus"] <- "M. abscessus"
  df$species_short[df$Species == "Morexalla_catarrhalis"] <- "M. catarrhalis"
  df$species_short[df$Species == "Mycobacterium_tuberculosis"] <- "M. tuberculosis"
  df$species_short[df$Species == "Neisseria_gonorrhoeae"] <- "N. gonorrhoeae"
  df$species_short[df$Species == "Neisseria_meningitidis"] <- "N. meningitidis"
  df$species_short[df$Species == "Pseudomonas_aeruginosa"] <- "P. aeruginosa"
  df$species_short[df$Species == "Salmonella_enterica"] <- "S. enterica"
  df$species_short[df$Species == "Streptococcus_dysgalactiae_subspecies_equisimilis"] <- "S. equisimilis"
  df$species_short[df$Species == "Stenotrophomonas_maltophilia"] <- "S. maltophilia"
  df$species_short[df$Species == "Staphylococcus_aureus"] <- "S. aureus"
  df$species_short[df$Species == "Streptococcus_pyogenes"] <- "S. pyogenes"
  
  df
}

#compare with Bobay & Ochman 2018 https://pmc.ncbi.nlm.nih.gov/articles/PMC6186134/
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
}

search_string <- function(filename, string){
  label = if (grepl(string, filename)) {
    TRUE
  } else {
    FALSE
  }
  return(label)
}

parse_speces <- function(filename) {
  species <- ""
  for (species_name in species_list) {
    if (grepl(species_name, filename)) {
      species <- species_name
    } else {
      next
    }
  }
  species
}

# v placehold stands for value
parse_filename <- function(filename) {
  list(
    species = parse_speces(filename),
    all_params = search_string(filename, "all_params"),
    rate_genes2 = search_string(filename, "rate_genes2"),
    core_genes = search_string(filename, "core_genes"),
    core_nuc = search_string(filename, "core_nuc"),
    pan_size = search_string(filename, "pan_size"),
    n_gen = search_string(filename, "n_gen"),
    avg_gene_freq = search_string(filename, "avg_gene_freq"),
    pop_size = search_string(filename, "pop_size"),
    prop_positive = search_string(filename, "prop_positive"),
    HR_rate = search_string(filename, "HR_rate"),
    HGT_rate = search_string(filename, "HGT_rate"),
    rate_genes1 = search_string(filename, "rate_genes1"),
    prop_genes2 = search_string(filename, "prop_genes2"),
    core_mu = search_string(filename, "core_mu"),
    competition_strength = search_string(filename, "competition_strength"),
    pos_lambda = search_string(filename, "pos_lambda"),
    neg_lambda = search_string(filename, "neg_lambda")
  )
}

get_summary <- function(df, experiment_id, species)
{
  idx_keep <- which(!grepl("_v$", colnames(df)))
  col.names.ori <- colnames(df)[idx_keep]
  
  # strip "_v" and compare
  names_keep <- match(col.names.ori, ori.col.names)
  col.names.new <- col.names[names_keep]
  
  # Generate marginal plots for each pair
  idx <- 4
  
  # here, need to generate a row table of one row, with the unmatched paramter, its value, and all subsequent parameters changed with their true value, median and 95% credible
  # interval, as well as whether the true parameter lies within
  summary_df <- sapply(idx_keep, function(idx) {
    xvar <- idx_keep[idx]
    
    # get original column names
    x.name.ori <- col.names.ori[idx]
    
    # get new column names
    x.name <- col.names.new[xvar]
    
    x.median <- median(df[[x.name.ori]])
    if (is.na(x.median)) {
      return(c(experiment_id, species, x.name.ori, x.name, NA, NA, NA, NA))
    }
    
    x.5quantile <- quantile(df[[x.name.ori]], 0.025)
    x.95quantile <- quantile(df[[x.name.ori]], 0.975)
    
    range.size <- (x.95quantile - x.5quantile) / abs(x.median)
    if (is.nan(range.size)) {
      return(c(experiment_id, species, x.name.ori, x.name, NA, NA, NA, NA))
    }
    
    return(c(experiment_id, species, x.name.ori, x.name, x.5quantile, x.median, x.95quantile, range.size))
  })
  
  summary_df <- as.data.frame(t(summary_df))
  
  colnames(summary_df) <- c("Experiment", "Species", "Fitted_param_short", "Fitted_param_long", "LCI", "Median", "UCI", "CI_size")
  summary_df <- na.omit(summary_df)
  
  # determine which parameters were fitted, assign to experiment id
  fitted_params <- unique(summary_df$Fitted_param_short)
  
  summary_df$Experiment <- paste(sort(fitted_params), collapse="_")
  
  summary_df
}

parse_results <- function(df_paths, outpref)
{
  j <- 18
  for (j in 1:length(df_paths))
  {
    filename <- df_paths[j]
    df <- read.table(filename, sep = ",", comment.char = "", header=TRUE)
    
    # add metadata only if not present
    parsed_values <- parse_filename(filename)
    if (parsed_values['all_params'] == TRUE) {
      for (name in names(parsed_values)) {
          df[[name]] <- TRUE
        }
      } else {
        for (name in names(parsed_values)) {
          if (parsed_values[[name]] != TRUE) {
            df[[name]] <- parsed_values[[name]]
          }
        }
    }

    
    # remove values that aren't present
    #df <- df[,colSums(is.na(df))<nrow(df)]
    
    base.name <- file_path_sans_ext(basename(filename))
    species <- parsed_values[["species"]]
    
    df <- get_summary(df, j, species)
    
    # create df_all if first iteration
    if (j == 1)
    {
      df_all <- df
    } else {
      df_all <- rbind(df_all, df)
    }
  }
  
  df_all
}

plot_experiments <- function(summary_df, outpref, for_grouping)
{
  # get unique experiments
  experiments <- unique(summary_df$Experiment)
  
  experiment <- "core_mu_prop_genes2_rate_genes1"
  for (experiment in experiments) {
    subset_df <- subset(summary_df, Experiment == experiment)
    
    fitted_params <- unique(subset_df$Fitted_param_long)
    fitted_param <- "Core mutation rate"
    for (fitted_param in fitted_params) {
      param_df <- subset(subset_df, Fitted_param_long == fitted_param)
      fitted_param_short <- unique(param_df$Fitted_param_short)
      param_df$Species <- factor(param_df$Species)
      param_df$Median <- as.numeric(param_df$Median)
      param_df$LCI <- as.numeric(param_df$LCI)
      param_df$UCI <- as.numeric(param_df$UCI)
      
      p <- ggplot(data = param_df, aes(x = reorder(species_short, Median), y=Median, colour = Life_Style)) + geom_point(size = 3) +
        geom_errorbar(data = param_df,
                      aes(ymin = LCI, ymax = UCI),
                      width = 0.2) +
        scale_y_log10() +
        theme_light() +
        theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
        theme(axis.text.x = element_text(size = 14), axis.text.y = element_text(size = 14), axis.title=element_text(size=16,face="bold"), strip.text.x = element_text(size = 12), 
              legend.text = element_text(size = 12), legend.title = element_text(size = 14, face="bold")) +
        xlab("Species") + ylab(fitted_param) + labs(color = "Life style")
      p
      
      if (grepl("all_params", experiment))
      {
        fname <- paste0(outpref, "exp_", "all_params", "_param_", fitted_param_short, ".png")
      } else if (grepl("prop_positive", experiment)) {
        fname <- paste0(outpref, "exp_", "selection", "_param_", fitted_param_short, ".png")
      } else {
        fname <- paste0(outpref, "exp_", experiment, "_param_", fitted_param_short, ".png")
      }
      
      ggsave(fname, p, width = 9, height = 7)
      
    }
    
    # plot grouped figure
    if (experiment == for_grouping) {
      
      
      # 1. Define your desired Life_Style order
      subset_df$Life_Style <- factor(
        subset_df$Life_Style,
        levels = c("Animal pathogen", "Commensal", "Free living")
      )
      
      # 2. Create the combined ordering key (Life_Style first, then species)
      subset_df <- subset_df %>%
        mutate(
          combo = paste(Life_Style, species_short, sep = "_")
        ) %>%
        arrange(Life_Style, species_short) %>%  # <-- this does the ordering logic
        mutate(
          Species_ordered = factor(combo, levels = unique(combo))
        )
      
      # 3. Apply that ordering to species_short (so you can plot species_short directly)
      species_levels <- subset_df %>%
        arrange(Life_Style, species_short) %>%
        distinct(species_short, .keep_all = TRUE) %>%
        pull(species_short)
      
      subset_df$species_short <- factor(subset_df$species_short, levels = species_levels)
      
      
      fitted_params <- unique(subset_df$Fitted_param_long)
      subset_df$Species <- factor(subset_df$Species)
      subset_df$Median <- as.numeric(subset_df$Median)
      subset_df$LCI <- as.numeric(subset_df$LCI)
      subset_df$UCI <- as.numeric(subset_df$UCI)
      
      p <- ggplot(data = subset_df, aes(x = species_short, y=Median, colour = Life_Style)) + geom_point(size = 3) +
        geom_errorbar(data = subset_df,
                      aes(ymin = LCI, ymax = UCI),
                      width = 0.2) +
        scale_y_log10() +
        facet_wrap( ~ Fitted_param_long, ncol=1, scales = "free_y", strip.position = "left") +
        scale_color_npg() +
        theme_light() +
        theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, face = "italic")) +
        theme(axis.text.x = element_text(size = 14), axis.text.y = element_text(size = 14), axis.title=element_text(size=16,face="bold"), strip.text.y = element_text(size = 16, colour = "black", face = "bold"), 
              legend.text = element_text(size = 12), legend.title = element_text(size = 14, face="bold"),
              strip.placement = "outside", strip.background = element_blank()) +
        xlab("Species") + ylab("") + labs(color = "Life style")
      p
      
      fname <- paste0(outpref, "exp_", experiment, "_merged")
      ggsave(paste0(fname, ".png"), p, width = 9, height = 10)
      ggsave(paste0(fname, ".svg"), p, width = 9, height = 10)
      
      write.csv(subset_df, file = paste0(fname, ".csv"), row.names = FALSE)
      
      # plot_list <- list()
      # 
      # num_params <- length(fitted_params)
      # counter <- 0
      # fitted_param <- "Core mutation rate"
      # for (fitted_param in fitted_params) {
      #   counter <- counter + 1
      #   param_df <- subset(subset_df, Fitted_param_long == fitted_param)
      #   fitted_param_short <- unique(param_df$Fitted_param_short)
      #   param_df$Species <- factor(param_df$Species)
      #   param_df$Median <- as.numeric(param_df$Median)
      #   param_df$LCI <- as.numeric(param_df$LCI)
      #   param_df$UCI <- as.numeric(param_df$UCI)
      #   
      #   p <- ggplot(data = param_df, aes(x = reorder(species_short, species_short), y=Median, colour = Life_Style)) + geom_point(size = 3) +
      #     geom_errorbar(data = param_df,
      #                   aes(ymin = LCI, ymax = UCI),
      #                   width = 0.2) +
      #     scale_y_log10() +
      #     theme_light() +
      #     theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
      #     theme(axis.text.x = element_text(size = 14), axis.text.y = element_text(size = 14), axis.title=element_text(size=16,face="bold"), strip.text.x = element_text(size = 12), 
      #           legend.text = element_text(size = 12), legend.title = element_text(size = 14, face="bold")) +
      #     xlab("Species") + ylab(fitted_param) + labs(color = "Life style")
      #   p
      #   
      #   if (counter < num_params)
      #   {
      #     p <- p + theme(axis.title.x=element_blank(),
      #                     axis.text.x=element_blank(),
      #                     axis.ticks.x=element_blank())
      #     p
      #   }
      #   
      #   plot_list[[counter]] <- p
      # }
      # 
      # #combine plots
      # group_plots <- ggarrange(
      #   plotlist=plot_list,
      #   ncol = 1,
      #   nrow = num_params,
      #   common.legend = TRUE,
      #   legend = "right",
      #   align = "v")
      # group_plots
    }
  }
}

indir <- "real_genome_fits/"
df_paths <- Sys.glob(paste0(indir,"data/*_mcmc_posterior_samples.csv"))
outpref <- paste0(indir, "figures/")
summary_df <- parse_results(df_paths, outpref)
summary_df <- rename_species(summary_df)
summary_df <- merge(summary_df, ne_dataframe, by.x = "species_short", by.y = "Species")

for_grouping = "core_mu_prop_genes2_rate_genes1"
plot_experiments(summary_df, outpref, for_grouping)

#merge results table with running parameters
fit.df <- read.csv(paste0(outpref, "exp_", for_grouping, "_merged.csv"))
params.df <- read.csv(paste0(outpref, "ppmod_real_input_with_coremu.txt"), sep = "\t")

fitted.params <- "core_mu"
subset.fit.df <- subset(fit.df, Fitted_param_short == fitted.params)
subset.fit.df <- subset.fit.df[c("Species", "Median")]
colnames(subset.fit.df) <- c("name", fitted.params)
merged.df <- merge(params.df, subset.fit.df)

write.table(merged.df, file = paste0(outpref, "ppmod_real_input_with_coremu.txt"), sep = "\t", quote = FALSE, row.names = FALSE)
