library(ggplot2)
library(dplyr)
library(ggsci)

search_string <- function(filename, string){
  search_rg <- paste0(".*", string, "_(-?[0-9.]+(?:[eE][-+]?[0-9]+)?)(?:_|\\.|$).*")
  label = if (grepl(string, filename)) {
    as.numeric(sub(search_rg, "\\1", filename))
  } else {
    0
  }
  return(label)
}

parse_filename <- function(filename) {
  list(
    prop_positive = search_string(filename, "prop_positive"),
    HR_rate = search_string(filename, "HR_rate"),
    HGT_rate = search_string(filename, "HGT_rate"),
    rate_genes1 = search_string(filename, "r1"),
    rate_genes2 = search_string(filename, "r2"),
    prop_genes2 = search_string(filename, "prop2"),
    core_mu = search_string(filename, "core_mu"),
    competition = search_string(filename, "competition"),
    pos_lambda = search_string(filename, "pos_lambda"),
    neg_lambda = search_string(filename, "neg_lambda")
  )
}

parse_results <- function(df_paths)
{
  df_all <- data.frame(Est_parameter = c(), Mean = c(), Median = c(), UCI = c(), LCI = c(), prop_positive = c(), HR_rate = c(), HGT_rate = c(), rate_genes1 = c(), prop_genes2 = c(), competition = c())
  params_all <- data.frame(prop_positive = c(), HR_rate = c(), HGT_rate = c(), rate_genes1 = c(), prop_genes2 = c(), competition = c(), baseline = c())
  
  j <- 1
  for (j in 1:length(df_paths))
  {
    filename <- df_paths[j]
    df <- read.table(filename, sep = ",", comment.char = "", header=TRUE)
    colnames(df) <- c("Est_parameter", "Mean", "Median", "LCI", "UCI")
    
    # add metadata
    parsed_values <- parse_filename(filename)
    for (name in names(parsed_values)) {
      df[[name]] <- parsed_values[[name]]
    }
    
    df_all <- rbind(df_all, df)
  }
  
  return(df_all)
}

df_paths <- Sys.glob("credible_intervals/RBF_KDE_smoothing_0.25/*.csv")
outpref <- "credible_intervals"
df_all <- parse_results(df_paths)
df_all$rate_genes1 <- factor(df_all$rate_genes1)
df_all$prop_genes2 <- factor(df_all$prop_genes2)
df_all$core_mu <- factor(df_all$core_mu)

# summary graphs
{
  df_subset <- subset(df_all, Est_parameter == "rate_genes1" & core_mu == 1e-5)
  df_subset$real_value <- as.numeric(as.character(df_subset$rate_genes1))
  df_subset$in_range <-  ifelse((df_subset$real_value >= df_subset$ LCI & df_subset$real_value <= df_subset$UCI),
                             "Within range", "Out of range")
  p <- ggplot(df_subset, aes(y = rate_genes1)) +
    # Horizontal range for Cred_2.5 to Cred_97.5
    geom_errorbarh(aes(xmin = LCI, xmax = UCI), height = 0.2) +
    # Point for Mean
    geom_point(aes(x = Median), color = "black", size = 3) +
    # Conditional point for mu or prop with color based on range inclusion
    geom_point(
      aes(x = real_value, color = ifelse((real_value >= LCI & real_value <= UCI),
                       "Within range", "Out of range")
      ),
      size = 3
    ) +
    # Define custom colors for the conditional points
    scale_color_manual(values = c("Within range" = "green", "Out of range" = "red")) +
    scale_x_log10() +
    # Labels and theme
    facet_grid(prop_genes2 ~ .) +
    labs(x = "Gene turnover rate 1", y = "Real value", color = "") +
    theme_light()
  p
  
  df_overall <- df_subset
  ggsave(paste(outpref, "gene_rate1.png", sep = "/"))
  
  df_subset <- subset(df_all, Est_parameter == "prop_genes2" & core_mu == 1e-5)
  df_subset$real_value <- as.numeric(as.character(df_subset$prop_genes2))
  df_subset$in_range <-  ifelse((df_subset$real_value >= df_subset$ LCI & df_subset$real_value <= df_subset$UCI),
                                "Within range", "Out of range")
  p <- ggplot(df_subset, aes(y = prop_genes2)) +
    # Horizontal range for Cred_2.5 to Cred_97.5
    geom_errorbarh(aes(xmin = LCI, xmax = UCI), height = 0.2) +
    # Point for Mean
    geom_point(aes(x = Median), color = "black", size = 3) +
    # Conditional point for mu or prop with color based on range inclusion
    geom_point(
      aes(x = real_value, color = ifelse((real_value >= LCI & real_value <= UCI),
                                         "Within range", "Out of range")
      ),
      size = 3
    ) +
    # Define custom colors for the conditional points
    scale_color_manual(values = c("Within range" = "green", "Out of range" = "red")) +
    # Labels and theme
    facet_grid(rate_genes1 ~ .) +
    labs(x = "Proportion of fast genes", y = "Real value", color = "") +
    theme_light()
  p
  
  ggsave(paste(outpref, "prop2_genes.png", sep = "/"))
  df_overall <- rbind(df_overall, df_subset)
  df_overall$range_size <- df_overall$UCI - df_overall$LCI
}

# bar charts
{
  # get sensitivity
  # overall
  {
    grouped_df <- df_overall %>% 
      group_by(Est_parameter) %>%
      summarize(
        Perc_correct = mean(in_range == "Within range"), 
        range_size = mean(range_size)
      )

    p <- ggplot(grouped_df, aes(x = Est_parameter, y = Perc_correct, fill = Est_parameter)) + 
      geom_col() +
      scale_y_continuous(limits = c(0, 1)) + xlab("Parameter") + ylab("Average Sensitivity") + 
      theme_light() + 
      scale_fill_npg() +
      theme(strip.text = element_text(color = "white", size = 12),
            #strip.text.y.right = element_text(angle = 0),
            axis.title=element_text(size=14,face="bold"),
            axis.text=element_text(size=10),
            legend.title=element_text(size=12, face="bold"), 
            legend.text=element_text(size=12),
            legend.position = "None") +
      labs(fill = "Population size")
    p
    ggsave(paste(outpref,"overall_sensitivity.png", sep="/"), width = 9, height = 6)
    
    p <- ggplot(grouped_df, aes(x = Est_parameter, y = range_size, fill = Est_parameter)) + 
      geom_col()  + xlab("Parameter") + ylab("Average Credible Interval Size") + 
      theme_light() + 
      scale_fill_npg() +
      theme(strip.text = element_text(color = "white", size = 12),
            #strip.text.y.right = element_text(angle = 0),
            axis.title=element_text(size=14,face="bold"),
            axis.text=element_text(size=10),
            legend.title=element_text(size=12, face="bold"), 
            legend.text=element_text(size=12),
            legend.position = "None")
    p
    ggsave(paste(outpref,"overall_cred_interval_size.png", sep="/"), width = 9, height = 6)
  }
  
  # by rate genes 1
  {
    grouped_df <- df_overall %>% 
      group_by(rate_genes1) %>%
      summarize(
        Perc_correct = mean(in_range == "Within range"), 
        range_size = mean(range_size)
      )
    
    p <- ggplot(grouped_df, aes(x = rate_genes1, y = Perc_correct, fill = rate_genes1)) + 
      geom_col() +
      scale_y_continuous(limits = c(0, 1)) + xlab("Gene turnover rate 1") + ylab("Average Sensitivity") + 
      theme_light() + 
      scale_fill_npg() +
      theme(strip.text = element_text(color = "white", size = 12),
            #strip.text.y.right = element_text(angle = 0),
            axis.title=element_text(size=14,face="bold"),
            axis.text=element_text(size=10),
            legend.title=element_text(size=12, face="bold"), 
            legend.text=element_text(size=12),
            legend.position = "None")
    p
    ggsave(paste(outpref,"rate_genes1_sensitivity.png", sep="/"), width = 9, height = 6)
    
    p <- ggplot(grouped_df, aes(x = rate_genes1, y = range_size, fill = rate_genes1)) + 
      geom_col()+ xlab("Gene turnover rate 1") + ylab("Average Credible Interval Size") + 
      theme_light() + 
      scale_fill_npg() +
      theme(strip.text = element_text(color = "white", size = 12),
            #strip.text.y.right = element_text(angle = 0),
            axis.title=element_text(size=14,face="bold"),
            axis.text=element_text(size=10),
            legend.title=element_text(size=12, face="bold"), 
            legend.text=element_text(size=12),
            legend.position = "None")
    p
    ggsave(paste(outpref,"rate_genes1_cred_interval_size.png", sep="/"), width = 9, height = 6)
  }
  
  # by prop 2
  {
    grouped_df <- df_overall %>% 
      group_by(prop_genes2) %>%
      summarize(
        Perc_correct = mean(in_range == "Within range"), 
        range_size = mean(range_size)
      )
    
    p <- ggplot(grouped_df, aes(x = prop_genes2, y = Perc_correct, fill = prop_genes2)) + 
      geom_col() +
      scale_y_continuous(limits = c(0, 1)) + xlab("Proportion of fast genes") + ylab("Average Sensitivity") + 
      theme_light() + 
      scale_fill_npg() +
      theme(strip.text = element_text(color = "white", size = 12),
            #strip.text.y.right = element_text(angle = 0),
            axis.title=element_text(size=14,face="bold"),
            axis.text=element_text(size=10),
            legend.title=element_text(size=12, face="bold"), 
            legend.text=element_text(size=12),
            legend.position = "None")
    p
    ggsave(paste(outpref,"prop_genes2_sensitivity.png", sep="/"), width = 9, height = 6)
    
    p <- ggplot(grouped_df, aes(x = prop_genes2, y = range_size, fill = prop_genes2)) + 
      geom_col() + xlab("Proportion of fast genes") + ylab("Average Credible Interval Size") + 
      theme_light() + 
      scale_fill_npg() +
      theme(strip.text = element_text(color = "white", size = 12),
            #strip.text.y.right = element_text(angle = 0),
            axis.title=element_text(size=14,face="bold"),
            axis.text=element_text(size=10),
            legend.title=element_text(size=12, face="bold"), 
            legend.text=element_text(size=12),
            legend.position = "None")
    p
    ggsave(paste(outpref,"prop_genes2_cred_interval_size.png", sep="/"), width = 9, height = 6)
  }
}
