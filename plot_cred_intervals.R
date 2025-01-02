library(ggplot2)
library(dplyr)

#df <- read.csv("results/predictions_negative_exp_3param_canberra_ngen200_hypparam_check_moderate_no_logit_BOLFI_parsed.txt", sep = "\t")
df <- read.csv("results/predictions_negative_exp_3param_canberra_ngen200_hypparam_check_moderate_no_logit_more_params_BOLFI_parsed.txt", sep = "\t")


df$core <- as.numeric(sub(".*core_(\\d+)_.*", "\\1", df$Name))
df$pan <- as.numeric(sub(".*pan_(\\d+)_.*", "\\1", df$Name))
df$freq <- as.numeric(sub(".*freq_([0-9.]+)_.*", "\\1", df$Name))
df$mu <- as.numeric(sub(".*mu_([0-9.]+)_.*", "\\1", df$Name))
df$prop <- as.numeric(sub(".*prop_([0-9.]+)_.*", "\\1", df$Name))
df$speed <- as.numeric(sub(".*speed_(\\d+)_.*", "\\1", df$Name))
df$initev <- as.numeric(sub(".*initev_(\\d+)_.*", "\\1", df$Name))
df$nevid <- as.numeric(sub(".*nevid_(\\d+)_.*", "\\1", df$Name))
df$updint <- as.numeric(sub(".*updint_(\\d+)_.*", "\\1", df$Name))
df$acqnoise <- as.numeric(sub(".*acqnoise_([0-9.]+)_.*", "\\1", df$Name))
df$thresh <- as.numeric(sub(".*thresh_([0-9.]+)_.*", "\\1", df$Name))

# split dfs
pan_mu_df <- subset(df, Param == "pan_mu")
prop_df <- subset(df, Param == "proportion")

pan_mu_df <- pan_mu_df %>%
  arrange(mu) 
pan_mu_df$Name <- factor(pan_mu_df$Name, levels = pan_mu_df$Name)
prop_df <- prop_df %>%
  arrange(prop)
prop_df$Name <- factor(prop_df$Name, levels = prop_df$Name)

# summary graphs
{
  p <- ggplot(pan_mu_df, aes(y = Name)) +
    # Horizontal range for Cred_2.5 to Cred_97.5
    geom_errorbarh(aes(xmin = Cred_2.5, xmax = Cred_97.5), height = 0.2) +
    # Point for Mean
    geom_point(aes(x = Mean), color = "black", size = 3) +
    # Conditional point for mu or prop with color based on range inclusion
    geom_point(
      aes(
        x = ifelse(Param == "pan_mu", mu, prop),
        color = ifelse((Param == "pan_mu" & mu >= Cred_2.5 & mu <= Cred_97.5) |
                         (Param == "proportion" & prop >= Cred_2.5 & prop <= Cred_97.5), 
                       "Within range", "Out of range")
      ),
      size = 3
    ) +
    # Define custom colors for the conditional points
    scale_color_manual(values = c("Within range" = "green", "Out of range" = "red")) +
    # Labels and theme
    labs(x = "Pan mu", y = "Simulation Name", color = "") +
    theme_light() +
    facet_grid()
  p
  
  ggsave("pan_mu_comparison.png")
  
  p <- ggplot(prop_df, aes(y = Name)) +
    # Horizontal range for Cred_2.5 to Cred_97.5
    geom_errorbarh(aes(xmin = Cred_2.5, xmax = Cred_97.5), height = 0.2) +
    # Point for Mean
    geom_point(aes(x = Mean), color = "black", size = 3) +
    # Conditional point for mu or prop with color based on range inclusion
    geom_point(
      aes(
        x = ifelse(Param == "pan_mu", mu, prop),
        color = ifelse((Param == "pan_mu" & mu >= Cred_2.5 & mu <= Cred_97.5) |
                         (Param == "proportion" & prop >= Cred_2.5 & prop <= Cred_97.5), 
                       "Within range", "Out of range")
      ),
      size = 3
    ) +
    # Define custom colors for the conditional points
    scale_color_manual(values = c("Within range" = "green", "Out of range" = "red")) +
    # Labels and theme
    labs(x = "Prop fast:slow genes", y = "Simulation Name", color = "") +
    theme_light()
  p
  ggsave("prop_comparison.png")
}


# param comp graphs
{
  p <- ggplot(pan_mu_df, aes(x = factor(acqnoise))) +
    # Horizontal range for Cred_2.5 to Cred_97.5
    geom_errorbar(aes(ymin = Cred_2.5, ymax = Cred_97.5), width = 0.2) +
    # Point for Mean
    geom_point(aes(y = Mean), color = "black", size = 3) +
    # Conditional point for mu or prop with color based on range inclusion
    geom_point(
      aes(
        y = ifelse(Param == "pan_mu", mu, prop),
        color = ifelse((Param == "pan_mu" & mu >= Cred_2.5 & mu <= Cred_97.5) |
                         (Param == "proportion" & prop >= Cred_2.5 & prop <= Cred_97.5), 
                       "Within range", "Out of range")
      ),
      size = 3
    ) +
    # Define custom colors for the conditional points
    scale_color_manual(values = c("Within range" = "green", "Out of range" = "red")) +
    # Labels and theme
    labs(y = "Pan mu", x = "Acquisition Noise", color = "") +
    theme_light() +
    facet_grid(mu~prop)
  p
  
  ggsave("pan_mu_comparison_acqnoise.png")
  
  p <- ggplot(prop_df, aes(x = factor(acqnoise))) +
    # Horizontal range for Cred_2.5 to Cred_97.5
    geom_errorbar(aes(ymin = Cred_2.5, ymax = Cred_97.5), width = 0.2) +
    # Point for Mean
    geom_point(aes(y = Mean), color = "black", size = 3) +
    # Conditional point for mu or prop with color based on range inclusion
    geom_point(
      aes(
        y = ifelse(Param == "pan_mu", mu, prop),
        color = ifelse((Param == "pan_mu" & mu >= Cred_2.5 & mu <= Cred_97.5) |
                         (Param == "proportion" & prop >= Cred_2.5 & prop <= Cred_97.5), 
                       "Within range", "Out of range")
      ),
      size = 3
    ) +
    # Define custom colors for the conditional points
    scale_color_manual(values = c("Within range" = "green", "Out of range" = "red")) +
    # Labels and theme
    labs(y = "Prop fast:slow genes", x = "Acquisition Noise", color = "") +
    theme_light() +
    facet_grid(mu~prop)
  p
  
  p <- ggplot(, aes(y = Name)) +
    # Horizontal range for Cred_2.5 to Cred_97.5
    geom_errorbarh(aes(xmin = Cred_2.5, xmax = Cred_97.5), height = 0.2) +
    # Point for Mean
    geom_point(aes(x = Mean), color = "black", size = 3) +
    # Conditional point for mu or prop with color based on range inclusion
    geom_point(
      aes(
        x = ifelse(Param == "pan_mu", mu, prop),
        color = ifelse((Param == "pan_mu" & mu >= Cred_2.5 & mu <= Cred_97.5) |
                         (Param == "proportion" & prop >= Cred_2.5 & prop <= Cred_97.5), 
                       "Within range", "Out of range")
      ),
      size = 3
    ) +
    # Define custom colors for the conditional points
    scale_color_manual(values = c("Within range" = "green", "Out of range" = "red")) +
    # Labels and theme
    labs(x = "Prop fast:slow genes", y = "Simulation Name", color = "") +
    theme_light()
  p
  ggsave("prop_comparison_acqnoise.png")
}

# for presentation
{
  sub_df <- subset(df, acqnoise == 0.01)
  sub_df$Mean <- round(sub_df$Mean, digits = 2)
  sub_df$Cred_2.5 <- round(sub_df$Cred_2.5, digits = 2)
  sub_df$Cred_97.5 <- round(sub_df$Cred_97.5, digits = 2)
  
  # convert 0.0 to 1.0
  sub_df$prop[sub_df$prop == 0.0] <- 1.0
  
  sub_df$Value = ifelse((sub_df$Param == "pan_mu"),
                        sub_df$mu, sub_df$prop)
  # for dodging
  sub_df$NonValue = ifelse((sub_df$Param == "pan_mu"),
                           sub_df$prop, sub_df$mu)
  
  sub_df$Param[sub_df$Param == "pan_mu"] <- "Mutation Rate"
  sub_df$Param[sub_df$Param == "proportion"] <- "Proportion Fast:Slow Genes"
  
  sub_df$in_range <- ifelse((sub_df$Param == "Mutation Rate" & sub_df$mu >= sub_df$Cred_2.5 & sub_df$mu <= sub_df$Cred_97.5) |
           (sub_df$Param == "Proportion Fast:Slow Genes" & sub_df$prop >= sub_df$Cred_2.5 & sub_df$prop <= sub_df$Cred_97.5), 
         "Within range", "Out of range")
  
  # Set the dodge width
  dodge_width <- 0.5
  
  p <- ggplot(sub_df, aes(x = factor(Value), group = NonValue)) +
    # Horizontal range for Cred_2.5 to Cred_97.5 with dodge
    geom_errorbar(
      aes(ymin = Cred_2.5, ymax = Cred_97.5), 
      width = 0.2, 
      position = position_dodge(width = dodge_width)
    ) +
    # Point for Mean with dodge
    geom_point(
      aes(y = Mean), 
      color = "black", 
      size = 3, 
      position = position_dodge(width = dodge_width)
    ) +
    # Conditional point for mu or prop with color based on range inclusion with dodge
    geom_point(
      aes(
        y = ifelse(Param == "Mutation Rate", mu, prop),
        color = in_range
      ),
      size = 3,
      position = position_dodge(width = dodge_width)
    ) +
    # Define custom colors for the conditional points
    scale_color_manual(values = c("Within range" = "green", "Out of range" = "red"),
                       limits=c("Within range", "Out of range"), drop = F) +
    # Labels and theme
    labs(y = "Parameter Estimate", x = "Known Parameter Value", color = "") +
    theme_light() +
    facet_wrap(.~Param, scales = "free_x") +
    theme(axis.text.x = element_text(size = 14), axis.text.y = element_text(size = 12), axis.title=element_text(size=16,face="bold"), strip.text.x = element_text(size = 12), legend.title=element_text(size=14,face="bold"), legend.text=element_text(size=12), legend.position = "none")
  p
  
  ggsave("acqnoise_0.01_simulation_comparisons.png")
}

# for presentation, large param comparison
{
  sub_df <- df
  sub_df$Mean <- round(sub_df$Mean, digits = 2)
  sub_df$Cred_2.5 <- round(sub_df$Cred_2.5, digits = 2)
  sub_df$Cred_97.5 <- round(sub_df$Cred_97.5, digits = 2)
  
  # convert 0.0 to 1.0
  sub_df$prop[sub_df$prop == 0.0] <- 1.0
  
  sub_df$Value = ifelse((sub_df$Param == "pan_mu"),
                        sub_df$mu, sub_df$prop)
  # for dodging
  sub_df$NonValue = ifelse((sub_df$Param == "pan_mu"),
                           sub_df$prop, sub_df$mu)
  
  sub_df$Param[sub_df$Param == "pan_mu"] <- "Mutation Rate"
  sub_df$Param[sub_df$Param == "proportion"] <- "Proportion Fast:Slow Genes"
  
  sub_df$in_range <- ifelse((sub_df$Param == "Mutation Rate" & sub_df$mu >= sub_df$Cred_2.5 & sub_df$mu <= sub_df$Cred_97.5) |
                              (sub_df$Param == "Proportion Fast:Slow Genes" & sub_df$prop >= sub_df$Cred_2.5 & sub_df$prop <= sub_df$Cred_97.5), 
                            "Within range", "Out of range")
  
  
  sub_df$prop_factor <- as.factor(sub_df$prop)
  sub_df$facet_group <- paste(paste("Core: ", sub_df$core, sep = ""), paste("Pan: ", sub_df$pan, sep = ""), paste("Freq: ", sub_df$freq, sep = ""), sep = "\n")
  sub_df$mu_factor <- paste("Mu: ", sub_df$mu, sep = "")
  
  sub_df$range_size <- sub_df$Cred_97.5 - sub_df$Cred_2.5
  
  
  # split up results
  pan_mu_df <- subset(sub_df, Param == "Mutation Rate")
  prop_df <- subset(sub_df, Param == "Proportion Fast:Slow Genes")
  
  
  # should be 14 per group, if not rerun
  prop_df %>% 
    group_by(facet_group) %>%
    summarise(no_rows = length(facet_group))
  
  # Set the dodge width
  dodge_width <- 0.5
  
  p <- ggplot(prop_df, aes(y = prop_factor)) +
    # Horizontal range for Cred_2.5 to Cred_97.5
    geom_errorbarh(aes(xmin = Cred_2.5, xmax = Cred_97.5), height = 0.2) +
    # Point for Mean
    geom_point(aes(x = Mean), color = "black", size = 3) +
    # Conditional point for mu or prop with color based on range inclusion with dodge
    geom_point(
      aes(
        x = ifelse(Param == "Mutation Rate", mu, prop),
        color = in_range
      ),
      size = 3,
      position = position_dodge(width = dodge_width)
    ) +
    # Define custom colors for the conditional points
    scale_color_manual(values = c("Within range" = "green", "Out of range" = "red"),
                       limits=c("Within range", "Out of range"), drop = F) +
    # Labels and theme
    labs(x = "Estimated Proportion Fast:Slow Genes", y = "Actual Proportion Fast:Slow Genes", color = "Actual Value") +
    theme_light() +
    facet_grid(facet_group~mu_factor) +
    theme(strip.text = element_text(color = "white", size = 12),
          strip.text.y.right = element_text(angle = 0),
          axis.title=element_text(size=14,face="bold"),
          axis.text=element_text(size=10),
          legend.title=element_text(size=12, face="bold"), 
          legend.text=element_text(size=12)) +
    scale_x_continuous(limits = c(0, 1))
  p
  
  ggsave("prop_gridsearch.png", width = 10, height = 12)
  
  p <- ggplot(prop_df, aes(y = prop_factor)) +
    geom_point(
      aes(
        x = range_size,
        color = in_range
      ),
      size = 3,
      position = position_dodge(width = dodge_width)
    ) +
    # Define custom colors for the conditional points
    scale_color_manual(values = c("Within range" = "green", "Out of range" = "red"),
                       limits=c("Within range", "Out of range"), drop = F) +
    # Labels and theme
    labs(x = "Credible Interval Size", y = "Actual Proportion Fast:Slow Genes", color = "Actual Value") +
    theme_light() +
    facet_grid(facet_group~mu_factor) +
    theme(strip.text = element_text(color = "white", size = 12),
          strip.text.y.right = element_text(angle = 0),
          axis.title=element_text(size=14,face="bold"),
          axis.text=element_text(size=10),
          legend.title=element_text(size=12, face="bold"), 
          legend.text=element_text(size=12)) +
    scale_x_continuous(limits = c(0, 1))
  p
  
  ggsave("prop_cred_interval_gridsearch.png", width = 10, height = 12)
  
  p <- ggplot(pan_mu_df, aes(y = prop_factor)) +
    # Horizontal range for Cred_2.5 to Cred_97.5
    geom_errorbarh(aes(xmin = Cred_2.5, xmax = Cred_97.5), height = 0.2) +
    # Point for Mean
    geom_point(aes(x = Mean), color = "black", size = 3) +
    # Conditional point for mu or prop with color based on range inclusion with dodge
    geom_point(
      aes(
        x = ifelse(Param == "Mutation Rate", mu, prop),
        color = in_range
      ),
      size = 3,
      position = position_dodge(width = dodge_width)
    ) +
    # Define custom colors for the conditional points
    scale_color_manual(values = c("Within range" = "green", "Out of range" = "red"),
                       limits=c("Within range", "Out of range"), drop = F) +
    # Labels and theme
    labs(x = "Estimated Mu", y = "Actual Proportion Fast:Slow Genes", color = "Actual Value") +
    theme_light() +
    facet_grid(facet_group~mu_factor) +
    theme(strip.text = element_text(color = "white", size = 12),
          strip.text.y.right = element_text(angle = 0),
          axis.title=element_text(size=14,face="bold"),
          axis.text=element_text(size=10),
          legend.title=element_text(size=12, face="bold"), 
          legend.text=element_text(size=12)) +
    scale_x_continuous(limits = c(0, 1))
  p
  ggsave("pan_mu_gridsearch.png", width = 10, height = 12)
  
  p <- ggplot(pan_mu_df, aes(y = prop_factor)) +
    geom_point(
      aes(
        x = range_size,
        color = in_range
      ),
      size = 3,
      position = position_dodge(width = dodge_width)
    ) +
    # Define custom colors for the conditional points
    scale_color_manual(values = c("Within range" = "green", "Out of range" = "red"),
                       limits=c("Within range", "Out of range"), drop = F) +
    # Labels and theme
    labs(x = "Credible Interval Size", y = "Actual Proportion Fast:Slow Genes", color = "Actual Value") +
    theme_light() +
    facet_grid(facet_group~mu_factor) +
    theme(strip.text = element_text(color = "white", size = 12),
          strip.text.y.right = element_text(angle = 0),
          axis.title=element_text(size=14,face="bold"),
          axis.text=element_text(size=10),
          legend.title=element_text(size=12, face="bold"), 
          legend.text=element_text(size=12)) +
    scale_x_continuous(limits = c(0, 1))
  p
  
  ggsave("pan_mu_cred_interval_gridsearch.png", width = 10, height = 12)
}


