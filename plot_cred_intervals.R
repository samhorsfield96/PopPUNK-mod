library(ggplot2)
library(dplyr)

df <- read.csv("results/predictions_negative_exp_3param_canberra_ngen200_hypparam_check_moderate_no_logit_BOLFI_parsed.txt", sep = "\t")

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
    theme(legend.position="none")
  p
  
  ggsave("acqnoise_0.01_simulation_comparisons.png")
}

