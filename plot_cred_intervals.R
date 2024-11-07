library(ggplot2)
library(dplyr)

df <- read.csv("results/predictions_negative_exp_3param_canberra_ngen200_BOLFI_parsed.txt", sep = "\t")

df$core <- as.numeric(sub(".*core_(\\d+)_.*", "\\1", df$Name))
df$pan <- as.numeric(sub(".*pan_(\\d+)_.*", "\\1", df$Name))
df$freq <- as.numeric(sub(".*freq_([0-9.]+)_.*", "\\1", df$Name))
df$mu <- as.numeric(sub(".*mu_([0-9.]+)_.*", "\\1", df$Name))
df$prop <- as.numeric(sub(".*prop_([0-9.]+)_.*", "\\1", df$Name))
df$speed <- as.numeric(sub(".*speed_(\\d+)_.*", "\\1", df$teNamext))

# split dfs
pan_mu_df <- subset(df, Param == "pan_mu")
prop_df <- subset(df, Param == "proportion")

pan_mu_df <- pan_mu_df %>%
  arrange(mu) 
pan_mu_df$Name <- factor(pan_mu_df$Name, levels = pan_mu_df$Name)
prop_df <- prop_df %>%
  arrange(prop)
prop_df$Name <- factor(prop_df$Name, levels = prop_df$Name)

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
  theme_light()
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
