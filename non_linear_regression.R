library(ggplot2)
library(data.table)
library(purrr)
library(parallel)
library(clue)
library(stringr)
library(MASS)
library(olsrr)
library(dplyr)
library(ggformula)
library(dbscan)
library(hash)
library(ggtrendline)
library(nlraa)
library(ggbeeswarm)
library(ggsci)
library("viridis")
library(plotly)
library(ggrepel)
library(ggfortify)
#library(nls2)
library(drc)
#library(devtools)
#install_github("onofriandreapg/aomisc")
library(aomisc)


files <- Sys.glob("./distances/samples/*.txt")
plot.fit <- TRUE
plot.coeff <- TRUE
sample.replace <- FALSE
sample.num <- 30000
constrained <- TRUE

# only for sample with replacement
if (sample.replace == TRUE)
{
  all_pairs <- data.frame(Species = c(), Sample = c(), Core = c(), Accessory = c())
} else {
  all_pairs <- data.frame(Species = c(), Core = c(), Accessory = c())
}

all.predictions <- data.frame(Mode = c(), Species = c(), Core = c(), Accessory = c())

for (j in 1:length(files))
{
  #j <-11
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
all_pairs$Species[all_pairs$Species == "lp"] <- "L. pneumophila"
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

unique(all_pairs$Species)


#mode_vec <- c("polynomial", "sigmoid", "quadplat", "asymp")
mode_vec <- c("asymp")
AIC.df <- data.frame(Mode = c(), Species = c(), AIC = c())
coeff.df <- data.frame(Mode = c(), Species = c(), Coeff = c(), Est = c(), Confhigh = c(), Conflow = c(), Sig = c(), row.names = NULL)

for (mode in mode_vec)
{
  all_species <- unique(all_pairs$Species)
  
  model.AIC <- 0
  model <- c()
  
  for (i in 1:length(all_species))
  {
    #i <- 1
    entry <- all_species[i]
    sample_points <- all_pairs[all_pairs$Species == entry,]
    p <- ggplot(sample_points, aes(x = Core, y = Accessory)) + theme_light() +
      geom_point(alpha=0.3, color="#2980B9") + xlab("Core distance (π)") + ylab("Accessory distance (a)")
    
    drop <- c("Species","Sample")
    sample_points = sample_points[,!(names(sample_points) %in% drop)]
    # generate evenly spaced points to predict from
    pred_points <- as.data.frame(seq(from = 0, to = max(sample_points$Core), by = max(sample_points$Core)/100))
    names(pred_points)[1] <- "Core"
    
    if (mode == "polynomial")
    {
      sample_points_temp <- sample_points
      pred_points_temp <- pred_points
      num_coeffs <- 2
      
      # add coefficients
      for (c in 2:num_coeffs)
      {
        name <- paste("Core", as.character(c), sep = "")
        sample_points_temp[name] <- sample_points_temp$Core ^ c
        pred_points_temp[name] <- pred_points_temp$Core ^ c
      }
      
      model <- lm(Accessory ~ ., data = sample_points_temp)
      
      
      # BIC
      #model.step <- stepAIC(model, k = log(length(model$coefficients)))
      # AIC
      model <- stepAIC(model, trace = FALSE)
      
      predicted_df <- data.frame(Accessory = predict(model, pred_points_temp), Core=pred_points_temp$Core)
      model.AIC <- AIC(model)
      
      #summary(model)
    } else if (mode == "sigmoid") {
      nlc <- nls.control(maxiter = 1000)
      model <- nls(Accessory ~ SSlogis(Core, Asym, xmid, scal), control=nlc, data = sample_points)
      
      model.AIC <- AIC(model)
      
      predicted_df <- data.frame(Accessory = predict(model, pred_points), Core=pred_points$Core)
    } else if (mode == "exp3p") {
      
      model <- nls(Accessory ~  SSexp3P(Core, a, b, c), data = sample_points, control = nls.control(maxiter = 1000))
      # tryCatch(
      #   expr = {
      #     model <- nls(Accessory ~  SSexp3P(Core, a, b, c), data = sample_points, control = nls.control(maxiter = 1000))
      #   },
      #   error = function(e) {
      #     message('Caught an error! Refitting with brute force')
      #     f <- function(x, a, b, c) {
      #       a * exp(b * x) + c
      #     }
      #     model <- nls2(Accessory ~ f(Core, a, b, c), start=list(a=0, b=-10000, c=0), data = sample_points, control = nls.control(maxiter = 1000), algorithm = "brute-force")
      #   }
      # )
      
      model.AIC <- AIC(model)
      
      predicted_df <- data.frame(Accessory = predict(model, pred_points), Core=pred_points$Core)
    } else if (mode == "asymp") {
      if (constrained == TRUE)
      {
        model <- drm(Accessory ~ Core, data = sample_points, fct = DRC.asymReg(), lowerl = c(0, 0, 0), upperl = c(1, Inf, 1))
      } else
      {
        model <- drm(Accessory ~ Core, data = sample_points, fct = DRC.asymReg())
      }
      
      model.AIC <- AIC(model)
      
      predicted_df <- data.frame(Accessory = predict(model, pred_points), Core=pred_points$Core)
    }
    else if (mode == "exp2p") {
      # remove all zero-values
      sample_points_temp <- sample_points[sample_points$Core > 0,]
      sample_points[sample_points$Core == 0,]
      sample_points_temp[sample_points_temp$Core == 0,]
      model <- nls(Accessory ~  a*exp(b*Core), data = sample_points_temp, start = list(a=0.01,b=0), control = nls.control(maxiter = 1000))
      
      model.AIC <- AIC(model)
      
      predicted_df <- data.frame(Accessory = predict(model, pred_points), Core=pred_points$Core)
    } else if (mode == "power3p") {
      # remove all zero-values
      sample_points_temp <- sample_points[sample_points$Core > 0,]
      sample_points[sample_points$Core == 0,]
      sample_points_temp[sample_points_temp$Core == 0,]
      model <- nls(Accessory ~  SSpower3P(Core, a, b, c), data = sample_points_temp, control = nls.control(maxiter = 10000, minFactor = 1/10240))

      model.AIC <- AIC(model)
      
      predicted_df <- data.frame(Accessory = predict(model, pred_points), Core=pred_points$Core)
    } else if (mode == "quadplat") {
      # remove all zero-values
      model <- nls(Accessory ~  SSquadp3xs(Core, a, b, jp), data = sample_points, control = nls.control(maxiter = 10000))

      model.AIC <- AIC(model)
      
      predicted_df <- data.frame(Accessory = predict(model, pred_points), Core=pred_points$Core)
    }
    
    sum.coeffs <- summary(model)$coefficients
    
    if (mode == "asymp")
    {
      rown <- str_split(rownames(sum.coeffs), ":")
      rown <- unlist(rown)[ c(TRUE,FALSE) ]
    } else {
      rown <- rownames(sum.coeffs)
    }
    
    coeff.df.temp <- data.frame(Mode = mode, Species = entry, Coeff = rown, Est = sum.coeffs[,1], Confhigh = sum.coeffs[,1] + 1.96*sum.coeffs[,2], Conflow = sum.coeffs[,1] - 1.96*sum.coeffs[,2], Sig = sum.coeffs[,4], row.names = NULL)
    
    
    if (mode == "asymp")
    {
      ceoeff.init.rate <- data.frame(Mode = mode, Species = coeff.df.temp[coeff.df.temp$Coeff == "m",]$Species, Coeff = "init_rate", Est = coeff.df.temp[coeff.df.temp$Coeff == "m",]$Est * (coeff.df.temp[coeff.df.temp$Coeff == "plateau",]$Est - coeff.df.temp[coeff.df.temp$Coeff == "init",]$Est), Confhigh = NA, Conflow = NA, Sig = 0)
      coeff.df.temp <- rbind(coeff.df.temp, ceoeff.init.rate)
    } else {
      rown <- rownames(sum.coeffs)
    }
    
    
    if (mode == "exp3p")
    {
      coeff.intercept <- data.frame(Mode = mode, Species = coeff.df.temp[coeff.df.temp$Coeff == "a",]$Species, Coeff = "Intercept", Est = coeff.df.temp[coeff.df.temp$Coeff == "a",]$Est + coeff.df.temp[coeff.df.temp$Coeff == "c",]$Est, Confhigh = coeff.df.temp[coeff.df.temp$Coeff == "a",]$Confhigh + coeff.df.temp[coeff.df.temp$Coeff == "c",]$Confhigh, Conflow = coeff.df.temp[coeff.df.temp$Coeff == "a",]$Conflow + coeff.df.temp[coeff.df.temp$Coeff == "c",]$Conflow, Sig = 0)
      coeff.df.temp <- rbind(coeff.df.temp, coeff.intercept)
      ceoeff.init.rate <- data.frame(Mode = mode, Species = coeff.df.temp[coeff.df.temp$Coeff == "a",]$Species, Coeff = "init_rate", Est = -coeff.df.temp[coeff.df.temp$Coeff == "b",]$Est * (coeff.df.temp[coeff.df.temp$Coeff == "c",]$Est - coeff.df.temp[coeff.df.temp$Coeff == "Intercept",]$Est), Confhigh = NA, Conflow = NA, Sig = 0)
      coeff.df.temp <- rbind(coeff.df.temp, ceoeff.init.rate)
    }
    
    
    coeff.df <- rbind(coeff.df, coeff.df.temp)
    
    AIC.df.temp <- data.frame(Mode = mode, Species = entry, AIC = model.AIC)
    AIC.df <- rbind(AIC.df, AIC.df.temp)
    
    predicted.append <- data.frame(Mode = mode, Species = entry, Core = predicted_df$Core, Accessory = predicted_df$Accessory)
    
    all.predictions <- rbind(all.predictions, predicted.append)
    
    if (plot.fit == TRUE)
    {
      p <- p + geom_line(color='black', linewidth=1.2, data = predicted_df, aes(x=Core, y=Accessory)) + ggtitle(paste(entry, "_", mode, sep=""))
      #p
      ggsave(paste(entry, "_", mode, "_fit.svg", sep=""), width=10, height = 6)
    }
  }
  
  #plot coeffs
  if (plot.coeff == TRUE)
  {
    sample_coeffs <- subset(coeff.df, Mode == mode)
    all_coeffs <- unique(sample_coeffs$Coeff)
    for (i in 1:length(all_coeffs))
    {
      c <- all_coeffs[[i]]
      sample_coeffs_temp <- subset(sample_coeffs, Coeff == c)
      # plot exp3p b parameter on -log scale
      if (mode == "exp3p" & c == "b") {
        sample_coeffs_temp$Est <- sample_coeffs_temp$Est * -1
        #p <- ggplot(sample_coeffs_temp, aes(x=reorder(Species, Est), y = Est, color=Species)) + theme_light() + geom_point() + xlab("Species") + ylab(paste(c,"(-log10)")) + theme(axis.text.x = element_text(angle = 45, size = 14, hjust = 1, vjust = 1), axis.text.y = element_text(size = 16), axis.title=element_text(size=16,face="bold"), strip.text.x = element_text(size = 12), legend.position="none") + scale_y_log10() # + geom_errorbar(aes(ymin=Confhigh, ymax=Conflow, color=Species, width=0.2))
        p <- ggplot(sample_coeffs_temp, aes(x=reorder(Species, Est), y = Est, color=Species)) + theme_light() + geom_point(size=4) + xlab("Species") + ylab(paste(c,"(-log10)")) + theme(axis.text.x = element_text(size = 18), axis.text.y = element_text(size = 18), axis.title=element_text(size=16,face="bold"), strip.text.x = element_text(size = 12), legend.position="none") + scale_y_log10() + coord_flip() # + geom_errorbar(aes(ymin=Confhigh, ymax=Conflow, color=Species, width=0.2))
      } else {
        #p <- ggplot(sample_coeffs_temp, aes(x=reorder(Species, Est), y = Est, color=Species)) + theme_light() + geom_point() + xlab("Species") + ylab(c) + theme(axis.text.x = element_text(angle = 45, size = 14, hjust = 1, vjust = 1), axis.text.y = element_text(size = 16), axis.title=element_text(size=16,face="bold"), strip.text.x = element_text(size = 12), legend.position="none") # + geom_errorbar(aes(ymin=Confhigh, ymax=Conflow, color=Species, width=0.2))
        p <- ggplot(sample_coeffs_temp, aes(x=reorder(Species, Est), y = Est, color=Species)) + theme_light() + geom_point(size=4) + xlab("Species") + ylab(c) + theme(axis.text.x = element_text(size = 18), axis.text.y = element_text(size = 18), axis.title=element_text(size=20,face="bold"), strip.text.x = element_text(size = 12), legend.position="none") + coord_flip() # + geom_errorbar(aes(ymin=Confhigh, ymax=Conflow, color=Species, width=0.2))
        if (mode == "asymp")
        {
          if (c == "m"){
            p <- p + ylab("Rate ratio (k)")
          } else if (c == "plateau"){
            p <- p + ylab("Plateau (i)")
          } else if (c == "init"){
            p <- p + ylab("Y-Intercept (j)")
          }
            
        }
      }
      p
      #plot.list[[i]] <- p
      ggsave(paste(mode, "_", c, "_comp.svg", sep=""), plot = p, height = 10, width = 8)
    }
  }
}


# plot interative views of predictions
for (mode in mode_vec){
  mode_sample <- subset(all.predictions, Mode == mode)
  p2 <- ggplot(mode_sample, aes(x = Core, y = Accessory, colour=Species)) +
    geom_line()+ theme_light() + guides(colour=guide_legend(title="Species")) + ylab("Accessory") + xlab("Core") + ggtitle(paste("Fit:", mode)) + scale_color_viridis(discrete = TRUE, option = "D") + scale_fill_viridis(discrete = TRUE) +
    theme(axis.text.x = element_text(size = 14), axis.text.y = element_text(size = 16), axis.title=element_text(size=16,face="bold"), strip.text.x = element_text(size = 12), legend.title=element_text(size=16), legend.text=element_text(size=14))
  q <- ggplotly(p2)
  show(q)
}

# plot paired plots of parameter estimates
for (mode in mode_vec){
  mode_sample <- subset(coeff.df, Mode == mode)
  
  # if exp3p, only look at intercept, initial rate and plataeu
  if (mode == "exp3p") {
    mode_sample <- subset(mode_sample, Coeff == "c" |  Coeff == "Intercept" |  Coeff == "init_rate")
  }
  
  all_coeffs <- unique(mode_sample$Coeff)
  merged_sample_coeffs <- data.frame(Species = as.factor(all_species))
  for (i1 in 1:length(all_coeffs))
  {
    c1 <- all_coeffs[[i1]]
    x <- data.frame(Species = subset(mode_sample, Coeff == c1)$Species, Est = subset(mode_sample, Coeff == c1)$Est)
    merged_sample_coeffs <- merge(merged_sample_coeffs, x, by = "Species", all.x=TRUE, all.y=TRUE)
    for (i2 in 1:length(all_coeffs))
    {
      if (i1 < i2) {
        
        c2 <- all_coeffs[[i2]]
        
        y <- data.frame(Species = subset(mode_sample, Coeff == c2)$Species, Est = subset(mode_sample, Coeff == c2)$Est)
        
        merged <- merge(x, y, by = "Species", all.x=TRUE, all.y=TRUE)
        
        if (plot.coeff == TRUE) {
          p <- ggplot(merged, aes(x=Est.x, y = Est.y, color=Species, label=Species)) + theme_light() + geom_point(size=2) + geom_text_repel(box.padding = 0.5) + ggtitle(mode) + xlab(c1) + ylab(c2) + theme(axis.text.x = element_text(size = 14), axis.text.y = element_text(size = 14), axis.title=element_text(size=16,face="bold"), strip.text.x = element_text(size = 12), legend.position="none")
          ggsave(paste(mode, "_", c1, "_vs_", c2, ".png", sep=""))
        }
      }
    }
  }
  # replace any NA with 0
  merged_sample_coeffs[is.na(merged_sample_coeffs)] <- 0
  pca_res <- prcomp(merged_sample_coeffs[-1], scale. = TRUE)
  p <- autoplot(pca_res, data = merged_sample_coeffs) + ggtitle(mode) + geom_point(aes(color = Species)) + geom_text_repel(aes(label = Species, color = Species))
  ggsave(paste(mode, "_PCA.png", sep=""))
}

if (plot.fit == TRUE)
{
  output_file = "all_params.txt"
  write.csv(coeff.df, output_file, row.names=FALSE, quote=FALSE) 
  output_file = "all_AIC.txt"
  write.csv(AIC.df, output_file, row.names=FALSE, quote=FALSE)
}

AIC.df.grouped <- AIC.df %>%
  group_by(Species) %>%
  summarise(prop = AIC/min(AIC), Mode = Mode, AIC = AIC) #%>%
mutate(prop = prop, Mode = Mode) %>% 
  arrange(desc(prop))

# plot comparisons of AIC between methods
p <- ggplot(AIC.df, aes(x = Species, y = prop, color = Mode)) + facet_wrap(~Species, scale="free") + theme_light() + ylab("Proportion of min. AIC") + geom_beeswarm(size=2, cex = 10) + theme(axis.text.x = element_blank(), axis.title.x=element_blank(), legend.title=element_text(size=16), legend.text=element_text(size=14)) + scale_color_npg()
p

# get number of top models by species
top.AIC <- AIC.df %>% group_by(Species) %>% top_n(1, -1 * AIC)
table(top.AIC$Mode)

rank.AIC <- AIC.df %>% group_by(Species) %>% mutate(rank = order(-1 * AIC, decreasing=TRUE))
rank.AIC <- rank.AIC %>% group_by(Mode) %>% summarise(avg_pos = mean(rank))
rank.AIC

# plot random species
species1 <- "M. abscessus"
species2 <- "M. tuberculosis"
plot.df <- subset(all_pairs, Species == species1 | Species == species2)

colour_vec <- c("M. abscessus" = "#2980B9", "M. tuberculosis" = "#A93226")

p <- ggplot(plot.df, aes(x = Core, y = Accessory, colour = Species)) + theme_light() + geom_point(alpha=0.3) + scale_colour_manual(values=colour_vec) + xlab("π") + ylab("a") + theme(axis.text.x = element_text(size = 12), axis.text.y = element_text(size = 12), axis.title=element_text(size=16,face="bold"), strip.text.x = element_text(size = 12), legend.title=element_text(size=14,face="bold"), legend.text=element_text(size=12)) + scale_y_continuous(limits = c(0,0.4)) + scale_x_continuous(limits = c(0, 0.035))
p
ggsave(file="Mtb_Ma_comparison.png", plot=p)

species1 <- "S. pneumoniae"
species2 <- "S. pyogenes"
plot.df <- subset(all_pairs, Species == species1 | Species == species2)

colour_vec <- c("S. pneumoniae" = "#28B463", "S. pyogenes" = "#D35400")

p <- ggplot(plot.df, aes(x = Core, y = Accessory, colour = Species)) + theme_light() + geom_point(alpha=0.3) + scale_colour_manual(values=colour_vec) + xlab("π") + ylab("a") + theme(axis.text.x = element_text(size = 12), axis.text.y = element_text(size = 12), axis.title=element_text(size=16,face="bold"), strip.text.x = element_text(size = 12), legend.title=element_text(size=14,face="bold"), legend.text=element_text(size=12)) + scale_y_continuous(limits = c(0,0.4)) + scale_x_continuous(limits = c(0, 0.035))
p
ggsave(file="Spneumo_GAS_comparison.png", plot=p)
