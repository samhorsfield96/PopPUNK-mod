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
#library(aomisc)
library(statforbiology)
library(readxl)
library(scales)
library(ggpubr)
library(rstatix)


files <- Sys.glob("./real_genome_distances/data/*.txt")
outdir <- "./real_genome_distances/figures/"
plot.fit <- TRUE
plot.coeff <- TRUE
sample.replace <- FALSE
sample.num <- 1000000000
constrained <- TRUE
resample <- FALSE


if (resample == FALSE)
{
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
  all_pairs$Abbreviation <- all_pairs$Species
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
  #write.csv(all_pairs, "popPUNK_random_30K_sample.csv", row.names=FALSE)
} else
{
  all_pairs <- read.csv("popPUNK_random_30K_sample.csv")
  all.predictions <- data.frame(Mode = c(), Species = c(), Core = c(), Accessory = c())
}

# get effective populatiom size from universal genes
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
ne_dataframe$dNdS <- as.numeric(ne_dataframe$dNdS)
ne_dataframe$Phi <- as.numeric(ne_dataframe$Phi)
ne_dataframe$Adjusted_pangenome <- as.numeric(ne_dataframe$Adjusted_pangenome)

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
  
  i <- 1
  for (i in 1:length(all_species))
  {
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
      p <- p + geom_line(color='black', linewidth=1.2, data = predicted_df, aes(x=Core, y=Accessory))# + ggtitle(paste(entry, "_", mode, sep=""))
      #p
      ggsave(paste(outdir, entry, "_", mode, "_fit.jpg", sep=""), width=10, height = 6)
    }
  }
  
  #plot coeffs
  if (plot.coeff == TRUE)
  {
    spearmans.rho.bobay <- data.frame(Coeff = c(), Comparator = c(), Rho = c(), P.value = c())
    sample_coeffs <- subset(coeff.df, Mode == mode)
    sample_coeffs <- merge(sample_coeffs, ne_dataframe, by="Species")
    my_comparisons <- list( c("Animal pathogen", "Commensal"), c("Animal pathogen", "Free living"), c("Commensal", "Free living") )
    
    all_coeffs <- unique(sample_coeffs$Coeff)
    i <- 2
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
        p <- ggplot(sample_coeffs_temp, aes(x=reorder(Species, Est), y = Est, color=Life_Style)) + theme_light() + geom_point(size=4) + xlab("Species") + ylab(c) + theme(axis.text.x = element_text(size = 18), axis.text.y = element_text(size = 18), axis.title=element_text(size=20,face="bold"), strip.text.x = element_text(size = 12), legend.title=element_text(size=20,face="bold"), legend.text=element_text(size=16)) + guides(colour=guide_legend(title="Lifestyle")) + coord_flip() # + geom_errorbar(aes(ymin=Confhigh, ymax=Conflow, color=Species, width=0.2))
        pNe <- ggplot(sample_coeffs_temp, aes(x=Ne, y = Est, color=Life_Style)) + theme_light() + geom_point(size=4) + xlab("Effective population size (Ne)") + ylab(c) + theme(axis.text.x = element_text(size = 18), axis.text.y = element_text(size = 18), axis.title=element_text(size=20,face="bold"), strip.text.x = element_text(size = 12), legend.title=element_text(size=20,face="bold"), legend.text=element_text(size=16)) + guides(colour=guide_legend(title="Lifestyle")) + scale_x_continuous(labels = label_scientific(digits = 2))# + coord_flip() # + geom_errorbar(aes(ymin=Confhigh, ymax=Conflow, color=Species, width=0.2))
        prm <- ggplot(sample_coeffs_temp, aes(x=rm, y = Est, color=Life_Style)) + theme_light() + geom_point(size=4) + xlab("r/m") + ylab(c) + theme(axis.text.x = element_text(size = 18), axis.text.y = element_text(size = 18), axis.title=element_text(size=20,face="bold"), strip.text.x = element_text(size = 12), legend.title=element_text(size=20,face="bold"), legend.text=element_text(size=16)) + guides(colour=guide_legend(title="Lifestyle"))# + coord_flip() # + geom_errorbar(aes(ymin=Confhigh, ymax=Conflow, color=Species, width=0.2))
        pdnds <- ggplot(sample_coeffs_temp, aes(x=dNdS, y = Est, color=Life_Style)) + theme_light() + geom_point(size=4) + xlab("dN/dS") + ylab(c) + theme(axis.text.x = element_text(size = 18), axis.text.y = element_text(size = 18), axis.title=element_text(size=20,face="bold"), strip.text.x = element_text(size = 12), legend.title=element_text(size=20,face="bold"), legend.text=element_text(size=16)) + guides(colour=guide_legend(title="Lifestyle"))# + coord_flip() # + geom_errorbar(aes(ymin=Confhigh, ymax=Conflow, color=Species, width=0.2))
        ppan <- ggplot(sample_coeffs_temp, aes(x=Adjusted_pangenome, y = Est, color=Life_Style)) + theme_light() + geom_point(size=4) + xlab("Pangenome size (no. genes)") + ylab(c) + theme(axis.text.x = element_text(size = 18), axis.text.y = element_text(size = 18), axis.title=element_text(size=20,face="bold"), strip.text.x = element_text(size = 12), legend.title=element_text(size=20,face="bold"), legend.text=element_text(size=16)) + guides(colour=guide_legend(title="Lifestyle"))# + coord_flip() # + geom_errorbar(aes(ymin=Confhigh, ymax=Conflow, color=Species, width=0.2))
        pbox <- ggplot(sample_coeffs_temp, aes(x=Life_Style, y = Est, color=Life_Style)) + theme_light() + geom_boxplot(outlier.shape = NA) + geom_jitter(size=2, alpha=0.4) + xlab("Lifestyle") + ylab(c) + theme(axis.text.x = element_text(size = 18), axis.text.y = element_text(size = 18), axis.title=element_text(size=20,face="bold"), strip.text.x = element_text(size = 12), legend.position = "none") + stat_compare_means(comparisons = my_comparisons, method = "wilcox.test", aes(label = ..p.adj..), p.adjust.method = "bonferroni")
        
        if (mode == "asymp")
        {
          if (c == "m"){
            p <- p + ylab("Rate ratio (k)") + scale_y_log10()
            pNe <- pNe + ylab("Rate ratio (k)") + geom_smooth(method = "lm", se=FALSE, color="black", formula = y ~ x)
            prm <- prm + ylab("Rate ratio (k)") + geom_smooth(method = "lm", se=FALSE, color="black", formula = y ~ x)
            pdnds <- pdnds + ylab("Rate ratio (k)") + geom_smooth(method = "lm", se=FALSE, color="black", formula = y ~ x)
            ppan <- ppan + ylab("Rate ratio (k)") + geom_smooth(method = "lm", se=FALSE, color="black", formula = y ~ x)
            pbox <- pbox + ylab("Rate ratio (k)")
          } else if (c == "plateau"){
            p <- p + ylab("Plateau (i)")
            pNe <- pNe + ylab("Plateau (i)") + geom_smooth(method = "lm", se=FALSE, color="black", formula = y ~ x)
            prm <- prm + ylab("Plateau (i)") + geom_smooth(method = "lm", se=FALSE, color="black", formula = y ~ x)
            pdnds <- pdnds + ylab("Plateau (i)") + geom_smooth(method = "lm", se=FALSE, color="black", formula = y ~ x)
            ppan <- ppan + ylab("Plateau (i)") + geom_smooth(method = "lm", se=FALSE, color="black", formula = y ~ x)
            pbox <- pbox + ylab("Plateau (i)")
          } else if (c == "init"){
            p <- p + ylab("Y-Intercept (j)")
            pNe <- pNe + ylab("Y-Intercept (j)") + geom_smooth(method = "lm", se=FALSE, color="black", formula = y ~ x)
            prm <- prm + ylab("Y-Intercept (j)") + geom_smooth(method = "lm", se=FALSE, color="black", formula = y ~ x)
            pdnds <- pdnds + ylab("Y-Intercept (j)") + geom_smooth(method = "lm", se=FALSE, color="black", formula = y ~ x)
            ppan <- ppan + ylab("Y-Intercept (j)") + geom_smooth(method = "lm", se=FALSE, color="black", formula = y ~ x)
            pbox <- pbox + ylab("Y-Intercept (j)")
          }
          #ggsave(paste(outdir, mode, "_", c, "_comp_Ne.svg", sep=""), plot = pNe, height = 10, width = 12)
          #ggsave(paste(outdir, mode, "_", c, "_comp_rm.svg", sep=""), plot = prm, height = 10, width = 12)
          #ggsave(paste(outdir, mode, "_", c, "_comp_dNdS.svg", sep=""), plot = pdnds, height = 10, width = 12)
          #ggsave(paste(outdir, mode, "_", c, "_comp_pan_size.svg", sep=""), plot = ppan, height = 10, width = 12)
          
          #combined.plot <- ggarrange(pNe, prm, ppan, labels = "AUTO", align = "hv", common.legend = TRUE, nrow = 3, legend = "right", font.label = list(size = 20))
          #ggsave(paste(outdir, mode, "_", c, "_comp_combined.svg", sep=""), plot = combined.plot, height = 16, width = 10)
          ggsave(paste(outdir, mode, "_", c, "_boxplot.svg", sep=""), plot = pbox, height = 6, width = 8)
          #ggsave(paste(outdir, mode, "_", c, "_comp_combined.png", sep=""), plot = combined.plot, height = 16, width = 10)
          ggsave(paste(outdir, mode, "_", c, "_boxplot.png", sep=""), plot = pbox, height = 6, width = 8)
          
          
          # calculate spearman's rho, used to test monotonic relationship
          
          corr <- cor.test(x=sample_coeffs_temp$Est, y=sample_coeffs_temp$Ne, method = 'spearman')
          spearmans.rho.bobay <- rbind(spearmans.rho.bobay, data.frame(Coeff = c, Comparator = "Ne", Rho = corr$statistic, P.value = corr$p.value))
          corr <- cor.test(x=sample_coeffs_temp$Est, y=sample_coeffs_temp$rm, method = 'spearman')
          spearmans.rho.bobay <- rbind(spearmans.rho.bobay, data.frame(Coeff = c, Comparator = "rm", Rho = corr$statistic, P.value = corr$p.value))
          corr <- cor.test(x=sample_coeffs_temp$Est, y=sample_coeffs_temp$Adjusted_pangenome, method = 'spearman')
          spearmans.rho.bobay <- rbind(spearmans.rho.bobay, data.frame(Coeff = c, Comparator = "Pan_size", Rho = corr$statistic, P.value = corr$p.value))
          corr <- cor.test(x=sample_coeffs_temp$Est, y=sample_coeffs_temp$Phi, method = 'spearman')
          spearmans.rho.bobay <- rbind(spearmans.rho.bobay, data.frame(Coeff = c, Comparator = "Phi", Rho = corr$statistic, P.value = corr$p.value))
          
        }
      }
      p
      #plot.list[[i]] <- p
      ggsave(paste(outdir, mode, "_", c, "_comp.svg", sep=""), plot = p, height = 10, width = 8)
      ggsave(paste(outdir, mode, "_", c, "_comp.png", sep=""), plot = p, height = 10, width = 8)
    }
  }
  spearmans.rho.bobay$P.value.adj <- p.adjust(spearmans.rho.bobay$P.value, method="bonferroni")
  write.csv(spearmans.rho.bobay, "spearman_rho_bobay_comp.csv", row.names=FALSE)
}

