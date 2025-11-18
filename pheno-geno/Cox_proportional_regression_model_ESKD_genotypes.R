# ESKD survival analysis
# ESKD data depoists at folder test_data

library(ggfortify)
library(survminer)
library(ggpubr)
library(ggplot2)
library(survival)
library(xlsx)
library(dplyr)
library(cowplot)

# load data 
eskd <- read.table("./test_data/ESKD_data.xls", header = T, sep = "\t")[,1:3]
# columns: time status genotype
# status: 1 > normal; 2 > ESKD

# remove genotype with too few cases
eskd %>% group_by(genotype) %>% mutate(freq = n()) %>% 
                       filter(freq > 10) %>% as.data.frame -> eskd



# Cox proportional hazards regression models
cox1 <- coxph(Surv(time, status) ~ 1 + genotype, data = eskd)
cox1 <- summary(cox1)
print(cox1_sum$sctest[3]) # log-rank p
# P = 1.599728e-08
print(cox1_sum$coefficients[,1]) # HR
"""
genotypeCOL4A3_collagen_hom  genotypeCOL4A4_collagen_het           genotypeCOL4A4_hom 
                  0.93121161                   0.33043814                   0.03444387 
genotypeCOL4A5_collagen_hemi  genotypeCOL4A5_collagen_het      genotypeCOL4A5_NC1_hemi 
                  1.92766118                  -0.23594279                   2.06459783 
      genotypeCOL4A5_NC1_het 
                  1.02955536 
"""
cox1_newdata <- data.frame(genotype = cox1$xlevels, time = rep(mean(eskd$time), 2))
rownames(cox1_newdata) = cox1_newdata$genotype
cox1_fit <- survfit(cox1, newdata = cox1_newdata)

# strata is following the rownames of cox4_newdata

# plotting
ggsurvplot(
  cox1_fit,             
  data = eskd,         
  #risk.table = TRUE,     
  pval = TRUE,          
  conf.int = TRUE,     
  palette = "npg",
  xlim = c(0,65),         
  xlab = "Ages",  
  break.time.by = 10,    
  ggtheme = theme_cowplot(),
  #risk.table.y.text.col = T,
  #risk.table.height = 0.25, 
  #risk.table.y.text = FALSE,
  #ncensor.plot = TRUE,     
  #ncensor.plot.height = 0.25,
  conf.int.style = "step", 
  surv.median.line = "hv"
)









