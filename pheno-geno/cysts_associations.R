# Glm model for genotype association analysis with renal cysts phenotype
# Chinese AS cohort
# Author: Zhen Y

library(MASS)
library(vcdExtra)
library(tidyverse)
library(car)
library(RColorBrewer)
library(mosaic)
library(gridExtra)
library(ggpubr)
library(grid)
library(gtable)
library(ComplexHeatmap)
library(boot)

# load cysts genotype and phenotype data from test_data folder
cysts = read.table("./test_data/cyst_data.xls", sep = "\t", header = T)
cysts <- cysts[rep(seq_len(nrow(cysts)), cysts$patient.counts), , drop = FALSE]

# remove GT with few cases that lower the statitics power
cysts %>% group_by(GT) %>% mutate(freq = n()) %>% filter(freq > 10) %>% as.data.frame -> cysts

# reform phenotypic data; 1 for normal; 2 for 1-2 cysts; 3 for greater than 3 cysts
cysts$cysts = as.numeric(factor(cysts$cysts, levels = c("Normal", "1-2 cysts", "gt 3 cysts")))

# poisson linked glm model
model <- glm(cysts ~ GT - 1, data = cysts, family = "poisson")

summ = summary(model)

# table drawing
g <- tableGrob(signif(summ$coefficients[order(summ$coefficients[,4]),c(1,4)],2), theme = ttheme_minimal())
separators <- replicate(nrow(g),
                        segmentsGrob(x0 = unit(0, "npc"), x1 = unit(1, "npc"),
                                     y0 = unit(0, "npc"), y1 = unit(0, "npc"),
                                     gp=gpar(lwd=2)),
                        simplify=FALSE)

ggg = signif(summ$coefficients[order(summ$coefficients[,4]),c(1,4)],2)
rownames(ggg) = gsub("^GT", "", rownames(ggg))
g <- tableGrob(ggg, theme = ttheme_minimal())
separators <- replicate(nrow(g),
                        segmentsGrob(x0 = unit(0, "npc"), x1 = unit(1, "npc"),
                                     y0 = unit(0, "npc"), y1 = unit(0, "npc"),
                                     gp=gpar(lwd=2)),
                        simplify=FALSE)

g <- gtable::gtable_add_grob(g, grobs = separators, t = 1, r = ncol(g), b = seq_len(nrow(g)), l = 1)
grid.newpage()
grid.draw(g)

# plotting
plot1 = cysts %>% mutate(cysts=as.factor(cysts)) %>% group_by(cysts,GT) %>%
  tally %>% group_by(GT) %>% mutate(freq = n/sum(n)) %>% 
  ggplot(aes(x = GT, y=freq)) + geom_bar(aes(fill=cysts), position="stack", stat="identity", width = 0.5) +
  theme_bw() + geom_text(aes(label=GT), y=0, angle=90, hjust=0) +
  theme(legend.title = element_blank(), axis.text.x = element_blank()) +
  theme(axis.text.y = element_text(size = 12, color="black")) + xlab("") + ylab("") +
  theme(plot.margin = margin(0,1,1,1, "cm")) + theme(axis.ticks.x = element_blank()) +
  scale_fill_brewer(palette = "Blues")

# AIC and statistics for GLM model
a1 = anova(model, test = "LRT")
a1
stepAIC(model)

# coefficients of each genotype
coeff = coefficients(model); names(coeff) = gsub("GT", "", names(coeff))
cysts$coef = coeff[match(cysts$GT, names(coeff))]

plot2 = cysts %>% ggplot(aes(x = 1, y = GT, fill = coef)) + geom_tile(color="black") + geom_text(aes(label=round(coef,2)), size=4) +
  scale_fill_gradient2(low="white", mid = "#dbc5ca", high = "#dc143c", midpoint = 0.2) +
  theme_void() + theme(axis.title = element_blank(), axis.text = element_blank(), 
                       axis.ticks = element_blank(), plot.margin = margin(0.5,2.6,0,2.4, "cm")) + theme(legend.position = "none") +
  ggtitle(paste("P =", formatC(a1[5][2,1], format = "e", digits = 2))) + theme(plot.title = element_text(size = 12, face = "bold.italic")) + coord_flip()

ppp1 = summ$coefficients[order(summ$coefficients[,4]),c(1,4)][,2]
names(ppp1) = gsub("GT", "", names(ppp1))
cysts$pval = ppp1[match(cysts$GT, names(ppp1))]
plot3 = cysts %>% ggplot(aes(x = 1, y = GT, fill = -log10(pval))) + geom_tile(color="black") + geom_text(aes(label=signif(pval,2)), size=2) +
  scale_fill_fermenter(palette = "PRGn") +
  theme_void() + theme(axis.title = element_blank(), axis.text = element_blank(), 
                       axis.ticks = element_blank(), plot.margin = margin(0.5,2.6,0,2.4, "cm")) + theme(legend.position = "none") +
  ggtitle(paste("P =", formatC(a1[5][2,1], format = "e", digits = 2))) + theme(plot.title = element_text(size = 12, face = "bold.italic")) + coord_flip()

ggarrange(plot3, plot2, plot1, ncol=1, heights = c(1,1,6))




