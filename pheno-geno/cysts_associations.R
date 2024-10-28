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

# load "merged_phe_geno.rdata" from local which contains genotype and phenotype cyst information for each patient
load("merged_phe_geno.rdata")

# merged_phe_geno is the data frame from merged_phe_geno.rdata

# shrink COL4A4 genotypes cuz the small sample size
makehom = c("COL4A4_C-NC_A-NC_hom", "COL4A4_C-NC_collagen_hom", "COL4A4_C-NC_hom", 
            "COL4A4_collagen_A-NC_hom", "COL4A4_collagen_hom")

merged_phe_geno$GT[merged_phe_geno$GT %in% makehom] = "COL4A4_hom"

# subset cyst and GT columns
test1 = merged_phe_geno[,c("cyst", "GT")]
# remove no-record lines
test1 = test1[which(test1$cysts != "\\"),]

test1$cysts = as.numeric(factor(test1$cysts, levels = c("normal", "type1", "type2")))
test1 = na.omit(test1)
test1 = test1[which(test1$GT != ""), ]
# remove genotypes with recorded patients fewer than 10
test1 %>% group_by(GT) %>% mutate(freq = n()) %>% filter(freq > 10) %>% as.data.frame -> test1
test1$GT = factor(test1$GT)

# poisson linked glm model
model1 <- glm(cysts ~ GT - 1, data = test1, family = "poisson")

summ1 = summary(model1)

# table drawing
g <- tableGrob(signif(summ1$coefficients[order(summ1$coefficients[,4]),c(1,4)],2), theme = ttheme_minimal())
separators <- replicate(nrow(g),
                        segmentsGrob(x0 = unit(0, "npc"), x1 = unit(1, "npc"),
                                     y0 = unit(0, "npc"), y1 = unit(0, "npc"),
                                     gp=gpar(lwd=2)),
                        simplify=FALSE)

ggg = signif(summ1$coefficients[order(summ1$coefficients[,4]),c(1,4)],2)
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
plot1 = test1 %>% mutate(cysts=as.factor(cysts)) %>% group_by(cysts,GT) %>%
  tally %>% group_by(GT) %>% mutate(freq = n/sum(n)) %>% 
  ggplot(aes(x = GT, y=freq)) + geom_bar(aes(fill=cysts), position="stack", stat="identity", width = 0.5) +
  theme_bw() + geom_text(aes(label=GT), y=0, angle=90, hjust=0) +
  theme(legend.title = element_blank(), axis.text.x = element_blank()) +
  theme(axis.text.y = element_text(size = 12, color="black")) + xlab("") + ylab("") +
  theme(plot.margin = margin(0,1,1,1, "cm")) + theme(axis.ticks.x = element_blank()) +
  scale_fill_brewer(palette = "Blues")

# AIC and statistics for GLM model
a1 = anova(model1, test = "LRT")
a1
stepAIC(model1)

# coefficients of each genotype
coeff = coefficients(model1); names(coeff) = gsub("GT", "", names(coeff))
test1$coef = coeff[match(test1$GT, names(coeff))]

plot2 = test1 %>% ggplot(aes(x = 1, y = GT, fill = coef)) + geom_tile(color="black") + geom_text(aes(label=round(coef,2)), size=4) +
  scale_fill_gradient2(low="white", mid = "#dbc5ca", high = "#dc143c", midpoint = 0.2) +
  theme_void() + theme(axis.title = element_blank(), axis.text = element_blank(), 
                       axis.ticks = element_blank(), plot.margin = margin(0.5,2.6,0,2.4, "cm")) + theme(legend.position = "none") +
  ggtitle(paste("P =", formatC(a1[5][2,1], format = "e", digits = 2))) + theme(plot.title = element_text(size = 12, face = "bold.italic")) + coord_flip()

ppp1 = summ1$coefficients[order(summ1$coefficients[,4]),c(1,4)][,2]
names(ppp1) = gsub("GT", "", names(ppp1))
test1$pval = ppp1[match(test1$GT, names(ppp1))]
plot3 = test1 %>% ggplot(aes(x = 1, y = GT, fill = -log10(pval))) + geom_tile(color="black") + geom_text(aes(label=signif(pval,2)), size=2) +
  scale_fill_fermenter(palette = "PRGn") +
  theme_void() + theme(axis.title = element_blank(), axis.text = element_blank(), 
                       axis.ticks = element_blank(), plot.margin = margin(0.5,2.6,0,2.4, "cm")) + theme(legend.position = "none") +
  ggtitle(paste("P =", formatC(a1[5][2,1], format = "e", digits = 2))) + theme(plot.title = element_text(size = 12, face = "bold.italic")) + coord_flip()

ggarrange(plot3, plot2, plot1, ncol=1, heights = c(1,1,6))



# mosaic plot for residuals
test2 = xtabs(~GT + cysts, data = test1)
mosaic(model1, ~ GT + cysts, residuals = residuals(model1), 
       gp = shading_hsv(observed = test2, h = c(0.98, 0.625), s=c(1, 0), 
                        v=c(1,0.5),lty = c("solid", "dashed"), 
                        interpolate = c(0,2,4,6,8,10), p.value = NA), 
       rot_labels=c(0,90,0,0), residuals_type = "rstandard", 
       labeling =  labeling_border(gp_labels = gpar(fontsize = 6), 
                                   rot_labels = c(0,90, 0, 0), 
                                   just_labels=c("center","right"), 
                                   tl_varnames = FALSE))




