library(rstan)
library(brms)
library(bayesrules)
library(ggsci)
library(cowplot)
library(ggdist)
library(tidybayes)
# hearing

# generate prior
moment_estimation <- function(params) {
  alpha <- params[1]
  beta <- params[2]
  
  # Theoretical moments
  theoretical_mean <- alpha / (alpha + beta)
  theoretical_var <- (alpha * beta) / ((alpha + beta)^2 * (alpha + beta + 1))
  
  # Sample moments
  sample_mean <- m / N
  sample_var <- (m * (N - m)) / (N^2 * (N - 1))
  
  # Objective function: sum of squared differences between theoretical and sample moments
  obj <- (theoretical_mean - sample_mean)^2 + (theoretical_var - sample_var)^2
  
  return(obj)
}



hear_df = merged_phe_geno[,c("hearing", "GT", "Gene", "Age")]

hear_df = na.omit(hear_df)
hear_df = hear_df[which(hear_df$GT != ""), ]
hear_df %>% group_by(GT) %>% mutate(freq = n()) %>% filter(freq > 10) %>% as.data.frame -> hear_df
hear_df$GT = factor(hear_df$GT)
hear_df$Gene = factor(hear_df$Gene)

N <- 274
m <- 104


formula <- bf(hearing | trials(274) ~ Age + (1 | GT))
#get_prior(formula = formula,
#          data = hear_df,
#          family = binomial(link = "logit"))

# non-informative prior
model <- brm(
    formula,
    family = binomial(link = "logit"),
    data = hear_df,
    prior = prior(beta(1, 1), class = b, lb = 0, ub = 1),
    chains = 4, warmup = 1000, iter = 2000, seed = 123,
  refresh = 0
)



model |> tidy_draws() |> pivot_longer(cols = starts_with("r_GT")) |> select(name, value) |> 
mutate(name = str_remove(str_remove(name, "r_GT\\["), ",Intercept]")) |> ggplot(aes(x = value, y = name, fill=name, color=name)) + 
stat_halfeye(point_interval = median_qi, .width = .95, color="black") + 
labs(title = "Posterior Distributions of Mutation Effects on Hearing Phenotype",
     x = "Effect on Hearing Phenotype",
     y = "Density") + 
theme_bw(base_size = 16, base_rect_size = 1, base_line_size = 1) + scale_fill_flatui() + scale_color_flatui() + theme(legend.position ="none") + 
theme(plot.margin = unit(c(1,1,1,1), "cm")) + xlim(c(-1.5, 1.5))


get11 = posterior_samples(model)
mutations <- colnames(get11)
results <- lapply(mutations, function(mut) {
  hypothesis(model, paste0(mut, " != 0"))
})

# Print p-values
p_values <- sapply(results, function(res) {
  res$p.value
})
print(p_values)

plotdf = stack(get11[,5:12])
plotdf$ind = str_remove(str_remove(plotdf$ind, "r_GT\\["), ",Intercept]")
posterior_df = plotdf

draws_fit <- as_draws_array(model)
posterior_df <- as.data.frame(draws_fit)
ggplot(posterior_df, aes(x = values, fill = ind, color = ind)) +
    #geom_density(size=0, alpha = 0.5) +
    stat_halfeye(point_interval = mean_qi, .width = .95) + 
    #facet_wrap(~parameter, scales = "free") +
    labs(title = "Posterior Distributions of Mutation Effects on Hearing Phenotype",
         x = "Effect on Hearing Phenotype",
         y = "Density") + 
    theme_cowplot(font_size = 16, line_size = 1) + scale_fill_flatui() + scale_color_flatui()


barplot(posterior_summary(model)[,1][5:12], col = pal_flatui()(8), xaxt="n")
axis(side = 1, at = 1:8, labels = rownames(posterior_summary(model))[5:12], las=2)



saveRDS(model, file = "hearing_bayesian_regression_model.rds")



# eyesight is applied as the same method
