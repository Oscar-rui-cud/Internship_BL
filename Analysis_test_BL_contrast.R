library(readr)
library(brms)
library(ggplot2)
library(tidyr)
library(tidybayes)
library(modelr)
library(ggdist)
library(rstan)
library(ggrepel)
library(RColorBrewer)
library(posterior)
library(distributional)
library(HDInterval)
library(bayesplot)
library("see")
library(bayestestR)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
BL_DATA_CONTR <- read_csv("BL_DATA_CONTR.csv")

#Filtering the 3 first duration because of a clear lacks of data points in theses durations
df <- subset(BL_DATA_CONTR, BL_DATA_CONTR$duration_test > 2000) 


#We rescale Log(E) contrast and Log(D) so their value are contained between 0 and 1
df$Ldu <- df$Ldu - min(df$Ldu)
df$Ldu <- df$Ldu/max(df$Ldu)

df$Len <- df$Len - min(df$Len)
df$Len <- df$Len/max(df$Len)

df$contrast_test <- df$contrast_test - min(df$contrast_test)
df$contrast_test <- df$contrast_test/max(df$contrast_test)


# We plot the data to have an idea of what we are doing/dealing with
one <- ggplot(data = df, aes(x=Ldu,y=Len)) + geom_point(aes(color=task))
one
two <- ggplot(data = df, aes(x=Ldu,y=contrast_test)) + geom_point(aes(color=task))
two

# First model

#formula1
f1 <- bf(
  Len ~ b0 +
    b1 * (Ldu - tC) * step(tC - Ldu) +
    b2 * (Ldu - tC) * step(Ldu - tC),
  b0 + b1 + b2 + tC ~ 1 + (1 |task |subj),  # ← correlated random effects here
  nl = TRUE
)
#priors1 (pretty flat but mostly contained within 0 and 1)
p1 <- c(
  set_prior('normal(0, 0.5)', nlpar = 'b0'),
  set_prior('normal(0, 1)', nlpar = 'b1'),
  set_prior('normal(0, 1)', nlpar = 'b2'),
  set_prior('normal(0.5, 0.25)', nlpar = 'tC'),
  set_prior('exponential(6)', class = 'sd', nlpar = 'b0'),
  set_prior('exponential(3)', class = 'sd', nlpar = 'b1'),
  set_prior('exponential(3)', class = 'sd', nlpar = 'b2'),
  set_prior('exponential(6)', class = 'sd', nlpar = 'tC'),
  set_prior('exponential(6)', class = 'sigma'),
  set_prior('lkj(1)', class = 'cor')
)

#run the model 

m1 <- brm(
  f1,
  df,
  prior = p1,
  cores = 4,
  chains = 4,
  iter = 10000,
  warmup = 3000,
  sample_prior = 'yes', 
  backend = "cmdstanr",
  save_pars = save_pars(all = TRUE),
  control = list(adapt_delta = 0.95, max_treedepth = 15)
)
pp_check(m1,ndraws = 100)
plot(conditional_effects(m1), points = TRUE)

fixef(m1)
hypothesis(m1, "b1_Intercept = 0")
#So it would seem that their is 2 times more chance for b1_Intercept to be different
#from 0 than it is to be equal to 0







#Broken stick model where the slope is null until a breaking point tC.
fsimple <- bf(
  Len ~ b0 +
    b2 * (Ldu - tC) * step(Ldu - tC),
  b0  + b2 + tC ~ 1 + (1|task|subj),  # ← correlated random effects here
  nl = TRUE
)
default_prior(fsimple, data = df)

m1simple <-brm(
  fsimple,
  data = df,
  prior = default_prior(fsimple,data = df),
  cores = 4,
  chains = 4,
  iter = 10000,
  warmup = 3000,
  sample_prior = 'yes', 
  backend = "cmdstanr",
  save_pars = save_pars(all = TRUE),
  control = list(adapt_delta = 0.95, max_treedepth = 15)
)
pp_check(m1simple,ndraws = 100)
plot(conditional_effects(m1simple), points = TRUE)


bayes_factor(m1,m1simple, log = T)




get_variables(m1)
get_variables(m2)

loo1 <- loo(m1, reloo=T)
loo2 <- loo(m1simple,reloo=T)
loo_compare(loo1,loo2)

#         elpd_diff se_diff
#m1simple  0.0       0.0   
#m1       -6.0       8.4  
summary(m1)
plot_pars(m1)



parnames(m1)
fixef(m1)
# Focus on describing posterior
hdi_range <- bayestestR::hdi(m1, ci = c(0.65, 0.70, 0.80, 0.89, 0.95))
hdi_range

plot(hdi_range, show_intercept = T)

hdi_range_simple <- bayestestR::hdi(m1simple, ci = c(0.65, 0.70, 0.80, 0.89, 0.95))
hdi_range_simple

plot(hdi_range_simple, show_intercept = T)

formula7 <- bf(  Len ~ b0 +    
                   b1 * (Ldu - tC) * step(tC - Ldu) +    
                   b2 * (Ldu - tC) * step(Ldu - tC),  
                 b0 + b1 + b2 + tC ~ 1 + (1 | gr | subj),  # ← correlated random effects here  
                 sigma ~ 1 + Ldu,  
                 nl = TRUE)
formula99 <- bf(  contrast_test ~ b0 +    
                   b1 * (Ldu - tC) * step(tC - Ldu) +    
                   b2 * (Ldu - tC) * step(Ldu - tC),  
                 b0 + b1 + b2 + tC ~ 1 + (1 | gr | subj),  # ← correlated random effects here  
                 sigma ~ 1 + Ldu,  
                 nl = TRUE)

priors7 <- c(  
  set_prior('normal(0, 0.5)', nlpar = 'b0'),  
  set_prior('normal(0, 1)', nlpar = 'b1'),  
  set_prior('normal(0, 1)', nlpar = 'b2'),  
  set_prior('normal(0.5, 0.25)', nlpar = 'tC'),  
  set_prior('exponential(6)', class = 'sd', nlpar = 'b0'),  
  set_prior('exponential(3)', class = 'sd', nlpar = 'b1'),  
  set_prior('exponential(3)', class = 'sd', nlpar = 'b2'),  
  set_prior('exponential(6)', class = 'sd', nlpar = 'tC'),  
  set_prior('normal(0, 1)', class = 'b', dpar = 'sigma'),         # prior for slope of sigma model
  set_prior('normal(0, 1)', class = 'Intercept', dpar = 'sigma'), # prior for intercept of sigma model
  set_prior('lkj(1)', class = 'cor')
)

bmm8 <- brm(
  formula7,
  df,
  family = "student",
  prior = priors7,
  cores = 4,
  chains = 4,
  iter = 10000,
  warmup = 3000,
  sample_prior = 'yes', 
  backend = "cmdstanr",
  save_pars = save_pars(all = TRUE),
  control = list(adapt_delta = 0.95, max_treedepth = 15)
)
pp_check(bmm8, ndraws = 50)
plot(conditional_effects(bmm8), points = T)
loo8 <- loo(bmm8, save_psis = T)
plot(loo8,label_points =T)
psis8 <- loo8$psis_object
psis8W <- weights(psis8)
yrep8 <- posterior_predict(bmm8)
ppc_loo_pit_overlay(df$Len,yrep8, lw = psis8W)

fixef(bmm8)
hypothesis(bmm8, "b1_Intercept = 0")


posterior <- as.matrix(bmm8)
posterior_plot <- mcmc_areas(posterior,
           pars = c("b_b1_Intercept","b_b0_Intercept","b_b2_Intercept"),
           prob = 0.8) + vline_0() + theme_classic()
posterior_plot + labs(title = "posterior distribution of bmm8", y = "parameters") + scale_color_material()


post = predicted_draws(bmm8, newdata=df, ndraws=100)

bmm9<- brm(
  formula99,
  df,
  family = "student",
  prior = priors7,
  cores = 4,
  chains = 4,
  iter = 10000,
  warmup = 3000,
  sample_prior = 'yes', 
  backend = "cmdstanr",
  save_pars = save_pars(all = TRUE),
  control = list(adapt_delta = 0.95, max_treedepth = 15)
)
pp_check(bmm9,ndraws = 100)
plot(conditional_effects(bmm9),points = T)
loo9 <- loo(bmm9, save_psis = T)
plot(loo9)
pipi <- loo9$psis_object
PIPI <- weights(pipi)
yrep9 <- posterior_predict(bmm9)
ppc_loo_pit_overlay(df$contrast_test,yrep9, lw = PIPI)

posterior <- as.matrix(bmm9)
mcmc_areas(posterior,
           pars = c("b_b1_Intercept","b_b0_Intercept","b_b2_Intercept"),
           prob = 0.8) + vline_0()


ppc_pit_ecdf(df$Len, yrep8, prob = 0.99, plot_diff = TRUE)

result <- hdi(bmm8, effects = "fixed", component = "distributional", ci =  c(0.5, 0.75, 0.89, 0.95))
plot(result, priors = T, n_columns = 2)

#to check priors
result_bf <- bayesfactor_parameters(bmm8, verbose = FALSE, component = "distributional")

result_bf

plot(result_bf) +
  scale_color_material() +
  scale_fill_material()





