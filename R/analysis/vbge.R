#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 2020.06.16: Max Lindmark
#
# Fit VBGE models that vary with respect to common or area-specific parameters by
# adding a dummy-variable (0, 1) and use WAIC to compare them.
# Because we use back-calculated data (few repetitions within the same individual),
# we fit models with catch-age-only.
# The models have sigma varying with age and cohort-varying K and L_inf.
# 
# A. Load libraries
# 
# B. Read data
# 
# C. Fit models
#
# D. Compare models
# 
# E. Produce figures
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# A. LOAD LIBRARIES ================================================================
rm(list = ls())

# Load libraries (install first if needed)
library(tidyverse); theme_set(theme_classic(base_size = 12))
library(brms)
library(nlstools)
library(viridis)
library(bayesplot)
library(tidylog)
library(tidybayes)
library(RColorBrewer)
library(patchwork)
library(modelr)

# Print package versions for versions
sessionInfo() 

# For parallel processing
options(mc.cores = parallel::detectCores()) 


# B. READ DATA =====================================================================
# df <- read.csv("data/Growth_data_BT_FM_1970-2004.csv", sep = ";") # This is the 
# original data from Huss et al (2019)
df <- read.csv("data/size_at_age_BT_FM_1970-2004.csv", sep = ";")

# Add an area-specific ID
df$ID2 <- paste(df$ID, df$area, sep = ".")

# Filter years to match the rest of the data sets and analyses
df <- df %>% filter(birth_year > 1980 & birth_year < 1998,
                    catch_year < 2004 & catch_age == age)

# Plot data
p1 <- ggplot(df, aes(age, length, color = area)) +
  geom_point() +
  facet_wrap(~ birth_year)
p1

# Now we need to reshape the data frame a little bit to add in the dummy variable
# for area, and we'll use this dataframe for the models
bt <- filter(df, area == "BT")
fm <- filter(df, area == "FM")

dfm <- data.frame(rbind(cbind(bt, areaBT=1, areaFM=0), cbind(fm, areaBT=0, areaFM=1)))

dfm$age <- as.numeric(dfm$age)

ggplot(filter(dfm, birth_year == 1988), aes(age, length, color = area)) +
  geom_point()

# They seem to come from a different gear
filter(dfm, area == "FM" & birth_year == 1988 & length > 310)
ggplot(filter(dfm, area == "FM"), aes(age, length, color = factor(gear))) +
  geom_point(size = 0.4, alpha = 0.5) +
  facet_wrap(~ birth_year)

# Remove the odd gear
sort(unique(dfm$gear))
dfm <- dfm %>% filter(!gear == 32)
sort(unique(dfm$gear))

# Now it looks better
ggplot(filter(dfm, area == "FM"), aes(age, length, color = factor(gear))) +
  geom_point(size = 0.4, alpha = 0.5) +
  facet_wrap(~ birth_year)

# And in BT more or less only gear 9 is used
ggplot(filter(dfm, area == "BT"), aes(age, length, color = factor(gear))) +
  geom_point(size = 0.4, alpha = 0.5) +
  facet_wrap(~ birth_year)

# Lastly, fit the model using cm not mm
dfm$length_cm <- dfm$length / 10
dfm$log_length_cm <- log(dfm$length_cm)


# C. FIT MODELS ====================================================================
# Here are some guides I followed
# https://rstudio-pubs-static.s3.amazonaws.com/57692_215e844f73e949ada4854dc688677dc1.html
# vignette("brms_nonlinear", package = "brms"); https://cran.r-project.org/web/packages/brms/vignettes/brms_nonlinear.html
# https://cran.r-project.org/web/packages/brms/vignettes/brms_distreg.html

# All models have year as a ID within year as random factors for K & L_inf & a model on sigma:
# https://cran.r-project.org/web/packages/brms/vignettes/brms_distreg.html

# Preliminary analysis result in bimodal posterior distributions, so therefore I 
# will put informative priors on the models.

# First selected informative priors through a prior predictive check on all data combined
# hist(rnorm(10000, mean = 45, sd = 10))    # L_inf
# hist(rnorm(10000, mean = -0.3, sd = 0.2)) # t0 
# hist(rnorm(10000, mean = 0.2, sd = 0.1))  # K

# However, this did not work with the full model. I achieved convergence by setting 
# much tighter priors.
# prior(normal(-0.5, 0.1), nlpar = "t0FM") +
# prior(normal(-0.5, 0.1), nlpar = "t0BT") +
# prior(normal(0.17, 0.02), nlpar = "KFM") +
# prior(normal(0.17, 0.02), nlpar = "KBT") +
# prior(normal(40, 1), nlpar = "LinfFM") +
# prior(normal(40, 1), nlpar = "LinfBT")

# I then relaxed them step-wise. That resulted in the following priors: 
# prior(normal(-0.5, 0.4), nlpar = "t0FM") +
# prior(normal(-0.5, 0.4), nlpar = "t0BT") +
# prior(normal(0.17, 0.05), nlpar = "KFM") +
# prior(normal(0.17, 0.05), nlpar = "KBT") +
# prior(normal(40, 5), nlpar = "LinfFM") +
# prior(normal(40, 5), nlpar = "LinfBT")

# However, they seemed to narrow.
# Therefore I tried to constrain parameters by using uniform distributions. I landed 
# on using it only for K, which was sufficient in order to be able to put broad normal
# priors on L_inf and t_0
# https://rdrr.io/cran/brms/man/set_prior.html


# M0: Prior predictive check: BT+FM ================================================
hist(rnorm(100000, mean = 45, sd = 20)) 

M0fmbt <- brm(
  bf(length_cm ~ Linf*(1-exp(-K*(age-t0))),
     Linf ~ 1, t0 ~ 1, K ~ 1, nl = TRUE),
  data = dfm, family = gaussian(),
  prior = c(prior(normal(45, 20), nlpar = "Linf"),
            prior(normal(-0.5, 1), nlpar = "t0"),
            prior("uniform(0, 0.6)", nlpar = "K", lb = 0, ub = 0.6)),
            sample_prior = "only", 
  iter = 2000, thin = 1, cores = 3, chains = 3, seed = 1)

plot(conditional_effects(M0fmbt), points = TRUE)
ggsave("figures/supp/vbge_prior_pred_check.png", width = 6.5, height = 6.5, dpi = 600)

pal <- rev(brewer.pal(n = 6, name = "Paired")[c(2, 6)])

p0 <- dfm %>% 
  data_grid(age = seq_range(age, by = 1)) %>%
  add_predicted_draws(M0fmbt, re_formula = NA) %>%
  ggplot(aes(x = factor(age), y = length_cm)) +
  stat_lineribbon(aes(y = .prediction), .width = c(.8), alpha = 0.2) +
  geom_jitter(data = dfm, aes(age, length_cm, color = area),
              alpha = 0.2, width = 0.3, height = 0, size = 0.8) +
  stat_lineribbon(aes(y = .prediction), .width = c(.8), alpha = 0.8, fill = NA,
                  color = "black") +
  guides(fill = FALSE, color = guide_legend(override.aes = list(alpha = 1))) +
  scale_color_manual(values = pal, labels = c("Warm", "Cold")) +
  labs(y = "Length (cm)", x = "Age (yrs)", fill = "Area", colour = "Area") +
  NULL

pWord0 <- p0 + theme(text = element_text(size = 12), 
                     legend.position = c(0.1, 0.9), 
                     legend.title = element_text(size = 10),
                     legend.text = element_text(size = 10))

ggsave("figures/supp/vbge_prior_pred_check.png", width = 6.5, height = 6.5, dpi = 600)

prior_summary(M0fmbt)


# M1: All parameters specific by area ==============================================
prior <-
  prior(normal(-0.5, 1), nlpar = "t0FM") +
  prior(normal(-0.5, 1), nlpar = "t0BT") +
  prior("uniform(0, 0.6)", nlpar = "KFM", lb = 0, ub = 0.6) +
  prior("uniform(0, 0.6)", nlpar = "KBT", lb = 0, ub = 0.6) +
  prior(normal(45, 20), nlpar = "LinfFM") +
  prior(normal(45, 20), nlpar = "LinfBT")
  
start_time <- Sys.time()
m1 <- 
  brm(
    bf(length_cm ~ areaBT*LinfBT*(1-exp(-KBT*(age-t0BT))) + areaFM*LinfFM*(1-exp(-KFM*(age-t0FM))),
       sigma ~ age,
       t0FM ~ 1,
       t0BT ~ 1,
       KFM ~ 1 + (1|birth_year),    # parameter varying by birth_year
       KBT ~ 1 + (1|birth_year),    # parameter varying by birth_year
       LinfFM ~ 1 + (1|birth_year), # parameter varying by birth_year
       LinfBT ~ 1 + (1|birth_year), # parameter varying by birth_year
       nl = TRUE),
    data = dfm,
    family = gaussian(),
    prior = prior,
    iter = 3000, thin = 1, cores = 3, chains = 3, inits = "0",
    control = list(max_treedepth = 13, adapt_delta = 0.9))
end_time <- Sys.time()
end_time - start_time

summary(m1)
plot(m1)

# Save model object to not have to rerun it...
saveRDS(m1, "output/vbge/m1.rds")
# m1 <- readRDS("output/vbge/m1.rds")


# M2: L_inf common parameter, t_0 & K specific =====================================
prior2 <-
  prior(normal(-0.5, 1), nlpar = "t0FM") +
  prior(normal(-0.5, 1), nlpar = "t0BT") +
  prior("uniform(0, 0.6)", nlpar = "KFM", lb = 0, ub = 0.6) +
  prior("uniform(0, 0.6)", nlpar = "KBT", lb = 0, ub = 0.6) +
  prior(normal(45, 20), nlpar = "Linf")

start_time <- Sys.time()
m2 <- 
  brm(
    bf(length_cm ~ areaBT*Linf*(1-exp(-KBT*(age-t0BT))) + areaFM*Linf*(1-exp(-KFM*(age-t0FM))),
       sigma ~ age,
       t0FM ~ 1,
       t0BT ~ 1,
       KFM ~ 1 + (1|birth_year),    # parameter varying by birth_year
       KBT ~ 1 + (1|birth_year),    # parameter varying by birth_year
       Linf ~ 1 + (1|birth_year), # parameter varying by birth_year
       nl = TRUE),
    data = dfm,
    family = gaussian(),
    prior = prior2,
    iter = 3000, thin = 1, cores = 3, chains = 3, inits = "0",
    control = list(max_treedepth = 13, adapt_delta = 0.9))
end_time <- Sys.time()
end_time - start_time

summary(m2)
plot(m2)

# Save model object to not have to rerun it...
saveRDS(m2, "output/vbge/m2.rds")
# m2 <- readRDS("output/vbge/m2.rds")


# M3: K common parameter, t_0 & L_inf specific =====================================
prior3 <-
  prior(normal(-0.5, 1), nlpar = "t0FM") +
  prior(normal(-0.5, 1), nlpar = "t0BT") +
  prior("uniform(0, 0.6)", nlpar = "K", lb = 0, ub = 0.6) +
  prior(normal(45, 20), nlpar = "LinfFM") +
  prior(normal(45, 20), nlpar = "LinfBT")

start_time <- Sys.time()
m3 <- 
  brm(
    bf(length_cm ~ areaBT*LinfBT*(1-exp(-K*(age-t0BT))) + areaFM*LinfFM*(1-exp(-K*(age-t0FM))),
       sigma ~ age,
       t0FM ~ 1,
       t0BT ~ 1,
       K ~ 1 + (1|birth_year),      # parameter varying by birth_year
       LinfFM ~ 1 + (1|birth_year), # parameter varying by birth_year
       LinfBT ~ 1 + (1|birth_year), # parameter varying by birth_year
       nl = TRUE),
    data = dfm,
    family = gaussian(),
    prior = prior3,
    iter = 3000, thin = 1, cores = 3, chains = 3, inits = "0",
    control = list(max_treedepth = 13, adapt_delta = 0.9))
end_time <- Sys.time()
end_time - start_time

summary(m3)
plot(m3)

# Save model object to not have to rerun it...
saveRDS(m3, "output/vbge/m3.rds")
# m3 <- readRDS("output/vbge/m3.rds")


# M4: t_0 common parameter, K & L_inf specific =====================================
prior4 <-
  prior(normal(-0.5, 1), nlpar = "t0") +
  prior("uniform(0, 0.6)", nlpar = "KFM", lb = 0, ub = 0.6) +
  prior("uniform(0, 0.6)", nlpar = "KBT", lb = 0, ub = 0.6) +
  prior(normal(45, 20), nlpar = "LinfFM") +
  prior(normal(45, 20), nlpar = "LinfBT")
  
start_time <- Sys.time()
m4 <- 
  brm(
    bf(length_cm ~ areaBT*LinfBT*(1-exp(-KBT*(age-t0))) + areaFM*LinfFM*(1-exp(-KFM*(age-t0))),
       sigma ~ age,
       t0 ~ 1,
       KFM ~ 1 + (1|birth_year),    # parameter varying by birth_year
       KBT ~ 1 + (1|birth_year),    # parameter varying by birth_year
       LinfFM ~ 1 + (1|birth_year), # parameter varying by birth_year
       LinfBT ~ 1 + (1|birth_year), # parameter varying by birth_year
       nl = TRUE),
    data = dfm,
    family = gaussian(),
    prior = prior4,
    iter = 3000, thin = 1, cores = 3, chains = 3, inits = "0",
    control = list(max_treedepth = 13, adapt_delta = 0.9))
end_time <- Sys.time()
end_time - start_time

summary(m4)
plot(m4)

# Save model object to not have to rerun it...
saveRDS(m4, "output/vbge/m4.rds")
# m4 <- readRDS("output/vbge/m4.rds")


# M5: L_inf & K common parameter, t_0 specific =====================================
prior5 <-
  prior(normal(-0.5, 1), nlpar = "t0FM") +
  prior(normal(-0.5, 1), nlpar = "t0BT") +
  prior("uniform(0, 0.6)", nlpar = "K", lb = 0, ub = 0.6) +
  prior(normal(45, 20), nlpar = "Linf")
  
start_time <- Sys.time()
m5 <- 
  brm(
    bf(length_cm ~ areaBT*Linf*(1-exp(-K*(age-t0BT))) + areaFM*Linf*(1-exp(-K*(age-t0FM))),
       sigma ~ age,
       t0FM ~ 1,
       t0BT ~ 1,
       K ~ 1 + (1|birth_year),    # parameter varying by birth_year
       Linf ~ 1 + (1|birth_year), # parameter varying by birth_year
       nl = TRUE),
    data = dfm,
    family = gaussian(),
    prior = prior5,
    iter = 3000, thin = 1, cores = 3, chains = 3, inits = "0",
    control = list(max_treedepth = 13, adapt_delta = 0.9))
end_time <- Sys.time()
end_time - start_time

summary(m5)
plot(m5)

# Save model object to not have to rerun it...
saveRDS(m5, "output/vbge/m5.rds")
# m5 <- readRDS("output/vbge/m5.rds")


# M6: L_inf & t_0 common parameter, K specific =====================================
prior6 <-
  prior(normal(-0.5, 1), nlpar = "t0") +
  prior("uniform(0, 0.6)", nlpar = "KFM", lb = 0, ub = 0.6) +
  prior("uniform(0, 0.6)", nlpar = "KBT", lb = 0, ub = 0.6) +
  prior(normal(45, 20), nlpar = "Linf")

start_time <- Sys.time()
m6 <- 
  brm(
    bf(length_cm ~ areaBT*Linf*(1-exp(-KBT*(age-t0))) + areaFM*Linf*(1-exp(-KFM*(age-t0))),
       sigma ~ age,
       t0 ~ 1,
       KFM ~ 1 + (1|birth_year),  # parameter varying by birth_year
       KBT ~ 1 + (1|birth_year),  # parameter varying by birth_year
       Linf ~ 1 + (1|birth_year), # parameter varying by birth_year
       nl = TRUE),
    data = dfm,
    family = gaussian(),
    prior = prior6,
    iter = 3000, thin = 1, cores = 3, chains = 3, inits = "0",
    control = list(max_treedepth = 13, adapt_delta = 0.9))
end_time <- Sys.time()
end_time - start_time

summary(m6)
plot(m6)

# Save model object to not have to rerun it...
saveRDS(m6, "output/vbge/m6.rds")
# m6 <- readRDS("output/vbge/m6.rds")


# M7: K & t_0 common parameter, L_inf specific =====================================
prior7 <-
  prior(normal(-0.5, 1), nlpar = "t0") +
  prior("uniform(0, 0.6)", nlpar = "K", lb = 0, ub = 0.6) +
  prior(normal(45, 20), nlpar = "LinfFM") +
  prior(normal(45, 20), nlpar = "LinfBT")

start_time <- Sys.time()
m7 <- 
  brm(
    bf(length_cm ~ areaBT*LinfBT*(1-exp(-K*(age-t0))) + areaFM*LinfFM*(1-exp(-K*(age-t0))),
       sigma ~ age,
       t0 ~ 1,
       K ~ 1 + (1|birth_year),      # parameter varying by birth_year
       LinfFM ~ 1 + (1|birth_year), # parameter varying by birth_year
       LinfBT ~ 1 + (1|birth_year), # parameter varying by birth_year
       nl = TRUE),
    data = dfm,
    family = gaussian(),
    prior = prior7,
    iter = 3000, thin = 1, cores = 3, chains = 3, inits = "0",
    control = list(max_treedepth = 13, adapt_delta = 0.9))
end_time <- Sys.time()
end_time - start_time

summary(m7)
plot(m7)

# Save model object to not have to rerun it...
saveRDS(m7, "output/vbge/m7.rds")
# m7 <- readRDS("output/vbge/m7.rds")


# M8: All parameters common ========================================================
prior8 <-
  prior(normal(-0.5, 1), nlpar = "t0") +
  prior("uniform(0, 0.6)", nlpar = "K", lb = 0, ub = 0.6) +
  prior(normal(45, 20), nlpar = "Linf")

start_time <- Sys.time()
m8 <- 
  brm(
    bf(length_cm ~ areaBT*Linf*(1-exp(-K*(age-t0))) + areaFM*Linf*(1-exp(-K*(age-t0))),
       sigma ~ age,
       t0 ~ 1,
       K ~ 1 + (1|birth_year),    # parameter varying by birth_year
       Linf ~ 1 + (1|birth_year), # parameter varying by birth_year
       nl = TRUE),
    data = dfm,
    family = gaussian(),
    prior = prior8,
    iter = 3000, thin = 1, cores = 3, chains = 3, inits = "0",
    control = list(max_treedepth = 13, adapt_delta = 0.9))
end_time <- Sys.time()
end_time - start_time

summary(m8)
plot(m8)

# Save model object to not have to rerun it...
saveRDS(m8, "output/vbge/m8.rds")
# m8 <- readRDS("output/vbge/m8.rds")



# D. COMPARE MODELS ================================================================
loo_m1 <- loo(m1)
loo_m2 <- loo(m2)
loo_m3 <- loo(m3)
loo_m4 <- loo(m4)
loo_m5 <- loo(m5)
loo_m6 <- loo(m6)
loo_m7 <- loo(m7)
loo_m8 <- loo(m8)

# Compare models
loo_compare(loo_m1, loo_m2, loo_m3, loo_m4, loo_m5, loo_m6, loo_m7, loo_m8)
#     elpd_diff   se_diff
# m1     0.0       0.0
# m4    -6.7       3.3
# m2   -95.6      15.2
# m3  -142.2      19.5
# m7  -149.4      19.6
# m6  -159.5      20.4
# m5  -881.5      52.7
# m8 -1831.3      70.4


# Using model 1


# E. PRODUCE FIGURES ===============================================================
# https://mjskay.github.io/tidybayes/articles/tidy-brms.html

# summary(m1)
# Population-Level Effects: 
# Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# sigma_Intercept      0.21      0.02     0.17     0.25 1.00     5769     2397
# t0FM_Intercept      -0.44      0.07    -0.58    -0.31 1.00     3311     2316
# t0BT_Intercept      -0.14      0.03    -0.21    -0.08 1.00     4034     2328
# KFM_Intercept        0.14      0.01     0.12     0.17 1.00      975     1624
# KBT_Intercept        0.19      0.02     0.15     0.24 1.00      919     1349
# LinfFM_Intercept    39.68      2.27    35.53    44.33 1.00     1280     1549
# LinfBT_Intercept    46.46      4.97    37.40    57.08 1.00      705     1164
# sigma_age            0.15      0.00     0.14     0.16 1.00     5592     2381

# Plot main predictions
p1 <- dfm %>% 
  data_grid(age = seq_range(age, by = 1),
            area = c("FM", "BT")) %>%
  mutate(areaFM = ifelse(area == "FM", 1, 0),
         areaBT = ifelse(area == "BT", 1, 0)) %>% 
  add_predicted_draws(m1, re_formula = NA) %>%
  ggplot(aes(x = factor(age), y = length_cm, color = area, fill = area)) +
    stat_lineribbon(aes(y = .prediction), .width = c(.8), alpha = 0.2) +
    geom_jitter(data = dfm, alpha = 0.2, width = 0.3,
                height = 0, size = 0.8) +
    stat_lineribbon(aes(y = .prediction), .width = c(.8), alpha = 0.8, fill = NA) +
    scale_fill_manual(values = pal, labels = c("Warm", "Cold")) +
    scale_color_manual(values = pal, labels = c("Warm", "Cold")) +
    labs(y = "Length (cm)", x = "Age (yrs)", fill = "Area", colour = "Area") +
  NULL

pWord1 <- p1 + theme(text = element_text(size = 12), # 12 for word doc
                     legend.position = c(0, 1), 
                     legend.justification = c(0, 1),
                     legend.title = element_text(size = 20),
                     legend.text = element_text(size = 20))

ggsave("figures/vbge/vbge_pred.png", width = 6.5, height = 6.5, dpi = 600)
  

# Plot predictions by year:
p2 <- dfm %>% 
  data_grid(age = seq_range(age, by = 1),
            birth_year = seq_range(birth_year, by = 1),
            area = c("FM", "BT")) %>%
  mutate(areaFM = ifelse(area == "FM", 1, 0),
         areaBT = ifelse(area == "BT", 1, 0)) %>% 
  add_predicted_draws(m1) %>%
  ggplot(aes(x = factor(age), y = length_cm, color = area, fill = area)) +
    stat_lineribbon(aes(y = .prediction), .width = c(.8), alpha = 0.2) +
    geom_jitter(data = dfm, alpha = 0.2, width = 0.3,
                height = 0, size = 0.6) +
    stat_lineribbon(aes(y = .prediction), .width = c(.8), alpha = 0.8,
                    fill = NA, size = 0.8) +
    facet_wrap(~birth_year) +
    scale_fill_manual(values = pal, labels = c("Warm", "Cold")) +
    scale_color_manual(values = pal, labels = c("Warm", "Cold")) +
    labs(y = "Length (cm)", x = "Age (yrs)", fill = "Area", colour = "Area") +
  NULL

pWord2 <- p2 + theme(text = element_text(size = 12), # 12 for word doc
                     legend.position = c(0.7, 0.1), 
                     legend.title = element_text(size = 20),
                     legend.text = element_text(size = 20))

ggsave("figures/supp/vbge_pred_year.png", width = 6.5, height = 6.5, dpi = 600)


# Plot posterior predictive checks 
pp_check(m1, nsamples = 50) + 
  theme(text = element_text(size = 12),
        legend.position = c(0.8, 0.8), 
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 12))

ggsave("figures/supp/vbge_ppc.png", width = 6.5, height = 6.5, dpi = 600)


# Plotting mcmc_dens and use patchwork to plot them together. Note I add the vertical
# lines manually simply by extracting the fixed effects
m1_fe <- fixef(m1, probs = c(0.1, 0.9)) %>% as.data.frame()
posterior <- as.array(m1)

# Define matching palette
pal2 <- alpha(pal, alpha = 0.5)

color_scheme_set(rep("white", 6)) # This is to be able to have a fill color with alpha

Linf_BT <- mcmc_dens(posterior, pars = c("b_LinfBT_Intercept"),
          facet_args = list(nrow = 2)) + 
  geom_density(fill = pal2[1]) + 
  geom_vline(xintercept = m1_fe$Estimate[7], linetype = 1, color = "white") +
  geom_vline(xintercept = m1_fe$Q10[7], linetype = 2, color = "white") +
  geom_vline(xintercept = m1_fe$Q90[7], linetype = 2, color = "white") +
  coord_cartesian(xlim = c(32, 65)) +
  labs(x = expression(paste(italic(L[inf]), " [cm]")), y = "Warm") + 
  #annotate("text", -Inf, Inf, label = "Warm", size = 5, hjust = -0.5, vjust = 1.3) +
  theme(text = element_text(size = 12)) 

Linf_FM <- mcmc_dens(posterior, pars = c("b_LinfFM_Intercept"),
                     facet_args = list(nrow = 2)) + 
  geom_density(fill = pal2[2]) + 
  geom_vline(xintercept = m1_fe$Estimate[6], linetype = 1, color = "white") +
  geom_vline(xintercept = m1_fe$Q10[6], linetype = 2, color = "white") +
  geom_vline(xintercept = m1_fe$Q90[6], linetype = 2, color = "white") +
  coord_cartesian(xlim = c(32, 65)) +
  labs(x = expression(paste(italic(L[inf]), " [cm]")), y = "Cold") + 
  #annotate("text", -Inf, Inf, label = "Cold", size = 5, hjust = -0.5, vjust = 1.3) +
  theme(text = element_text(size = 12))

K_BT <- mcmc_dens(posterior, pars = c("b_KBT_Intercept"),
                     facet_args = list(nrow = 2)) + 
  geom_density(fill = pal2[1]) + 
  geom_vline(xintercept = m1_fe$Estimate[5], linetype = 1, color = "white") +
  geom_vline(xintercept = m1_fe$Q10[5], linetype = 2, color = "white") +
  geom_vline(xintercept = m1_fe$Q90[5], linetype = 2, color = "white") +
  coord_cartesian(xlim = c(0.1, 0.27)) +
  xlab(expression(paste(italic(K)))) + 
  theme(text = element_text(size = 12)) 

K_FM <- mcmc_dens(posterior, pars = c("b_KFM_Intercept"),
                     facet_args = list(nrow = 2)) + 
  geom_density(fill = pal2[2]) + 
  geom_vline(xintercept = m1_fe$Estimate[4], linetype = 1, color = "white") +
  geom_vline(xintercept = m1_fe$Q10[4], linetype = 2, color = "white") +
  geom_vline(xintercept = m1_fe$Q90[4], linetype = 2, color = "white") +
  coord_cartesian(xlim = c(0.1, 0.27)) +
  xlab(expression(paste(italic(K)))) + 
  theme(text = element_text(size = 12))

Linf_BT + K_BT + Linf_FM + K_FM

ggsave("figures/vbge/K_Linf_posterior.png", width = 6.5, height = 6.5, dpi = 600)


# Plot random effects (as offsets, check this is true!)
# FM_Linf
color_scheme_set("blue")
random_FM_Linf <- mcmc_areas_ridges(posterior, regex_pars = "r_birth_year__LinfFM")
word_random_FM_Linf <- random_FM_Linf + theme(text = element_text(size = 12))
ggsave("figures/vbge/random_FM_Linf.png", width = 6.5, height = 6.5, dpi = 600)

# BT_Linf
color_scheme_set("red")
random_BT_Linf <- mcmc_areas_ridges(posterior, regex_pars = "r_birth_year__LinfBT")
word_random_BT_Linf <- random_BT_Linf + theme(text = element_text(size = 12))
ggsave("figures/vbge/random_BT_Linf.png", width = 6.5, height = 6.5, dpi = 600)

# FM_K
color_scheme_set("blue")
random_FM_K <- mcmc_areas_ridges(posterior, regex_pars = "r_birth_year__KFM")
word_random_FM_K <- random_FM_K + theme(text = element_text(size = 12))
ggsave("figures/vbge/random_FM_K.png", width = 6.5, height = 6.5, dpi = 600)

# BT_K
color_scheme_set("red")
random_BT_K <- mcmc_areas_ridges(posterior, regex_pars = "r_birth_year__KBT")
word_random_BT_K <- random_BT_K + theme(text = element_text(size = 12))
ggsave("figures/vbge/random_BT_K.png", width = 6.5, height = 6.5, dpi = 600)


# More plotting options from Bayesplot
# posterior <- as.array(m1)

# Check parameter names
# names(posterior[1,1,])
#
# # Plot MCMC draws
# color_scheme_set("pink")
#
# # L_inf
# p_inf <- mcmc_areas(posterior, pars = c("b_LinfFM_Intercept", "b_LinfBT_Intercept"),
#   prob = 0.8, point_est = "mean") +
#   scale_y_discrete(labels = c("Cold L_infinity", "Warm L_infinity")) +
#   xlab("Length (cm)")
#
# pWordp_inf <- p_inf + theme(text = element_text(size = 12))
#
# ggsave("figures/vbge/L_inf_posterior.png", width = 6.5, height = 6.5, dpi = 600)
#
#
# # K
# p_K <- mcmc_areas(posterior, pars = c("b_KBT_Intercept", "b_KFM_Intercept"),
#   prob = 0.8, prob_outer = 0.99, point_est = "mean") +
#   scale_y_discrete(labels = c("Warm K", "Cold K")) +
#   xlab("K (t^-1)")
#
# pWordp_K <- p_K + theme(text = element_text(size = 12))
#
# ggsave("figures/vbge/K_posterior.png", width = 6.5, height = 6.5, dpi = 600)