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

dfm <- data.frame(rbind(cbind(bt, areaW=1, areaC=0), cbind(fm, areaW=0, areaC=1)))

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

# Change age to integer
dfm$age <- as.integer(dfm$age)


# C. FIT MODELS ====================================================================
# Here are some guides I followed
# https://rstudio-pubs-static.s3.amazonaws.com/57692_215e844f73e949ada4854dc688677dc1.html
# vignette("brms_nonlinear", package = "brms"); https://cran.r-project.org/web/packages/brms/vignettes/brms_nonlinear.html
# https://cran.r-project.org/web/packages/brms/vignettes/brms_distreg.html

# All models have year as a ID within year as random factors for K & L_inf & a model on sigma:
# https://cran.r-project.org/web/packages/brms/vignettes/brms_distreg.html

# Preliminary analysis result in bimodal posterior distributions, so therefore I 
# will put informative priors on the models.

# Below BT is the warm area, FM is the cold. After this code I use warm/cold instead

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


# M0: Prior predictive check: Warm+Cold merged =====================================
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
  labs(y = "Length [cm]", x = "Age [yrs]", fill = "Area", colour = "Area") +
  NULL

pWord0 <- p0 + theme(text = element_text(size = 12), 
                     legend.position = c(0.1, 0.9), 
                     legend.title = element_text(size = 10),
                     legend.text = element_text(size = 10))

ggsave("figures/supp/vbge_prior_pred_check.png", width = 6.5, height = 6.5, dpi = 600)

prior_summary(M0fmbt)


# M1: All parameters specific by area ==============================================
prior <-
  prior(normal(-0.5, 1), nlpar = "t0C") +
  prior(normal(-0.5, 1), nlpar = "t0W") +
  prior("uniform(0, 0.6)", nlpar = "KC", lb = 0, ub = 0.6) +
  prior("uniform(0, 0.6)", nlpar = "KW", lb = 0, ub = 0.6) +
  prior(normal(45, 20), nlpar = "LinfC") +
  prior(normal(45, 20), nlpar = "LinfW")
  
start_time <- Sys.time()
m1 <- 
  brm(
    bf(length_cm ~ areaW*LinfW*(1-exp(-KW*(age-t0W))) + areaC*LinfC*(1-exp(-KC*(age-t0C))),
       t0C ~ 1,
       t0W ~ 1,
       KC ~ 1 + (1|birth_year),    # parameter varying by birth_year
       KW ~ 1 + (1|birth_year),    # parameter varying by birth_year
       LinfC ~ 1 + (1|birth_year), # parameter varying by birth_year
       LinfW ~ 1 + (1|birth_year), # parameter varying by birth_year
       nl = TRUE),
    data = dfm,
    family = student(),
    prior = prior,
    iter = 3000, thin = 1, cores = 3, chains = 3, inits = "0",
    control = list(max_treedepth = 13, adapt_delta = 0.9))
end_time <- Sys.time()
end_time - start_time
# between 0.5-2 hours on macbook pro, forgot to record this

summary(m1)
plot(m1)

# Save model object to not have to rerun it...
#saveRDS(m1, "output/vbge/m1.rds")
#m1 <- readRDS("output/vbge/m1.rds")

# > prior_summary(m1)
#               prior class      coef      group resp dpar nlpar               bound
# 1       uniform(0, 0.6)     b                                   KC <lower=0,upper=0.6>
# 2                           b Intercept                         KC                    
# 3       uniform(0, 0.6)     b                                   KW <lower=0,upper=0.6>
# 4                           b Intercept                         KW                    
# 5        normal(45, 20)     b                                LinfC                    
# 6                           b Intercept                      LinfC                    
# 7        normal(45, 20)     b                                LinfW                    
# 8                           b Intercept                      LinfW                    
# 9       normal(-0.5, 1)     b                                  t0C                    
# 10                          b Intercept                        t0C                    
# 11      normal(-0.5, 1)     b                                  t0W                    
# 12                          b Intercept                        t0W                    
# 13        gamma(2, 0.1)    nu                                                         
# 14 student_t(3, 0, 5.8)    sd                                   KC                    
# 15 student_t(3, 0, 5.8)    sd                                   KW                    
# 16 student_t(3, 0, 5.8)    sd                                LinfC                    
# 17 student_t(3, 0, 5.8)    sd                                LinfW                    
# 18                         sd           birth_year              KC                    
# 19                         sd Intercept birth_year              KC                    
# 20                         sd           birth_year              KW                    
# 21                         sd Intercept birth_year              KW                    
# 22                         sd           birth_year           LinfC                    
# 23                         sd Intercept birth_year           LinfC                    
# 24                         sd           birth_year           LinfW                    
# 25                         sd Intercept birth_year           LinfW                    
# 26 student_t(3, 0, 5.8) sigma  


# M2: L_inf common parameter, t_0 & K specific =====================================
prior2 <-
  prior(normal(-0.5, 1), nlpar = "t0C") +
  prior(normal(-0.5, 1), nlpar = "t0W") +
  prior("uniform(0, 0.6)", nlpar = "KC", lb = 0, ub = 0.6) +
  prior("uniform(0, 0.6)", nlpar = "KW", lb = 0, ub = 0.6) +
  prior(normal(45, 20), nlpar = "Linf")

start_time <- Sys.time()
m2 <- 
  brm(
    bf(length_cm ~ areaW*Linf*(1-exp(-KW*(age-t0W))) + areaC*Linf*(1-exp(-KC*(age-t0C))),
       t0C ~ 1,
       t0W ~ 1,
       KC ~ 1 + (1|birth_year),    # parameter varying by birth_year
       KW ~ 1 + (1|birth_year),    # parameter varying by birth_year
       Linf ~ 1 + (1|birth_year),   # parameter varying by birth_year
       nl = TRUE),
    data = dfm,
    family = student(),
    prior = prior2,
    iter = 3000, thin = 1, cores = 3, chains = 3, inits = "0",
    control = list(max_treedepth = 13, adapt_delta = 0.9))
end_time <- Sys.time()
end_time - start_time
# Time difference of 1.643709 hours

summary(m2)
plot(m2)

# Save model object to not have to rerun it...
saveRDS(m2, "output/vbge/m2.rds")
# m2 <- readRDS("output/vbge/m2.rds")


# M3: K common parameter, t_0 & L_inf specific =====================================
prior3 <-
  prior(normal(-0.5, 1), nlpar = "t0C") +
  prior(normal(-0.5, 1), nlpar = "t0W") +
  prior("uniform(0, 0.6)", nlpar = "K", lb = 0, ub = 0.6) +
  prior(normal(45, 20), nlpar = "LinfC") +
  prior(normal(45, 20), nlpar = "LinfW")

start_time <- Sys.time()
m3 <- 
  brm(
    bf(length_cm ~ areaW*LinfW*(1-exp(-K*(age-t0W))) + areaC*LinfC*(1-exp(-K*(age-t0C))),
       t0C ~ 1,
       t0W ~ 1,
       K ~ 1 + (1|birth_year),      # parameter varying by birth_year
       LinfC ~ 1 + (1|birth_year), # parameter varying by birth_year
       LinfW ~ 1 + (1|birth_year), # parameter varying by birth_year
       nl = TRUE),
    data = dfm,
    family = student(),
    prior = prior3,
    iter = 3000, thin = 1, cores = 3, chains = 3, inits = "0",
    control = list(max_treedepth = 13, adapt_delta = 0.9))
end_time <- Sys.time()
end_time - start_time
# Time difference of 1.179841 hours

summary(m3)
plot(m3)

# Save model object to not have to rerun it...
saveRDS(m3, "output/vbge/m3.rds")
# m3 <- readRDS("output/vbge/m3.rds")


# M4: t_0 common parameter, K & L_inf specific =====================================
prior4 <-
  prior(normal(-0.5, 1), nlpar = "t0") +
  prior("uniform(0, 0.6)", nlpar = "KC", lb = 0, ub = 0.6) +
  prior("uniform(0, 0.6)", nlpar = "KW", lb = 0, ub = 0.6) +
  prior(normal(45, 20), nlpar = "LinfC") +
  prior(normal(45, 20), nlpar = "LinfW")
  
start_time <- Sys.time()
m4 <- 
  brm(
    bf(length_cm ~ areaW*LinfW*(1-exp(-KW*(age-t0))) + areaC*LinfC*(1-exp(-KC*(age-t0))),
       t0 ~ 1,
       KC ~ 1 + (1|birth_year),    # parameter varying by birth_year
       KW ~ 1 + (1|birth_year),    # parameter varying by birth_year
       LinfC ~ 1 + (1|birth_year), # parameter varying by birth_year
       LinfW ~ 1 + (1|birth_year), # parameter varying by birth_year
       nl = TRUE),
    data = dfm,
    family = student(),
    prior = prior4,
    iter = 3000, thin = 1, cores = 3, chains = 3, inits = "0",
    control = list(max_treedepth = 13, adapt_delta = 0.9))
end_time <- Sys.time()
end_time - start_time
# Time difference of 2.186266 hours

summary(m4)
plot(m4)

# Save model object to not have to rerun it...
saveRDS(m4, "output/vbge/m4.rds")
# m4 <- readRDS("output/vbge/m4.rds")


# M5: L_inf & K common parameter, t_0 specific =====================================
prior5 <-
  prior(normal(-0.5, 1), nlpar = "t0C") +
  prior(normal(-0.5, 1), nlpar = "t0W") +
  prior("uniform(0, 0.6)", nlpar = "K", lb = 0, ub = 0.6) +
  prior(normal(45, 20), nlpar = "Linf")
  
start_time <- Sys.time()
m5 <- 
  brm(
    bf(length_cm ~ areaW*Linf*(1-exp(-K*(age-t0W))) + areaC*Linf*(1-exp(-K*(age-t0C))),
       t0C ~ 1,
       t0W ~ 1,
       K ~ 1 + (1|birth_year),    # parameter varying by birth_year
       Linf ~ 1 + (1|birth_year), # parameter varying by birth_year
       nl = TRUE),
    data = dfm,
    family = student(),
    prior = prior5,
    iter = 3000, thin = 1, cores = 3, chains = 3, inits = "0",
    control = list(max_treedepth = 13, adapt_delta = 0.9))
end_time <- Sys.time()
end_time - start_time
# Time difference of 37.9494 mins

summary(m5)
plot(m5)

# Save model object to not have to rerun it...
saveRDS(m5, "output/vbge/m5.rds")
# m5 <- readRDS("output/vbge/m5.rds")


# M6: L_inf & t_0 common parameter, K specific =====================================
prior6 <-
  prior(normal(-0.5, 1), nlpar = "t0") +
  prior("uniform(0, 0.6)", nlpar = "KC", lb = 0, ub = 0.6) +
  prior("uniform(0, 0.6)", nlpar = "KW", lb = 0, ub = 0.6) +
  prior(normal(45, 20), nlpar = "Linf")

start_time <- Sys.time()
m6 <- 
  brm(
    bf(length_cm ~ areaW*Linf*(1-exp(-KW*(age-t0))) + areaC*Linf*(1-exp(-KC*(age-t0))),
       t0 ~ 1,
       KC ~ 1 + (1|birth_year),  # parameter varying by birth_year
       KW ~ 1 + (1|birth_year),  # parameter varying by birth_year
       Linf ~ 1 + (1|birth_year), # parameter varying by birth_year
       nl = TRUE),
    data = dfm,
    family = student(),
    prior = prior6,
    iter = 3000, thin = 1, cores = 3, chains = 3, inits = "0",
    control = list(max_treedepth = 13, adapt_delta = 0.9))
end_time <- Sys.time()
end_time - start_time
# Time difference of 1.077574 hours

summary(m6)
plot(m6)

# Save model object to not have to rerun it...
saveRDS(m6, "output/vbge/m6.rds")
# m6 <- readRDS("output/vbge/m6.rds")


# M7: K & t_0 common parameter, L_inf specific =====================================
prior7 <-
  prior(normal(-0.5, 1), nlpar = "t0") +
  prior("uniform(0, 0.6)", nlpar = "K", lb = 0, ub = 0.6) +
  prior(normal(45, 20), nlpar = "LinfC") +
  prior(normal(45, 20), nlpar = "LinfW")

start_time <- Sys.time()
m7 <- 
  brm(
    bf(length_cm ~ areaW*LinfW*(1-exp(-K*(age-t0))) + areaC*LinfC*(1-exp(-K*(age-t0))),
       t0 ~ 1,
       K ~ 1 + (1|birth_year),      # parameter varying by birth_year
       LinfC ~ 1 + (1|birth_year), # parameter varying by birth_year
       LinfW ~ 1 + (1|birth_year), # parameter varying by birth_year
       nl = TRUE),
    data = dfm,
    family = student(),
    prior = prior7,
    iter = 3000, thin = 1, cores = 3, chains = 3, inits = "0",
    control = list(max_treedepth = 13, adapt_delta = 0.9))
end_time <- Sys.time()
end_time - start_time
# Time difference of 58.70017 mins

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
    bf(length_cm ~ areaW*Linf*(1-exp(-K*(age-t0))) + areaC*Linf*(1-exp(-K*(age-t0))),
       t0 ~ 1,
       K ~ 1 + (1|birth_year),    # parameter varying by birth_year
       Linf ~ 1 + (1|birth_year), # parameter varying by birth_year
       nl = TRUE),
    data = dfm,
    family = student(),
    prior = prior8,
    iter = 3000, thin = 1, cores = 3, chains = 3, inits = "0",
    control = list(max_treedepth = 13, adapt_delta = 0.9))
end_time <- Sys.time()
end_time - start_time
# Time difference of 32.73364 mins

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
# elpd_diff se_diff
# m1     0.0       0.0
# m4    -8.9       5.1
# m2  -104.3      21.5
# m3  -127.9      24.8
# m7  -140.3      24.6
# m6  -181.1      26.0
# m5  -789.1      50.2
# m8 -1617.7      68.4

# Using model 1


# E. PRODUCE FIGURES ===============================================================
# https://mjskay.github.io/tidybayes/articles/tidy-brms.html

pal <- rev(brewer.pal(n = 6, name = "Paired")[c(2, 6)])

##### Plot predictions ============================================================
# summary(m1)
# Population-Level Effects: 
#   Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS

# Plot main predictions
p1 <- dfm %>% 
  data_grid(age = seq_range(age, by = 1),
            area = c("FM", "BT")) %>%
  mutate(areaC = ifelse(area == "FM", 1, 0),
         areaW = ifelse(area == "BT", 1, 0)) %>% 
  add_predicted_draws(m1, re_formula = NA) %>%
  ggplot(aes(x = factor(age), y = length_cm, color = area, fill = area)) +
  stat_lineribbon(aes(y = .prediction), .width = c(.5, 0.9), alpha = 0.2, size = 0.8) +
  geom_jitter(data = dfm, alpha = 0.1, width = 0.3, height = 0, size = 0.8) +
  stat_lineribbon(aes(y = .prediction), .width = 0, alpha = 0.8, size = 0.8) +
  guides(fill = FALSE,
         color = guide_legend(override.aes = list(linetype = 0, size = 3, shape = 16, alpha = 0.5))) +
  scale_fill_manual(values = pal, labels = c("Warm", "Cold")) +
  scale_color_manual(values = pal, labels = c("Warm", "Cold")) +
  labs(y = "Length [cm]", x = "Age [yrs]", fill = "Area", colour = "Area") +
  annotate("text", 8, 10, label = paste("n=", nrow(dfm), sep = ""), size = 3) +
  NULL

pWord1 <- p1 + theme(text = element_text(size = 12), # 12 for word doc
                     legend.position = c(0.1, 0.9), 
                     legend.spacing.y = unit(0, 'cm'),
                     legend.key.size = unit(0, "cm"),
                     legend.title = element_text(size = 10),
                     legend.text = element_text(size = 10))

#ggsave("figures/supp/vbge_pred.png", width = 6.5, height = 6.5, dpi = 600)
  

# Plotting mcmc_dens and use patchwork to plot them together. Note I add the vertical
# lines manually simply by extracting the fixed effects
m1_fe <- fixef(m1, probs = c(0.1, 0.9)) %>% as.data.frame()
posterior <- as.array(m1)
# 
# # Define matching palette
# pal2 <- alpha(pal, alpha = 0.8)
# 
# color_scheme_set(rep("white", 6)) # This is to be able to have a fill color with alpha
# 
# Linf_warm <- mcmc_dens(posterior, pars = c("b_LinfW_Intercept"),
#           facet_args = list(nrow = 2)) + 
#   geom_density(fill = pal2[1], color = NA) + 
#   geom_vline(xintercept = m1_fe$Estimate[6], linetype = 1, color = "white") +
#   geom_vline(xintercept = m1_fe$Q10[6], linetype = 2, color = "white") +
#   geom_vline(xintercept = m1_fe$Q90[6], linetype = 2, color = "white") +
#   coord_cartesian(xlim = c(27, 75)) +
#   #labs(x = expression(paste(italic(L[inf]), " [cm]")), y = "") +
#   labs(x = "", y = "") + 
#   #annotate("text", -Inf, Inf, label = "Warm", size = 5, hjust = -0.5, vjust = 1.3) +
#   theme(text = element_text(size = 12)) 
# 
# Linf_cold <- mcmc_dens(posterior, pars = c("b_LinfC_Intercept"),
#                      facet_args = list(nrow = 2)) + 
#   geom_density(fill = pal2[2], color = NA) + 
#   geom_vline(xintercept = m1_fe$Estimate[5], linetype = 1, color = "white") +
#   geom_vline(xintercept = m1_fe$Q10[5], linetype = 2, color = "white") +
#   geom_vline(xintercept = m1_fe$Q90[5], linetype = 2, color = "white") +
#   coord_cartesian(xlim = c(27, 75)) +
#   labs(x = expression(paste(italic(L[inf]), " [cm]")), y = "") + 
#   #annotate("text", -Inf, Inf, label = "Cold", size = 5, hjust = -0.5, vjust = 1.3) +
#   theme(text = element_text(size = 12))
# 
# K_warm <- mcmc_dens(posterior, pars = c("b_KW_Intercept"),
#                      facet_args = list(nrow = 2)) + 
#   geom_density(fill = pal2[1], color = NA) + 
#   geom_vline(xintercept = m1_fe$Estimate[4], linetype = 1, color = "white") +
#   geom_vline(xintercept = m1_fe$Q10[4], linetype = 2, color = "white") +
#   geom_vline(xintercept = m1_fe$Q90[4], linetype = 2, color = "white") +
#   coord_cartesian(xlim = c(0.08, 0.28)) +
#   #xlab(expression(paste(italic(K), " [", yr^-1,"]", sep = ""))) + 
#   xlab("") + 
#   theme(text = element_text(size = 12)) 
# 
# K_cold <- mcmc_dens(posterior, pars = c("b_KC_Intercept"),
#                      facet_args = list(nrow = 2)) + 
#   geom_density(fill = pal2[2], color = NA) + 
#   geom_vline(xintercept = m1_fe$Estimate[3], linetype = 1, color = "white") +
#   geom_vline(xintercept = m1_fe$Q10[3], linetype = 2, color = "white") +
#   geom_vline(xintercept = m1_fe$Q90[3], linetype = 2, color = "white") +
#   coord_cartesian(xlim = c(0.08, 0.28)) +
#   xlab(expression(paste(italic(K), " [", yr^-1,"]", sep = ""))) + 
#   theme(text = element_text(size = 12))
# Linf_warm + K_warm + Linf_cold + K_cold
# ggsave("figures/supp/K_Linf_posterior.png", width = 6.5, height = 6.5, dpi = 600)

# http://mjskay.github.io/tidybayes/articles/tidy-brms.html
post_K <- 
  m1 %>%
  gather_draws(b_KC_Intercept, b_KW_Intercept) %>%
  ggplot(aes(x = .value, fill = .variable, color = .variable)) +
  stat_halfeye(alpha = 0.5, size = 5, .width = c(0.7)) +
  guides(fill = guide_legend(override.aes = list(size = 1, shape = NA, linetype = 0)),
         color = FALSE) + 
  scale_fill_manual(values = rev(pal), labels = c("Cold", "Warm")) +
  scale_color_manual(values = rev(pal)) +
  labs(x = expression(paste(italic(K), " [", yr^-1,"]", sep = "")), fill = "") +
  theme(legend.position = c(0.9, 0.9),
        legend.key.size = unit(0.2, "cm"),
        legend.background = element_blank())

post_L_inf <- 
  m1 %>%
  gather_draws(b_LinfC_Intercept, b_LinfW_Intercept) %>%
  ggplot(aes(x = .value, fill = .variable, color = .variable)) +
  stat_halfeye(alpha = 0.5, size = 5, .width = c(0.7)) +
  # guides(fill = guide_legend(override.aes = list(size = 1, shape = NA, linetype = 0)),
  #        color = FALSE) +
  guides(fill = FALSE, color = FALSE) + 
  scale_fill_manual(values = rev(pal), labels = c("Cold", "Warm")) +
  scale_color_manual(values = rev(pal)) +
  labs(x = expression(paste(italic(L[infinity]), " [cm]")), fill = "") +
  theme(legend.position = c(0.9, 0.9),
        legend.key.size = unit(0.2, "cm"),
        legend.background = element_blank())


# Plot distribution of differences
# http://mjskay.github.io/tidybayes/articles/tidy-brms.html
diff <- m1 %>%
  spread_draws(b_LinfC_Intercept, b_LinfW_Intercept, b_KC_Intercept, b_KW_Intercept) %>%
  mutate(diff_K = b_KW_Intercept - b_KC_Intercept,
         diff_L_inf = b_LinfW_Intercept - b_LinfC_Intercept) 

prop_diff_K <- summarise(diff, Proportion_of_the_difference_below_0 = sum(diff_K < 0) / length(diff_K))
prop_diff_L_inf <- summarise(diff, Proportion_of_the_difference_below_0 = sum(diff_L_inf < 0) / length(diff_L_inf))

# https://bookdown.org/content/3890/interactions.html
post_diff_K <- ggplot(diff, aes(x = diff_K, fill = stat(x > 0))) +
  stat_halfeye(alpha = 0.5, size = 5, .width = 0) +
  guides(fill = guide_legend(override.aes = list(size = 1, shape = NA, linetype = 0)), color = FALSE) + 
  scale_fill_manual(values = c("grey10", "grey70")) +
  annotate("text", 0.11, 0.95, size = 3, label = paste("prop. diff<0=", round(prop_diff_K, 2), sep = "")) +
  labs(x = expression(~italic(K[warm])~-~italic(K[cold]))) +
  theme(legend.position = c(0.2, 0.8),
        legend.key.size = unit(0.2, "cm"),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 10),
        legend.background = element_blank())

post_diff_K

post_diff_L_inf <- ggplot(diff, aes(x = diff_L_inf, fill = stat(x > 0))) +
  stat_halfeye(alpha = 0.5, size = 5, .width = 0) +
  guides(fill = guide_legend(override.aes = list(size = 1, shape = NA, linetype = 0)), color = FALSE) + 
  scale_fill_manual(values = c("grey10", "grey70")) +
  annotate("text", 22, 0.95, size = 3, label = paste("prop. diff<0=", round(prop_diff_L_inf, 2), sep = "")) +
  labs(x = expression(paste(~italic(L[infinity][warm])~-~italic(L[infinity][cold])))) +
  theme(legend.position = c(0.15, 0.8),
        legend.key.size = unit(0.2, "cm"),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 10),
        legend.background = element_blank())

post_diff_L_inf

# Plotting all together
# pWord1 / (Linf_warm + K_warm + Linf_cold + K_cold) +
#   plot_layout(heights = c(2, 1)) +
#   plot_annotation(tag_levels = 'A')

pWord1 / ((post_K/post_diff_K) | (post_L_inf/post_diff_L_inf)) +
  plot_layout(heights = c(1.2, 1)) +
  plot_annotation(tag_levels = 'A')

ggsave("figures/vbge_pred_K_Linf_post.png", width = 6.5, height = 6.5, dpi = 600)


##### Random year effects ==========================================================
# http://mjskay.github.io/tidybayes/articles/tidy-brms.html

# Plot predictions by cohort:
p2 <- dfm %>%
  data_grid(age = seq_range(age, by = 1),
            birth_year = seq_range(birth_year, by = 1),
            area = c("FM", "BT")) %>%
  mutate(areaC = ifelse(area == "FM", 1, 0),
         areaW = ifelse(area == "BT", 1, 0)) %>%
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
  labs(y = "Length [cm]", x = "Age [yrs]", fill = "Area", colour = "Area") +
  NULL

pWord2 <- p2 + theme(text = element_text(size = 12), # 12 for word doc
                     legend.position = c(0.7, 0.1),
                     legend.title = element_text(size = 12),
                     legend.text = element_text(size = 12))

ggsave("figures/supp/vbge_pred_year.png", width = 6.5, height = 6.5, dpi = 600)

# Cohort-specific VBGE parameters
get_variables(m1)

# Warm K
pKW <- m1 %>%
  spread_draws(b_KW_Intercept,
               r_birth_year__KW[birth_year, Intercept]) %>%
  mutate(year_mean_KW = b_KW_Intercept + r_birth_year__KW) %>% # The random effects are offsets
  ggplot(aes(y = factor(birth_year), x = year_mean_KW)) +
  stat_halfeye(fill = pal2[1], alpha = 0.8) + 
  labs(y = "Year", x = expression(paste(italic(K), " [", yr^-1,"]", sep = ""))) + 
  ggtitle("Warm")

# Cold K
pKC <- m1 %>%
  spread_draws(b_KC_Intercept,
               r_birth_year__KC[birth_year, Intercept]) %>%
  mutate(year_mean_KC = b_KC_Intercept + r_birth_year__KC) %>% # The random effects are offsets
  ggplot(aes(y = factor(birth_year), x = year_mean_KC)) +
  stat_halfeye(fill = pal2[2], alpha = 0.8) + 
  labs(y = "Year", x = expression(paste(italic(K), " [", yr^-1,"]", sep = ""))) + 
  ggtitle("Cold")

pKW + pKC

ggsave("figures/supp/vbge_random_K.png", width = 6.5, height = 6.5, dpi = 600)

# Warm L_inf
pLinfW <- m1 %>%
  spread_draws(b_LinfW_Intercept,
               r_birth_year__LinfW[birth_year, Intercept]) %>%
  mutate(year_mean_LinfW = b_LinfW_Intercept + r_birth_year__LinfW) %>% # The random effects are offsets
  ggplot(aes(y = factor(birth_year), x = year_mean_LinfW)) +
  stat_halfeye(fill = pal2[1], alpha = 0.8) + 
  labs(y = "Year", x = expression(paste(italic(L[inf]), " [cm]"))) + 
  coord_cartesian(xlim = c(29, 140)) +
  ggtitle("Warm")

# Cold L_inf
pLinfC <- m1 %>%
  spread_draws(b_LinfC_Intercept,
               r_birth_year__LinfC[birth_year, Intercept]) %>%
  mutate(year_mean_LinfC = b_LinfC_Intercept + r_birth_year__LinfC) %>% # The random effects are offsets
  ggplot(aes(y = factor(birth_year), x = year_mean_LinfC)) +
  stat_halfeye(fill = pal2[2], alpha = 0.8) + 
  labs(y = "Year", x = expression(paste(italic(L[inf]), " [cm]"))) + 
  coord_cartesian(xlim = c(29, 140)) +
  ggtitle("Cold")

pLinfW + pLinfC

ggsave("figures/supp/vbge_random_Linf.png", width = 6.5, height = 6.5, dpi = 600)


##### Model diagnostics & fit ======================================================
pal_diag <- rev(brewer.pal(n = 3, name = "Dark2"))

# Chain convergence
posterior <- as.array(m1)
dimnames(posterior)

d1 <- mcmc_trace(posterior,
                 pars = c("b_t0C_Intercept", "b_t0W_Intercept", "b_KC_Intercept", 
                          "b_KW_Intercept", "b_LinfC_Intercept", "b_LinfW_Intercept",
                          "sd_birth_year__KC_Intercept", "sd_birth_year__KW_Intercept",
                          "sd_birth_year__LinfC_Intercept", "sd_birth_year__LinfW_Intercept",
                          "sigma"),#, "nu"),
                 facet_args = list(ncol = 3, strip.position = "left")) + 
  theme(text = element_text(size = 12),
        strip.text = element_text(size = 4),
        legend.position = "top") + 
  scale_color_manual(values = alpha(pal_diag, alpha = 0.8))

# Resid vs fitted
d2 <- dfm %>%
  add_residual_draws(m1) %>%
  ggplot(aes(x = .row, y = .residual)) +
  stat_pointinterval(alpha = 0.5, size = 0.7) + 
  theme(text = element_text(size = 12))

# qq-plot
# d3 <- dfm %>%
#   add_residual_draws(m1) %>%
#   median_qi() %>%
#   ggplot(aes(sample = .residual)) +
#   geom_qq_line() +
#   geom_qq(alpha = 0.8) +
#   theme(text = element_text(size = 12))

# Student QQ plot
# https://stackoverflow.com/questions/42493048/computation-failed-for-stat-summary-what-must-be-a-character-string-or-a-func
# https://www.seascapemodels.org/rstats/2017/10/06/qqplot-non-normal-glm.html

summary(m1) # Extract "fixed" effects from m2 for plotting the equation 
nu <- 4.34

# "Base" version
# t <- dfm_dummy %>%
#  add_residual_draws(m3s) %>%
#  median_qi()
# resids <- t$.residual
# n <- nrow(dfm_dummy)
# qqplot(qt(ppoints(n), df = nu), resids,
# xlab = "Theoretical quantile", ylab = "residuals")
# qqline(resids, lty = 2)

# Below ggplot version (check they are the same!)
#?geom_qq_line. Does not take a df argument but dparams, a bit strange
d3 <- dfm %>%
  add_residual_draws(m1) %>%
  median_qi() %>%
  ggplot(aes(sample = .residual)) +
  geom_qq_line(distribution = qt, dparams = nu) +
  geom_qq(alpha = 0.8, distribution = qt, dparams = nu) +
  theme(text = element_text(size = 12))

# Posterior predictive
d4 <- pp_check(m1) + 
  theme(text = element_text(size = 12),
        legend.position = c(0.15, 0.95),
        legend.background = element_rect(fill = NA)) + 
  scale_color_manual(values = rev(pal_diag)) +
  labs(color = "")

d1 / (d2 / (d3 + d4)) + 
  plot_annotation(tag_levels = 'A')

ggsave("figures/supp/vbge_diag_fit2.png", width = 6.5, height = 10.5, dpi = 600)


