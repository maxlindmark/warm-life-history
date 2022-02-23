#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 2020.09.02: Max Lindmark
#
# Fit allometric growth models of the form: log(G) ~ log(L) * Area
# with individual and year-varying intercepts, and use WAIC to compare them.
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
library(janitor)

# Print package versions for versions
sessionInfo() 

# For parallel processing
options(mc.cores = parallel::detectCores()) 


# B. READ DATA =====================================================================
#** Rearrange and calculate specific growth ========================================
# df <- read.csv("data/Growth_data_BT_FM_1970-2004.csv", sep = ";") # This is the 
# original data from Huss et al (2019)
df <- read.csv("data/size_at_age_BT_FM_1970-2004.csv", sep = ";")

# test
df %>% arrange(length)
df %>% arrange(desc(length))

max(df$length)

# Check out how many ("length") unique ID observations there are
length(unique(df$ID)) 

# Remove gear 32, see vbge script
df <- df %>% filter(!gear == 32)

# I need to make the ID completely unique (since I should include ID information in my
# mixed model), otherwise they may have the same ID in different areas.
# We solve this by writing ID with area + ID. The paste function pastes things together
df$ID <- paste(df$ID, df$area, sep = "")

# Check catch age corresponds to number of read length-at-ages!
df <- df %>%
  group_by(ID) %>%
  mutate(n = n()) %>%
  mutate(length = length/10) %>%
  rename("back_calc_age" = "age") %>% 
  ungroup()

ggplot(df, aes(factor(catch_age), factor(n))) + geom_point()

# Remove individuals with the same ID (bug)
df <- df %>% filter(n <= catch_age)

# Now every individual has the same number of rows as it's age in years
ggplot(df, aes(factor(catch_age), factor(n))) + geom_point()

# Now I create a "wide" data frame, so that I can easily create new columns by taking
# one column minus another
df_wide <- df %>%
  filter(birth_year > 1980) %>%
  spread(back_calc_age, length) %>% # Here I convert values in the "age" column to my new columns!
  # The value in these columns is taken from the "length" column
  arrange(ID)

# Check a single ID
df %>% filter(ID == "1983127BT")
df_wide %>% filter(ID == "1983127BT")

# Calculation of growth and geometric length for all ages
df_growth <- data.frame(df_wide %>%
                          mutate(G1 = 100 * (log(`2`) - log(`1`)),
                                 G2 = 100 * (log(`3`) - log(`2`)),
                                 G3 = 100 * (log(`4`) - log(`3`)),
                                 G4 = 100 * (log(`5`) - log(`4`)),
                                 G5 = 100 * (log(`6`) - log(`5`)),
                                 G6 = 100 * (log(`7`) - log(`6`)),
                                 G7 = 100 * (log(`8`) - log(`7`)),
                                 G8 = 100 * (log(`9`) - log(`8`)), # 9 is the max age in the data!
                                 L1 = `1`,
                                 L2 = `2`,
                                 L3 = `3`,
                                 L4 = `4`,
                                 L5 = `5`,
                                 L6 = `6`,
                                 L7 = `7`,
                                 L8 = `8`))

# Check it went OK
df_wide %>% filter(ID == "1983127BT")
df_growth %>% filter(ID == "1983127BT")

# Now make the data "long" again (each observation is a row),
# so that we can plot growth ~ length
# We need to do it separately (length and growth) to avoid duplicates

# First subset the important columns... length and growth data separately
df_g <- df_growth %>% select(c(area, ID, catch_age, G1, G2, G3, G4, G5, G6, G7, G8))

df_l <- df_growth %>% select(c(area, ID, catch_age, L1, L2, L3, L4, L5, L6, L7, L8))

# Now we have two wide data frames... let's make them long using the gather function
# - separately - and then merge them
df_g_l <- df_g %>% gather(g_age, growth, 4:11) # columns 4-11 are gathered
df_l_l <- df_l %>% gather(g_length, length, 4:11)

# Now join datasets. Before I do that I need to add a new common column ("age"), so that
# R knows which length to go with each growth. Check data frames again...
head(arrange(subset(df_g_l, catch_age == 3), ID), 5)
head(arrange(subset(df_l_l, catch_age == 3), ID), 5)

# Here I'm splitting the L1, G1 stuff so that I get G and 1, because then the numbers
# will match in the two datasets, and I make sure the correct length is matched with growth
df_g_l <- df_g_l %>% separate(g_age, c("g", "back_calc_age"), sep = 1) %>% arrange(ID)
df_l_l <- df_l_l %>% separate(g_length, c("g", "back_calc_age"), sep = 1) %>% arrange(ID)

# Now I'll do a left_join to add in length to growth data. Match by ID and age
df_l_l_subset <- df_l_l %>% select(ID, back_calc_age, length)
df_all <- left_join(df_g_l, df_l_l_subset, by = c("ID", "back_calc_age"))

# Check it went OK
df_growth %>% filter(ID == "1983127BT")
df_all %>% filter(ID == "1983127BT")

# Remove NA growth and length individuals from this analysis. Calculate other things as well
# centering and squaring etc
dfm <- df_all %>%
  drop_na(growth) %>%
  filter(growth > 0) %>%
  drop_na(length) %>% 
  mutate(log_length = log(length),
         log_growth = log(growth),
         log_length_ct = log_length - mean(log_length),
         log_length_ct_sq = log_length_ct*log_length_ct) %>% 
  separate(ID, c("catch_year", "ID2"), sep = 4) %>% 
  mutate_at(c("catch_year", "catch_age", "back_calc_age"), as.numeric) %>% 
  mutate(birth_year = catch_year - catch_age,
         ID = paste(catch_year, ID2, sep = "")) %>% 
  select(-g) %>% 
  filter(birth_year < 1998) %>%
  filter(catch_year < 2003) %>% 
  mutate_at(c("ID", "ID2", "area"), as.factor)


#** Plot data ======================================================================
# Plot distribution of catch years
dfm %>% 
  group_by(catch_age, catch_year, area) %>% 
  summarise(n = n()) %>% 
  ungroup() %>% 
  ggplot(., aes(factor(catch_age), n, fill = area)) +
  geom_bar(stat = "identity") + 
  facet_wrap(~catch_year, scales = "free")

# Filter to have at least 4 data points per individual
dfm_3 <- dfm %>% filter(catch_age > 3)
dfm_4 <- dfm %>% filter(catch_age > 4)

# Plot distribution of data
p3 <- ggplot(dfm_3, aes(x = log_growth)) + geom_density()
p4 <- ggplot(dfm_4, aes(x = log_growth)) + geom_density()

p3 + p4

dfm <- dfm %>% filter(catch_age > 4)

dfm %>% 
  group_by(catch_age, catch_year, area) %>% 
  summarise(n = n()) %>% 
  ungroup() %>% 
  ggplot(., aes(factor(catch_age), n, fill = area)) +
  geom_bar(stat = "identity") + 
  facet_wrap(~catch_year, scales = "free")

# Plot sample size per individual
min(dfm$catch_age)

dfm %>% 
  group_by(ID) %>% 
  mutate(n = n()) %>% 
  ungroup() %>% 
  ggplot(., aes(n)) +
  geom_histogram() 

# Check relationship between catch age and # of back-calculated ages
dfm %>% 
  group_by(ID) %>%
  mutate(n = n()) %>%
  ggplot(., aes(n, catch_age)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1) +
  geom_abline(intercept = 1, slope = 1, col = "red")

# Some data without growth, hence log 0
dfm <- dfm %>% group_by(ID) %>% mutate(n = n()) %>% mutate(test = catch_age - n) %>% ungroup()

# Remove these
dfm <- dfm %>% filter(test == 1)

# Create area 2 variable
dfm <- dfm %>%
  mutate(area2 = ifelse(area == "BT", "Warm", "Cold"))

# Plot again
dfm %>% 
  group_by(catch_age, catch_year, area) %>% 
  summarise(n = n()) %>% 
  ungroup() %>% 
  ggplot(., aes(factor(catch_age), n, fill = area)) +
  geom_bar(stat = "identity") + 
  facet_wrap(~catch_year, scales = "free")

# And now again but with sample size per age by catch year
# No catch data in FM in 1986
dfm %>% 
  group_by(back_calc_age, catch_year, area) %>% 
  summarise(n = n()) %>% 
  ungroup() %>% 
  ggplot(., aes(factor(back_calc_age), n, fill = area)) +
  geom_bar(stat = "identity") + 
  facet_wrap(~catch_year)

# Same but with birth year
dfm %>% 
  group_by(back_calc_age, birth_year, area) %>% 
  summarise(n = n()) %>% 
  ungroup() %>% 
  ggplot(., aes(factor(back_calc_age), n, fill = area)) +
  geom_bar(stat = "identity") + 
  facet_wrap(~birth_year)

# Calculate samples sizes
dfm %>% group_by(area) %>% summarise(n = n())
dfm %>% group_by(area) %>% distinct(ID) %>% summarise(n = n())

# Average number of data points per individual?
dfm %>% group_by(ID) %>% summarise(n = n()) %>% ungroup() %>% summarize(mean_n = mean(n))

# Now plot full data.
dfm %>% 
  ggplot(., aes(x = log_length, y = log_growth, color = ID)) +
  facet_wrap(~area) +
  geom_point(size = 1, alpha = 0.2) +
  geom_line(size = 1, alpha = 0.2) +
  scale_color_viridis(discrete = T, direction = -1) + 
  guides(color = FALSE) +
  NULL

# Check unique years
dfm %>% group_by(area, catch_year) %>% summarise(catch_year = (unique(catch_year))) %>% data.frame()

# Check sample size:
nrow(dfm)
length(unique(dfm$ID))

nrow(dfm) / length(unique(dfm$ID))

max(dfm$catch_year)


# C. FIT MODELS ====================================================================
# I'm following the multilevel model vignette here: https://cran.r-project.org/web/packages/brms/index.html
# And this one for specifying random structure: https://discourse.mc-stan.org/t/levels-within-levels/8814/3
# https://stats.stackexchange.com/questions/400700/why-do-we-do-crossed-vs-nested-vs-other-random-effects
# https://stackoverflow.com/questions/29717308/how-do-i-code-the-individual-in-to-an-lme4-nested-model

# https://cran.r-project.org/web/packages/insight/vignettes/insight.html
# See the insight package for verifying model structure

# Non-linear models ================================================================
# Due to non-optimal QQ-plots with the log-linear, despite adding quad terms etc. I
# here instead fit a non-linear model, dummy coded as the VBGE model

# Follow the method in VBGE and use dummy coding
bt <- filter(dfm, area == "BT")
fm <- filter(dfm, area == "FM")

dfm_dummy <- data.frame(rbind(cbind(bt, areaW=1, areaC=0), cbind(fm, areaW=0, areaC=1)))

# Plot non-linear relationship
ggplot(dfm_dummy, aes(length, growth, color = area2)) +
  geom_point()

# Simulate from the prior predictive distribution
# M0: Prior predictive check: Warm+Cold merged =====================================
# Define priors
prior0 <-
  prior(normal(500, 100), nlpar = "alpha") +
  prior(normal(-1.2, 0.3), nlpar = "theta")
  
M0fmbt <- brm(
  bf(growth ~ alpha*length^theta, alpha ~ 1, theta ~ 1, nl = TRUE),
  data = dfm_dummy, family = student(),
  prior = prior0,
  sample_prior = "only", 
  iter = 4000, thin = 1, cores = 3, chains = 3, seed = 9)

# From add_fitted_draws {tidybayes}	which I use for the general predictions
# add_predicted_draws adds draws from posterior predictions to the data. It corresponds to ... or brms::predict.brmsfit() in brms.
pp <- conditional_effects(M0fmbt, method = "posterior_predict")

plot(pp, plot = FALSE)[[1]] +
  labs(y = expression(paste("Growth [%", yr^-1, "]")), x = "Length [cm]")

ggsave("figures/supp/growth_prior_pred_check.png", width = 6.5, height = 6.5, dpi = 600)


# m1 ===============================================================================
# all parameters area-specific 
prior <-
  prior(normal(500, 100), nlpar = "alphaW") +
  prior(normal(500, 100), nlpar = "alphaC") +
  prior(normal(-1.2, 0.3), nlpar = "thetaW") +
  prior(normal(-1.2, 0.3), nlpar = "thetaC")

ptm <- proc.time()
m1 <- brm(bf(growth ~ areaW*alphaW*length^thetaW + areaC*alphaC*length^thetaC, 
             alphaW ~ 1 + (1|birth_year/ID), # parameter varying by ID within birth_year
             alphaC ~ 1 + (1|birth_year/ID), # parameter varying by ID within birth_year
             thetaW + thetaC ~ 1, nl = TRUE),
           family = student(),
           data = dfm_dummy, prior = prior, iter = 4000, cores = 3, chains = 3,
           seed = 9, save_pars = save_pars(all = TRUE),
           control = list(adapt_delta = 0.99))
proc.time() - ptm
# user   system  elapsed 
# 95.962   16.307 6181.922 

# summary(m1)
# plot(m1)
# prior_summary(m1)

# Save model object to not have to rerun it...
# saveRDS(m1, "output/growth_scaling/m1.rds")
# m1 <- readRDS("output/growth_scaling/m1.rds")

prior_summary(m1)
# prior class      coef         group resp dpar  nlpar bound       source
# normal(500, 100)     b                                   alphaC               user
# normal(500, 100)     b Intercept                         alphaC       (vectorized)
# normal(500, 100)     b                                   alphaW               user
# normal(500, 100)     b Intercept                         alphaW       (vectorized)
# normal(-1.2, 0.3)     b                                   thetaC               user
# normal(-1.2, 0.3)     b Intercept                         thetaC       (vectorized)
# normal(-1.2, 0.3)     b                                   thetaW               user
# normal(-1.2, 0.3)     b Intercept                         thetaW       (vectorized)
# gamma(2, 0.1)    nu                                                     default
# student_t(3, 0, 13.3)    sd                                   alphaC            default
# student_t(3, 0, 13.3)    sd                                   alphaW            default
# student_t(3, 0, 13.3)    sd              birth_year           alphaC       (vectorized)
# student_t(3, 0, 13.3)    sd Intercept    birth_year           alphaC       (vectorized)
# student_t(3, 0, 13.3)    sd              birth_year           alphaW       (vectorized)
# student_t(3, 0, 13.3)    sd Intercept    birth_year           alphaW       (vectorized)
# student_t(3, 0, 13.3)    sd           birth_year:ID           alphaC       (vectorized)
# student_t(3, 0, 13.3)    sd Intercept birth_year:ID           alphaC       (vectorized)
# student_t(3, 0, 13.3)    sd           birth_year:ID           alphaW       (vectorized)
# student_t(3, 0, 13.3)    sd Intercept birth_year:ID           alphaW       (vectorized)
# student_t(3, 0, 13.3) sigma                                                     default


# m2 ===============================================================================
# m2 has a common theta

# Define priors
prior2 <-
  prior(normal(500, 100), nlpar = "alphaW") +
  prior(normal(500, 100), nlpar = "alphaC") +
  prior(normal(-1.2, 0.3), nlpar = "theta")
 
ptm <- proc.time()
m2 <- brm(bf(growth ~ areaW*alphaW*length^theta + areaC*alphaC*length^theta, 
             alphaW ~ 1 + (1|birth_year/ID), # parameter varying by ID within birth_year
             alphaC ~ 1 + (1|birth_year/ID), # parameter varying by ID within birth_year
             theta ~ 1, nl = TRUE),
           family = student(),
           data = dfm_dummy, prior = prior2, iter = 4000, cores = 3, chains = 3,
           seed = 9, save_pars = save_pars(all = TRUE),
           control = list(adapt_delta = 0.99))
proc.time() - ptm


# summary(m2)
# plot(m2)
# prior_summary(m2)

# Save model object to not have to rerun it...
# saveRDS(m2, "output/growth_scaling/m2.rds")
# m2 <- readRDS("output/growth_scaling/m2.rds")


# D. COMPARE MODELS ================================================================
# Compare models: https://mc-stan.org/loo/articles/loo2-example.html
# Expected log pointwise predictive density

# loo_m1 <- loo(m1, moment_match = TRUE)
# loo_m2 <- loo(m2, moment_match = TRUE)
loo_m1 <- loo(m1, moment_match = TRUE)
loo_m2 <- loo(m2, moment_match = TRUE)

plot(loo_m1)
abline(a = 0.7, b = 0)

plot(loo_m2)
abline(a = 0.7, b = 0)

loo_compare(loo_m1, loo_m2)
# elpd_diff se_diff
# m1  0.0       0.0   
# m2 -2.7       4.4 


# E. PRODUCE FIGURES ===============================================================
##### Plot predictions =============================================================
pal <- brewer.pal(n = 6, name = "Paired")[c(2, 6)]

as.data.frame(fixef(m1)) # Extract "fixed" effects from m1 for plotting the equation 

dfm_dummy <- dfm_dummy %>% 
  mutate(area2_plot = ifelse(area2 == "Cold", "Ref", "Heat"))
  
pscatter <- dfm_dummy %>%
  ungroup() %>%
  data_grid(length = seq_range(length, n = 101),
            area2 = c("Warm", "Cold")) %>%
  mutate(areaC = ifelse(area2 == "Cold", 1, 0),
         areaW = ifelse(area2 == "Warm", 1, 0)) %>%
  add_predicted_draws(m1, re_formula = NA) %>%
  mutate(area2_plot = ifelse(area2 == "Cold", "Ref", "Heat")) %>%
  ggplot(aes(x = length, y = growth, color = area2_plot, fill = area2_plot)) +
  stat_lineribbon(aes(y = .prediction), .width = c(.90), alpha = 1/4) +
  geom_point(data = dfm_dummy, alpha = 0.05, size = 0.8) +
  stat_lineribbon(aes(y = .prediction), .width = 0, alpha = 0.8) +
  scale_fill_manual(values = rev(pal)) +
  scale_color_manual(values = rev(pal)) +
  guides(fill = "none",
         color = guide_legend(override.aes = list(linetype = 0, size = 3, fill = NA,
                                                  shape = 16, alpha = 0.5))) +
  labs(y = expression(paste("Growth [%", yr^-1, "]")),
       x = "Length [cm]", fill = "Area", colour = "Area") +
  annotate("text", 35, 42, label = paste("n=", nrow(dfm), sep = ""), size = 3.5) +
  annotate("text", 35, 36, size = 3.5, color = pal[1],
           label = expression(italic("y=433.45×"~length^-1.18))) + # Cold
  annotate("text", 35, 30, size = 3.5, color = pal[2],
           label = expression(italic("y=509.69×"~length^-1.13))) + # Warm
  theme(text = element_text(size = 12), 
        legend.position = c(0.9, 0.9), 
        legend.spacing.y = unit(0, 'cm'),
        legend.key.size = unit(0, "cm"),
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 10))

# Plot posteriors of parameters
# http://mjskay.github.io/tidybayes/articles/tidy-brms.html

post_alpha <- 
  m1 %>%
  gather_draws(b_alphaC_Intercept, b_alphaW_Intercept) %>%
  ggplot(aes(x = .value, fill = .variable, color = .variable)) +
  stat_halfeye(alpha = 0.5, size = 5, .width = c(0.9)) +
  guides(color = FALSE, fill = FALSE) +
  scale_fill_manual(values = pal, labels = c("Cold", "Warm")) +
  scale_color_manual(values = pal) +
  labs(x = expression(italic(alpha)), fill = "") +
  theme(legend.position = c(0.9, 0.9),
        legend.key.size = unit(0.2, "cm"),
        legend.background = element_blank())

post_theta <- 
  m1 %>%
  gather_draws(b_thetaC_Intercept, b_thetaW_Intercept) %>%
  ggplot(aes(x = .value, fill = .variable, color = .variable)) +
  stat_halfeye(alpha = 0.5, size = 5, .width = c(0.9)) +
  guides(fill = FALSE, color = FALSE) + 
  scale_fill_manual(values = pal, labels = c("Cold", "Warm")) +
  scale_color_manual(values = pal) +
  labs(x = expression(italic(theta)), fill = "") +
  theme(legend.position = c(0.9, 0.9),
        legend.key.size = unit(0.2, "cm"),
        legend.background = element_blank())

# Plot distribution of differences: see statistical rethinging 2 p.157 and:
# https://bookdown.org/content/3890/interactions.html 
# http://mjskay.github.io/tidybayes/articles/tidy-brms.html
diff <- m1 %>%
  spread_draws(b_alphaC_Intercept, b_alphaW_Intercept, b_thetaC_Intercept, b_thetaW_Intercept) %>%
  mutate(diff_alpha = b_alphaW_Intercept - b_alphaC_Intercept,
         diff_theta = b_thetaW_Intercept - b_thetaC_Intercept) 

prop_diff_alpha <- diff %>% 
  summarise(Proportion_of_the_difference_below_0 = sum(diff_alpha < 0) / length(diff_alpha))

prop_diff_theta <- diff %>% 
  summarise(Proportion_of_the_difference_below_0 = sum(diff_theta < 0) / length(diff_theta))

# https://bookdown.org/content/3890/interactions.html
post_diff_alpha <- ggplot(diff, aes(x = diff_alpha, fill = stat(x > 0))) +
  stat_halfeye(alpha = 0.5, size = 5, .width = 0) +
  guides(fill = guide_legend(override.aes = list(size = 1, shape = NA, linetype = 0)), color = "none") + 
  scale_fill_manual(values = c("grey10", "grey70")) +
  annotate("text", 130, 0.95, size = 3, label = paste("prop. diff<0=", round(prop_diff_alpha, 3), sep = "")) +
  labs(x = expression(~italic(alpha[heat])~-~italic(alpha[ref]))) +
  theme(legend.position = c(0.2, 0.7),
        legend.key.size = unit(0.2, "cm"),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 10),
        legend.background = element_blank())

post_diff_theta <- ggplot(diff, aes(x = diff_theta, fill = stat(x > 0))) +
  stat_halfeye(alpha = 0.5, size = 5, .width = 0) +
  guides(fill = guide_legend(override.aes = list(size = 1, shape = NA, linetype = 0)), color = "none") + 
  scale_fill_manual(values = c("grey10", "grey70")) +
  annotate("text", 0.07, 0.95, size = 3, label = paste("prop. diff<0=", round(prop_diff_theta, 3), sep = "")) +
  labs(x = expression(~italic(theta[heat])~-~italic(theta[ref]))) +
  theme(legend.position = c(0.2, 0.7),
        legend.key.size = unit(0.2, "cm"),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 10),
        legend.background = element_blank())

# Plotting all together
pscatter / ((post_alpha/post_diff_alpha) | (post_theta/post_diff_theta)) +
  plot_layout(heights = c(1.2, 1)) +
  plot_annotation(tag_levels = 'A')

ggsave("figures/growth_pred.png", width = 6.5, height = 6.5, dpi = 600)


##### Prior vs posterior ===========================================================
# https://discourse.mc-stan.org/t/presenting-influence-of-different-priors/23393
# Refit model and sample prior
prior <-
  prior(normal(500, 100), nlpar = "alphaW") +
  prior(normal(500, 100), nlpar = "alphaC") +
  prior(normal(-1.2, 0.3), nlpar = "thetaW") +
  prior(normal(-1.2, 0.3), nlpar = "thetaC")

ptm <- proc.time()
m1_w_prior <- brm(bf(growth ~ areaW*alphaW*length^thetaW + areaC*alphaC*length^thetaC, 
                     alphaW ~ 1 + (1|birth_year/ID), # parameter varying by ID within birth_year
                     alphaC ~ 1 + (1|birth_year/ID), # parameter varying by ID within birth_year
                     thetaW + thetaC ~ 1, nl = TRUE),
                  family = student(), data = dfm_dummy, prior = prior,
                  iter = 4000, cores = 3, chains = 3, seed = 9,
                  save_pars = save_pars(all = TRUE), sample_prior = "yes",
                  control = list(adapt_delta = 0.99))

# saveRDS(m1_w_prior, "output/growth_scaling/m1_w_prior.rds")
# m1_w_prior <- readRDS("output/growth_scaling/m1_w_prior.rds")

# This is just to check variables names in the samples...
# test <- brm(bf(growth ~ areaW*alphaW*length^thetaW + areaC*alphaC*length^thetaC, 
#                alphaW ~ 1,
#                alphaC ~ 1,
#                thetaW + thetaC ~ 1, nl = TRUE),
#             family = student(), data = dfm_dummy, prior = prior,
#             iter = 1, cores = 1, chains = 1, seed = 9,
#             save_pars = save_pars(all = TRUE), sample_prior = "yes",
#             control = list(adapt_delta = 0.99))
# test %>% posterior_samples() %>% clean_names()

post <- m1_w_prior %>%
  posterior_samples() %>%
  clean_names() %>% 
  dplyr::select(b_alpha_w_intercept, b_alpha_c_intercept, prior_b_alpha_w, prior_b_alpha_c,
                b_theta_w_intercept, b_theta_c_intercept, prior_b_theta_w, prior_b_theta_c)

post_long <- post %>% pivot_longer(cols = c(1:8), names_to = "Parameter", values_to = "value")

# parameter "alpha"
prior_post_alpha <- post_long %>%
  filter(Parameter %in% c("b_alpha_w_intercept", "b_alpha_c_intercept", "prior_b_alpha_w", "prior_b_alpha_c")) %>% 
  ggplot(., aes(value, fill = Parameter, color = Parameter))+
  geom_density(alpha = 0.4) +
  labs(x = expression(alpha)) +
  coord_cartesian(expand = 0) +
  scale_color_manual(values = c(NA, NA, "gray50", "gray50")) +
  scale_fill_manual(values = c(pal[1], pal[2], NA, NA),
                    labels = c(expression(alpha[cold]), expression(alpha[warm]), 
                               expression(paste(Prior~italic(alpha[cold]))), 
                               expression(paste(Prior~italic(alpha[warm]))))) + 
  guides(color = FALSE,
         fill = guide_legend(override.aes = list(color = c(NA, NA, "gray50", "gray50")))) +
  theme(legend.position = c(0.2, 0.8),
        legend.text.align = 0)

# parameter "theta"
prior_post_theta <- post_long %>%
  filter(Parameter %in% c("b_theta_w_intercept", "b_theta_c_intercept", "prior_b_theta_w", "prior_b_theta_c")) %>% 
  ggplot(., aes(value, fill = Parameter, color = Parameter))+
  geom_density(alpha = 0.4) +
  labs(x = expression(theta)) +
  coord_cartesian(expand = 0) +
  scale_color_manual(values = c(NA, NA, "gray50", "gray50")) +
  scale_fill_manual(values = c(pal[1], pal[2], NA, NA),
                    labels = c(expression(theta[cold]), expression(theta[warm]), 
                               expression(paste(Prior~italic(theta[cold]))), 
                               expression(paste(Prior~italic(theta[warm]))))) + 
  guides(color = FALSE,
         fill = guide_legend(override.aes = list(color = c(NA, NA, "gray50", "gray50")))) +
  theme(legend.position = c(0.2, 0.8),
        legend.text.align = 0)

prior_post_theta_ins <- post_long %>%
  filter(Parameter %in% c("b_theta_w_intercept", "b_theta_c_intercept", "prior_b_theta_w", "prior_b_theta_c")) %>% 
  ggplot(., aes(value, fill = Parameter, color = Parameter))+
  geom_density(alpha = 0.4) +
  labs(x = expression(theta)) +
  coord_cartesian(expand = 0, ylim = c(0, 1.8), xlim = c(-1.25, -1.05)) +
  scale_color_manual(values = c(NA, NA, "gray50", "gray50")) +
  scale_fill_manual(values = c(pal[1], pal[2], NA, NA)) + 
  guides(color = FALSE, fill = FALSE)
# prior_post_theta_ins

top <- prior_post_alpha
bottom <- prior_post_theta + inset_element(prior_post_theta_ins, left = 0.55, bottom = 0.5, right = 0.99, top = 0.99)

top / bottom + plot_annotation(tag_levels = list(c('A', 'B'), ''))

ggsave("figures/supp/growth_prior_post.png", width = 6.5, height = 8.5, dpi = 600)


##### Model diagnostics & fit ======================================================
pal_diag <- rev(brewer.pal(n = 3, name = "Dark2"))

# Chain convergence
posterior <- as.array(m1)
dimnames(posterior)

d1 <- mcmc_trace(posterior,
                 pars = c("b_alphaW_Intercept", "b_alphaC_Intercept", "b_thetaW_Intercept", "b_thetaC_Intercept",
                          "sd_birth_year__alphaW_Intercept", "sd_birth_year:ID__alphaW_Intercept",
                          "sd_birth_year__alphaC_Intercept", "sd_birth_year:ID__alphaC_Intercept", 
                          "sigma", "nu"),
                 facet_args = list(ncol = 3, strip.position = "left")) + 
  theme(text = element_text(size = 12),
        strip.text = element_text(size = 3),
        legend.position = "top") + 
  scale_color_manual(values = alpha(pal_diag, alpha = 0.8))

# Resid vs fitted
# The following two plots exhaust the memory, following this helps:
# https://stackoverflow.com/questions/51248293/error-vector-memory-exhausted-limit-reached-r-3-5-0-macos
d2 <- dfm_dummy %>%
  add_residual_draws(m1) %>%
  ggplot(aes(x = .row, y = .residual)) +
  stat_pointinterval(alpha = 0.5, size = 0.7) + 
  theme(text = element_text(size = 12))

# qq-plot
# d3 <- dfm_dummy %>%
#   add_residual_draws(m1) %>%
#   median_qi() %>%
#   ggplot(aes(sample = .residual)) +
#   geom_qq_line() +
#   geom_qq(alpha = 0.8) +
#   theme(text = element_text(size = 12))

# Student QQ plot
# https://stackoverflow.com/questions/42493048/computation-failed-for-stat-summary-what-must-be-a-character-string-or-a-func
# https://www.seascapemodels.org/rstats/2017/10/06/qqplot-non-normal-glm.html

summary(m1)$spec_pars # Extract "fixed" effects from m1 for plotting the equation 
nu <- summary(m1)$spec_pars[2, 1]
nu

# "Base" version
# t <- dfm_dummy %>%
#  add_residual_draws(m1) %>%
#  median_qi()
# resids <- t$.residual
# n <- nrow(dfm_dummy)
# qqplot(qt(ppoints(n), df = nu), resids,
# xlab = "Theoretical quantile", ylab = "residuals")
# qqline(resids, lty = 2)

# Below ggplot version (check they are the same!)
#?geom_qq_line. Does not take a df argument but dparams, a bit strange
d3 <- dfm_dummy %>%
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

ggsave("figures/supp/growth_diag_fit.png", width = 6.5, height = 8.5, dpi = 600)


