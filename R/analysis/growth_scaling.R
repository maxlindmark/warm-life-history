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

# Now plot full data.
dfm %>% 
  ggplot(., aes(x = log_length, y = log_growth, color = ID)) +
  facet_wrap(~area) +
  geom_point(size = 1, alpha = 0.2) +
  geom_line(size = 1, alpha = 0.2) +
  scale_color_viridis(discrete = T, direction = -1) + 
  guides(color = FALSE) +
  NULL


# C. FIT MODELS ====================================================================
# I'm following the multilevel model vignette here: https://cran.r-project.org/web/packages/brms/index.html
# And this one for specifying random structure: https://discourse.mc-stan.org/t/levels-within-levels/8814/3
# https://stats.stackexchange.com/questions/400700/why-do-we-do-crossed-vs-nested-vs-other-random-effects
# https://stackoverflow.com/questions/29717308/how-do-i-code-the-individual-in-to-an-lme4-nested-model

# https://cran.r-project.org/web/packages/insight/vignettes/insight.html
# See the insight package for verifying model structure

# m1 ==============================================================================
# With quadratic interaction

m1a <- brm(bf(log_growth ~ log_length_ct*area + log_length_ct_sq*area + (1|birth_year/ID),
              sigma ~ log_length_ct),
          family = gaussian(), data = dfm, iter = 4000, cores = 3, chains = 3,
          save_all_pars = TRUE)

summary(m1a)
plot(m1a)
prior_summary(m1a)

# Save model object to not have to rerun it...
#saveRDS(m1a, "output/growth_scaling/m1a.rds")
#m1a <- readRDS("output/growth_scaling/m1a.rds")


# m1b ==============================================================================
# No quadratic interaction
m1b <- brm(bf(log_growth ~ log_length_ct*area + log_length_ct_sq + (1|birth_year/ID),
              sigma ~ log_length_ct),
          family = gaussian(), data = dfm, iter = 4000, cores = 3, chains = 3,
          save_all_pars = TRUE)

summary(m1b)
plot(m1b)
prior_summary(m1b)

# Save model object to not have to rerun it...
#saveRDS(m1b, "output/growth_scaling/m1b.rds")
#m1b <- readRDS("output/growth_scaling/m1b.rds")


# m2a ==============================================================================
# Student model
# With quadratic interaction

m2a <- brm(bf(log_growth ~ log_length_ct*area + log_length_ct_sq*area + (1|birth_year/ID),
              sigma ~ log_length_ct),
          family = student(), data = dfm, iter = 4000, cores = 3, chains = 3,
          save_all_pars = TRUE)

summary(m2a)
plot(m2a)
prior_summary(m2a)

# Save model object to not have to rerun it...
#saveRDS(m2a, "output/growth_scaling/m2a.rds")
#m2a <- readRDS("output/growth_scaling/m2a.rds")


# m2b ==============================================================================
# Student model
# No quadratic interaction

m2b <- brm(bf(log_growth ~ log_length_ct*area + log_length_ct_sq + (1|birth_year/ID),
              sigma ~ log_length_ct),
           family = student(), data = dfm, iter = 4000, cores = 3, chains = 3,
           save_all_pars = TRUE)

summary(m2b)
plot(m2b)
prior_summary(m2b)

# Save model object to not have to rerun it...
#saveRDS(m2b, "output/growth_scaling/m2b.rds")
#m2b <- readRDS("output/growth_scaling/m2b.rds")


# D. COMPARE MODELS ================================================================
# Compare models: https://mc-stan.org/loo/articles/loo2-example.html
# Expected log pointwise predictive density

loo_m1a <- loo(m1a)
loo_m1b <- loo(m1b)
loo_m2a <- loo(m2a)
loo_m2b <- loo(m2b)

loo_compare(loo_m1a, loo_m1b, loo_m2a, loo_m2b)
# > loo_compare(loo_m1a, loo_m1b, loo_m2a, loo_m2b)
# elpd_diff se_diff
# m2b  0.0       0.0   
# m2a -1.0       0.7   
# m1b -6.0       6.7   
# m1a -9.2       6.7   

summary(m2b)


# E. PRODUCE FIGURES ===============================================================
# Plot prediction and data
pal <- rev(brewer.pal(n = 6, name = "Paired")[c(2, 6)])

p1 <- dfm %>%
  ungroup() %>%
  data_grid(log_length_ct = seq_range(log_length_ct, n = 101),
            area = c("FM", "BT")) %>%
  mutate(log_length_ct_sq = log_length_ct*log_length_ct) %>% 
  add_predicted_draws(m2b, re_formula = NA) %>%
  mutate(log_length = log_length_ct + mean(dfm$log_length)) %>% 
  ggplot(aes(x = log_length, y = log_growth, color = area, fill = area)) +
  stat_lineribbon(aes(y = .prediction), .width = c(.90, .50), alpha = 1/4) +
  geom_point(data = dfm, alpha = 0.1, size = 0.8) +
  scale_fill_manual(values = pal, labels = c("Warm", "Cold")) +
  scale_color_manual(values = pal, labels = c("Warm", "Cold")) +
  labs(y = "log(growth)", x = "log(length)", fill = "Area", colour = "Area") +
  NULL

pWord1 <- p1 + theme(text = element_text(size = 12), 
                     legend.position = c(0.9, 0.9), 
                     legend.title = element_text(size = 10),
                     legend.text = element_text(size = 10))

ggsave("figures/growth_scaling/growth_pred.png", width = 6.5, height = 6.5, dpi = 600)

# Posterior predictive checks
pp_check(m2b) + 
  theme(text = element_text(size = 12),
        legend.position = c(0.9, 0.9), 
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 12))
ggsave("figures/supp/growth_ppc.png", width = 6.5, height = 6.5, dpi = 600)

# Posterior predictive checks: summary statistics median
pp_check(m2b, type = "stat", stat = 'median', nsamples = NULL) + 
  theme(text = element_text(size = 12),
        legend.position = c(0.9, 0.9), 
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 12))
ggsave("figures/supp/growth_sumstat_median_ppc.png", width = 6.5, height = 6.5, dpi = 600)

# Posterior predictive checks: summary statistics mean
pp_check(m2b, type = "stat", stat = 'mean', nsamples = NULL) + 
  theme(text = element_text(size = 12),
        legend.position = c(0.9, 0.9), 
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 12))
ggsave("figures/supp/growth_sumstat_mean_ppc.png", width = 6.5, height = 6.5, dpi = 600)

# Chain convergence
posterior <- as.array(m2b)
dimnames(posterior)

color_scheme_set("mix-blue-red")
mcmc_trace(posterior,
           pars = c("b_Intercept", "b_sigma_Intercept", "b_log_length_ct", "b_areaFM", 
                    "b_log_length_ct_sq", "b_log_length_ct:areaFM", "b_sigma_log_length_ct"),
           facet_args = list(ncol = 2, strip.position = "left")) + 
  theme(text = element_text(size = 12),
        #legend.position = c(0.7, 0.1), 
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 12))

ggsave("figures/supp/growth_chain_convergence.png", width = 7.5, height = 7.5, dpi = 600)


