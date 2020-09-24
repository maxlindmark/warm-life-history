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
# Rearrange and calculate specific growth ==========================================
# df <- read.csv("data/Growth_data_BT_FM_1970-2004.csv", sep = ";") # This is the 
# original data from Huss et al (2019)
df <- read.csv("data/size_at_age_BT_FM_1970-2004.csv", sep = ";")

# Check out how many ("length") unique ID observations there are
length(unique(df$ID)) 

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
                                 L1 = (`2`*`1`)^0.5,
                                 L2 = (`3`*`2`)^0.5,
                                 L3 = (`4`*`3`)^0.5,
                                 L4 = (`5`*`4`)^0.5,
                                 L5 = (`6`*`5`)^0.5,
                                 L6 = (`7`*`6`)^0.5,
                                 L7 = (`8`*`7`)^0.5,
                                 L8 = (`9`*`8`)^0.5))

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
dfm <- df_all %>%
  drop_na(growth) %>%
  filter(growth > 0) %>%
  drop_na(length) %>% 
  mutate(log_length = log(length),
         log_growth = log(growth)) %>% 
  separate(ID, c("catch_year", "ID2"), sep = 4) %>% 
  mutate_at(c("catch_year", "catch_age", "back_calc_age"), as.numeric) %>% 
  mutate(birth_year = catch_year - catch_age,
         ID = paste(catch_year, ID2, sep = "")) %>% 
  select(-g) %>% 
  filter(birth_year < 1998) %>%
  filter(catch_year < 2003) %>% # Change in ageing formula
  mutate_at(c("ID", "ID2", "area"), as.factor) %>% 
  mutate(log_length_sq = log_length*log_length)


# Plot data ========================================================================
# Plot distribution of catch years
dfm %>% 
  group_by(catch_age, catch_year, area) %>% 
  summarise(n = n()) %>% 
  ungroup() %>% 
  ggplot(., aes(factor(catch_age), n, fill = area)) +
  geom_bar(stat = "identity") + 
  facet_wrap(~catch_year, scales = "free")

# Filter to have at least 3 data points per individual
dfm <- dfm %>% filter(catch_age > 3)

dfm %>% 
  group_by(catch_age, catch_year, area) %>% 
  summarise(n = n()) %>% 
  ungroup() %>% 
  ggplot(., aes(factor(catch_age), n, fill = area)) +
  geom_bar(stat = "identity") + 
  facet_wrap(~catch_year, scales = "free")

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

summary(dfm)

dfm_test <- dfm %>% filter(birth_year > 1992)

dfm_test %>% 
  ggplot(., aes(x = log_length, y = log_growth, color = ID)) +
  facet_grid(birth_year ~ area) +
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


# Prior predictive check for the log-linear model
# https://discourse.mc-stan.org/t/help-understanding-and-setting-informative-priors-in-brms/9574/12
# First set priors
# To set a prior on the fixed intercept, we need to write 0+intercept, see
# https://www.rensvandeschoot.com/tutorials/brms-priors/ and https://rdrr.io/cran/brms/man/brmsformula.html

# Basic log-log regression
#plot(density(rnorm(n = 100000, mean = 0, sd = 5)), main = "prior b")
#plot(density(rnorm(n = 100000, mean = 5, sd = 5)), main = "prior b") 

priors <- c(set_prior("normal(0, 5)", class = "b", coef = "areaFM"),
            set_prior("normal(-1, 5)", class = "b", coef = "log_length"),
            set_prior("normal(0, 5)", class = "b", coef = "log_length:areaFM"),
            set_prior("normal(5, 5)", class = "b", coef = "Intercept"))

m0a <- brm(bf(log_growth ~ 0 + Intercept + log_length*area + (1 + log_length|birth_year/ID)),
           family = gaussian(), data = dfm, inits = "0", iter = 500, cores = 2, chains = 2,
           prior = priors, sample_prior = "only") 

prior_summary(m0a)

pp_check(m0a, nsamples = 50)

conditional_effects(m0a)
# See this for plotting conditional effects: https://bookdown.org/content/3890/interactions.html

# log-log regression with quadratic term
priors_sq <- c(set_prior("normal(0, 5)", class = "b", coef = "areaFM"),
               set_prior("normal(0, 5)", class = "b", coef = "areaFM:log_length_sq"),
               set_prior("normal(0, 5)", class = "b", coef = "Intercept"),
               set_prior("normal(-1, 5)", class = "b", coef = "log_length"), 
               set_prior("normal(0, 5)", class = "b", coef = "log_length_sq"), 
               set_prior("normal(5, 5)", class = "b", coef = "log_length:areaFM"))

#priors_sq <- c(set_prior("normal(0, 10)", class = "b"))

m0b <- brm(bf(log_growth ~ 0 + Intercept + log_length*area + log_length_sq*area +
                (1 + log_length|birth_year/ID)), 
           family = gaussian(), data = dfm, inits = "0", iter = 500, cores = 2, chains = 2,
           prior = priors_sq,
           sample_prior = "only") 

prior_summary(m0b)

pp_check(m0b, nsamples = 50)

conditional_effects(m0b)


# m1 ==============================================================================
# log_growth ~ log_length*area
# https://cran.r-project.org/web/packages/insight/vignettes/insight.html
# See the insight package for verifying model structure
m1 <- brm(bf(log_growth ~ 0 + Intercept + log_length*area + (1 + log_length|birth_year/ID)),
          family = gaussian(), prior = priors,
          data = dfm, iter = 3000, cores = 2, chains = 2, control = list(adapt_delta = 0.9)) 

summary(m1)
plot(m1)

# Save model object to not have to rerun it...
saveRDS(m1, "output/growth_scaling/m1.rds")
#m1 <- readRDS("output/growth_scaling/m1.rds")

pal <- rev(brewer.pal(n = 6, name = "Paired")[c(2, 6)])

# Plot prediction and data
dfm %>% 
  ungroup() %>% 
  data_grid(log_length = seq_range(log_length, n = 101),
           area = c("FM", "BT")) %>%
  add_predicted_draws(m1, re_formula = NA) %>%
  ggplot(aes(x = log_length, y = log_growth, color = area, fill = area)) +
  stat_lineribbon(aes(y = .prediction), .width = c(.90, .50), alpha = 1/4) +
  geom_jitter(data = dfm, alpha = 0.1, width = 0.4,
              height = 0, size = 0.8) +
  scale_fill_manual(values = pal, labels = c("Warm", "Cold")) +
  scale_color_manual(values = pal, labels = c("Warm", "Cold")) +
  labs(y = "log(growth)", x = "log(length)", fill = "Area", colour = "Area") +
  NULL


# m2 ==============================================================================
# log_growth ~ log_length*area; sigma ~ log_length
# https://cran.r-project.org/web/packages/insight/vignettes/insight.html
# See the insight package for verifying model structure
m2 <- brm(bf(log_growth ~ 0 + Intercept + log_length*area + (1 + log_length|birth_year/ID),
             sigma ~ log_length), family = gaussian(), prior = priors,
          data = dfm, iter = 3000, cores = 2, chains = 2, control = list(adapt_delta = 0.9)) 

summary(m2)
plot(m2)

# Save model object to not have to rerun it...
saveRDS(m2, "output/growth_scaling/m2.rds")
#m2 <- readRDS("output/growth_scaling/m2.rds")

# Plot prediction and data
dfm %>% 
  ungroup() %>% 
  data_grid(log_length = seq_range(log_length, n = 101),
            area = c("FM", "BT")) %>%
  add_predicted_draws(m2, re_formula = NA) %>%
  ggplot(aes(x = log_length, y = log_growth, color = area, fill = area)) +
  stat_lineribbon(aes(y = .prediction), .width = c(.90, .50), alpha = 1/4) +
  geom_jitter(data = dfm, alpha = 0.1, width = 0.4,
              height = 0, size = 0.8) +
  scale_fill_manual(values = pal, labels = c("Warm", "Cold")) +
  scale_color_manual(values = pal, labels = c("Warm", "Cold")) +
  labs(y = "log(growth)", x = "log(length)", fill = "Area", colour = "Area") +
  NULL



# Posterior predictive checks ======================================================
p1 <- pp_check(m1) + ggtitle("log-linear") + theme(aspect.ratio = 1)
p2 <- pp_check(m2) + ggtitle("log-linear + sigma") + theme(aspect.ratio = 1)

p1 + p2  + plot_layout(ncol = 2) 



# E. COMPARE MODELS ================================================================
# Compare models: https://mc-stan.org/loo/articles/loo2-example.html
# Expected log pointwise predictive density

# Log linear models
loo_m1 <- loo(m1)
loo_m2 <- loo(m2)

loo_compare(loo_m1, loo_m2)

