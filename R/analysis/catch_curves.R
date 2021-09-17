#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 2021.08.04: Max Lindmark
#
# Fit catch curves to catch-at-age for both areas with year as random effect
# 
# A. Load libraries
# 
# B. Read data
# 
# C. Fit models
# 
# D. Produce figures
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
df <- read.csv("data/aged_catch_BT_FM_1987-2003.csv")

head(df)

# df %>%
#   filter(Area == "BT") %>%
#   group_by(year, netID) %>%
#   summarise(n = n()) %>%
#   ggplot(., aes(year, n)) + geom_bar(stat="identity")

# df %>% group_by(Area) %>% summarise(mean_age = mean(age), sd_age = sd(age))
# df %>% filter(age > 2) %>% group_by(Area) %>% summarise(mean_age = mean(age), sd_age = sd(age))
# ggplot(df, aes(Area, age)) + geom_boxplot()
# df %>% filter(age > 2) %>% ggplot(., aes(Area, age)) + geom_boxplot()
# summary(lm(age~Area, data = df))
# summary(lm(age~Area, data = filter(df, age > 2)))

# Currently the data is in format 1 row, 1 ind. We want a column that has the CPUE
# by length class.

# How many nets in total per area and year? (For scaling with effort later)
df <- df %>% group_by(Area, year) %>%
  mutate(n_nets_year = length(unique(netID))) %>%
  ungroup() %>%
  as.data.frame()

# Test
df %>% filter(year == 1995 & Area == "BT")
df %>% filter(year == 1995 & Area == "BT") %>% distinct(netID)

# We want the data as follows:
# Year, LngtClass, Age, Number

# Group by year, Area and age, summarize and get n() per length, year and area...
df2 <- df %>%
  group_by(year, Area, age) %>% 
  summarise(catch_n = n()) %>% # Get catch as a column, this is now the number of rows a given length occurs each year and area
  ungroup() %>% 
  mutate(effort_id = paste(year, Area, sep = ".")) %>% # For calculating the total effort that gave that catch (total for year and area) 
  as.data.frame()

# Now we need to get the effort back in there
df_effort <- df %>%
  mutate(effort_id = paste(year, Area, sep = ".")) %>% 
  select(effort_id, n_nets_year) %>% 
  distinct(effort_id, .keep_all = TRUE) %>% 
  as.data.frame()

# Now do a left_join
df3 <- left_join(df2, df_effort, by = "effort_id") %>% as.data.frame()

# Looks Ok!
df %>% filter(year == 1987 & Area == "FM")
df %>% filter(year == 1987 & Area == "FM") %>% distinct(n_nets_year, .keep_all = TRUE)
df2 %>% filter(effort_id == "1987.FM")
df3 %>% filter(effort_id == "1987.FM")

# Go from total catch per year to catch by size-class / number of nets that year 
df4 <- df3 %>% 
  ungroup() %>% 
  rename("area" = "Area") %>% 
  mutate(cpue_numbers = catch_n/n_nets_year) # Get numbers CPUE, divide by the previously create n_nets, which is # of unique net ID's in each area and year
         
# Test I get 1 unique row per age, year and area
df4 %>%
  group_by(year, area, age) %>%
  summarise(n = n()) %>% 
  filter(!n == 1)

# Plot the log cpue as function of age to find the ages that correspond to the descending limb:
ggplot(df4, aes(factor(age), log(cpue_numbers), color = factor(year))) + 
  geom_point() +
  facet_wrap(~ area) +
  scale_color_viridis(discrete = TRUE) + 
  coord_cartesian(expand = 0) + 
  theme(aspect.ratio = 1)

# Let's use age 3 and older
ggplot(df4, aes(factor(age), log(cpue_numbers), color = factor(year))) + 
  geom_point() +
  facet_wrap(~ area) +
  scale_color_viridis(discrete = TRUE) + 
  coord_cartesian(expand = 0) + 
  theme(aspect.ratio = 1)

# Plot total catch
df2 %>%
  group_by(Area, year) %>%
  summarize(tot_catch = sum(catch_n)) %>% 
  ggplot(., aes(year, tot_catch)) + geom_bar(stat = "identity") +
  facet_wrap(~ Area)

# Plot cpue
df4 %>%
  group_by(area, year) %>%
  ggplot(., aes(year, cpue_numbers)) + geom_bar(stat = "identity") +
  facet_wrap(~ area)

# Filter and rename data
# d <- df4 %>% filter(age > 2)

# Something strange in BT 1996, station is wrong, hence netID, total nets and thus CPUE
# Removing it for now!
d <- df4 %>%
  filter(age > 2) %>% 
  mutate(keep = ifelse(area == "BT" & year == 1996, "N", "Y")) %>% 
  filter(keep == "Y")


# C. FIT MODELS ====================================================================
##### Catch curves =================================================================
# See this for notation: https://solomonkurz.netlify.app/post/2020-12-09-multilevel-models-and-the-index-variable-approach/ 
# Fitting models of log catch ~ age with interactive or additive effects of area, 
# using catch year as a random effect

# Edit variables
d <- d %>%
  mutate(log_cpue = log(cpue_numbers),
         area2 = ifelse(area == "BT", "Warm", "Cold"))

# Random intercepts only because of better convergence (posteriors of random year slopes looking funny)
# m1 <- brm(log_cpue ~ -1 + age * area2 + (0 + area2|year),

prior0 <-
  prior(normal(0, 10), class = "b") +
  prior(gamma(2, 0.1), class = "nu") +
  prior(student_t(3, 0, 2.5), class = "sd") +
  prior(student_t(3, 0, 2.5), class = "sigma")

# No random slopes (same slopes for all years, but still area-specific)
m0 <- brm(
  log_cpue ~ -1 + age * area2 + (0 + area2||year),          
  family = student(), data = d, iter = 4000, cores = 3, chains = 3,
  save_all_pars = TRUE,
  prior = prior0,
  )

summary(m0)
loo_m0 <- loo(m0, moment_match = TRUE)


# Area-specific slopes that also vary by year (uncorrelated random effects)
m1 <- brm(
  log_cpue ~ -1 + age * area2 + (0 + area2*age||year),
  family = student(), data = d, iter = 4000, cores = 3, chains = 3,
  save_all_pars = TRUE,
  prior = prior0,
  control = list(adapt_delta = 0.99)
  )

summary(m1)
loo_m1 <- loo(m1, moment_match = TRUE)

# Check how bad they are. I think ok for now, below 0.7 at least
plot(loo_m1, diagnostic = c("k", "n_eff"), label_points = TRUE)

# Area-specific slopes that also vary by year (correlated random effects)
m1b <- brm(log_cpue ~ -1 + age * area2 + (0 + area2*age||year),
           family = student(), data = d, iter = 4000, cores = 3, chains = 3,
           save_all_pars = TRUE,
           control = list(adapt_delta = 0.99),
           prior = prior0)

summary(m1b)
loo_m1b <- loo(m1b, moment_match = TRUE)

# Check how bad they are. I think ok for now, below 0.7 at least
plot(loo_m1b, diagnostic = c("k", "n_eff"), label_points = TRUE)

loo_compare(loo_m0, loo_m1, loo_m1b)


# D. PRODUCE FIGURES ===============================================================
##### Plot Predictions =============================================================
pal <- rev(brewer.pal(n = 6, name = "Paired")[c(2, 6)])

as.data.frame(fixef(m1)) # Extract "fixed" effects from m2 for plotting the equation 

p1 <- d %>%
  ungroup() %>%
  data_grid(age = seq_range(age, by = 1),
            area2 = c("Warm", "Cold")) %>%
  add_predicted_draws(m1, re_formula = NA) %>%
  ggplot(aes(factor(age), y = log_cpue, color = area2, fill = area2)) +
  stat_lineribbon(aes(y = .prediction), .width = c(.90, .50), alpha = 0.25, size = 0.5) +
  geom_jitter(data = d, alpha = 0.9, size = 1.3, width = 0.2, height = 0, shape = 21, color = "white") +
  stat_lineribbon(aes(y = .prediction), .width = 0, alpha = 0.8, size = 0.5) +
  scale_fill_manual(values = rev(pal)) +
  scale_color_manual(values = rev(pal)) +
  labs(color = "Area", fill = "Area", x = "Age [yrs]", y = "Log(CPUE)") +
  guides(color = guide_legend(override.aes = list(linetype = 0, size = 2, shape = 16, alpha = 0.5,
                                                  color = rev(pal), fill = NA))) +
  annotate("text", 2.4, -0.45, label = paste("n=", nrow(d), sep = ""), size = 2.5) +
  annotate("text", 2.4, -0.9, size = 2.5, color = pal[2],
           label = expression(paste("y=6.55-0.64×age; ", italic(Z), "=0.64 [0.60, 0.69]", sep = ""))) + # Cold
  annotate("text", 2.4, -1.35, size = 2.5, color = pal[1],
           label = expression(paste("y=5.56-0.75×age; ", italic(Z), "=0.75 [0.65, 0.85]", sep = ""))) + # Warm
  NULL

pWord1 <- p1 + theme(text = element_text(size = 12), 
                     legend.position = c(0.9, 0.9),
                     legend.spacing.y = unit(0, 'cm'),
                     legend.key.size = unit(0, "cm"),
                     legend.title = element_text(size = 10),
                     legend.text = element_text(size = 8))

# ggsave("figures/catch_curve.png", width = 6.5, height = 6.5, dpi = 600)


##### Random year effects ==========================================================
# http://mjskay.github.io/tidybayes/articles/tidy-brms.html
get_variables(m1)

# Random intercepts
# Warm
warm_intercept <- m1 %>%
  spread_draws(b_area2Warm, r_year[year, area2Warm]) %>%
  filter(area2Warm == "area2Warm") %>% 
  mutate(warm_intercepts = b_area2Warm + r_year) %>% # The random effects are offsets
  ggplot(aes(y = factor(year), x = warm_intercepts)) +
  stat_halfeye(fill = pal[1], alpha = 0.8) + 
  labs(y = "Year", x = "Intercept") + 
  ggtitle("Warm")

# Cold
cold_intercept <- m1 %>%
  spread_draws(b_area2Cold, r_year[year, area2Cold]) %>%
  filter(area2Cold == "area2Cold") %>% 
  mutate(cold_intercepts = b_area2Cold + r_year) %>% # The random effects are offsets
  ggplot(aes(y = factor(year), x = cold_intercepts)) +
  stat_halfeye(fill = pal[2], alpha = 0.8) + 
  labs(y = "", x = "Intercept") + 
  ggtitle("Cold")

warm_intercept + cold_intercept

ggsave("figures/supp/catch_curves_random_intercepts.png", width = 6.5, height = 6.5, dpi = 600)

# Warm slopes
#https://stackoverflow.com/questions/57379091/how-to-extract-tidy-draws-from-brms-models-tidybayes-for-interaction-terms
# Check it corresponds to the summary
# m1 %>% spread_draws(b_age) %>% summarize(mean_warm = mean(b_age))
# m1 %>%
#   spread_draws(b_age, `b_age:area2Warm`, r_year[year, age]) %>% 
#   filter(age %in% c("age", "area2Warm:age")) %>% 
#   mutate(warm_mean_slopes = (b_age + `b_age:area2Warm`)) %>% 
#   ungroup() %>% 
#   summarise(mean_warm_slope = mean(warm_mean_slopes))
# 
# warm_slope <- m1 %>%
#   spread_draws(b_age, `b_age:area2Warm`, r_year[year, age]) %>% # Extract random effects and global effects
#   filter(age %in% c("age", "area2Warm:age")) %>% # Filter the correct parameters. Because the warm slope is an offset, I need also the cold slope. For some reason, I get all levels of r_year when I include multiple parameters, not when I do only r_year[year, age]
#   mutate(warm_mean_slope = (b_age + `b_age:area2Warm`)) %>% # Calculate the average warm slope (not it's a difference to the b_age, i.e. the cold slope)
#   filter(age == "area2Warm:age") %>% # Filter only draws corresponding to the random effect that is the warm slope (not the cold slope)
#   mutate(warm_slopes = warm_mean_slope + r_year) %>% # The random effects are offsets, so calculate the mean slope, we don't want the difference
#   ggplot(aes(y = factor(year), x = warm_slopes)) +
#   stat_halfeye(fill = pal[1], alpha = 0.8) + 
#   labs(y = "", x = "Slope") + 
#   coord_flip() +
#   coord_cartesian(x = c(-0.4, -1))
# 
# # Cold slopes
# cold_slope <- m1 %>%
#   spread_draws(b_age, r_year[year, age]) %>% # Extract random effects and global effects
#   filter(age %in% c("age")) %>% # Filter the correct parameters. For some reason, I get all levels of r_year when I include multiple parameters, not when I do only r_year[year, age]
#   mutate(cold_mean_slope = b_age) %>%
#   mutate(cold_slopes = cold_mean_slope + r_year) %>% # The random effects are offsets, so calculate the mean slope, we don't want the difference
#   ggplot(aes(y = factor(year), x = cold_slopes)) +
#   stat_halfeye(fill = pal[2], alpha = 0.8) + 
#   labs(y = "", x = "Slope") + 
#   coord_flip() +
#   coord_cartesian(x = c(-0.4, -1))
# 
# cold_slope / warm_slope + plot_layout(widths = 0.5)

# Plotting together...
cold_slope_df <- m1 %>%
  spread_draws(b_age, r_year[year, age]) %>% # Extract random effects and global effects
  filter(age %in% c("age")) %>% # Filter the correct parameters. For some reason, I get all levels of r_year when I include multiple parameters, not when I do only r_year[year, age]
  mutate(mean_slope = b_age,
         Area = "Cold") %>%
  mutate(slopes = mean_slope + r_year) 

warm_slope_df <- m1 %>%
  spread_draws(b_age, `b_age:area2Warm`, r_year[year, age]) %>% # Extract random effects and global effects
  filter(age %in% c("age", "area2Warm:age")) %>% # Filter the correct parameters. Because the warm slope is an offset, I need also the cold slope. For some reason, I get all levels of r_year when I include multiple parameters, not when I do only r_year[year, age]
  mutate(mean_slope = (b_age + `b_age:area2Warm`),
         Area = "Warm") %>% # Calculate the average warm slope (not it's a difference to the b_age, i.e. the cold slope)
  filter(age == "area2Warm:age") %>% # Filter only draws corresponding to the random effect that is the warm slope (not the cold slope)
  mutate(slopes = mean_slope + r_year)
  
full_df <- bind_rows(cold_slope_df, warm_slope_df)

p_random <- full_df %>% 
  mutate(Z = slopes*-1) %>% # Convert from slopes to Z
  ggplot(., aes(y = factor(year), x = Z, fill = factor(Area), color = factor(Area))) +
  geom_vline(xintercept = c(0.64, 0.75), color = rev(pal), linetype = 2, alpha = 0.4, size = 0.5) +
  stat_slab(alpha = 0.13, position = position_nudge(y = 0.05), color = NA) + 
  stat_pointinterval(.width = c(.95), position = position_dodge(width = 0.25),
                     size = 0.1, alpha = 0.8) +
  scale_fill_manual(values = rev(pal)) +
  scale_color_manual(values = rev(pal)) +
  labs(y = "Year", x = expression(italic(Z))) + 
  guides(fill = F, color = F) +
  theme(text = element_text(size = 12),
        axis.text.x = element_text(angle = 30, size = 9)) + 
  coord_flip(xlim = c(0.48, 0.9))
  
pWord1 / p_random + plot_annotation(tag_levels = 'A')

ggsave("figures/catch_curve_random.png", width = 5.5, height = 6.5, dpi = 600)


# Testing stat_pointinterval plots median 
# data(RankCorr_u_tau, package = "ggdist")
# t <- RankCorr_u_tau %>% group_by(i) %>% mutate(median = median(u_tau),
#                                                mean = mean(u_tau))
# t %>%
#   ggplot(aes(y = factor(i), x = u_tau)) +
#   stat_pointinterval(.width = c(.66, .95)) +
#   geom_point(aes(y = factor(i), x = mean), color = "red") +
#   geom_point(aes(y = factor(i), x = median), color = "blue")



##### Model diagnostics & fit ======================================================
pal_diag <- rev(brewer.pal(n = 3, name = "Dark2"))

# Chain convergence
posterior <- as.array(m1)
dimnames(posterior)

d1 <- mcmc_trace(posterior,
                 pars = c("b_age", "b_area2Cold", "b_area2Warm", 
                          "b_age:area2Warm", "sd_year__area2Cold", "sd_year__area2Warm",
                          "sd_year__age", "sd_year__area2Warm:age", "nu", "sigma"),
                 facet_args = list(ncol = 3, strip.position = "left")) + 
  theme(text = element_text(size = 12),
        strip.text = element_text(size = 6),
        legend.position = "top") + 
  scale_color_manual(values = alpha(pal_diag, alpha = 0.8))

# Resid vs fitted
d2 <- d %>%
  add_residual_draws(m1) %>%
  ggplot(aes(x = .row, y = .residual)) +
  stat_pointinterval(alpha = 0.5, size = 0.7) + 
  theme(text = element_text(size = 12))

# qq-plot
# d3 <- d %>%
#   add_residual_draws(m1) %>%
#   median_qi() %>%
#   ggplot(aes(sample = .residual)) +
#   geom_qq_line() +
#   geom_qq(alpha = 0.8) +
#   theme(text = element_text(size = 12))

summary(m1)$spec_pars
nu <- summary(m1)$spec_pars[2, 1]
nu

d3 <- d %>%
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

ggsave("figures/supp/catch_curve_diag_fit.png", width = 6.5, height = 8.5, dpi = 600)
