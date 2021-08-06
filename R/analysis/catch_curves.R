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
# Fitting models of log catch ~ age with interactive or additive effects of area, 
# using catch year as a random effect

# Center year variable and rename area
d <- d %>%
  mutate(year_ct = year - min(year),
         year_ct = as.integer(year_ct),
         log_cpue = log(cpue_numbers),
         area2 = ifelse(area == "BT", "Warm", "Cold"))

# Without interaction
m1 <- brm(log_cpue ~ age + area2 + (1 + age|year_ct),
          family = gaussian(), data = d, iter = 4000, cores = 3, chains = 3,
          save_all_pars = TRUE,
          control = list(adapt_delta = 0.99))

# With interaction
m2 <- brm(log_cpue ~ age * area2 + (1 + age|year_ct),
          family = gaussian(), data = d, iter = 4000, cores = 3, chains = 3,
          save_all_pars = TRUE,
          control = list(adapt_delta = 0.99))

loo_m1 <- loo(m1)
loo_m2 <- loo(m2)

loo_compare(loo_m1, loo_m2)

summary(m2)
plot(m2)


# D. PRODUCE FIGURES ===============================================================
##### Plot Predictions =============================================================

pal <- rev(brewer.pal(n = 6, name = "Paired")[c(2, 6)])

#as.data.frame(fixef(m1))
as.data.frame(fixef(m2)) # Extract "fixed" effects from m2 for plotting the equation 

p1 <- d %>%
  ungroup() %>%
  data_grid(age = seq_range(age, by = 1),
            area2 = c("Warm", "Cold")) %>%
  add_predicted_draws(m2, re_formula = NA) %>%
  ggplot(aes(factor(age), y = log_cpue, color = area2, fill = area2)) +
  stat_lineribbon(aes(y = .prediction), .width = c(.90, .50), alpha = 0.25) +
  geom_jitter(data = d, alpha = 0.9, size = 3, width = 0.2, height = 0, shape = 21, color = "white") +
  stat_lineribbon(aes(y = .prediction), .width = 0, alpha = 0.8) +
  scale_fill_manual(values = rev(pal)) +
  scale_color_manual(values = rev(pal)) +
  labs(color = "Area", fill = "Area",
       x = "Age",
       y = "Log(CPUE)") +
  guides(fill = guide_legend(override.aes = list(alpha = 0.5))) + 
  # Using m2
  annotate("text", 2.2, -1, size = 3.5, color = pal[2],
           label = expression(paste("y=6.56-0.64×age; ", italic(Z), "=0.65 [0.59, 0.70]", sep = ""))) + # Cold
  annotate("text", 2.2, -1.25, size = 3.5, color = pal[1],
           label = expression(paste("y=5.60-0.75×age; ", italic(Z), "=0.75 [0.63, 0.86]", sep = ""))) + # Warm
  NULL

pWord1 <- p1 + theme(text = element_text(size = 12), 
                     legend.position = c(0.9, 0.9), 
                     legend.title = element_text(size = 10),
                     legend.text = element_text(size = 10))

ggsave("figures/catch_curve.png", width = 6.5, height = 6.5, dpi = 600)


##### Model diagnostics ============================================================
pal_diag <- rev(brewer.pal(n = 3, name = "Dark2"))

# Chain convergence
posterior <- as.array(m2)
dimnames(posterior)

d1 <- mcmc_trace(posterior,
                 pars = c("b_Intercept", "b_age", "b_area2Warm", 
                          "b_age:area2Warm", "sd_year_ct__Intercept", "sd_year_ct__age",
                          "cor_year_ct__Intercept__age", "sigma", "Intercept"),
                 facet_args = list(ncol = 3, strip.position = "left")) + 
  theme(text = element_text(size = 12),
        strip.text = element_text(size = 6),
        legend.position = "top") + 
  scale_color_manual(values = alpha(pal_diag, alpha = 0.8))

# Resid vs fitted
d2 <- d %>%
  add_residual_draws(m2) %>%
  ggplot(aes(x = .row, y = .residual)) +
  stat_pointinterval(alpha = 0.5, size = 0.7) + 
  theme(text = element_text(size = 12))

# qq-plot
d3 <- d %>%
  add_residual_draws(m2) %>%
  median_qi() %>%
  ggplot(aes(sample = .residual)) +
  geom_qq_line() +
  geom_qq(alpha = 0.8) +
  theme(text = element_text(size = 12))

# Posterior predictive
d4 <- pp_check(m2) + 
  theme(text = element_text(size = 12),
        legend.position = c(0.15, 0.95),
        legend.background = element_rect(fill = NA)) + 
  scale_color_manual(values = rev(pal_diag)) +
  labs(color = "")

d1 / (d2 / (d3 + d4)) + 
  plot_annotation(tag_levels = 'A')

ggsave("figures/supp/catch_curve_diag.png", width = 6.5, height = 8.5, dpi = 600)
