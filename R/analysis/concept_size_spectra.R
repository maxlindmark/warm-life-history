# Install new mizer to a specific location

#install.packages("mizer")

# install.packages("mizer", lib = "/Users/maxlindmark/Desktop/R_STUDIO_PROJECTS/warm_life_history")

library(mizer, lib.loc = "/Users/maxlindmark/Desktop/R_STUDIO_PROJECTS/warm_life_history")  

rm(list = ls())

# Using this as a guide
# https://sizespectrum.org/mizer/articles/single_species_size-spectrum_dynamics.html

sessionInfo()

library(mizer)
library(tidyverse)
library(ggplot2)
library(patchwork)
library(ggpmisc)
library(RColorBrewer)

params_high_growth <- newSingleSpeciesParams(lambda = 2.05)
params_high_growth@species_params$w_mat <- 20
params_high_growth <- steady(params_high_growth)

plotSpectra(params_high_growth, power = 0, resource = FALSE)

######## DEFINE MODELS ########
# Lower growth (this will be my default actually)
params_def <- params_high_growth # changing growth via gamma
species_params(params_def)$gamma <- species_params(params_high_growth)$gamma*0.8
params_def <- steady(params_def)

# Lower L_inf but high growth (TSR)
params_TSR <- params_high_growth # changing growth via gamma (from the high growth model) and w_inf
species_params(params_TSR)$w_inf <- species_params(params_high_growth)$w_inf*0.5
params_TSR <- steady(params_TSR)

# Lower L_inf and same growth!
params_TSR_no_growth <- params_def # changing growth via gamma (from the high growth model) and w_inf
species_params(params_TSR_no_growth)$w_inf <- species_params(params_def)$w_inf*0.5
params_TSR_no_growth <- steady(params_TSR_no_growth)

####### Testing growth models
g1 <- getGrowthCurves(params_high_growth)
plot_dat1 <- reshape2::melt(g1)
plot_dat1$type <- "high growth"
plot_dat1

g2 <- getGrowthCurves(params_def)
plot_dat2 <- reshape2::melt(g2)
plot_dat2$type <- "def"
plot_dat2

g3 <- getGrowthCurves(params_TSR)
plot_dat3 <- reshape2::melt(g3)
plot_dat3$type <- "TSR"
plot_dat3

g4 <- getGrowthCurves(params_TSR_no_growth)
plot_dat4 <- reshape2::melt(g4)
plot_dat4$type <- "TSR no growth"
plot_dat4

dat <- rbind(plot_dat1, plot_dat2, plot_dat3, plot_dat4)
ggplot(dat, aes(Age, value, color = type, linetype = type)) + geom_line(size = 1.5) + 
  theme_light() + coord_cartesian(xlim = c(0, 12))
#######

# Increase mortality (currently only for large fish)
params_fishing <- params_def
species_params(params_fishing)$sel_func <- "knife_edge"
species_params(params_fishing)$knife_edge_size <- 30
# We also need to specify the fishing effort and can then plot the resulting fishing mortality.
params_fishing <- setFishing(params_fishing, initial_effort = 1)
plotFMort(params_fishing) + geom_vline(xintercept = 30, linetype = "dashed")
params_fishing <- steady(params_fishing)
params_fishing@selectivity
plot(params_fishing)

######## Compare their spectras ########
# High growth
nf_high_growth <- melt(initialN(params_high_growth))
nf_high_growth$type <- "high growth"
nf_high_growth$value2 <- nf_high_growth$value * w(params_high_growth)

# Def
nf_def <- melt(initialN(params_def))
nf_def$type <- "default"
nf_def$value2 <- nf_def$value * w(params_def)

# TSR
nf_TSR <- melt(initialN(params_TSR))
nf_TSR$type <- "TSR"
nf_TSR$value2 <- nf_TSR$value * w(params_TSR)

# TSR no growth
nf_TSR_ng <- melt(initialN(params_TSR_no_growth))
nf_TSR_ng$type <- "TSR no growth"
nf_TSR_ng$value2 <- nf_TSR_ng$value * w(params_TSR_no_growth)

# Fishing mort
nf_fishing <- melt(initialN(params_fishing))
nf_fishing$type <- "fishing mort"
nf_fishing$value2 <- nf_fishing$value * w(params_fishing)

nf_combined <- rbind(nf_high_growth, nf_def, nf_TSR, nf_TSR_ng, nf_fishing)

# Plot and check slopes and intercepts
# give a name to a formula
formula <- value2 ~ w

# no weights
ggplot(my.data, aes(x, y)) +
  geom_point() +
  geom_smooth(method = "lm", formula = formula) +
  stat_poly_eq(formula = formula, parse = TRUE)

# Test plot to see lines are good fits
#ggplot(filter(nf_combined, w < 99 & value2 > 1e-09),
ggplot(filter(nf_combined, value2 > 1e-10),
#ggplot(filter(nf_combined, value2 > 1e-13),
#ggplot(nf_combined,
       aes(x = w, y = value2, colour = type)) + #y = value, # in the vignette they use value, I want biomass spectra so I use value 2
  scale_x_log10() +
  scale_y_log10() +
  geom_point(size = 1.5) +
  facet_wrap(~type) +
  coord_cartesian(xlim = c(1e-03, 100), ylim = c(1e-10, 1e-01)) +
  stat_ma_line(se = FALSE) +
  stat_ma_eq(aes(label = stat(eq.label)), label.x = "left", label.y = "bottom",) +
  scale_color_brewer(palette = "Dark2") +
  theme_light() + 
  NULL

# Real plot
ggplot(filter(nf_combined, value2 > 1e-10),
       aes(x = w, y = value2, colour = type, linetype = type)) + #y = value, # in the vignette they use value, I want biomass spectra so I use value 2
  scale_x_log10() +
  scale_y_log10() +
  coord_cartesian(xlim = c(1e-03, 100), ylim = c(1e-07, 1e-01)) +
  stat_ma_line(se = FALSE, size = 1.5) +
  stat_ma_eq(aes(label = stat(eq.label)), label.x = "left", label.y = "bottom", size = 6) +
  scale_color_brewer(palette = "Dark2", name = "Scenario") +
  theme_light(base_size = 18) + 
  guides(linetype = FALSE) +
  theme(legend.position = c(0.8, 0.8)) +
  NULL

# raw spectra no equation, focus on TSR and fishing mort
ggplot(filter(nf_combined, type %in% c("default", "fishing mort", "TSR")),
       aes(x = w, y = value2, colour = type, linetype = type)) + #y = value, # in the vignette they use value, I want biomass spectra so I use value 2
  scale_x_log10() +
  scale_y_log10() +
  geom_point(size = 1.5, alpha = 0.2) +
  stat_smooth(method = "lm", se = FALSE) +
  coord_cartesian(xlim = c(1e-03, 100), ylim = c(1e-06, 1e-01)) +
  scale_color_brewer(palette = "Dark2", name = "Scenario") +
  theme_light(base_size = 18) + 
  guides(linetype = FALSE) +
  theme(legend.position = c(0.8, 0.8)) +
  NULL


# Maybe, I should make facet on these plots with different mortality effects, i.e. No effects, all sizes mortality and size-dep mortality
# Because all scenarios really need mort