#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 2020.06.16: Max Lindmark
#
# Fit models to explore differences in size spectra slopes
# 
# A. Load libraries
# 
# B. Read data
# 
# C. Fit models of spectra by year following Edwards (MEPS 2020)
#
# D. Fit models to size spectra slopes
#
# E. Produce figures
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# A. LOAD LIBRARIES ================================================================
rm(list = ls())

# Load libraries (install first if needed)
library(tidyverse); theme_set(theme_classic(base_size = 12))
library(tidylog)
library(RColorBrewer)
library(patchwork)
library(sizeSpectra)
library(brms)
library(nlstools)
library(bayesplot)
library(tidylog)
library(tidybayes)
library(RColorBrewer)
library(modelr)

# other attached packages:
# [1] sizeSpectra_1.0.0.0 patchwork_1.0.1     RColorBrewer_1.1-2  tidylog_1.0.2       forcats_0.5.0       stringr_1.4.0      
# [7] dplyr_1.0.2         purrr_0.3.4         readr_1.3.1         tidyr_1.1.0         tibble_3.0.3        ggplot2_3.3.2      
# [13] tidyverse_1.3.0  

# For parallel processing
options(mc.cores = parallel::detectCores()) 


# B. READ DATA =====================================================================
# df <- read.csv("data/Growth_data_BT_FM_1970-2004.csv", sep = ";") # This is the 
# original data from Huss et al (2019)
df <- read.csv("data/catch_BT_FM_1987-2003.csv")

# How many nets in total? (For scaling with effort later)
df %>%
  group_by(netID) %>%
  mutate(n = n()) %>%
  ggplot(., aes(netID, n)) +
  geom_histogram(stat ="identity")

df <- df %>% group_by(Area, year) %>% mutate(n_nets_year = length(unique(netID))) %>% ungroup()

# Test
df %>% filter(year == 1995 & Area == "BT") %>% distinct(netID)

# Now we need to find representative catch sizes
# Plot n per length group (cm)
p1 <- ggplot(df, aes(factor(length_group))) +   
  geom_histogram(stat = "count") +
  facet_wrap(~Area, scales = "free") + # Note the differences in # data
  theme_classic() + 
  coord_cartesian(expand = 0)

# Filter out fish > 13 cm
df <- df %>% filter(length_group > 13) %>% as.data.frame()

# Plot again
p2 <- ggplot(df, aes(factor(length_group))) +   
  geom_histogram(stat = "count") +
  facet_wrap(~Area, scales = "free") + # Note the differences in # data
  theme_classic() + 
  coord_cartesian(expand = 0)

p1/p2

# Before we process data any further, save it because we will use it in this format
# for analysing mean size (i.e. we need 1 row = 1 ind)

size_df <- df

# Now we want to process the data a bit further... Following Edwards sizeSpectra package
# we want the data as follows:
# Year 	SpecCode 	LngtClass 	Number 	LWa 	LWb 	bodyMass 	Biomass
# Hence, for each year, we need to calculate the CPUE
# https://htmlpreview.github.io/?https://raw.githubusercontent.com/andrew-edwards/sizeSpectra/master/doc/MEPS_IBTS_MLEbins.html

# Group by year, Area and length group, summarize and get n()
df2 <- df %>%
  group_by(year, Area, length_group) %>% 
  summarise(catch_n = n()) %>% # Get catch as a column, this is now the number of rows a given length occurs each year and area
  ungroup() %>% 
  mutate(effort_id = paste(year, Area, sep = ".")) %>% 
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
df %>% filter(year == 1987 & Area == "FM") %>% distinct(n_nets_year, .keep_all = TRUE)
df2 %>% filter(effort_id == "1987.FM")
df3 %>% filter(effort_id == "1987.FM")

# Now do some additional calculations (weights etc, to match data structure in sizeSpectra package)
# LW parameters from FishBase (2020.09.24)
a <- 0.01 
b <- 3.08

df4 <- df3 %>% 
  ungroup() %>% 
  rename("area" = "Area",
         "min_length_group_cm" = "length_group") %>% # What we have now is the minimum length in each class
  mutate(max_length_group_cm = min_length_group_cm + 2.4, # Get max length in bin (2.5 cm length-classes)
         wmin = a*min_length_group_cm^b,     # Get min mass in bin
         wmax = a*max_length_group_cm^b) %>% # Get max mass in bin
  mutate(cpue_numbers = catch_n/n_nets_year, # Get numbers CPUE, divide by the previously create n_nets, which is # of unique net ID's in each area and year
         cpue_biom = (catch_n*((wmin + wmax)/2))/n_nets_year) %>% # Get biomass CPUE, use mean of mass in size range
  mutate(SpecCode = "Perch")

# Test I get 1 unique row per size class, year and area
df4 %>%
  group_by(year, area, min_length_group_cm) %>%
  summarise(n = n()) %>% 
  filter(n == 1)

# Plot the log biomass cpue as function of log weight
ggplot(df4, aes(log(wmax), log(cpue_biom), color = factor(year))) + 
  stat_smooth(se = FALSE) +
  geom_point() +
  facet_wrap(~ area) +
  scale_color_viridis(discrete = TRUE) + 
  coord_cartesian(expand = 0) + 
  theme(aspect.ratio = 1)

# Plot the total number of catches by length-group in this processed data
# Compare with p2, should be the same plot...
ggplot(df4, aes(min_length_group_cm, catch_n)) + 
  geom_bar(stat = "identity") + 
  facet_wrap(~ area, scales = "free") +
  theme_classic() + 
  coord_cartesian(expand = 0)

p2

# Now that we have CPUE and not COUNT, we can see if the filter to maximum catch size is
# appropriate
p3 <- ggplot(df4, aes(factor(round(wmin, digits = 3)), cpue_biom)) + 
  geom_bar(stat = "identity") + 
  facet_wrap(~ area, scales = "free") +
  theme_classic() + 
  coord_cartesian(expand = 0)

# Filter out fish > 13 cm
df4 <- df4 %>% filter(wmin > 50) %>% as.data.frame()

# Plot again
p4 <- ggplot(df4, aes(factor(round(wmin, digits = 3)), cpue_biom)) + 
  geom_bar(stat = "identity") + 
  facet_wrap(~ area, scales = "free") +
  theme_classic() + 
  coord_cartesian(expand = 0)

p3/p4


# C. FIT MODELS ====================================================================
# Test with my data
dataBintest <- df4 %>% rename("Year" = "year",
                              "Number" = "cpue_numbers")


# Forsmark =========================================================================
# Following and modifying this vignette:
# https://htmlpreview.github.io/?https://raw.githubusercontent.com/andrew-edwards/sizeSpectra/master/doc/MEPS_IBTS_MLEbins.html

dataBin_FM <- dataBintest %>%
  filter(area == "FM") %>%
  select(SpecCode, wmin, wmax, Number, Year)

#** Loop through all years and model ===============================================
# Get a vector of all years to loop over
fullYears = sort(unique(dataBin_FM$Year))

for(iii in 1:length(fullYears))
  {
    dataBin_FM_loop = dplyr::filter(dataBin_FM, Year == fullYears[iii])
    
    n_FM = sum(dataBin_FM_loop$Number)
    xmin_FM = min(dataBin_FM_loop$wmin)
    xmax_FM = max(dataBin_FM_loop$wmax)

    FM_spectra_one_year  = calcLike(negLL.fn = negLL.PLB.bins.species,
                                  p = -1.5,
                                  suppress.warnings = FALSE,
                                  vecDiff = 2, # Default is 0.5, I get warning thogh
                                  dataBinForLike = dataBin_FM_loop,
                                  n = n_FM,
                                  xmin = xmin_FM,
                                  xmax = xmax_FM)
    
    if(iii == 1)
    {
      FM_spectra = data.frame(Year = fullYears[iii],
                            xmin = xmin_FM,
                            xmax = xmax_FM,
                            n = n_FM,
                            b = FM_spectra_one_year$MLE,
                            confMin = FM_spectra_one_year$conf[1],
                            confMax = FM_spectra_one_year$conf[2])
    } else {
      FM_spectra = rbind(FM_spectra,
                       c(fullYears[iii],
                         xmin_FM,
                         xmax_FM,
                         n_FM,
                         FM_spectra_one_year$MLE,
                         FM_spectra_one_year$conf[1],
                         FM_spectra_one_year$conf[2]))
    }
}

FM_spectra


#** Loop and extract data for plotting =============================================
# https://htmlpreview.github.io/?https://raw.githubusercontent.com/andrew-edwards/sizeSpectra/master/doc/MEPS_IBTS_recommend.html
# These calculations are to get the required input for the recommended plot
# (see ?sizeSpectra::ISD_bin_plot for the structure). This could maybe be functionalised
# like the plotting function ISD_bin_plot().

dataRecommend.isd_FM = dplyr::select(dataBin_FM, Year, wmin, wmax, Number)

data.year.list_FM = list() # To save results for each year

fullYears = sort(unique(dataBin_FM$Year))

for(i in 1:length(fullYears))
{
  data.year = dplyr::filter(dataRecommend.isd_FM, Year == fullYears[i])
  data.year = dplyr::arrange(data.year, desc(wmin))
  sumNumber = sum(data.year$Number)
  
  # data.year = dplyr::mutate(data.year,
  #                          cumSum = cumsum(Number))
  # This is wrong when we have two species with the same
  #  length-weight coefficients in the same year, so need countGTEwmin below
  # Have to do not with dplyr:
  wmin.vec = data.year$wmin
  wmax.vec = data.year$wmax
  num.vec  = data.year$Number
  
  countGTEwmin = rep(NA, length(num.vec)) # to do a manual count
  lowCount = countGTEwmin
  highCount = countGTEwmin
  
  for(iii in 1:length(countGTEwmin))
  {
    countGTEwmin[iii] = sum( (wmin.vec >= wmin.vec[iii]) * num.vec)
    lowCount[iii]  = sum( (wmin.vec >= wmax.vec[iii]) * num.vec)
    highCount[iii] = sum( (wmax.vec >  wmin.vec[iii]) * num.vec)
  }
  data.year = cbind(data.year,
                    "countGTEwmin" = countGTEwmin,
                    "lowCount" = lowCount,
                    "highCount" = highCount)
  data.year = dplyr::tbl_df(data.year) # This is one of the desired input for the plotting function below
  
  data.year.list_FM[[i]] = data.year
}

xlim.global = c(min(dataRecommend.isd_FM$wmin),
                max(dataRecommend.isd_FM$wmax))   # x-axis limits to be common for all plots


#** Loop plots and save by year ====================================================
# https://htmlpreview.github.io/?https://raw.githubusercontent.com/andrew-edwards/sizeSpectra/master/doc/MEPS_IBTS_MLEbins.html
# In the vignette for plotting they don't define MLEbins.res, so I checked the other
# vignettes and see that it's the dataframe created in the fitting loop

for(i in 1:length(fullYears))
{
  postscript(paste0("figures/supp/annual_size_spectra/FM/", fullYears[i], ".eps"),
             height = 8, width = 5.36,
             horizontal=FALSE, paper="special")
  
  ISD_bin_plot(data.year = data.year.list_FM[[i]],
               b.MLE = dplyr::filter(FM_spectra, Year == fullYears[i])$b,
               b.confMin = dplyr::filter(FM_spectra, Year == fullYears[i])$confMin,
               b.confMax = dplyr::filter(FM_spectra, Year == fullYears[i])$confMax,
               year = fullYears[i],
               xlim = xlim.global,
               xmin = dplyr::filter(FM_spectra, Year == fullYears[i])$xmin,
               xmax = dplyr::filter(FM_spectra, Year == fullYears[i])$xmax
  )
  dev.off()
}


# Biotest ==========================================================================
# Following and modifying this vignette:
# https://htmlpreview.github.io/?https://raw.githubusercontent.com/andrew-edwards/sizeSpectra/master/doc/MEPS_IBTS_MLEbins.html

dataBin_BT <- dataBintest %>%
  filter(area == "BT") %>%
  select(SpecCode, wmin, wmax, Number, Year)


#** Loop through all years and model ===============================================
# Get a vector of all years to loop over
fullYears = sort(unique(dataBin_BT$Year))

for(iii in 1:length(fullYears))
{
  dataBin_BT_loop = dplyr::filter(dataBin_BT, Year == fullYears[iii])
  
  n_BT = sum(dataBin_BT_loop$Number)
  xmin_BT = min(dataBin_BT_loop$wmin)
  xmax_BT = max(dataBin_BT_loop$wmax)
  
  BT_spectra_one_year  = calcLike(negLL.fn = negLL.PLB.bins.species,
                                  p = -1.5,
                                  suppress.warnings = FALSE,
                                  vecDiff = 2, # Default is 0.5, I get warning though
                                  dataBinForLike = dataBin_BT_loop,
                                  n = n_BT,
                                  xmin = xmin_BT,
                                  xmax = xmax_BT)
  
  if(iii == 1)
  {
    BT_spectra = data.frame(Year = fullYears[iii],
                            xmin = xmin_BT,
                            xmax = xmax_BT,
                            n = n_BT,
                            b = BT_spectra_one_year$MLE,
                            confMin = BT_spectra_one_year$conf[1],
                            confMax = BT_spectra_one_year$conf[2])
  } else {
    BT_spectra = rbind(BT_spectra,
                       c(fullYears[iii],
                         xmin_BT,
                         xmax_BT,
                         n_BT,
                         BT_spectra_one_year$MLE,
                         BT_spectra_one_year$conf[1],
                         BT_spectra_one_year$conf[2]))
  }
}

BT_spectra


#** Loop and extract data for plotting =============================================
# https://htmlpreview.github.io/?https://raw.githubusercontent.com/andrew-edwards/sizeSpectra/master/doc/MEPS_IBTS_recommend.html
# These calculations are to get the required input for the recommended plot
# (see ?sizeSpectra::ISD_bin_plot for the structure). This could maybe be functionalised
# like the plotting function ISD_bin_plot().

dataRecommend.isd_BT = dplyr::select(dataBin_BT, Year, wmin, wmax, Number)

data.year.list_BT = list() # To save results for each year

fullYears = sort(unique(dataBin_BT$Year))

for(i in 1:length(fullYears))
{
  data.year = dplyr::filter(dataRecommend.isd_BT, Year == fullYears[i])
  data.year = dplyr::arrange(data.year, desc(wmin))
  sumNumber = sum(data.year$Number)
  
  # data.year = dplyr::mutate(data.year,
  #                          cumSum = cumsum(Number))
  # This is wrong when we have two species with the same
  #  length-weight coefficients in the same year, so need countGTEwmin below
  # Have to do not with dplyr:
  wmin.vec = data.year$wmin
  wmax.vec = data.year$wmax
  num.vec  = data.year$Number
  
  countGTEwmin = rep(NA, length(num.vec)) # to do a manual count
  lowCount = countGTEwmin
  highCount = countGTEwmin
  
  for(iii in 1:length(countGTEwmin))
  {
    countGTEwmin[iii] = sum( (wmin.vec >= wmin.vec[iii]) * num.vec)
    lowCount[iii]  = sum( (wmin.vec >= wmax.vec[iii]) * num.vec)
    highCount[iii] = sum( (wmax.vec >  wmin.vec[iii]) * num.vec)
  }
  data.year = cbind(data.year,
                    "countGTEwmin" = countGTEwmin,
                    "lowCount" = lowCount,
                    "highCount" = highCount)
  data.year = dplyr::tbl_df(data.year) # This is one of the desired input for the plotting function below
  
  data.year.list_BT[[i]] = data.year
}

xlim.global = c(min(dataRecommend.isd_BT$wmin),
                max(dataRecommend.isd_BT$wmax)) # x-axis limits to be common for all plots


#** Loop plots and save by year ====================================================
# https://htmlpreview.github.io/?https://raw.githubusercontent.com/andrew-edwards/sizeSpectra/master/doc/MEPS_IBTS_MLEbins.html
# In the vignette for plotting they don't define MLEbins.res, so I checked the other
# vignettes and see that it's the dataframe created in the fitting loop

for(i in 1:length(fullYears))
{
  postscript(paste0("figures/supp/annual_size_spectra/BT/", fullYears[i], ".eps"),
             height = 8, width = 5.36,
             horizontal=FALSE, paper="special")
  
  ISD_bin_plot(data.year = data.year.list_BT[[i]],
               b.MLE = dplyr::filter(BT_spectra, Year == fullYears[i])$b,
               b.confMin = dplyr::filter(BT_spectra, Year == fullYears[i])$confMin,
               b.confMax = dplyr::filter(BT_spectra, Year == fullYears[i])$confMax,
               year = fullYears[i],
               xlim = xlim.global,
               xmin = dplyr::filter(BT_spectra, Year == fullYears[i])$xmin,
               xmax = dplyr::filter(BT_spectra, Year == fullYears[i])$xmax
  )
  dev.off()
}


# Combine areas and plot slopes across years =======================================================
# https://htmlpreview.github.io/?https://raw.githubusercontent.com/andrew-edwards/sizeSpectra/master/doc/MEPS_IBTS_MLEbins.html

FM_spectra = dplyr::tbl_df(FM_spectra)
FM_spectra = dplyr::mutate(FM_spectra, stdErr = (abs(confMin-b) + abs(confMax-b))/(2*1.96) )
FM_spectra = dplyr::mutate(FM_spectra, area = "FM")

BT_spectra = dplyr::tbl_df(BT_spectra)
BT_spectra = dplyr::mutate(BT_spectra, stdErr = (abs(confMin-b) + abs(confMax-b))/(2*1.96) )
BT_spectra = dplyr::mutate(BT_spectra, area = "BT")

spectra <- rbind(BT_spectra, FM_spectra)

# This is the default vignette plot
res = timeSerPlot(FM_spectra,
                  legName = "(a) MLEbins",
#                  yLim = c(-2.2, -0.9),
                  xLab = "Year",
                  method = "MLEbins",
                  legPos = "bottomleft",
                  weightReg = TRUE,
                  xTicksSmallInc = 1,
                  yTicksSmallInc = 0.05)

res = timeSerPlot(BT_spectra,
                  legName = "(a) MLEbins",
                  #                  yLim = c(-2.2, -0.9),
                  xLab = "Year",
                  method = "MLEbins",
                  legPos = "bottomleft",
                  weightReg = TRUE,
                  xTicksSmallInc = 1,
                  yTicksSmallInc = 0.05)


# D. FIT MODELS ====================================================================
# Mean center year variable
spectra <- spectra %>%
  mutate(Year_ct = Year - min(Year),
         Year_ct = as.integer(Year_ct),
         area2 = ifelse(area == "BT", "Warm", "Cold"))

# Difference in mean size spectrum slopes between area, year as random factor
m1 <- brm(b ~ area2 + (1|Year_ct),
          family = gaussian(), data = spectra, iter = 4000, cores = 3, chains = 3,
          save_all_pars = TRUE,
          control = list(adapt_delta = 0.99))

summary(m1)
plot(m1)


# E. PRODUCE FIGURES ===============================================================
pal <- brewer.pal(n = 6, name = "Paired")[c(2, 6)]
pal2 <- brewer.pal(n = 6, name = "Paired")[c(2, 6)]
pal2 <- alpha(pal2, alpha = 0.2)

##### Plot predictions ============================================================
p1 <- m1 %>%
  spread_draws(b_Intercept, b_area2Warm) %>%
  mutate("Warm" = b_Intercept + b_area2Warm) %>%
  rename("Cold" = "b_Intercept") %>% 
  pivot_longer(c("Warm", "Cold"), names_to = "area", values_to = c("b")) %>% 
  ggplot(aes(x = area, y = b, fill = area)) +
  coord_flip() +
  stat_halfeye(position = position_nudge(x = .1, y = 0), alpha = 0.8) +
  scale_fill_manual(values = pal2) +
  geom_jitter(data = spectra, aes(x = area2, y = b, color = Year), width = 0.05,  
              alpha = 0.8, size = 3, inherit.aes = FALSE) +
  scale_color_viridis() +
  labs(x = "Area",
       y = expression(paste("Size-spectrum slope ", italic((b)))),
       color = "Year") +
  guides(fill = FALSE) +
  NULL
  
pWord1 <- p1 + theme(text = element_text(size = 12), # 12 for word doc
                    legend.title = element_text(size = 10),
                    legend.text = element_text(size = 8), 
                    legend.position = c(1, 0.001), 
                    legend.justification = c(1, 0),
                    legend.direction = "horizontal")

ggsave("figures/size_spec.png", width = 6.5, height = 6.5, dpi = 600)


##### Model diagnostics ============================================================
pal_diag <- rev(brewer.pal(n = 3, name = "Dark2"))

# Chain convergence
posterior <- as.array(m1)
dimnames(posterior)

d1 <- mcmc_trace(posterior,
                 pars = c("b_Intercept", "b_area2Warm", "sd_Year_ct__Intercept", 
                          "sigma", "Intercept"),
                 facet_args = list(ncol = 3, strip.position = "left")) + 
  theme(text = element_text(size = 12),
        strip.text = element_text(size = 6),
        legend.position = "top") + 
  scale_color_manual(values = alpha(pal_diag, alpha = 0.8))

# Resid vs fitted
d2 <- spectra %>%
  add_residual_draws(m1) %>%
  ggplot(aes(x = .row, y = .residual)) +
  stat_pointinterval(alpha = 0.5, size = 0.7) + 
  theme(text = element_text(size = 12))

# qq-plot
d3 <- spectra %>%
  add_residual_draws(m1) %>%
  median_qi() %>%
  ggplot(aes(sample = .residual)) +
  geom_qq_line() +
  geom_qq(alpha = 0.8) +
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

ggsave("figures/supp/size_spec_diag.png", width = 6.5, height = 8.5, dpi = 600)
