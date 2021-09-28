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
library(viridis)
library(ggridges)

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


# C. FIT SIZE SPECTRUM SLOPE MODELS ================================================
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
  mutate(area2 = ifelse(area == "BT", "Warm", "Cold"),
         Year_ct = as.integer(Year - min(Year)), 
         area2 = as.factor(area2))

# https://stats.stackexchange.com/questions/352772/what-do-blank-cells-mean-in-the-output-of-prior-summary-in-the-brms-package
# > prior_summary(m1)

# https://vasishth.github.io/Freq_CogSci/from-the-paired-t-test-to-the-linear-mixed-model.html
# No year effect
prior0 <-
  prior(normal(-2, 10), class = "b", coef = area2Cold) +
  prior(normal(-2, 10), class = "b", coef = area2Warm) +
  prior(student_t(3, 0, 2.5), class = "sigma")

m0 <- brm(b ~ -1 + area2,
          family = gaussian(), data = spectra, iter = 4000, cores = 3, chains = 3,
          seed = 9,
          save_pars = save_pars(all = TRUE),
          prior = prior0
          )

summary(m0)
plot(m0)

# Add year as predictor
prior1 <-
  prior(normal(-2, 10), class = "b", coef = area2Cold) +
  prior(normal(-2, 10), class = "b", coef = area2Warm) +
  prior(normal(0, 10), class = "b", coef = Year_ct) +
  prior(student_t(3, 0, 2.5), class = "sigma")

m1 <- brm(
  b ~ - 1 + area2 + Year_ct,
  family = gaussian(), data = spectra, iter = 4000, cores = 3, chains = 3,
  seed = 9,
  save_pars = save_pars(all = TRUE),
  prior = prior1
  )

summary(m1)
plot(m1)

# Interaction between year and area
prior2 <-
  prior(normal(-2, 10), class = "b", coef = area2Cold) +
  prior(normal(-2, 10), class = "b", coef = area2Warm) +
  prior(normal(0, 10), class = "b", coef = area2Warm:Year_ct) +
  prior(normal(0, 10), class = "b", coef = Year_ct) +
  prior(student_t(3, 0, 2.5), class = "sigma")

m2 <- brm(
  b ~ - 1 + area2*Year_ct,
  family = gaussian(), data = spectra, iter = 4000, cores = 3, chains = 3,
  seed = 9,
  save_pars = save_pars(all = TRUE),
  prior = prior2
  )

summary(m2)
plot(m2)

loo_m0 <- loo(m0, moment_match = TRUE)
loo_m1 <- loo(m1, moment_match = TRUE)
loo_m2 <- loo(m2, moment_match = TRUE)

loo_compare(loo_m0, loo_m1, loo_m2)
# elpd_diff se_diff
# m1  0.0       0.0   
# m2 -1.1       0.4   
# m0 -8.7       4.2


# E. PRODUCE FIGURES ===============================================================
pal <- rev(brewer.pal(n = 6, name = "Paired")[c(2, 6)])

##### Plot predictions ============================================================
# Ridge plots of size distribution 
# ggplot(size_df, aes(x = length_group, y = factor(year))) +
#   geom_density_ridges(stat = "binline", scale = 0.95, draw_baseline = FALSE)
# 
# sort(unique(size_df$length_group))
# 
# hist(size_df$length_group, breaks = length(unique(size_df$length_group)))
# 
# t1 <- ggplot(filter(size_df, Area == "BT"), aes(x = length_group, y = after_stat(count/sum(count)))) +   
#   geom_bar() + 
#   scale_y_continuous(labels = scales::percent) +
#   theme_classic() + 
#   ggtitle("BT") +
#   coord_cartesian(expand = 0)
# 
# t2 <- ggplot(filter(size_df, Area == "FM"), aes(x = length_group, y = after_stat(count/sum(count)))) +   
#   geom_bar() + 
#   scale_y_continuous(labels = scales::percent) +
#   theme_classic() + 
#   ggtitle("FM") +
#   coord_cartesian(expand = 0)

# Plot proportion histograms to illustrtate the size-distribution
# By year
# size_df %>%
#   mutate(Area2 = ifelse(Area == "BT", "Warm", "Cold")) %>% 
#   ggplot(., aes(x = length_group, fill = Area2, group = Area2)) +
#   stat_count(mapping = aes(x = length_group, y = ..prop.., group = Area2),
#              position = position_dodge(), alpha = 0.8) +
#   scale_fill_manual(values = rev(pal)) +
#   coord_cartesian(expand = 0) +
#   labs(x = "Length group [cm]") +
#   facet_wrap(~ year, ncol = 2) +
#   guides(fill = FALSE)

p0 <- 
  size_df %>%
  mutate(Area2 = ifelse(Area == "BT", "Warm", "Cold")) %>% 
  ggplot(., aes(x = length_group, fill = Area2, group = Area2)) +
  stat_count(mapping = aes(x = length_group, y = ..prop.., group = Area2),
             position = position_dodge(), alpha = 0.8) +
  scale_fill_manual(values = rev(pal)) +
  theme_light() + 
  coord_cartesian(expand = 0) +
  labs(x = "Length group [cm]") +
  guides(fill = FALSE)

pWord0 <- p0 + theme(text = element_text(size = 12), 
                     axis.text = element_text(angle = 30)#,
                     # legend.spacing.y = unit(0, 'cm'),
                     # legend.key.size = unit(0, "cm"),
                     # legend.title = element_text(size = 10),
                     # legend.text = element_text(size = 10)
)


# Compare 
# (t1 + t2) / pWord0

as.data.frame(fixef(m1)) # Extract "fixed" effects from m2 for plotting the equation 

p1 <- spectra %>%
  ungroup() %>%
  data_grid(Year_ct = seq_range(Year_ct, by = 1),
            area2 = c("Warm", "Cold")) %>%
  mutate(Year = Year_ct + min(spectra$Year)) %>% 
  add_predicted_draws(m1, re_formula = NA) %>%
  ggplot(aes(factor(Year), y = b, color = area2, fill = area2, group = area2)) +
  stat_lineribbon(aes(y = .prediction), .width = c(.90), alpha = 0.2) +
  geom_point(data = spectra, alpha = 0.6, size = 3, shape = 21, color = "white",
             position = position_dodge(width = 0.5)) +
  geom_errorbar(data = spectra, aes(x = factor(Year), ymin = confMin, ymax = confMax,
                                    color = area2, group = area2),
                alpha = 0.25, size = 1, position = position_dodge(width = 0.5), width = 0) +
  stat_lineribbon(aes(y = .prediction), .width = c(.0), alpha = 0.8) +
  scale_fill_manual(values = rev(pal)) +
  scale_color_manual(values = rev(pal)) +
  labs(color = "Area", fill = "Area", x = "Year",
       y = expression(paste("Size-spectrum slope ", italic((b))))) +
  guides(color = guide_legend(override.aes = list(linetype = 0, size = 2, shape = 16, alpha = 0.5,
                                                  color = rev(pal), fill = NA))) +
  annotate("text", 15, -4.4, size = 3.5, color = pal[2],
           label = "y=-3.50 + 0.08×year") + # Cold
  annotate("text", 15, -4.6, size = 3.5, color = pal[1],
           label = "y=-3.13 + 0.08×year") + # Warm
  NULL
  
pWord1 <- p1 + theme(text = element_text(size = 12), 
                     axis.text = element_text(angle = 30),
                     legend.spacing.y = unit(0, 'cm'),
                     legend.key.size = unit(0, "cm"),
                     legend.position = c(0.11, 0.9), 
                     legend.title = element_text(size = 10),
                     legend.text = element_text(size = 10),
                     legend.key = element_blank())

# Now plot the posterior intercepts
# lines manually simply by extracting the fixed effects
m1_fe <- fixef(m1, probs = seq(0, 1, by = 0.05)) %>% as.data.frame()
posterior <- as.array(m1)
dimnames(posterior)

# Define matching palette
# pal2 <- rev(alpha(pal, alpha = 0.8))
# 
# color_scheme_set(rep("white", 6)) # This is to be able to have a fill color with alpha
# 
# intercept_warm <- mcmc_dens(posterior, pars = c("b_area2Warm"),
#                             facet_args = list(nrow = 2)) + 
#   geom_density(fill = pal2[1], color = NA) + 
#   geom_vline(xintercept = m1_fe$Estimate[2], linetype = 1, color = "white") +
#   geom_vline(xintercept = m1_fe$Q10[2], linetype = 2, color = "white") +
#   geom_vline(xintercept = m1_fe$Q90[2], linetype = 2, color = "white") +
#   coord_cartesian(xlim = c(-4.2, -2.5)) + 
#   labs(x = "", y = "") + 
#   theme(axis.text.x = element_text(size = 6))
# 
# intercept_cold <- mcmc_dens(posterior, pars = c("b_area2Cold"),
#                             facet_args = list(nrow = 2)) + 
#   geom_density(fill = pal2[2], color = NA) + 
#   geom_vline(xintercept = m1_fe$Estimate[1], linetype = 1, color = "white") +
#   geom_vline(xintercept = m1_fe$Q10[1], linetype = 2, color = "white") +
#   geom_vline(xintercept = m1_fe$Q90[1], linetype = 2, color = "white") +
#   coord_cartesian(xlim = c(-4.2, -2.5)) + 
#   labs(x = expression(italic(alpha)), y = "") +
#   theme(axis.text.x = element_text(size = 6))

# http://mjskay.github.io/tidybayes/articles/tidy-brms.html
diff <- m1 %>%
  spread_draws(b_area2Cold, b_area2Warm) %>%
  mutate(diff = b_area2Warm - b_area2Cold) 

prop_diff <- summarise(diff, Proportion_of_the_difference_below_0 = sum(diff < 0) / length(diff))

# https://bookdown.org/content/3890/interactions.html
post_diff <- m1 %>%
  spread_draws(b_area2Cold, b_area2Warm) %>%
  mutate(diff = b_area2Warm - b_area2Cold) %>% 
  ggplot(aes(x = diff, fill = stat(x > 0))) +
  stat_halfeye(alpha = 0.5, size = 5, .width = 0) +
  guides(fill = guide_legend(override.aes = list(size = 1, shape = NA, linetype = 0)), color = FALSE) + 
  scale_fill_manual(values = c("grey10", "grey70")) +
  annotate("text", 0.72, 0.95, size = 3.5, label = paste("prop. diff<0=", round(prop_diff, 2), sep = "")) +
  labs(x = expression(~italic(alpha[warm])~-~italic(alpha[cold]))) +
  theme(legend.position = c(0.2, 0.8),
        legend.key.size = unit(0.2, "cm"),
        legend.background = element_blank())

post_diff

# This is the same as in here: https://bookdown.org/content/3890/interactions.html
# Which is also the same as in stat rethink v1
# TEST
# test <- posterior_samples(m1, pars = c("b_area2Cold", "b_area2Warm"))
# t1 <- test %>%
#   gather(key, value) %>%
#   ggplot(aes(x = value, group = key, color = key, fill = key)) +
#   geom_density(alpha = 1/4) +
#   theme(text = element_text(family = "Times"),
#         legend.position = "none") +
#   coord_cartesian(expand = 0)
# 
# test %>%
#   mutate(diff = b_area2Warm -b_area2Cold) %>%
#   summarise(Proportion_of_the_difference_below_0 = sum(diff < 0) / length(diff))
# 
# t2 <- test %>%
#   mutate(diff = b_area2Warm -b_area2Cold) %>%
#   ggplot(aes(x = diff)) +
#   geom_density(alpha = 1/4) +
#   geom_vline(xintercept = 0) +
#   coord_cartesian(expand = 0)
# 
# t1 / t2

post_inter <- m1 %>%
  gather_draws(b_area2Cold, b_area2Warm) %>%
  ggplot(aes(x = .value, fill = .variable, color = .variable)) +
  stat_halfeye(alpha = 0.5, size = 5, .width = c(0.7)) +
  guides(fill = guide_legend(override.aes = list(size = 1, shape = NA, linetype = 0)),
         color = FALSE) + 
  scale_fill_manual(values = rev(pal), labels = c("Cold", "Warm")) +
  scale_color_manual(values = rev(pal)) +
  labs(x = expression(italic(alpha)), fill = "") +
  theme(legend.position = c(0.15, 0.9),
        legend.key.size = unit(0.2, "cm"),
        legend.background = element_blank())

post_inter

pWord1 / ((post_inter/ post_diff) | pWord0) +
  plot_annotation(tag_levels = "A")

ggsave("figures/size_spec.png", width = 6.5, height = 6.5, dpi = 600)


##### Model diagnostics & fit ======================================================
pal_diag <- rev(brewer.pal(n = 3, name = "Dark2"))

# Chain convergence
posterior <- as.array(m1)
dimnames(posterior)

d1 <- mcmc_trace(posterior,
                 pars = c("b_area2Cold", "b_area2Warm", "b_Year_ct", 
                          "sigma"),
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

ggsave("figures/supp/size_spec_diag_fit.png", width = 6.5, height = 8.5, dpi = 600)
