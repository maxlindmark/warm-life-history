#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 2020.06.16: Max Lindmark
#
# Fit models of size spectra following Edwards (MEPS 2020)
# 
# A. Load libraries
# 
# B. Read data
# 
# C. Fit models of spectra by year using Andrew's model
#
# D. Fit models of slopes and percentiles
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
         cpue_biom = (catch_n*((wmin + wmax)/2))/n_nets) %>% # Get biomass CPUE, use mean of mass in size range
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



# D. FIT MODELS OF SLOPES ==========================================================
# Size spectrum slopes =============================================================

# Mean centre year variable
spectra <- spectra %>%
  mutate(Year_ct = Year - min(Year),
         Year_ct = as.integer(Year_ct))

# Without interaction
m1 <- brm(b ~ Year_ct + area + (1|Year_ct),
          family = gaussian(), data = spectra, iter = 4000, cores = 3, chains = 3,
          save_all_pars = TRUE,
          control = list(adapt_delta = 0.99))

# With interaction
m2 <- brm(b ~ Year_ct * area + (1|Year_ct),
          family = gaussian(), data = spectra, iter = 4000, cores = 3, chains = 3,
          save_all_pars = TRUE,
          control = list(adapt_delta = 0.99))

loo_m1 <- loo(m1)
loo_m2 <- loo(m2, moment_match = TRUE)

loo_compare(loo_m1, loo_m2)

summary(m1)
plot(m1)

# Plot
pal <- rev(brewer.pal(n = 6, name = "Paired")[c(2, 6)])

p1 <- spectra %>%
  ungroup() %>%
  data_grid(Year_ct = seq_range(Year_ct, by = 1),
            area = c("FM", "BT")) %>%
  add_predicted_draws(m1, re_formula = NA) %>%
  mutate(Year = Year_ct + min(spectra$Year)) %>% 
  ggplot(aes(x = Year, y = b, color = area, fill = area)) +
  stat_lineribbon(aes(y = .prediction), .width = c(.90, .50), alpha = 1/4) +
  geom_point(data = spectra, alpha = 1, size = 3) +
  scale_fill_manual(values = pal, labels = c("Warm", "Cold")) +
  scale_color_manual(values = pal, labels = c("Warm", "Cold")) +
  labs(color = "Area", fill = "Area",
       x = "Year",
       y = expression(paste("Size-spectrum slope ", italic((b))))) +
  NULL

pWord1 <- p1 + theme(text = element_text(size = 12), 
                     legend.position = c(0.1, 0.9), 
                     legend.title = element_text(size = 10),
                     legend.text = element_text(size = 10))

ggsave("figures/size_spectra_slopes.png", width = 6.5, height = 6.5, dpi = 600)


# Size percentiles =================================================================
#perc_dat <- 
  
# df %>%
#   group_by(Area, year) %>% 
#   summarise(q95 = quantile(length_group, 0.95)) %>% 
#   as.data.frame() %>% 
#   ggplot(., aes(year, q95, color = Area)) + geom_point() + stat_smooth(method = "lm")
# 
# quantile(filter(df, Area == "BT")$length_group, c(0.5, 0.95, .999))
# quantile(filter(df, Area == "FM")$length_group, c(0.5, 0.95, .999))

# Skip for now, doesn't work well with size-groups...
