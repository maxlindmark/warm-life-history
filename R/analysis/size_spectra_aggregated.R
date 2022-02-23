#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 2020.06.16: Max Lindmark
#
# Fit models to explore differences in size spectra slopes
# 
# A. Load libraries
# 
# B. Read data
# 
# C. Fit models of spectra all years pooled following Edwards (MEPS 2020)
#
# D. Produce figures
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

# Do NOT group by year here... this is the aggregated
# Group by Area and length group, summarize and get n()
df2 <- df %>%
  group_by(Area, length_group) %>% 
  summarise(catch_n = n()) %>% # Get catch as a column, this is now the number of rows a given length occurs each year and area
  ungroup() %>% 
  mutate(effort_id = paste(Area, sep = ".")) %>% 
  as.data.frame()

# Now we need to get the effort back in there
df_effort <- df %>%
  mutate(effort_id = paste(Area, sep = ".")) %>% 
  select(effort_id, n_nets_year) %>% 
  distinct(effort_id, .keep_all = TRUE) %>% 
  as.data.frame()

# Now do a left_join
df3 <- left_join(df2, df_effort, by = "effort_id") %>% as.data.frame()

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

# Test I get 1 unique row per size class and area
df4 %>%
  group_by(area, min_length_group_cm) %>%
  summarise(n = n()) %>% 
  filter(n == 1) %>% 
  as.data.frame()

# Plot the log biomass cpue as function of log weight
ggplot(df4, aes(log(wmax), log(cpue_biom))) + 
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
res <- species_bins_plots(df4)

# Test with my data
dataBintest <- df4 %>% rename("Number" = "cpue_numbers")

# Following and modifying this vignette:
# https://htmlpreview.github.io/?https://raw.githubusercontent.com/andrew-edwards/sizeSpectra/master/doc/MEPS_IBTS_MLEbins.html

# Forsmark
dataBin_FM <- dataBintest %>%
  filter(area == "FM") %>%
  select(SpecCode, wmin, wmax, Number)

  n_FM = sum(dataBin_FM$Number)
  xmin_FM = min(dataBin_FM$wmin)
  xmax_FM = max(dataBin_FM$wmax)
  
  FM_spectra  = calcLike(negLL.fn = negLL.PLB.bins.species,
                         p = -1.5,
                         suppress.warnings = FALSE,
                         vecDiff = 0.5, # Default is 0.5, I get warning though
                         dataBinForLike = dataBin_FM,
                         n = n_FM,
                         xmin = xmin_FM,
                         xmax = xmax_FM)
  
  FM_spectra_df = data.frame(Year = 2000,
                             xmin = xmin_FM,
                             xmax = xmax_FM,
                             n = n_FM,
                             b = FM_spectra$MLE,
                             confMin = FM_spectra$conf[1],
                             confMax = FM_spectra$conf[2])
  
  
# Biotest
dataBin_BT <- dataBintest %>%
  filter(area == "BT") %>%
  select(SpecCode, wmin, wmax, Number)

n_BT = sum(dataBin_BT$Number)
xmin_BT = min(dataBin_BT$wmin)
xmax_BT = max(dataBin_BT$wmax)

BT_spectra  = calcLike(negLL.fn = negLL.PLB.bins.species,
                       p = -1.5,
                       suppress.warnings = FALSE,
                       vecDiff = 0.5, # Default is 0.5, I get warning though
                       dataBinForLike = dataBin_BT,
                       n = n_BT,
                       xmin = xmin_BT,
                       xmax = xmax_BT)

BT_spectra_df = data.frame(Year = 2000,
                           xmin = xmin_BT,
                           xmax = xmax_BT,
                           n = n_BT,
                           b = BT_spectra$MLE,
                           confMin = BT_spectra$conf[1],
                           confMax = BT_spectra$conf[2])

FM_spectra_df
BT_spectra_df


#** Loop and extract data for plotting =============================================
# https://htmlpreview.github.io/?https://raw.githubusercontent.com/andrew-edwards/sizeSpectra/master/doc/MEPS_IBTS_recommend.html
# These calculations are to get the required input for the recommended plot
# (see ?sizeSpectra::ISD_bin_plot for the structure). This could maybe be functionalised
# like the plotting function ISD_bin_plot().

# Forsmark ==========================================================================
dataBin_FM$Year <- 2000
dataRecommend.isd_FM = dplyr::select(dataBin_FM, Year, wmin, wmax, Number)

data.year.fm <- dplyr::arrange(dataRecommend.isd_FM, desc(wmin))

sumNumber = sum(data.year.fm$Number)

wmin.vec = data.year.fm$wmin
wmax.vec = data.year.fm$wmax
num.vec  = data.year.fm$Number
  
countGTEwmin = rep(NA, length(num.vec)) # to do a manual count
lowCount = countGTEwmin
highCount = countGTEwmin
  
for(iii in 1:length(countGTEwmin))
  {
    countGTEwmin[iii] = sum( (wmin.vec >= wmin.vec[iii]) * num.vec)
    lowCount[iii]  = sum( (wmin.vec >= wmax.vec[iii]) * num.vec)
    highCount[iii] = sum( (wmax.vec >  wmin.vec[iii]) * num.vec)
  }

data.year.fm = cbind(data.year.fm,
                 "countGTEwmin" = countGTEwmin,
                 "lowCount" = lowCount,
                 "highCount" = highCount)

data.year.fm = dplyr::tbl_df(data.year.fm)
  
xlim.global = c(min(data.year.fm$wmin),
                max(data.year.fm$wmax))   # x-axis limits to be common for all plots

ISD_bin_plot(data.year = data.year.fm,
             b.MLE = FM_spectra_df$b,
             b.confMin = FM_spectra_df$confMin,
             b.confMax = FM_spectra_df$confMax,
             year = 2000,
             xlim = xlim.global,
             xmin = FM_spectra_df$xmin,
             xmax = FM_spectra_df$xmax
)


# Biotest ==========================================================================
dataBin_BT$Year <- 2000

dataRecommend.isd_BT = dplyr::select(dataBin_BT, Year, wmin, wmax, Number)

data.year.bt <- dplyr::arrange(dataRecommend.isd_BT, desc(wmin))

sumNumber = sum(data.year$Number)

wmin.vec = data.year.bt$wmin
wmax.vec = data.year.bt$wmax
num.vec  = data.year.bt$Number

countGTEwmin = rep(NA, length(num.vec)) # to do a manual count
lowCount = countGTEwmin
highCount = countGTEwmin

for(iii in 1:length(countGTEwmin))
{
  countGTEwmin[iii] = sum( (wmin.vec >= wmin.vec[iii]) * num.vec)
  lowCount[iii]  = sum( (wmin.vec >= wmax.vec[iii]) * num.vec)
  highCount[iii] = sum( (wmax.vec >  wmin.vec[iii]) * num.vec)
}

data.year.bt = cbind(data.year.bt,
                  "countGTEwmin" = countGTEwmin,
                  "lowCount" = lowCount,
                  "highCount" = highCount)

data.year.bt = dplyr::tbl_df(data.year.bt)

xlim.global = c(min(data.year.bt$wmin),
                max(data.year.bt$wmax))   # x-axis limits to be common for all plots

ISD_bin_plot(data.year = data.year.bt,
             b.MLE = BT_spectra_df$b,
             b.confMin = BT_spectra_df$confMin,
             b.confMax = BT_spectra_df$confMax,
             year = 2000,
             xlim = xlim.global,
             xmin = BT_spectra_df$xmin,
             xmax = BT_spectra_df$xmax
)


# D. PRODUCE FIGURES ===============================================================

# Combine areas and plot ===========================================================
# https://htmlpreview.github.io/?https://raw.githubusercontent.com/andrew-edwards/sizeSpectra/master/doc/MEPS_IBTS_MLEbins.html

FM_spectra = dplyr::tbl_df(FM_spectra_df)
FM_spectra = dplyr::mutate(FM_spectra, stdErr = (abs(confMin-b) + abs(confMax-b))/(2*1.96) )
FM_spectra = dplyr::mutate(FM_spectra, area = "FM")

BT_spectra = dplyr::tbl_df(BT_spectra_df)
BT_spectra = dplyr::mutate(BT_spectra, stdErr = (abs(confMin-b) + abs(confMax-b))/(2*1.96) )
BT_spectra = dplyr::mutate(BT_spectra, area = "BT")

spectra <- rbind(BT_spectra, FM_spectra)

pal <- rev(brewer.pal(n = 6, name = "Paired")[c(2, 6)])

# Now make a ggplot-version to put in the paper as an example
# First put in in a dataframe instead of adding them as lines in the base plot
# BT (warm)
data_year_bt <- data.year.bt
data_year_bt$area <- "Warm"
sumNumber = sum(data_year_bt$Number)

xmin = BT_spectra$xmin
xmax = BT_spectra$xmax
x.PLB = seq(xmin_BT, xmax_BT, length = 10001)
b.MLE = BT_spectra$b
b.confMin = BT_spectra$confMin
b.confMax = BT_spectra$confMax

data_year_bt2 <- data.frame(x.PLB = x.PLB,
                            y.PLB = (1 - pPLB(x = x.PLB, b = b.MLE, xmin = min(x.PLB), xmax = max(x.PLB))) * sumNumber,
                            y.PLB.confMin = (1 - pPLB(x = x.PLB, b = b.confMin, xmin = min(x.PLB), xmax = max(x.PLB))) * sumNumber,
                            y.PLB.confMax = (1 - pPLB(x = x.PLB, b = b.confMax, xmin = min(x.PLB), xmax = max(x.PLB))) * sumNumber,
                            area = "Warm")

# FM (cold)
data_year_fm <- data.year.fm
data_year_fm$area <- "Cold"
sumNumber = sum(data_year_fm$Number)

xmin = FM_spectra$xmin
xmax = FM_spectra$xmax
x.PLB = seq(xmin_FM, xmax_FM, length = 10001)
b.MLE = FM_spectra$b
b.confMin = FM_spectra$confMin
b.confMax = FM_spectra$confMax

data_year_fm2 <- data.frame(x.PLB = x.PLB,
                            y.PLB = (1 - pPLB(x = x.PLB, b = b.MLE, xmin = min(x.PLB), xmax = max(x.PLB))) * sumNumber,
                            y.PLB.confMin = (1 - pPLB(x = x.PLB, b = b.confMin, xmin = min(x.PLB), xmax = max(x.PLB))) * sumNumber,
                            y.PLB.confMax = (1 - pPLB(x = x.PLB, b = b.confMax, xmin = min(x.PLB), xmax = max(x.PLB))) * sumNumber,
                            area = "Cold")

data_year <- rbind(data_year_bt, data_year_fm)
data_year2 <- rbind(data_year_bt2, data_year_fm2)

b_warm <- round(BT_spectra$b, 2) 
b_cold <- round(FM_spectra$b, 2)
b_warm_upr <- round(BT_spectra$confMax, 2)
b_cold_upr <- round(FM_spectra$confMax, 2)
b_warm_lwr <- round(BT_spectra$confMin, 2)
b_cold_lwr <- round(FM_spectra$confMin, 2)

data_year2 <- data_year2 %>% 
  mutate(area_plot = ifelse(area == "Warm", "Heat", "Ref"))

pmle <- ggplot(data_year2) +
  geom_rect(data = data_year,
            aes(xmin = wmin, xmax = wmax, ymin = lowCount, ymax = highCount, fill = area),
            alpha = 0.2) +
  geom_line(aes(x.PLB, y.PLB, color = area_plot), size = 1.3, alpha = 0.8) +
  geom_line(aes(x.PLB, y.PLB.confMin, color = area_plot), linetype = 2, alpha = 0.8) +
  geom_line(aes(x.PLB, y.PLB.confMax, color = area_plot), linetype = 2, alpha = 0.8) +
  scale_color_manual(values = pal, name = "Area") +
  scale_fill_manual(values = rev(pal)) + 
  scale_x_log10() +
  scale_y_log10() +
  guides(fill = FALSE) + 
  coord_cartesian(ylim = c(0.1, 1000)) +
  labs(x = "Body mass, x [g]", y = "Number of values ≥ x") +
  theme(legend.position = c(0.8, 0.9)) +
  guides(color = guide_legend(override.aes = list(linetype = 1, size = 2, alpha = 0.3, color = pal))) +
  theme(text = element_text(size = 12), 
        legend.position = c(0.9, 0.9),
        legend.spacing.y = unit(0, 'cm'),
        legend.key.size = unit(0, "cm"),
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 8))

pmle

# Now plot the CI's
BT_spectra
spec_tot <- rbind(BT_spectra, FM_spectra) %>% 
  mutate(Area = ifelse(area == "BT", "Heated", "Reference"))

est_plot <- ggplot(spec_tot, aes(Area, b, color = Area)) + 
  geom_point(size = 3) +
  scale_color_manual(values = pal) +
  guides(color = "none") +
  labs(y = expression(paste("Size-spectrum exponent ", italic((gamma))))) +
  geom_errorbar(aes(x = Area, ymin = confMin, ymax = confMax), width = 0.05, size = 0.75)

phist <- 
  size_df %>%
  mutate(Area2 = ifelse(Area == "BT", "Warm", "Cold")) %>% 
  ggplot(., aes(x = length_group, fill = Area2, group = Area2)) +
  stat_count(mapping = aes(x = length_group, y = ..prop.., group = Area2),
             position = position_dodge(), alpha = 0.8) +
  scale_fill_manual(values = rev(pal)) +
  theme_light() + 
  coord_cartesian(expand = 0) +
  labs(x = "Length group [cm]", y = "Proportion") +
  guides(fill = FALSE) +
  theme(text = element_text(size = 12))

# Combine plots
pmle / (est_plot | phist) + plot_annotation(tag_levels = "A") #+ plot_layout(widths = c(5, 1, 1))

ggsave("figures/size_spec_v2.png", width = 6.5, height = 6.5, dpi = 600)


