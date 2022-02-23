#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 2022.02.2: Max Lindmark
#
# Check sensitivity to sample individual sample size for estimating b
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

# This is the last step where we have 1 row = 1 ind
# I will use 1995 as an example...

df_bt <- df %>% filter(year == 1995 & Area == "BT")

# Now make 10 random subsets of these data, starting from 50% and up
df_bt_50 <- df_bt[sample(nrow(df_bt), 0.50*nrow(df_bt)), ]
df_bt_55 <- df_bt[sample(nrow(df_bt), 0.55*nrow(df_bt)), ]
df_bt_60 <- df_bt[sample(nrow(df_bt), 0.60*nrow(df_bt)), ]
df_bt_65 <- df_bt[sample(nrow(df_bt), 0.65*nrow(df_bt)), ]
df_bt_70 <- df_bt[sample(nrow(df_bt), 0.70*nrow(df_bt)), ]
df_bt_75 <- df_bt[sample(nrow(df_bt), 0.75*nrow(df_bt)), ]
df_bt_80 <- df_bt[sample(nrow(df_bt), 0.80*nrow(df_bt)), ]
df_bt_85 <- df_bt[sample(nrow(df_bt), 0.85*nrow(df_bt)), ]
df_bt_90 <- df_bt[sample(nrow(df_bt), 0.90*nrow(df_bt)), ]
df_bt_95 <- df_bt[sample(nrow(df_bt), 0.95*nrow(df_bt)), ]


# Now we want to process the data a bit further... Following Edwards sizeSpectra package
# we want the data as follows:
# Year 	SpecCode 	LngtClass 	Number 	LWa 	LWb 	bodyMass 	Biomass
# Hence, for each year, we need to calculate the CPUE
# https://htmlpreview.github.io/?https://raw.githubusercontent.com/andrew-edwards/sizeSpectra/master/doc/MEPS_IBTS_MLEbins.html

# Group by year, Area and length group, summarize and get n()
df2_50 <- df_bt_50 %>%
  group_by(year, Area, length_group) %>% 
  summarise(catch_n = n()) %>% # Get catch as a column, this is now the number of rows a given length occurs each year and area
  ungroup() %>% 
  mutate(effort_id = paste(year, Area, sep = ".")) %>% 
  as.data.frame()

df2_55 <- df_bt_55 %>%
  group_by(year, Area, length_group) %>% 
  summarise(catch_n = n()) %>% # Get catch as a column, this is now the number of rows a given length occurs each year and area
  ungroup() %>% 
  mutate(effort_id = paste(year, Area, sep = ".")) %>% 
  as.data.frame()

df2_60 <- df_bt_60 %>%
  group_by(year, Area, length_group) %>% 
  summarise(catch_n = n()) %>% # Get catch as a column, this is now the number of rows a given length occurs each year and area
  ungroup() %>% 
  mutate(effort_id = paste(year, Area, sep = ".")) %>% 
  as.data.frame()

df2_65 <- df_bt_65 %>%
  group_by(year, Area, length_group) %>% 
  summarise(catch_n = n()) %>% # Get catch as a column, this is now the number of rows a given length occurs each year and area
  ungroup() %>% 
  mutate(effort_id = paste(year, Area, sep = ".")) %>% 
  as.data.frame()

df2_70 <- df_bt_70 %>%
  group_by(year, Area, length_group) %>% 
  summarise(catch_n = n()) %>% # Get catch as a column, this is now the number of rows a given length occurs each year and area
  ungroup() %>% 
  mutate(effort_id = paste(year, Area, sep = ".")) %>% 
  as.data.frame()

df2_75 <- df_bt_75 %>%
  group_by(year, Area, length_group) %>% 
  summarise(catch_n = n()) %>% # Get catch as a column, this is now the number of rows a given length occurs each year and area
  ungroup() %>% 
  mutate(effort_id = paste(year, Area, sep = ".")) %>% 
  as.data.frame()

df2_80 <- df_bt_80 %>%
  group_by(year, Area, length_group) %>% 
  summarise(catch_n = n()) %>% # Get catch as a column, this is now the number of rows a given length occurs each year and area
  ungroup() %>% 
  mutate(effort_id = paste(year, Area, sep = ".")) %>% 
  as.data.frame()

df2_85 <- df_bt_85 %>%
  group_by(year, Area, length_group) %>% 
  summarise(catch_n = n()) %>% # Get catch as a column, this is now the number of rows a given length occurs each year and area
  ungroup() %>% 
  mutate(effort_id = paste(year, Area, sep = ".")) %>% 
  as.data.frame()

df2_90 <- df_bt_90 %>%
  group_by(year, Area, length_group) %>% 
  summarise(catch_n = n()) %>% # Get catch as a column, this is now the number of rows a given length occurs each year and area
  ungroup() %>% 
  mutate(effort_id = paste(year, Area, sep = ".")) %>% 
  as.data.frame()

df2_95 <- df_bt_95 %>%
  group_by(year, Area, length_group) %>% 
  summarise(catch_n = n()) %>% # Get catch as a column, this is now the number of rows a given length occurs each year and area
  ungroup() %>% 
  mutate(effort_id = paste(year, Area, sep = ".")) %>% 
  as.data.frame()


# Now we need to get the effort back in there
df_effort_50 <- df_bt_50 %>%
  mutate(effort_id = paste(year, Area, sep = ".")) %>% 
  select(effort_id, n_nets_year) %>% 
  distinct(effort_id, .keep_all = TRUE) %>% 
  as.data.frame()

df_effort_55 <- df_bt_55 %>%
  mutate(effort_id = paste(year, Area, sep = ".")) %>% 
  select(effort_id, n_nets_year) %>% 
  distinct(effort_id, .keep_all = TRUE) %>% 
  as.data.frame()

df_effort_60 <- df_bt_60 %>%
  mutate(effort_id = paste(year, Area, sep = ".")) %>% 
  select(effort_id, n_nets_year) %>% 
  distinct(effort_id, .keep_all = TRUE) %>% 
  as.data.frame()

df_effort_65 <- df_bt_65 %>%
  mutate(effort_id = paste(year, Area, sep = ".")) %>% 
  select(effort_id, n_nets_year) %>% 
  distinct(effort_id, .keep_all = TRUE) %>% 
  as.data.frame()

df_effort_70 <- df_bt_70 %>%
  mutate(effort_id = paste(year, Area, sep = ".")) %>% 
  select(effort_id, n_nets_year) %>% 
  distinct(effort_id, .keep_all = TRUE) %>% 
  as.data.frame()

df_effort_75 <- df_bt_75 %>%
  mutate(effort_id = paste(year, Area, sep = ".")) %>% 
  select(effort_id, n_nets_year) %>% 
  distinct(effort_id, .keep_all = TRUE) %>% 
  as.data.frame()

df_effort_80 <- df_bt_80 %>%
  mutate(effort_id = paste(year, Area, sep = ".")) %>% 
  select(effort_id, n_nets_year) %>% 
  distinct(effort_id, .keep_all = TRUE) %>% 
  as.data.frame()

df_effort_85 <- df_bt_85 %>%
  mutate(effort_id = paste(year, Area, sep = ".")) %>% 
  select(effort_id, n_nets_year) %>% 
  distinct(effort_id, .keep_all = TRUE) %>% 
  as.data.frame()

df_effort_90 <- df_bt_90 %>%
  mutate(effort_id = paste(year, Area, sep = ".")) %>% 
  select(effort_id, n_nets_year) %>% 
  distinct(effort_id, .keep_all = TRUE) %>% 
  as.data.frame()

df_effort_95 <- df_bt_95 %>%
  mutate(effort_id = paste(year, Area, sep = ".")) %>% 
  select(effort_id, n_nets_year) %>% 
  distinct(effort_id, .keep_all = TRUE) %>% 
  as.data.frame()

# Now do a left_join
df3_50 <- left_join(df2_50, df_effort_50, by = "effort_id") %>% as.data.frame() %>% mutate(sample = "50")
df3_55 <- left_join(df2_55, df_effort_55, by = "effort_id") %>% as.data.frame() %>% mutate(sample = "55")
df3_60 <- left_join(df2_60, df_effort_60, by = "effort_id") %>% as.data.frame() %>% mutate(sample = "60")
df3_65 <- left_join(df2_65, df_effort_65, by = "effort_id") %>% as.data.frame() %>% mutate(sample = "65")
df3_70 <- left_join(df2_70, df_effort_70, by = "effort_id") %>% as.data.frame() %>% mutate(sample = "70")
df3_75 <- left_join(df2_75, df_effort_75, by = "effort_id") %>% as.data.frame() %>% mutate(sample = "75")
df3_80 <- left_join(df2_80, df_effort_80, by = "effort_id") %>% as.data.frame() %>% mutate(sample = "80")
df3_85 <- left_join(df2_85, df_effort_85, by = "effort_id") %>% as.data.frame() %>% mutate(sample = "85")
df3_90 <- left_join(df2_90, df_effort_90, by = "effort_id") %>% as.data.frame() %>% mutate(sample = "90")
df3_95 <- left_join(df2_95, df_effort_95, by = "effort_id") %>% as.data.frame() %>% mutate(sample = "95")

# Merge data...

df3 <- rbind(df3_50, df3_55, df3_60, df3_65, df3_70, df3_75, df3_80, df3_85, df3_90, df3_95)


# Now do some additional calculations (weights etc, to match data structure in sizeSpectra package)
# LW parameters from FishBase (2020.09.24)
a <- 0.01 
b <- 3.08

df4 <- df3 %>% 
  ungroup() %>%
  group_by(sample) %>% 
  rename("area" = "Area",
         "min_length_group_cm" = "length_group") %>% # What we have now is the minimum length in each class
  mutate(max_length_group_cm = min_length_group_cm + 2.4, # Get max length in bin (2.5 cm length-classes)
         wmin = a*min_length_group_cm^b,     # Get min mass in bin
         wmax = a*max_length_group_cm^b) %>% # Get max mass in bin
  mutate(cpue_numbers = catch_n/n_nets_year, # Get numbers CPUE, divide by the previously create n_nets, which is # of unique net ID's in each area and year
         cpue_biom = (catch_n*((wmin + wmax)/2))/n_nets_year) %>% # Get biomass CPUE, use mean of mass in size range
  mutate(SpecCode = "Perch")

# Test I get 1 unique row per size class, year, area and sample
df4 %>%
  group_by(year, area, min_length_group_cm, sample) %>%
  summarise(n = n()) %>% 
  ungroup() %>% 
  distinct(n)

# Plot the log biomass cpue as function of log weight
ggplot(df4, aes(log(wmax), log(cpue_biom), color = factor(year))) + 
  stat_smooth(se = FALSE) +
  geom_point() +
  facet_wrap(~ sample) +
  scale_color_viridis(discrete = TRUE) + 
  coord_cartesian(expand = 0) + 
  theme(aspect.ratio = 1)

# Plot the total number of catches by length-group in this processed data
# Compare with p2, should be the same plot...
ggplot(df4, aes(min_length_group_cm, catch_n)) + 
  geom_bar(stat = "identity") + 
  facet_wrap(~ sample, scales = "free") +
  theme_classic() + 
  coord_cartesian(expand = 0)

# Now that we have CPUE and not COUNT, we can see if the filter to maximum catch size is
# appropriate
p3 <- ggplot(df4, aes(factor(round(wmin, digits = 3)), cpue_biom)) + 
  geom_bar(stat = "identity") + 
  facet_wrap(~ sample, scales = "free") +
  theme_classic() + 
  coord_cartesian(expand = 0)

# Filter out small fish
sort(unique(df4$wmin))

df4 <- df4 %>% filter(wmin > 105) %>% as.data.frame()

# Plot again
p4 <- ggplot(df4, aes(factor(round(wmin, digits = 3)), cpue_biom)) + 
  geom_bar(stat = "identity") + 
  facet_wrap(~ sample, scales = "free") +
  theme_classic() + 
  coord_cartesian(expand = 0)

p3/p4




# C. FIT SIZE SPECTRUM SLOPE MODELS ================================================
# Test with my data
dataBintest <- df4 %>% rename("Year" = "year",
                              "Number" = "cpue_numbers")

colnames(dataBintest)
head(dataBintest)
sort(unique(dataBintest$effort_id))

# Biotest ==========================================================================
# Following and modifying this vignette:
# https://htmlpreview.github.io/?https://raw.githubusercontent.com/andrew-edwards/sizeSpectra/master/doc/MEPS_IBTS_MLEbins.html

dataBin_BT <- dataBintest %>%
  filter(area == "BT") %>%
  select(SpecCode, wmin, wmax, Number, Year, sample)


#** Loop through all years and model ===============================================
# Get a vector of all years to loop over
samples = sort(unique(dataBin_BT$sample))

for(iii in sort(unique(samples)))
{
  dataBin_BT_loop = dplyr::filter(dataBin_BT, sample == iii)
  
  n_BT = sum(dataBin_BT_loop$Number)
  xmin_BT = min(dataBin_BT_loop$wmin)
  xmax_BT = max(dataBin_BT_loop$wmax)
  
  BT_spectra_one_sample  = calcLike(negLL.fn = negLL.PLB.bins.species,
                                    p = -1.5,
                                    suppress.warnings = FALSE,
                                    vecDiff = 2, # Default is 0.5, I get warning though
                                    dataBinForLike = dataBin_BT_loop,
                                    n = n_BT,
                                    xmin = xmin_BT,
                                    xmax = xmax_BT)
  
  if(iii == 50)
  {
    BT_spectra = data.frame(sample = iii,
                            xmin = xmin_BT,
                            xmax = xmax_BT,
                            n = n_BT,
                            b = BT_spectra_one_sample$MLE,
                            confMin = BT_spectra_one_sample$conf[1],
                            confMax = BT_spectra_one_sample$conf[2])
  } else {
    BT_spectra = rbind(BT_spectra,
                       c(iii,
                         xmin_BT,
                         xmax_BT,
                         n_BT,
                         BT_spectra_one_sample$MLE,
                         BT_spectra_one_sample$conf[1],
                         BT_spectra_one_sample$conf[2]))
  }
}

BT_spectra <- BT_spectra %>% mutate(b_num = as.numeric(b),
                                    confMin_num = as.numeric(confMin),
                                    confMax_num = as.numeric(confMax))

# Plot
str(BT_spectra)

ggplot(BT_spectra, aes(sample, b_num)) + 
  geom_point() + 
  geom_errorbar(aes(x = sample, ymin = confMin_num, ymax = confMax_num), width = 0.2)


