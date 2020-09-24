#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 2020.09.02: Malin Karlsson & Max Lindmark
#
# Read in catch data, clean, filter and process. FM data, unlike the BT data, 
# comes in three data files which correspond to the tabs in the spreadsheet:
# FM09_alla_grundutrdag_combinedbyAG_190916
# They all have different layout so I will clean them separately...
#
# A. Load libraries
# 
# B. Read & process FORSMARK (REFERENCE) data sets
#
# C. Read & process BIOTEST (WARMING) data set
# 
# D. Merge areas and standardize length codes!
# 
# TO DO: REMOVE HASHTAGS IF WE DECIDE TO USE FULL AREAS...
# 
# The data contains two different length group standards std 2 (2.5 cm intervals) and
# std 3 (1 cm intervals). I need to rewrite std 3 to std 2 to get comparable data. 
# See the instruction manual for how they are defined ("Standardisering av längdgrupper")
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# A. LOAD LIBRARIES ================================================================
rm(list = ls())

# Load libraries (install first if needed)
library(tidyverse)
library(tidylog)
library(RColorBrewer)
library(patchwork)

# Print package versions for versions
sessionInfo() 

# other attached packages:
# [1] patchwork_1.0.1     RColorBrewer_1.1-2  tidylog_1.0.2       forcats_0.5.0       stringr_1.4.0       dplyr_1.0.2        
# [7] purrr_0.3.4         readr_1.3.1         tidyr_1.1.0         tibble_3.0.3        ggplot2_3.3.2       tidyverse_1.3.0    
# [13] sizeSpectra_1.0.0.0


# B. READ AND PROCESS FORSMARK DATA ================================================
#** FM 83-86 =======================================================================
# Need to set fileEncoding here, else error: "invalid multibyte string 18"
df83 <- read.csv("data/raw/Catch_data_FM09__1983-86_190916.csv", sep = ";", fileEncoding = "latin1")

# Tidy data. Remove unnecessary columns
df83 <- df83 %>%
  filter(Art == "ABBO", Årtal < 2004) %>%
  rename(year = Årtal, week = Vecka, day = Dag, effort = Ansträngning, species = Art,
         weight = Vikt, n = Antal) %>%
  select(-c(Vtn.stånd, VindRiktn.I, VindSt.I, Vind_upp_rikn, VindSt.Upp, Ström_I_rikn,
            Ström_upp_rikn, Salthalt_I_yta, Salthalt_I_botten, Salthalt_upp_yta, 
            Salthalt_upp_botten, Drift_i, Drift_u, Drift_dim, Siktdjup, Lufttryck_i,
            Lufttryck_upp, Sjuk_kontroll, X..))

# We would much rather prefer to have a column for each length, so that 1 row = observation.
dat83 <- df83 %>%
  gather(length, n2, c(X1, X2, X3, X4, X5, X6, X7, X8, X9, X10, X11, X12, X13, X14), na.rm = T)

head(dat83) # Now the data is in a long, tidy format. 
head(df83)

# Test it is correct:
subset(df83, year == 1983 & week == 30 & day == 2 & Stations_namn == "Asphällan")
subset(dat83, year == 1983 & week == 30 & day == 2 & Stations_namn == "Asphällan")

# Create a new empty column for numeric "length"
dat83$length_group <- as.numeric(substring(dat83$length, 2))

head(dat83)                            
str(dat83)

# Our n2 column tells us how many fish are in that size-group. We want one row for each observation!
# So we need to repeat that observation by n2.

# Now repeat
dat83 <- dat83[rep(seq(nrow(dat83)), dat83$n2),]
head(dat83, 50)


#** FM 1987-1990 ===================================================================
df87 <- read.csv("data/raw/Catch_data_FM09__1987-90_190916.csv", sep = ";", fileEncoding = "latin1")

# Tidy data. Remove unnecessary columns
df87 <- df87 %>%
  filter(Art == "ABBO", Årtal < 2004) %>%
  rename(year = Årtal, week = Vecka, day = Dag, effort = Ansträngning, species = Art,
         weight = Vikt, n = Antal) %>%
  select(-c(Vtn.stånd, VindRiktn.I, VindSt.I, Vind_upp_rikn, VindSt.Upp, Ström_I_rikn,
            Ström_upp_rikn, Salthalt_I_yta, Salthalt_I_botten, Salthalt_upp_yta, 
            Salthalt_upp_botten, Drift_i, Drift_u, Drift_dim, Siktdjup, Lufttryck_i,
            Lufttryck_upp, Sjuk_kontroll, X..))

# Convert to long data frame
# Now the X-columns will be rows AND you will get a new column that takes the old "colmn" value and put that
# in the new column n2. DOUBLE CHECK!!
dat87 <- df87 %>%
  gather(length, n2, c(X6, X9, X11, X14, X16, X19, X21, X24, X26, X29, X31, X34, X36, X39,
                       X41, X44, X46, X49, X51, X56, X59, X61, X64, X71), na.rm = T)

# Test it is correct:
subset(df87, year == 1987 & week == 32 & day == 2 & Stations_namn == "Asphällan")
subset(dat87, year == 1987 & week == 32 & day == 2 & Stations_namn == "Asphällan")

# Create a new empty column for numeric "length"
dat87$length_group <- as.numeric(substring(dat87$length, 2))

# Our n2 column tells us how many fish are in that size-group. We want one row for each observation!
# So we need to repeat that observation by n2.
dat87 <- dat87[rep(seq(nrow(dat87)), dat87$n2),]
head(dat87, 50)

# Check its correct:
head(subset(dat87, n2 == 3), 50)


#** FM 1991-2000 ===================================================================
df91 <- read.csv("data/raw/Catch_data_FM09__1991-00_190916.csv", sep = ";", fileEncoding = "latin1")

# Tidy data. Remove unnecessary columns
df91 <- df91 %>%
  filter(Art == "ABBO", Årtal < 2004) %>%
  rename(year = Årtal, week = Vecka, day = Dag, effort = Ansträngning, species = Art,
         weight = Vikt, n = Antal) %>%
  select(-c(Vtn.stånd, VindRiktn.I, VindSt.I, Vind_upp_rikn, VindSt.Upp, Ström_I_rikn,
            Ström_upp_rikn, Salthalt_I_yta, Salthalt_I_botten, Salthalt_upp_yta, 
            Salthalt_upp_botten, Drift_i, Drift_u, Drift_dim, Siktdjup, Lufttryck_i,
            Lufttryck_upp, Sjuk_kontroll))

# Go from wide to long format
dat91 <- df91 %>%
  gather(length, n2, c(X9, X11, X14, X16, X19, X21, X24, X26, X29, X31, X34, X36, X39, X41, X44,
                       X46, X49, X51, X54, X56,X59, X61, X66, X69, X71, X76, X81, X101), na.rm = T)

# Test it is correct:
subset(df91, year == 1992 & week == 32 & day == 2 & Stations_namn == "Asphällan")
subset(dat91, year == 1992 & week == 32 & day == 2 & Stations_namn == "Asphällan")

# Create a new empty column for numeric "length"
dat91$length_group <- as.numeric(substring(dat91$length, 2))

# Our n2 column tells us how many fish are in that size-group. We want one row for each observation!
# So we need to repeat that observation by n2.
dat91 <- dat91[rep(seq(nrow(dat91)), dat91$n2),]
head(dat91, 50)

# Check its correct:
head(subset(dat91, n2 == 3), 50)


#** FM 2001-2006 ===================================================================
df01 <- read.csv("data/raw/Catch_data_FM09__2001-06_190916.csv", sep = ";", fileEncoding = "latin1")

# Tidy data. Remove unnecessary columns
df01 <- df01 %>%
  filter(Art == "ABBO", Årtal < 2004) %>%
  rename(year = Årtal, week = Vecka, day = Dag, effort = Ansträngning, species = Art,
         weight = Vikt, n = Antal) %>%
  select(-c(Vtn.stånd, VindRiktn.I, VindSt.I, Vind_upp_rikn, VindSt.Upp, Ström_I_rikn,
            Ström_upp_rikn, Salthalt_I_yta, Salthalt_I_botten, Salthalt_upp_yta, 
            Salthalt_upp_botten, Drift_i, Drift_u, Drift_dim, Siktdjup, Lufttryck_i,
            Lufttryck_upp, Sjuk_kontroll, X..))

# Go from wide to long format
dat01 <- df01 %>%
  gather(length, n2, c(X7, X9, X10, X11, X12, X13, X14, X15, X16, X17, X18, X19, X20, X21, X22, X23, X24, X25,
                       X26, X27, X28, X29, X30, X31, X32, X33, X34, X35, X36, X37, X38, X39, X40, X41, X42, X43,
                       X44, X45, X46, X47, X48, X49, X50, X51, X52, X55, X62, X65, X80, X83, X90), na.rm = T)

# Test it is correct:
subset(df01, year == 2002 & week == 32 & day == 2 & Stations_namn == "Asphällan")
subset(dat01, year == 2002 & week == 32 & day == 2 & Stations_namn == "Asphällan")

# Create a new empty column for numeric "length"
dat01$length_group <- as.numeric(substring(dat01$length, 2))

# Our n2 column tells us how many fish are in that size-group. We want one row for each observation!
# So we need to repeat that observation by n2.
dat01 <- dat01[rep(seq(nrow(dat01)), dat01$n2),]
head(dat01, 50)

# Check its correct:
head(subset(dat01, n2 == 3), 50)

#** Combine all together =============================================================
catch_FM <- rbind(dat83, dat87, dat91, dat01)

catch_FM <- catch_FM %>% filter(year > 1982 & year < 2004)

#** Apply further filters ============================================================
# Remove disturbance
catch_FM <- catch_FM %>% filter(Störning == 0)

# They put nets many days in a row to get an overfishing affect.
# Catches decline after a few days. I don´t want this effect!
# Here we plot which days are fished (fill) over the week that has been fished (x-axis). "Overfishing" effect
# could happen if a year is fished a lot in many consecutive days. Its not clear here that is the case

# Full data
ggplot(catch_FM, aes(factor(week), fill = factor(day))) +
  facet_wrap(~year, scales = "free_y") +
  geom_bar() +
  scale_fill_brewer(palette="Set1", name = "Weekday") +
  theme_classic(base_size = 12) +
  ggtitle("Reference area")

# Filter autumn fishing?
sort(unique(catch_FM$week))

# This is not necessary because the data only contains autumn fishing! 

# Since BT data is filter for year > 1986 this has to be done here as well.
catch_FM <- catch_FM %>% filter(year > 1986)

# All nets have the same effort
unique(catch_FM$effort)
unique(catch_FM$week)

# In some years they have been fishing almost every week.
# To not have an overfishing-effect, we use only the first part of the
# fishing that season (first day and first week)

# We can probably see the overwfishin effect here in that the first boxes
# in each week tends to be larger (more fish), and it also declines in following weeks.
# We therefore need to select the first days of fishing.


# # Fick hjälp av Johan för detta!
# # Nu vill jag endast ha kvar data från första dagen som man har fiskat på varje vecka.
# # Det kommer inte vara dag 1 alla år utan vissa år har man börjat fiska på t.ex. en onsdag (dag 3).
# # Gör en ny vektor med 0 för varje rad.
# # Sen gör jag en loop som går igenom alla åren, veckor och stationer (Bara en sektion finns i data) i min data och ger första fiskedagen 1 och alla andra dagar 0.
# # Den nya vektorn har nu massor med 1 och 0.
# # Sen gör jag om den från 0 och 1 till TRUE och FALSE
# # Sen ny data frame = datsub1_dag1
# # sen plot
# filter_vector = rep(0,nrow(datsub1))
# 
# for(yr in min(datsub1$year):max(datsub1$year)){
#   for(we in min(datsub1$week):max(datsub1$week)){
#     for(sta in min(datsub1$Station):max(datsub1$Station)){
#       if(length(datsub1$day[datsub1$year == yr &
#                             datsub1$week == we &
#                             datsub1$Station == sta])>0){
#         
#         da = min(datsub1$day[datsub1$year == yr &
#                                datsub1$week == we &
#                                datsub1$Station==sta])
#         
#         filter_vector = filter_vector +
#           (datsub1$year == yr &
#              datsub1$week == we &
#              datsub1$day == da &
#              datsub1$Station == sta)
#       }
#     }
#   }
# }
# 
# filter_vector = as.logical(filter_vector)
# 
# datsub1_dag1 <- data.frame(datsub1[filter_vector,])
# 
# # If this worked, we should now only have one fishing occacsion (the first) per week 
# # Plot --> verkar fungera! :)
# p1 <- ggplot(datsub1_dag1, aes(factor(week), fill = factor(day))) +
#   facet_wrap(~year) +
#   geom_bar()
# 
# p2 <- ggplot(datsub1, aes(factor(week), fill = factor(day))) +
#   facet_wrap(~year) +
#   geom_bar()
# 
# # Compare!
# grid.arrange(p1, p2, nrow = 1)
# 
# # Now we need to also filter only the first week within a year
# datsub1_dag1 <- datsub1_dag1 %>% 
#   group_by(year) %>% 
#   filter(week == min(week)) %>% 
#   ungroup()
# 
# p3 <- ggplot(datsub1_dag1, aes(factor(week), fill = factor(day))) +
#   facet_wrap(~year) +
#   geom_bar()
# 
# grid.arrange(p1, p3, nrow = 1)
# 
# # All good!
# 
# # 3. SECTIONS AND STATIONS
# # Which sections are used?
# unique(datsub1_dag1$Sektion) # 1 and 2
# # which stations?
# unique(datsub1_dag1$Station) # 28, 29, 31, 32, 33, 34, 35, 36, 37, 38, 
# # Probably will no. of stations vary when I have removed disturbances.
# 
# # no. of stations each year
# ggplot(datsub1_dag1, aes(week, fill = factor(Station))) +
#   facet_wrap(~year, scales = "free_y") +
#   geom_bar() +
#   scale_fill_brewer(palette="Spectral", name = "Stations") +
#   theme_classic(base_size = 12) +
#   ggtitle("Reference area")
# 
# ggsave("C:/Users/malin/OneDrive - Sveriges Lantbruksuniversitet/Masterarbete/R_studio project/Figurer ANALYS/appendix/Stations_FM.png",
#        plot = last_plot(), scale = 1, width = 30, height = 20, units = "cm", dpi = 300)
# 
# # Actually, they do not vary much!
# 
# # The effort in my data is always 1 (=one night/net)
# ggplot(datsub1_dag1, aes(week, fill = factor(effort))) +
#   facet_wrap(~year) +
#   geom_bar()
# 
# # 4. Redskap?
# ggplot(datsub1_dag1, aes(week, fill = factor(Redskap))) +
#   facet_wrap(~year) +
#   geom_bar()
# 
# # only 9
# 
# # 5. Size-categories 
# unique(datsub1_dag1$Längdgr_std) # 2 = 2,5 cm intervall och 3 = 1 cm intervall
# sort(unique(datsub1_dag1$length_group))
# 
# # How many length groups are there in each std?
# # 2
# std2 <- data.frame(datsub1_dag1) %>%
#   filter(Längdgr_std == 2)
# sort(unique(std2$length_group))
# 
# # 3
# std3 <- data.frame(datsub1_dag1) %>%
#   filter(Längdgr_std == 3)
# sort(unique(std3$length_group))
# 
# # Plot std over year
# ggplot(datsub1_dag1, aes(week, fill = factor(Längdgr_std))) +
#   facet_wrap(~year) +
#   geom_bar() +
#   scale_fill_brewer(palette="Set1", name = "length group std") +
#   theme_classic(base_size = 12) +
#   ggtitle("Reference area")
# 
# # ggsave("C:/Users/malin/OneDrive - Sveriges Lantbruksuniversitet/Masterarbete/R_studio project/Figurer ANALYS/appendix/length_stdFM.png",
# #          plot = last_plot(), scale = 1, width = 16, height = 16, units = "cm", dpi = 300)
# 
# 
# datsub1_dag1 <- data.frame(datsub1_dag1)
# 
# 


# C. READ AND PROCESS BIOTEST DATA =================================================
# Need to set fileEncoding here, else error: "invalid multibyte string 18"
df <- read.csv("data/raw/Catch_data_BT09_140507.csv", sep = ";", fileEncoding = "latin1")

# Tidy data. Remove unnecessary columns
df <- df %>%
  filter(Art == "ABBO", Årtal < 2004) %>%
  dplyr::rename(year = Årtal, week = Vecka, day = Dag, effort = Ansträngning,
                species = Art, weight = Vikt, n = Antal) %>%
  select(-c(Vtn.stånd,  VindRiktn.I, VindSt.I, Vind_upp_rikn, VindSt.Upp, Ström_I_rikn,
            Ström_upp_rikn, Salthalt_I_yta, Salthalt_I_botten, Salthalt_upp_yta,
            Salthalt_upp_botten, Drift_i, Drift_u, Drift_dim, Siktdjup, Lufttryck_i,
            Lufttryck_upp, Sjuk_kontroll))

# Slightly cleaner but we have a problem: X1, X2... X9 are length.
# We would much rather prefer to have a column for each length, so that 1 row = observation.
# Luckily there is a super neat function in tidyr (part of the tidyverse) called "gather".
# What I whant to do is:

# 1. Give the new column names I want to create, e.g. "length" 
# 2. Specify which columns are moved and shuffled in these columns.
# These will be the current columns for unique lengths. So check the str() and the column number using names() of the data again 

# I will move columns 13:121 = the ones starting at x8 to x116

# Now the X-columns will be rows AND you will get a new column that takes the old "column"
# value and put that in the new column n2. DOUBLE CHECK!!
dat <- df %>%
  gather(length, n2, c(X8, X9, X10, X11, X12, X13, X14, X15, X16, X17, X18, X19, X20,
                       X21, X22, X23, X24, X25, X26, X27, X28, X29, X30, X31, X32, X33,
                       X34, X35, X36, X37, X38, X39, X40, X41, X42, X43, X44, X45, X46,
                       X47, X48, X49, X50, X51, X52, X53, X54, X55, X56, X57, X58, X59,
                       X60, X61, X62, X63, X64, X65, X66, X67, X68, X69, X70, X71, X72,
                       X73, X74, X75, X76, X77, X78, X79, X80, X81, X82, X84, X85, X86,
                       X87, X88, X89, X90, X91, X93, X94, X96, X97, X98, X99, X101, X103,
                       X104, X106, X109, X111, X114, X116), na.rm = T)

head(dat) # Now the data is in a long, tidy format. 
head(df)

# Test it is correct:
subset(df, year == 2003 & week == 19 & day == 4 & Stations_namn == "Ände pir")

subset(dat, year == 2003 & week == 19 & day == 4 & Stations_namn == "Ände pir")

# But "length" is not numeric but a character. Create a new empty column for numeric "length"
dat$length_group <- as.numeric(substring(dat$length, 2))

# Our n2 column tells us how many fish are in that size-group.
# We want one row for each observation! So we need to repeat that observation by n2.
dat <- dat[rep(seq(nrow(dat)), dat$n2),]
head(dat, 50)

# Check its correct:
head(subset(dat, n2 == 3), 50)

# Remove disturbance
# Disturbance kodes:
# 2: seals damage 
# 3: Strong algal growth on the gears
# 4: Clogging by drifting algae.
# 9: Other reason. (Damage by boat traffic, other human inference etc.)

# Plot disturbance
dat %>% 
  filter(Störning > 1) %>% 
  ggplot(., aes(x = Störning)) +
  geom_histogram(bins = 10) +
  scale_x_continuous(breaks = c(2, 3, 4, 9))

dat <- dat %>% filter(Störning == 0)

# Overfishing effect:
# They put nets many days in a row to get an overfishing affect.
# Catches decline after a few days. I don't want this effect!
# Here we plot which days are fished (fill) over the week that has been fished (x-axis). "Overfishing" effect
# could happen if a year is fished a lot in many consecutive days. Its not clear here that is the case

# Full data
ggplot(dat, aes(factor(week), fill = factor(day))) +
  facet_wrap(~year, scales = "free") +
  geom_bar() 

# Filter autumn
dat %>% 
  filter(week > 35) %>% 
  ggplot(., aes(factor(week), fill = factor(day))) +
  facet_wrap(~year, scales = "free_y") +
  geom_bar()

# Ok, 1984 is heavily fished all year actually.
# We'll remove it since it's very different from the rest. When that is removed, we
# can use malins for loop to select fishing days after .v40 (or any other week).
# The other two years with fishing before v.40 are 1996 & 2003, but that's w. 30 so we
# don't belive it has a big effect.

dat <- dat %>% filter(year > 1986)

# In some years they have been fishing almost every week. To not have an
# overfishing-effect, we use only the first part of the fishing that season (first
# day and first week)

# Filter autumn fishing
# By filtering week > 40 & week < 49 I know that I include al Oktober fishing
dat_oct <- dat %>%
  filter(week > 40 & week < 49) # OCTOBER

sort(unique(dat_oct$year))
sort(unique(dat_oct$week)) 
# Since this gives me v. 41 42 43 44 45 46 47 
# year -99 -00 is lost due to disturbance

# How does the filtered data set look?
ggplot(dat_oct, aes(x = factor(week), fill = factor(day))) +
  facet_wrap(~year) +
  geom_bar() +
  scale_fill_brewer(palette="Set1", name = "Weekday") +
  theme_classic(base_size = 12) +
  ggtitle("Biotest lake")

# Rename data
catch_BT <- dat_oct

# We can probably see the overwfishin effect here in that the first boxes in each week
# tends to be larger (more fish), and it also declines in following weeks.
# We therefore  need to select the first days of fishing.

# # How do I do that??
# 
# # Fick hjÃ¤lp av Johan fÃ¶r detta! 
# # Nu vill jag endast ha kvar data frÃ¥n fÃ¶rsta dagen som man har fiskat pÃ¥ varje vecka.
# # Det kommer inte vara dag 1 alla Ã¥r utan vissa Ã¥r har man bÃ¶rjat fiska pÃ¥ t.ex. en onsdag (dag 3).
# # GÃ¶r en ny vektor med 0 fÃ¶r varje rad.
# # Sen gÃ¶r jag en loop som gÃ¥r igenom alla Ã¥ren, veckor och stationer (Bara en sektion finns i data) i min data och ger fÃ¶rsta fiskedagen 1 och alla andra dagar 0.
# # Den nya vektorn har nu massor med 1 och 0.
# # Sen gÃ¶r jag om den frÃ¥n 0 och 1 till TRUE och FALSE
# # Sen ny data frame = datsubOKT_dag1
# # sen plot
# filter_vector = rep(0,nrow(dat_oct))
# 
# for(yr in min(dat_oct$year):max(dat_oct$year)){
#   for(we in min(dat_oct$week):max(dat_oct$week)){
#       for(sta in min(dat_oct$Station):max(dat_oct$Station)){
#         if(length(dat_oct$day[dat_oct$year==yr & 
#                                 dat_oct$week==we & 
#                                 dat_oct$Station==sta])>0){
#           
#           da = min(dat_oct$day[dat_oct$year==yr & 
#                                  dat_oct$week==we &
#                                  dat_oct$Station==sta])
# 
#           filter_vector = filter_vector +
#             (dat_oct$year==yr &
#              dat_oct$week==we & 
#              dat_oct$day==da &
#              dat_oct$Station==sta)
#         }
#       }
#     }
# }
# 
# filter_vector = as.logical(filter_vector)
# 
# dat_oct_dag1 <- dat_oct[filter_vector,]
# 
# # If this worked, we should now only have one fishing occasion (the first) per week 
# # Plot --> verkar fungera! :)
# ggplot(dat_oct_dag1, aes(x = factor(week), fill = factor(day))) +
#   facet_wrap(~year) +
#   geom_bar()
# 
# # OBS!! NÃ¥got skumt hÃ¤nder med -96 och -98. Men det hÃ¤r kanske inte gÃ¶r nÃ¥got
# # eftersom jag kommer ha emd endast fÃ¶rsta dagen fÃ¶rsta veckan.
# 
# ggplot(datsubOKT, aes(x = factor(week), fill = factor(day))) +
#   facet_wrap(~year) +
#   geom_bar()
# 
# # Now we need to also filter only the first week within a year
# datsubOKT_dag1 <- datsubOKT_dag1 %>% 
#   group_by(year) %>% 
#   filter(week == min(week)) %>% 
#   ungroup()
# 
# p3 <- ggplot(datsubOKT_dag1, aes(x = factor(week), fill = factor(day))) +
#   facet_wrap(~year) +
#   geom_bar()
# 
# grid.arrange(p1, p3, nrow = 1)
# 
# # All good!
# 
# # Temperature this period?
# # Calcualte mean temp. for each year
# datsubOKT_temp <- datsubOKT_dag1 %>%
#   group_by(year) %>%        
#   summarize(mean_temp  = mean(YtTmp.I)) %>% 
#   drop_na() %>% 
#   ungroup
# 
# # Plot mean temp. over year
# ggplot(datsubOKT_temp, aes(x = year, y = mean_temp)) + 
#   geom_point(size = 3) +
#   geom_line()
# 
# unique(datsubOKT_dag1$YtTmp.I)
# unique(datsubOKT_dag1$YtTmp.Upp)
# 
# ggplot(datsubOKT_dag1, aes(x = year, y = YtTmp.I)) +
#   geom_point(size = 1) +
#   geom_line(data=datsubOKT_temp, aes(x = year, y = mean_temp))
# 
# 
# # Save temperature data
# write.csv(datsubOKT_temp,"temp_BT.csv", row.names = FALSE)
# 
# # 3. SECTIONS AND STATIONS
# # Which sections are used?
# unique(datsubOKT_dag1$Sektion) # only 1
# 
# # Which stations?
# unique(datsubOKT_dag1$Station)
# # Probably will no. of stations vary when I have removed disturbances. 
# 
# # no. of stations each year
# ggplot(datsubOKT_dag1, aes(week, fill = factor(Station))) +
#   facet_wrap(~year, scales = "free_y") +
#   geom_bar() +
#   scale_fill_brewer(palette="Spectral", name = "Stations") +
#   theme_classic(base_size = 12) +
#   ggtitle("Biotest lake")
# 
# ggsave("C:/Users/malin/OneDrive - Sveriges Lantbruksuniversitet/Masterarbete/R_studio project/Figurer ANALYS/appendix/Stations_BT.png",
#        plot = last_plot(), scale = 1, width = 30, height = 20, units = "cm", dpi = 300)
# 
# # Actually, they do not vary much!
# # The effort in my data is always 1 (=one night/net)
# ggplot(datsubOKT_dag1, aes(week, fill = factor(effort))) +
#   facet_wrap(~year) +
#   geom_bar()
# 
# # 4. Redskap?
# ggplot(datsubOKT_dag1, aes(week, fill = factor(Redskap))) +
#   facet_wrap(~year) +
#   geom_bar()
# 
# # 5. Size-categories 
# unique(datsubOKT_dag1$length_group)
# unique(datsubOKT_dag1$LÃ¤ngdgr_std) # 2 = 2,5 cm intervall och 3 = 1 cm intervall
# 
# ggplot(datsubOKT_dag1, aes(week, fill = factor(LÃ¤ngdgr_std))) +
#   facet_wrap(~year) +
#   geom_bar() +
#   scale_fill_brewer(palette="Set1", name = "length group std") +
#   theme_classic(base_size = 12) +
#   ggtitle("Biotest lake")
# 
# ggsave("C:/Users/malin/OneDrive - Sveriges Lantbruksuniversitet/Masterarbete/R_studio project/Figurer ANALYS/appendix/length_stdBT.png",
#        plot = last_plot(), scale = 1, width = 16, height = 16, units = "cm", dpi = 300)
# 
# sort(unique(datsubOKT_dag1$year))
# 
# 
# 


# D. Standardize length groups =====================================================
# Clean data to make more readable
catch_BT <- catch_BT %>%
  select(Area, Sektion, Station, year, week, day, length_group, Station,
         Längdgr_std, length_group) %>% 
  as.data.frame()

catch_FM <- catch_FM %>%
  select(Area, Sektion, Station, year, week, day, length_group, Station,
         Längdgr_std, length_group) %>% 
  as.data.frame()

catch <- rbind(catch_BT, catch_FM)

# First split up the data i two separate dataframes, one with std 2 and one with std 3.
# Call them dat_std2 och dat_std3
dat_std2 <- catch %>%
  filter(Längdgr_std == 2)

dat_std3 <- catch %>%
  filter(Längdgr_std == 3)

# Convert std 3 to std 2
sort(unique(dat_std3$length_group))
# 7 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35
# 36 37 38 39 40 41 42 43 44 45 46 47

# Insert new column 
dat_std3$length_group2 <- NA
dat_std3$length_group2[dat_std3$length_group==7] <- 6
dat_std3$length_group2[dat_std3$length_group==9] <- 9
dat_std3$length_group2[dat_std3$length_group==10] <- 9
dat_std3$length_group2[dat_std3$length_group==11] <- 11
dat_std3$length_group2[dat_std3$length_group==12] <- 11
dat_std3$length_group2[dat_std3$length_group==13] <- 11
dat_std3$length_group2[dat_std3$length_group==14] <- 14
dat_std3$length_group2[dat_std3$length_group==15] <- 14
dat_std3$length_group2[dat_std3$length_group==16] <- 16
dat_std3$length_group2[dat_std3$length_group==17] <- 16
dat_std3$length_group2[dat_std3$length_group==18] <- 16
dat_std3$length_group2[dat_std3$length_group==19] <- 19
dat_std3$length_group2[dat_std3$length_group==20] <- 19
dat_std3$length_group2[dat_std3$length_group==21] <- 21
dat_std3$length_group2[dat_std3$length_group==22] <- 21
dat_std3$length_group2[dat_std3$length_group==23] <- 21
dat_std3$length_group2[dat_std3$length_group==24] <- 24
dat_std3$length_group2[dat_std3$length_group==25] <- 24
dat_std3$length_group2[dat_std3$length_group==26] <- 26
dat_std3$length_group2[dat_std3$length_group==27] <- 26
dat_std3$length_group2[dat_std3$length_group==28] <- 26
dat_std3$length_group2[dat_std3$length_group==29] <- 29
dat_std3$length_group2[dat_std3$length_group==30] <- 29
dat_std3$length_group2[dat_std3$length_group==31] <- 31
dat_std3$length_group2[dat_std3$length_group==32] <- 31
dat_std3$length_group2[dat_std3$length_group==33] <- 31
dat_std3$length_group2[dat_std3$length_group==34] <- 34
dat_std3$length_group2[dat_std3$length_group==35] <- 34
dat_std3$length_group2[dat_std3$length_group==36] <- 36
dat_std3$length_group2[dat_std3$length_group==37] <- 36
dat_std3$length_group2[dat_std3$length_group==38] <- 36
dat_std3$length_group2[dat_std3$length_group==39] <- 39
dat_std3$length_group2[dat_std3$length_group==40] <- 39
dat_std3$length_group2[dat_std3$length_group==41] <- 41
dat_std3$length_group2[dat_std3$length_group==42] <- 41
dat_std3$length_group2[dat_std3$length_group==43] <- 41
dat_std3$length_group2[dat_std3$length_group==44] <- 44
dat_std3$length_group2[dat_std3$length_group==45] <- 44
dat_std3$length_group2[dat_std3$length_group==46] <- 46
dat_std3$length_group2[dat_std3$length_group==47] <- 46
dat_std3$length_group2[dat_std3$length_group==48] <- 46

# Now compare actual values - Looks ok!
ggplot(dat_std3, aes(factor(length_group), factor(length_group2))) +
  geom_point(size = 2) 

# Remove length_group so that only std 2 is included
dat_std3 <- dat_std3 %>%
  select(-c(length_group))

# Rename length_group2 to length_group
dat_std3 <- dat_std3 %>%
  dplyr::rename(length_group = length_group2)

# Merge data frames
catch_full <- rbind(dat_std2, dat_std3)

# In order to get even length classes (with respect to the code) I gave a new code to std 2
# which represents the starting length in each interval, as opposed to the original
# which was the integer in the middle of 2.5 cm classes...
# E.g. if the interval is 37.6 - 40, the new code becomes 37.6.

sort(unique(catch_full$length_group))
# 6  9 11 14 16 19 21 24 26 29 31 34 36 39 41 44 46

# Insert new column 
catch_full$new_length_group <- NA
catch_full$new_length_group[catch_full$length_group==6] <- 5.1
catch_full$new_length_group[catch_full$length_group==9] <- 7.6
catch_full$new_length_group[catch_full$length_group==11] <- 10.1
catch_full$new_length_group[catch_full$length_group==14] <- 12.6
catch_full$new_length_group[catch_full$length_group==16] <- 15.1
catch_full$new_length_group[catch_full$length_group==19] <- 17.6
catch_full$new_length_group[catch_full$length_group==21] <- 20.1
catch_full$new_length_group[catch_full$length_group==24] <- 22.6
catch_full$new_length_group[catch_full$length_group==26] <- 25.1
catch_full$new_length_group[catch_full$length_group==29] <- 27.6
catch_full$new_length_group[catch_full$length_group==31] <- 30.1
catch_full$new_length_group[catch_full$length_group==34] <- 32.6
catch_full$new_length_group[catch_full$length_group==36] <- 35.1
catch_full$new_length_group[catch_full$length_group==39] <- 37.6
catch_full$new_length_group[catch_full$length_group==41] <- 40.1
catch_full$new_length_group[catch_full$length_group==44] <- 42.6
catch_full$new_length_group[catch_full$length_group==46] <- 45.1

# Compare values
ggplot(catch_full, aes(length_group, new_length_group)) +
  geom_point(size = 2) +
  geom_abline(slope = 1, intercept = 0, color = "red")

# Remove length_group so that only new_length_group is included
catch_full <- catch_full %>%
  select(-c(length_group))

# Rename length_group2 to length_group
catch_full <- catch_full %>%
  dplyr::rename(length_group = new_length_group)

# Insert netID
catch_full$netID <- paste(catch_full$Area, catch_full$year, catch_full$Station, sep = ".")

# Save data frame 
write.csv(catch_full, "data/catch_FM_BT_1987-2003.csv", row.names = FALSE)

