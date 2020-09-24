####
## 2019.12.06 - permutation test of Z from different runs of the semi-random ALK
####


#---- Clear the workspace
rm(list=ls(all=TRUE))

#---- Load libraries (install first if needed)
#install.packages("devtools")
#install.packages("plyr")
#install.packages("dplyr")
#install.packages("ggplot2")
#install.packages("grid")
#install.packages("gridExtra")
#install.packages("mosaic")

# For ALK 
#install.packages("FSA")
#install.packages("FSAdata")
#install.packages("plotrix")

library(devtools)
library(dplyr)
library(plyr)
library(ggplot2)
library(grid)
library(gridExtra)
library(FSA)
library(FSAdata) # for datafile
library(plotrix) # for histStack
library(mosaic)
library(RColorBrewer)
library(tidyr)


#---- Read in the data
df.age <- read.csv("Growth_data_BT_FM_1970-2004.csv", sep = ";")
df.lenBT <- read.csv("Sizespectra_BT.csv", sep = ",") 
df.lenFM <- read.csv("Sizespectra_FM.csv", sep = ",")

sort(unique(df.lenFM$year))
sort(unique(df.lenBT$year))

sort(unique(df.age$length))
sort(unique(df.lenBT$length_group))

# ========== ALK ==========
# One ALK for each area

# Restructure data'
# Split df.age by area
df.ageBT <- df.age %>% 
  filter(area == "BT", catch_year > 1986)

summary(df.ageBT)
sort(unique(df.ageBT$catch_year))

df.ageFM <- df.age %>% 
  filter(area == "FM", catch_year > 1986)

summary(df.ageFM)
sort(unique(df.ageFM$catch_year))

# Convert lengths in df.lenBT/FM from cm to mm
df.lenBT <- df.lenBT %>%
  dplyr::mutate(length_mm = length_group*10)

df.lenFM <- df.lenFM %>%
  dplyr::mutate(length_mm = length_group*10)

sort(unique(df.lenBT$length_mm))

# Rename length_group to lengthGr_cm
df.lenBT <- rename(df.lenBT, replace = c("length_group" = "length_group_cm"))
df.lenFM <- rename(df.lenFM, replace = c("length_group" = "length_group_cm"))

# Construction of length categorys 
# sort(unique(df.ageBT$length))
# sort(unique(df.ageFM$length))

# The starting category for 10-mm length categories is determined by finding the
# minimum length in the age sample with: 
Summarize(~length, data=df.ageBT, digits=1) # 43.0
Summarize(~length, data=df.ageFM, digits=1) # 40.0

#---- calculate XXX 'nloops' number of times
nloops = 1000
ZBTvector = numeric(nloops)
ZFMvector = numeric(nloops)

for(i in 1:nloops){
  
# and then starting the categories with the even-number 25-mm interval just below this value.
# Construction of the length category variable is completed and the first six rows are viewed
# with:
rb.ageBT <- lencat(~length, data=df.ageBT, startcat = 42, w = 25)
head(rb.ageBT)
  
rb.ageFM <- lencat(~length, data=df.ageFM, startcat = 39, w = 25)
head(rb.ageFM)
  
  
# This creates a variable in the age sample that identifies the length category to which
# each fish belongs.
sort(unique(rb.ageBT$LCat))
sort(unique(rb.ageFM$LCat))
  
# # Compare values
# ggplot(rb.ageBT1, aes(LCat, length)) +
#   geom_point(size = 2) +
#   geom_abline(slope = 1, intercept = 0, color = "red")
# 
# ggplot(rb.ageFM1, aes(LCat, length)) +
#   geom_point(size = 2) +
#   geom_abline(slope = 1, intercept = 0, color = "red")
  
  
# The summary contingency table and the row-proportion table (i.e., the age-length key)
# are constructed with:
rb.raw.BT <- with(rb.ageBT, table(LCat,age))
rb.raw.FM <- with(rb.ageFM, table(LCat,age))
  
rb.key.BT <- prop.table(rb.raw.BT, margin=1)
rb.key.FM <- prop.table(rb.raw.FM, margin=1)
  
round(rb.key.BT, 2) # rounded for display purposes only
round(rb.key.FM, 2) 
  
# ageKey() — DEPRECATED (will be removed by v1.0.0).  See alkIndivAge().
rb.lenBT <- alkIndivAge(rb.key.BT, ~length_mm, data = df.lenBT)
head(rb.lenBT)
  
rb.lenFM <- alkIndivAge(rb.key.FM, ~length_mm, data = df.lenFM)
head(rb.lenBT)
  
# Summary
(rb.sum.BT <- Summarize(length_mm~age, data=rb.lenBT, digits = 2))
(rb.sum.FM <- Summarize(length_mm~age, data=rb.lenFM, digits = 2))
  
# Merge rb.lenBT1 and rb.lenFM1
rb.len <- rbind(rb.lenBT, rb.lenFM)
  
sort(unique(rb.len$year))
  
# ========== Total mortality ==========
# fishR Vignette - Catch Curve Estimates of Mortality. By: Dr. Derek Ogle, Northland College December 16, 2013
  
df <- rb.len
  
# Add catch column
  df_sub <- df %>%
    dplyr::group_by(age, Area) %>%
    dplyr::summarise(count = n()) %>%
    dplyr::ungroup()
  
  # count --> catch (ct)
  df_sub <- df_sub %>%
    dplyr::rename(ct = count)
  
  # and the natural log of the catch data
  df_sub$logct <- log(df_sub$ct)
  
  # ========== Basic regression method for estimating Z ==========
  dat <- df_sub %>%
    filter(age > 2)
  
  # lm (ANCOVA)
  lm <- lm(logct~age * Area, data = dat)
  
  coef(lm)
  confint(lm)
  summary(lm)
  
  
  ZBTvector[i] <- summary(lm)$coefficient[2,1]
  ZFMvector[i] <- summary(lm)$coefficient[2,1] + summary(lm)$coefficient[4,1]
  print(i)
}

hist(ZBTvector)
hist(ZFMvector)

Zdat_BT <- data.frame(Z = ZBTvector,
                      Area = "BT")


Zdat_FM <- data.frame(Z = ZFMvector,
                      Area = "FM")

Zdat <- rbind(Zdat_BT, Zdat_FM)

unique(Zdat$Area)


# Confidence interval 
summary(lm)
confint(lm)

ZBT = summary(lm)$coefficient[2,1]
ZFM = summary(lm)$coefficient[2,1] + summary(lm)$coefficient[4,1]

# confint.merMod(lmer1ALL)
LCIBT <- confint(lm)[2,1] # LCI BT
UCIBT <- confint(lm)[2,2] # UCI BT
LCIFM <- confint(lm)[2,1] + confint(lm)[4,1] # LCI FM
UCIFM <- confint(lm)[2,2] + confint(lm)[4,2] # UCI FM

# New data frame with estimate and CI
Zconf <- data.frame(Area = c("BT", "FM"),
                   Z = c(ZBT, ZFM),
                   LCI = c(LCIBT, LCIFM),
                   UCI = c(UCIBT, UCIFM))

# Plot Z estimate and confidence intervals
ggplot(Zconf, aes(x = Area, y = -Z, color = factor(Area))) +
  geom_point(size=3) +
  geom_errorbar(aes(x = Area, ymin = -LCI, ymax = -UCI), width = 0.1) +
  theme_classic(base_size = 18) +
  xlab("") +
  scale_color_brewer(palette="Set1") +
  theme(legend.position = "none")

#ggsave("C:/Users/malin/OneDrive - Sveriges Lantbruksuniversitet/Masterarbete/R_studio project/Figurer ANALYS/Mortality/Z_CI.png", plot = last_plot(), scale = 1, width = 16, height = 16, units = "cm", dpi = 300)


# pred
dat$pred <- predict(lm)

# plot pred. lines
ggplot(dat, aes(x = age, y = logct)) +
  geom_point(data = dat, aes(x = age, y = logct), size = 2, shape = 15) +
  geom_line(data = dat, aes(x = age, y = pred), size = 0.8) +
  theme_bw() +
  facet_wrap(~Area)

# split dataframes for plot
dfBT <- dat %>%
  filter(Area == "BT")

dfFM <- dat %>%
  filter(Area == "FM")

# lm 
lmBT <- lm(logct~age, data = dfBT)
lmFM <- lm(logct~age, data = dfFM) 

# pred
dfBT$pred1 <- predict(lmBT)
dfFM$pred1 <- predict(lmFM)

# plot pred. lines
ggplot(dat, aes(x = age, y = exp(logct), color = Area)) +
  geom_point(data = dfBT, aes(x = age, y = exp(logct)), size = 2, shape = 15) +
  geom_line(data = dfBT, aes(x = age, y = exp(pred1)), size = 0.8) +
  geom_point(data = dfFM, aes(x = age, y = exp(logct)), size = 2, shape =  16) +
  geom_line(data = dfFM, aes(x = age, y = exp(pred1)), size = 0.8) +
  theme_classic(base_size = 18) +
  scale_color_brewer(palette="Set1") +
  scale_y_continuous(trans='log10') +
  ylab("CPUE (n)") +
  xlab("Age (years)") +
  scale_x_continuous(breaks = 3:9) +
  theme(legend.position = "none")

ggsave("C:/Users/malin/OneDrive - Sveriges Lantbruksuniversitet/Masterarbete/R_studio project/Figurer ANALYS/Mortality/Z.png",
       plot = last_plot(), scale = 1, width = 30, height = 16, units = "cm", dpi = 300)

# ============ Weighted Catch-Curve Regression (Recommended) ============
# The function "catchCurve()" does not accept an interaction of area and there is
# no applicable method for 'predict' applied to an object of class "catchCurve"

WlmBT <- catchCurve(ct ~ age, data=dfBT, ages2use=3:9, use.weights=TRUE)

WlmFM <- catchCurve(ct ~ age, data=dfFM, ages2use=3:9, use.weights=TRUE)



# Confidence interval 
summary(WlmBT)
confint(WlmBT)
summary(WlmFM)
confint(WlmFM)

# ZWlmBT = summary(WlmBT)$coefficient[2,1]
# ZFM = summary(lm)$coefficient[2,1] + summary(lm)$coefficient[4,1]
# 
# # confint.merMod(lmer1ALL)
# LCIBT <- confint(lm)[2,1] # LCI BT
# UCIBT <- confint(lm)[2,2] # UCI BT
# LCIFM <- confint(lm)[2,1] + confint(lm)[4,1] # LCI FM
# UCIFM <- confint(lm)[2,2] + confint(lm)[4,2] # UCI FM
# 
# # New data frame with estimate and CI
# Zconf <- data.frame(Area = c("BT", "FM"),
#                     Z = c(ZBT, ZFM),
#                     LCI = c(LCIBT, LCIFM),
#                     UCI = c(UCIBT, UCIFM))
# 
# # Plot b estimate and confidence intervals
# ggplot(Zconf, aes(Area, Z, color = factor(Area))) +
#   geom_point(size=3) +
#   geom_errorbar(aes(Area, ymin = LCI, ymax = UCI), width = 0.1) +
#   theme_classic(base_size = 18) +
#   xlab("") +
#   scale_color_brewer(palette="Set1") +
#   theme(legend.position = "none")


# ========== Mean age ==========
dfmeanageBT <- df %>%
  filter(Area == "BT")

dfmeanageFM <- df %>%
  filter(Area == "FM")

mean(dfmeanageBT$age)
mean(dfmeanageFM$age)

# ========== Permutation test Z estimates (basic regression method) ==========
# Dataframe with 1000 estimates of Z in each area
#Zdat <- read.csv("Z from different runs of the semi-random ALK_one for each area_1987-2003.csv", sep = ";", header = FALSE)
head(Zdat)
summary(Zdat)
str(Zdat)

#colnames(Zdat) = c("area", "Z")

# Calculate the observed difference in means
mosaic::mean(Z ~ Area, data = Zdat)
# Store diff. between observed
observed <-
  mean(Z ~ Area, data = Zdat) %>% 
  diff()

# diff() for two items gives the seckond minus the first (FM - BT)
observed

# To simulate4 a single trial, we need to shuffle the treatment labels (BT and FM)
mean(Z ~ shuffle(Area), data = Zdat)
mean(Z ~ shuffle(Area), data = Zdat) %>% 
  diff()

# Create a randomization distibution
Z_null <- do(5000) * mean(Z ~ shuffle(Area), data = Zdat) %>% 
  diff()

head(Z_null)

# plotting the randomization distibution
ggplot(data = Z_null) +
  geom_histogram(mapping = aes(x = FM), bins = 30) +
  xlab("mean difference")

# Superimpose a line indicating the observation
ggplot(data = Z_null) +
  geom_histogram(mapping = aes(x = FM), bins = 80) +
  xlab("mean difference") +
  geom_vline(xintercept = observed, linetype = 2, color = "dodgerblue4", size = 1) +
  scale_color_brewer(palette="Set1") +
  theme_classic(base_size = 18) +
  ylab("Count") +
  xlab("Mean difference") +
  ggtitle("A")

# ggsave("C:/Users/malin/OneDrive - Sveriges Lantbruksuniversitet/Masterarbete/R_studio project/Figurer ANALYS/Mortality/Randomization_distribution_with_obs_line.png",
#        plot = last_plot(), scale = 1, width = 16, height = 16, units = "cm", dpi = 300)


display.brewer.all()


# Calculate the proportion of simulated diff. in means as or more extreme than the observed
prop(~ FM <= observed, data = Z_null)

# Results indicates that my observed difference in Z is probably not random.
# Z in FM is higher (more negative) than in BT

# Boxplot of Z estimates
ggplot(Zdat, aes(x = Area, y = Z, color = factor(Area))) +
  geom_jitter(position=position_jitter(width=0.2, height=NULL), alpha=0.1) +
  geom_boxplot(width = 0.4) +
  ylab("Z") +
  theme_classic(base_size = 18) +
  scale_color_brewer(palette="Set1") +
  theme(legend.position = "none") +
  ggtitle("B")

ggplot(Zdat, aes(x = Area, y = -Z, color = factor(Area))) +
  geom_jitter(position=position_jitter(width=0.2, height=NULL), alpha=0.1) +
  geom_boxplot(width = 0.4) +
  ylab("Z") +
  theme_classic(base_size = 18) +
  scale_color_brewer(palette="Set1") +
  theme(legend.position = "none") +
  ggtitle("B")

# ggsave("C:/Users/malin/OneDrive - Sveriges Lantbruksuniversitet/Masterarbete/R_studio project/Figurer ANALYS/Mortality/Est_of_Z.png",
#         plot = last_plot(), scale = 1, width = 16, height = 16, units = "cm", dpi = 300)
