####
## 2019.10.29 - Total Mortality
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
#install.packages("FSA")

library(devtools)
library(dplyr)
library(plyr)
library(ggplot2)
library(grid)
library(gridExtra)
library(FSA)
library(RColorBrewer)


# Color palette RColorBrewer
par(mfrow=c(1, 1))
display.brewer.all()

# ========== Read in the data ==========
df <- read.csv("Length-at-age.csv", sep = ",") # Created in Age-length key script

min(df$length_mm)
min(df$age)

# fishR Vignette - Catch Curve Estimates of Mortality. By: Dr. Derek Ogle, Northland College December 16, 2013

# add catch column
df_sub <- df %>%
  dplyr::group_by(age, Area) %>%
  dplyr::summarise(count = n()) %>%
  dplyr::ungroup()

# count --> catch (ct)
df_sub <- df_sub %>%
  dplyr::rename(ct = count)

# and the natural log of the catch data
df_sub$logct <- log(df_sub$ct)

# Plot catch-curve
ggplot(df_sub, aes(x = age, y = ct)) +
  geom_point(size = 2) +
  facet_wrap(~Area) +
  scale_x_continuous(limits = c(1, 9), breaks = c(1, 2, 3, 4, 5, 6, 7, 8, 9)) +
  theme_bw()

ggplot(df_sub, aes(x = age, y = ct)) +
  geom_line(size = 1) +
  facet_wrap(~Area) +
  scale_x_continuous(limits = c(1, 9), breaks = c(1, 2, 3, 4, 5, 6, 7, 8, 9)) +
  theme_bw() +
  xlab("age") +
  ylab("Catch")
  

# and log
ggplot(df_sub, aes(x = age, y = logct)) +
  geom_point(size = 2) +
  scale_x_continuous(limits = c(1, 9), breaks = c(1, 2, 3, 4, 5, 6, 7, 8, 9)) +
  facet_wrap(~Area) +
  theme_bw()

ggplot(df_sub, aes(x = age, y = logct, color = Area)) +
  geom_line(size = 1.3) +
  scale_x_continuous(limits = c(1, 9), breaks = c(1, 2, 3, 4, 5, 6, 7, 8, 9)) +
  facet_wrap(~Area) +
  xlab("Age") +
  ylab("log(CPUE)") +
  theme_classic(base_size = 14) +
  scale_color_brewer(palette="Set1") +
  theme(legend.position = "none")

ggsave("C:/Users/malin/OneDrive - Sveriges Lantbruksuniversitet/Masterarbete/R_studio project/Figurer ANALYS/Mortality/Catch_curve_plot.png", plot = last_plot(), scale = 1, width = 16, height = 16, units = "cm", dpi = 300)
  

# ========== Basic regression method for estimating Z ==========

# From the initial catch curve plot it was determined that the descending limb consisted of
# ages 3 through 9 (BT) and 2 - 9 (FM).
# I need to create a new data frame which contains only the information for ages (3-9). 
dat <- df_sub %>%
  filter(age > 2)


# lm (ANCOVA)
# Testing the effect of area (categorical factor) on catches (dependent, y-var)
# while controlling for the effect of age (continous co-variable, x-var)
# Regression lines are compared by studying the interaction of area with the continuous
# independent variable (age).
lm <- lm(logct~age * Area, data = dat)
coef(lm)
confint(lm)

summary(lm)

# Significant effect of age and area, but no significant interaction.
# These results suggest that the slope of the regression
# is similar for both areas, meaning no difference in Z.

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
ggplot(dat, aes(x = age, y = logct, color = Area)) +
  geom_point(data = dfBT, aes(x = age, y = logct), size = 2, shape = 15) +
  geom_line(data = dfBT, aes(x = age, y = pred1), size = 0.8) +
  geom_point(data = dfFM, aes(x = age, y = logct), size = 2, shape =  16) +
  geom_line(data = dfFM, aes(x = age, y = pred1), size = 0.8) +
  theme_classic(base_size = 18) +
  scale_color_brewer(palette="Set1") +
  ylab("Log(Catch)") +
  xlab("Age") +
  theme(legend.position = "none")

# Estimate of Z with CI

# Z estimate
summary(lm)$coefficients
ZBT <- summary(lm)$coefficients[2,1]
ZFM <- summary(lm)$coefficients[2,1] + summary(lm)$coefficients[4,1]

# CI
confint(lm)
LCIBT <- confint(lm)[2,1]
UCIBT <- confint(lm)[2,2]
LCIFM <- confint(lm)[2,1] + confint(lm)[4,1]
UCIFM <- confint(lm)[2,2] + confint(lm)[4,2]

# New data frame with estimate and CI
Zdat <- data.frame(area = c("BT", "FM"),
                   Z = c(ZBT, ZFM),
                   LCI = c(LCIBT, LCIFM),
                   UCI = c(UCIBT, UCIFM))

# Plot Z estimate and confidence intervals
ggplot(Zdat, aes(area, Z, color = factor(area))) +
  geom_point(size=3) +
  geom_errorbar(aes(area, ymin = LCI, ymax = UCI), width = 0.1) +
  theme_classic(base_size = 18) +
  xlab("") +
  scale_color_brewer(palette="Set1") +
  theme(legend.position = "none")

# ============ Weighted Catch-Curve Regression (Recommended) ============
# The function "catchCurve()" does not accept an interaction of area and there is
# no applicable method for 'predict' applied to an object of class "catchCurve"

WlmBT <- catchCurve(ct ~ age, data=dfBT, ages2use=3:9, use.weights=TRUE)
summary(WlmBT)
confint(WlmBT)

WlmFM <- catchCurve(ct ~ age, data=dfFM, ages2use=3:9, use.weights=TRUE)
summary(WlmFM)
confint(WlmFM)

# ========== Assumption check ==========
lm <- lm(logct~age * Area, data = dat)
par(mfrow = c(2,2))
plot(lm)

# ================================== Average age ==================================

mean(dfBT$age) # 6.5
max(dfBT$age) # 9
mean(dfFM$age) # 6.5
max(dfFM$age) # 9

# t-test av. age
lmage <- lm(dat$age ~ dat$Area)
summary(lmage)


# ============== Av. age in certain percentiles ===================
# Kolla fler percentiler, från 90% och uppåt och se när det blir en skilland
# Gör till en data frame med en kolumn med percentiler (alltså den som varje tillhör)
# gör en ggplot som jag delar upp på percentiler

# Create dataframes with only age
# Split dataframes and group by
ageBT <- df %>%
  dplyr::filter(Area == "BT")

ageFM <- df %>%
  filter(Area == "FM")
  
ageBT = ageBT$age
ageFM = ageFM$age

# BT
qBT <- quantile(ageBT, probs = c(0.90, 0.92, 0.94, 0.96, 0.98))
qBT[1]
qBT[2]

qBT90 <- ageBT %>%
  dplyr::filter(V1 >= qBT[1]) %>% 
  dplyr::select(V1)

qBT90$Q <- c(rep(90))
qBT90$area <- c(rep("BT"))

qBT92 <- dfBT12 %>%
  filter(length_group >= qBT[2]) %>% 
  select(length_group)

qBT92$Q <- c(rep(92))
qBT92$area <- c(rep("BT"))

qBT94 <- dfBT12 %>%
  filter(length_group >= qBT[3]) %>% 
  select(length_group)

qBT94$Q <- c(rep(94))
qBT94$area <- c(rep("BT"))

qBT96 <- dfBT12 %>%
  filter(length_group >= qBT[4]) %>% 
  select(length_group)

qBT96$Q <- c(rep(96))
qBT96$area <- c(rep("BT"))

qBT98 <- dfBT12 %>%
  filter(length_group >= qBT[5]) %>% 
  select(length_group)

qBT98$Q <- c(rep(98))
qBT98$area <- c(rep("BT"))





dfqBT = df_perc_BT$age
qBT <- quantile(dfqBT, probs = c(0.90, 0.92, 0.94, 0.96, 0.98))
qBT


# BT 90%
qBT90 <- df_perc_BT %>%
  filter(age >= qBT[1]) %>% 
  select(age)

qBT90$Q <- c(rep(90))
qBT90$area <- c(rep("BT"))

# BT 92%
qBT92 <- df_perc_BT %>%
  filter(age >= qBT[2]) %>% 
  select(age)

qBT92$Q <- c(rep(92))
qBT92$area <- c(rep("BT"))

# BT 94%
qBT94 <- df_perc_BT %>%
  filter(age >= qBT[3]) %>% 
  select(age)

qBT94$Q <- c(rep(94))
qBT94$area <- c(rep("BT"))

# BT 96%
qBT96 <- df_perc_BT %>%
  filter(age >= qBT[4]) %>% 
  select(age)

qBT96$Q <- c(rep(96))
qBT96$area <- c(rep("BT"))

# BT 98%
qBT98 <- df_perc_BT %>%
  filter(age >= qBT[5]) %>% 
  select(age)

qBT98$Q <- c(rep(98))
qBT98$area <- c(rep("BT"))


# FM 90%
dfqFM = df_perc_FM$age
qFM <- quantile(dfqFM, probs = c(0.90, 0.92, 0.94, 0.96, 0.98))
qFM

# FM 90%
qFM90 <- df_perc_FM %>%
  filter(age >= qFM[1]) %>% 
  select(age)

qFM90$Q <- c(rep(90))
qFM90$area <- c(rep("FM"))

# FM 92%
qFM92 <- df_perc_FM %>%
  filter(age >= qFM[2]) %>% 
  select(age)

qFM92$Q <- c(rep(92))
qFM92$area <- c(rep("FM"))

# FM 94%
qFM94 <- df_perc_FM %>%
  filter(age >= qFM[3]) %>% 
  select(age)

qFM94$Q <- c(rep(94))
qFM94$area <- c(rep("FM"))

# FM 96%
qFM96 <- df_perc_FM %>%
  filter(age >= qFM[4]) %>% 
  select(age)

qFM96$Q <- c(rep(96))
qFM96$area <- c(rep("FM"))

# FM 98%
qFM98 <- df_perc_FM %>%
  filter(age >= qFM[5]) %>% 
  select(age)

qFM98$Q <- c(rep(98))
qFM98$area <- c(rep("FM"))

# New data frame with percentiles
Qdat <- rbind(qBT90, qFM90,
              qBT92, qFM92,
              qBT94, qFM94,
              qBT96, qFM96,
              qBT98, qFM98)

unique(Qdat$area)              
unique(Qdat$Q)

# Plot
ggplot(Qdat, aes(x = area, y =age)) +
  geom_boxplot() +
  facet_wrap(~Q, nrow = 1) +
  theme_bw()

# plot n-at-age  
ggplot(dat, aes(x = age, y = ct, fill = Area)) +
  geom_col(position = "dodge") +
  ylab("CPUE") +
  xlab("Age") +
  theme_classic(base_size = 18) +
  scale_color_brewer(palette="Set1") +
  theme(legend.position = "none")

# ========= % of total ==========
# Here I look at % of the total populations for each age
# (i.e. how large are the population of the 4, 5, 6, 7, 8, and 9 year olds)

# Starting by adding a column with total n for each area to dat
dat$total_n <- NA
dat$total_n[dat$Area=="BT"] <- sum(dat$ct[dat$Area=="BT"]) 
dat$total_n[dat$Area=="FM"] <- sum(dat$ct[dat$Area=="FM"]) 

# calc. %
dat$percent_of_tot <- ((dat$ct)/(dat$total_n))*100

# Plot
ggplot(dat, aes(x = age, y = percent_of_tot, fill = Area)) +
  geom_col(position = "dodge") +
  scale_x_continuous(breaks = 3:9) +
  theme_classic(base_size = 16) +
  ylab("% of tot.") +
  xlab("Age")
