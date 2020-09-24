####
## 2019.12.19 - Mean age from semi-random ALK
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
  filter(area == "BT")

summary(df.ageBT)

df.ageFM <- df.age %>% 
  filter(area == "FM")

summary(df.ageFM)

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
sort(unique(df.ageBT$length))
sort(unique(df.ageFM$length))

# The starting category for 10-mm length categories is determined by finding the
# minimum length in the age sample with: 
Summarize(~length, data=df.ageBT, digits=1) # 43.0
Summarize(~length, data=df.ageFM, digits=1) # 40.0


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

# ========== Mean age ==========
df <- rb.len

# add catch column
df_sub <- df %>%
  dplyr::group_by(age, Area) %>%
  dplyr::summarise(count = n()) %>%
  dplyr::ungroup()

# count --> catch (ct)
df_sub <- df_sub %>%
  dplyr::rename(ct = count)


# From the initial catch curve plot it was determined that the descending limb consisted of
# ages 3 through 9 (BT) and 2 - 9 (FM).
# I need to create a new data frame which contains only the information for ages (3-9). 
dat1 <- df_sub %>%
  filter(age > 2)

dat2 <- df %>%
  filter(age > 2)

# ============== Av. age in certain percentiles ===================

# Split data frames
ageBT <- dat2 %>%
  dplyr::filter(Area == "BT")

ageFM <- dat2 %>%
  filter(Area == "FM")

# Create vector with only age
ageBT = data.frame(age=ageBT$age)
ageFM = data.frame(age=ageFM$age)

# BT
qBT <- quantile(ageBT$age, probs = c(0.90, 0.92, 0.94, 0.96, 0.98))
qBT[1]
qBT[2]

qBT90 <- ageBT %>%
  filter(age >= qBT[1]) %>% 
  select(age)

qBT90$Q <- c(rep(90))
qBT90$area <- c(rep("BT"))

qBT92 <- ageBT %>%
  filter(age >= qBT[2]) %>% 
  select(age)

qBT92$Q <- c(rep(92))
qBT92$area <- c(rep("BT"))

qBT94 <- ageBT %>%
  filter(age >= qBT[3]) %>% 
  select(age)

qBT94$Q <- c(rep(94))
qBT94$area <- c(rep("BT"))

qBT96 <- ageBT %>%
  filter(age >= qBT[4]) %>% 
  select(age)

qBT96$Q <- c(rep(96))
qBT96$area <- c(rep("BT"))

qBT98 <- ageBT %>%
  filter(age >= qBT[5]) %>% 
  select(age)

qBT98$Q <- c(rep(98))
qBT98$area <- c(rep("BT"))


# FM
qFM <- quantile(ageFM$age, probs = c(0.90, 0.92, 0.94, 0.96, 0.98))
qFM[1]
qFM[2]

qFM90 <- ageFM %>%
  filter(age >= qFM[1]) %>% 
  select(age)

qFM90$Q <- c(rep(90))
qFM90$area <- c(rep("FM"))

qFM92 <- ageFM %>%
  filter(age >= qFM[2]) %>% 
  select(age)

qFM92$Q <- c(rep(92))
qFM92$area <- c(rep("FM"))

qFM94 <- ageFM %>%
  filter(age >= qFM[3]) %>% 
  select(age)

qFM94$Q <- c(rep(94))
qFM94$area <- c(rep("FM"))

qFM96 <- ageFM %>%
  filter(age >= qFM[4]) %>% 
  select(age)

qFM96$Q <- c(rep(96))
qFM96$area <- c(rep("FM"))

qFM98 <- ageFM %>%
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

# Merge dataframes for CI plot 
q90 <- rbind(qBT90, qFM90)
q92 <- rbind(qBT92, qFM92)
q94 <- rbind(qBT94, qFM94)
q96 <- rbind(qBT96, qFM96)
q98 <- rbind(qBT98, qFM98)

# lm for separate percentiles 
lm90 <- lm(q90$age ~ q90$area)
summary(lm90)
lm92 <- lm(q92$age ~ q92$area)
summary(lm92)
lm94 <- lm(q94$age ~ q94$area)
summary(lm94)
lm96 <- lm(q96$age ~ q96$area)
summary(lm96)
lm98 <- lm(q98$age ~ q98$area)
summary(lm98)

# Estimate of mean lenght for each percentile
est90BT <- summary(lm90)$coefficients[1,1]
est90FM <- summary(lm90)$coefficients[1,1] + summary(lm90)$coefficients[2,1]

est92BT <- summary(lm92)$coefficients[1,1]
est92FM <- summary(lm92)$coefficients[1,1] + summary(lm92)$coefficients[2,1]

est94BT <- summary(lm94)$coefficients[1,1]
est94FM <- summary(lm94)$coefficients[1,1] + summary(lm94)$coefficients[2,1]

est96BT <- summary(lm96)$coefficients[1,1]
est96FM <- summary(lm96)$coefficients[1,1] + summary(lm96)$coefficients[2,1]

est98BT <- summary(lm98)$coefficients[1,1]
est98FM <- summary(lm98)$coefficients[1,1] + summary(lm98)$coefficients[2,1]

# CI
LCI90BT <- confint(lm90)[1,1]
UCI90BT <- confint(lm90)[1,2]
LCI90FM <- confint(lm90)[1,1] + confint(lm90)[2,1]
UCI90FM <- confint(lm90)[1,2] + confint(lm90)[2,2]

LCI92BT <- confint(lm92)[1,1]
UCI92BT <- confint(lm92)[1,2]
LCI92FM <- confint(lm92)[1,1] + confint(lm92)[2,1]
UCI92FM <- confint(lm92)[1,2] + confint(lm92)[2,2]

LCI94BT <- confint(lm94)[1,1]
UCI94BT <- confint(lm94)[1,2]
LCI94FM <- confint(lm94)[1,1] + confint(lm94)[2,1]
UCI94FM <- confint(lm94)[1,2] + confint(lm94)[2,2]

LCI96BT <- confint(lm96)[1,1]
UCI96BT <- confint(lm96)[1,2]
LCI96FM <- confint(lm96)[1,1] + confint(lm96)[2,1]
UCI96FM <- confint(lm96)[1,2] + confint(lm96)[2,2]

LCI98BT <- confint(lm98)[1,1]
UCI98BT <- confint(lm98)[1,2]
LCI98FM <- confint(lm98)[1,1] + confint(lm98)[2,1]
UCI98FM <- confint(lm98)[1,2] + confint(lm98)[2,2]

# New data frame with estimate and CI
Qdat2 <- data.frame(Percentile = c("90%", "90%",
                                   "92%", "92%",
                                   "94%", "94%",
                                   "96%", "96%",
                                   "98%", "98%"),
                    area = c("BT", "FM",
                             "BT", "FM",
                             "BT", "FM",
                             "BT", "FM",
                             "BT", "FM"),
                    mean = c(est90BT, est90FM,
                             est92BT, est92FM,
                             est94BT, est94FM,
                             est96BT, est96FM,
                             est98BT, est98FM),
                    LCI = c(LCI90BT, LCI90FM,
                            LCI92BT, LCI92FM,
                            LCI94BT, LCI94FM,
                            LCI96BT, LCI96FM,
                            LCI98BT, LCI98FM),
                    UCI = c(UCI90BT, UCI90FM,
                            UCI92BT, UCI92FM,
                            UCI94BT, UCI94FM,
                            UCI96BT, UCI96FM,
                            UCI98BT, UCI98FM))

# Plot estimate and CI
ggplot(Qdat2, aes(x = Percentile, y = mean, color = factor(area))) +
  geom_point(size=3) +
  geom_errorbar(aes(Percentile, ymin = LCI, ymax = UCI), width = 0.1) +
  theme_classic(base_size = 18) +
  ylab("Mean age (years)") +
  scale_color_brewer(palette="Set1") +
  theme(legend.position = "none")

# ggsave("C:/Users/malin/OneDrive - Sveriges Lantbruksuniversitet/Masterarbete/R_studio project/Figurer ANALYS/Mortality/av_age_percentiles.png",
#        plot = last_plot(), scale = 1, width = 16, height = 16, units = "cm", dpi = 300)


# plot n-at-age  
ggplot(dat1, aes(x = age, y = ct, fill = Area)) +
  geom_col(position = "dodge") +
  ylab("CPUE") +
  xlab("Age") +
  theme_classic(base_size = 18) +
  scale_fill_brewer(palette="Set1") +
  theme(legend.position = "none")

# ggsave("C:/Users/malin/OneDrive - Sveriges Lantbruksuniversitet/Masterarbete/R_studio project/Figurer ANALYS/Mortality/CPUE_at_age.png",
#        plot = last_plot(), scale = 1, width = 16, height = 16, units = "cm", dpi = 300)

  

# ========= % of total for each age group ==========


# Add a column with total_n for each age group
dat1$total_n <- NA
dat1$total_n[dat1$Area=="BT"] <- sum(dat1$ct[dat1$Area=="BT"]) 
dat1$total_n[dat1$Area=="FM"] <- sum(dat1$ct[dat1$Area=="FM"])  

# calc. %
dat1$percent_of_tot <- ((dat1$ct)/(dat1$total_n))*100

unique(dat1$age)


# Plot
ggplot(dat1, aes(x = age, y = percent_of_tot, fill = Area)) +
  geom_col(position = "dodge") +
  theme_classic(base_size = 16) +
  scale_x_continuous(breaks = c(30.1, 32.6, 35.1, 37.6, 40.1, 42.6, 45.1)) +
  ylab("% of tot.") +
  xlab("Age (years)") +
  scale_fill_brewer(palette="Set1") +
  theme(legend.position = "none")

# ggsave("C:/Users/malin/OneDrive - Sveriges Lantbruksuniversitet/Masterarbete/R_studio project/Figurer ANALYS/Mortality/percent_of_total.png",
#        plot = last_plot(), scale = 1, width = 16, height = 16, units = "cm", dpi = 300)

