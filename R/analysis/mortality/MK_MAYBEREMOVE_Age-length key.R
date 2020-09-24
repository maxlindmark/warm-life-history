####
## 2019.10.17 - Age-length key
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



#---- Read in the data
df.age <- read.csv("Growth_data_BT_FM_1970-2004.csv", sep = ";")
df.lenBT <- read.csv("Sizespectra_BT.csv", sep = ",") 
df.lenFM <- read.csv("Sizespectra_FM.csv", sep = ",")

# Stratified or random sampling - look for differences between years to find when/if change in method
df.age.BT <- df.age %>%
  filter(area=="BT")
df.age.FM <- df.age %>%
  filter(area=="FM")

ggplot(df.age.BT, aes(length)) +
  geom_histogram(binwidth = 20) +
  facet_wrap(~catch_year) +
  theme_classic()

ggplot(df.age.FM, aes(length)) +
  geom_histogram(binwidth = 20) +
  facet_wrap(~catch_year) +
  theme_classic()


# =============== Restructure data ===============

# merge df.lenBT and df.lenFM
df.len <- rbind(df.lenBT, df.lenFM)

# Convert lengths in df.len from cm to mm
df.len <- df.len %>%
  mutate(length_mm = length_group*10)

# Rename length_group to lengthGr_cm
df.len <- rename(df.len, replace = c("length_group" = "length_group_cm"))

min(df.len$length_mm)

# I will use all data from the length sample (i.e. not remove individivuals < 170 mm) since
# assigned ages for this sample will be based on the proportion of ages in a length group
# in the age sample and the proportions in the length sample will not matter.

min(df.age$length)

# In the age sample are lengths based on back-calculated length-at-age.
# Hence, no filter of small individuals is needed.

# =============== Construktion of ALK ===============
# ============= Construction of length categorys =============

# The length category to which each fish belongs is recorded.
# For example, if 5-mm length categories are created that begin on the “0” and “5” units,
# then a 117 mm fish will be recorded as being in the 115-119 mm length category
# Generally, all length categories are of the same width; thus, for simplicity,
# only the beginning length in the length category is recorded (e.g., “115” mm).

sort(unique(df.age$length))

# The starting category for 10-mm length categories is determined by finding the
# minimum length in the age sample with: 
Summarize(~length, data=df.age, digits=1)

# and then starting the categories with the even-number 10-mm interval just below this value.
# Construction of the length category variable is completed and the first six rows are viewed
# with:
rb.age1 <- lencat(~length, data=df.age, startcat = 39, w = 10)
head(rb.age1)

# This creates a variable in the age sample that identifies the length category to which
# each fish belongs.

sort(unique(rb.age1$LCat))

# Compare values
ggplot(rb.age1, aes(LCat, length)) +
   geom_point(size = 2) +
   geom_abline(slope = 1, intercept = 0, color = "red")

# Looks ok

# Once the length category variable has been added to the age sample data frame,
# table() is used to construct the summary contingency table of numbers of fish
# in each combined length and age category.

# The row variable (length category) is the first and the column variable (age)
# is the second argument to this function.
# The results of table() should be assigned to an object and then submitted as the
# first argument to prop.table() along with margin=1 as a second argument to construct
# a row-proportions table.
# The resulting row-proportions table is the actual age-length key determined from
# the age sample and is ready to be applied to the length sample.

# The summary contingency table and the row-proportion table (i.e., the age-length key)
# are constructed with:
rb.raw <- with(rb.age1,table(LCat,age))
rb.key <- prop.table(rb.raw,margin=1)
round(rb.key,2) # rounded for display purposes only

write.table(rb.key, file = "Two-way contingency table.csv", sep = ",", quote = FALSE, row.names = F)

# =============== Age assignment ===============
# =========== The semi-random age assignment ============
# There are two methods for using the age-length key to assign ages to specific fish in the length sample:
# semi-random and completely random methods. The semi-random method is the preferred method.
# Reed more in Ogel Age-length key vignette section 1.2.

# The semi-random method is implemented with ageKey()

# 1) This function requires the numeric matrix containing the age-length key
# (e.g., as constructed with prop.table()) as the first argument)

# 2) a formula of the form ~lengt or age~lengt as the second argument,

# 3) and the name of the age sample data frame in data=.

# If the age-sample data frame contains a variable to receive the ages then use age~len.
# Alternatively, if the age-sample data frame does not contain a variable to receive the ages
# then use ~len and the function will create a new variable called age by default.

# The ageKey() function will determine the length categories to construct based on the
# age-length key sent in the first argument. The ageKey() function returns a data frame that
# has the same number of rows as the original length sample, but with each fish assigned an age
# according to the semi-random method. The results of ageKey() should be assigned to an object,
# preferably with a name different from the original length sample. For example, semi-random
# ages were assigned to the un-aged fish in the length sample with

# ageKey() — DEPRECATED (will be removed by v1.0.0).  See alkIndivAge().

rb.len1 <- alkIndivAge(rb.key, ~length_mm, data = df.len)
head(rb.len1)

# Compare values
# Modified length sample in coral and the age sample in blue in the background.
ggplot(df.age, aes(x = age, y = length)) +
  geom_jitter(alpha = 0.07, color = "coral2") +
  geom_jitter(data = rb.len1, aes(x = age, y = length_mm), alpha = 0.05, color = "deepskyblue3") +
  scale_x_continuous(breaks = 1:9) +
  facet_wrap(~area)
  
# Summary
(rb.sum <- Summarize(length_mm~age,data=rb.len1,digits=2))

# Age frequency
p1 <- ggplot(rb.len1, aes(age)) +
  geom_histogram(fill = "deepskyblue3", binwidth = 0.5) +
  facet_wrap(~Area) +
  scale_x_continuous(breaks = 1:9) +
  theme_classic() +
  ylab("n") +
  ggtitle("A")


# Length-at-age for all individuals in the sample together with mean length-at-age 
p2 <- ggplot(rb.len1, aes(x = age, y = length_mm)) +
  geom_jitter(alpha= 0.05, color = "coral2") +
  geom_line(data = rb.sum, aes(x = age, y = mean), size = 1, color="deepskyblue4") +
  scale_x_continuous(breaks = 1:9) +
  facet_wrap(~Area) +
  theme_classic() +
  ggtitle("C")


# Length frequency
p3 <- ggplot(rb.len1, aes(length_mm)) +
  geom_histogram(fill = "deepskyblue3", binwidth = 30) +
  facet_wrap(~Area) +
  theme_classic() +
  ylab("n") +
  ggtitle("B")
  
grid.arrange(p1, p3, p2)

# Save data frame 
write.csv(rb.len1,"Length-at-age.csv", row.names = FALSE)










