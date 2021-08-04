#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 2021.08.04: Max Lindmark
#
# Create age-length keys and assign ages to catch data
#
# A. Load libraries
# 
# B. Read data
# 
# C. Create age-length key and assign ages to catch data 
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# A. LOAD LIBRARIES ================================================================
rm(list = ls())

# Load libraries (install first if needed)
library(tidyverse); theme_set(theme_classic(base_size = 12))
library(tidylog)
library(RColorBrewer)
library(patchwork)
library(FSA)


# B. LOAD LIBRARIES ================================================================
length_at_age <- read.csv("data/size_at_age_BT_FM_1970-2004.csv", sep = ";")

# Remove gear 32, see vbge script
length_at_age <- length_at_age %>% filter(!gear == 32)

# Filter years I have catch data for
length_at_age <- length_at_age %>% filter(catch_year > 1986 & catch_year < 2004)

# I need to make the ID completely unique (since I should include ID information in my
# mixed model), otherwise they may have the same ID in different areas.
# We solve this by writing ID with area + ID. The paste function pastes things together
length_at_age$ID <- paste(length_at_age$ID, length_at_age$area, sep = "")

# Check catch age corresponds to number of read length-at-ages!
length_at_age <- length_at_age %>%
  group_by(ID) %>%
  mutate(n = n()) %>%
  mutate(length = length/10) %>%
  rename("back_calc_age" = "age") %>% 
  ungroup()

ggplot(length_at_age, aes(factor(catch_age), factor(n))) + geom_point()

# Remove individuals with the same ID (bug)
length_at_age <- length_at_age %>% filter(n <= catch_age)

# Now every individual has the same number of rows as it's age in years
ggplot(length_at_age, aes(factor(catch_age), factor(n))) + geom_point()

# Plot size distributions in sample
length_at_age %>%
  filter(area == "FM") %>% 
  ggplot(., aes(length)) +
  geom_histogram() +
  facet_wrap(~ catch_year) +
  theme_classic()

length_at_age %>%
  filter(area == "BT") %>% 
  ggplot(., aes(length)) +
  geom_histogram() +
  facet_wrap(~ catch_year) +
  theme_classic()

# Check if NA age...

# Now we need to use the same length categories as in the catch data
catch_data <- read.csv("data/catch_BT_FM_1987-2003.csv")

length_at_age <- length_at_age %>% 
  filter(length > 5.1) %>% # This is the smallest category in the catch data
  mutate(len_cat = lencat(length, w = 0.5, startcat = 5.1))
  
len_cats <- sort(unique(length_at_age$len_cat))

ggplot(length_at_age, aes(len_cat, length)) + geom_point()

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

# We need to make one of these for each year in the catch data...

rb_raw_list <- list()
rb_key_list <- list()
temp_df <- data.frame()

for(i in sort(unique(length_at_age$catch_year))) {
  
  temp_df <- length_at_age %>% filter(catch_year == i)

  # rb_raw <- with(temp_df, table(len_cat, back_calc_age))
  # rb_key <- prop.table(rb_raw, margin = 1)
  
  # Create age-length key
  alk_raw <- xtabs(~len_cat + back_calc_age, data = temp_df)
  alk <- prop.table(alk_raw, margin = 1)
  
  index <- i - min(unique(length_at_age$catch_year)) + 1 
  
  rb_raw_list[[index]] <- alk_raw
  rb_key_list[[index]] <- alk
  
}

names(rb_raw_list) <- sort(unique(length_at_age$catch_year))
names(rb_key_list) <- sort(unique(length_at_age$catch_year))

str(rb_raw_list)
str(rb_key_list)


# Now write a for loop that ages the catch data by each year, using the year-specific 
# age length key

# Store annual data frames with aged catch data from the year-specific ALK
aged_catch_dat <- list()

for(i in sort(unique(catch_data$year))) {
  
  # Subset years
  temp_catch <- catch_data %>% filter(year == i)
  temp_rb_key <- rb_key_list[i - min(catch_data$year) + 1][[1]]
  
  temp_catch <- alkIndivAge(temp_rb_key, ~length_group, data = temp_catch)

  index <- i - min(unique(length_at_age$catch_year)) + 1 
  
  aged_catch_dat[[index]] <- temp_catch
  
}

aged_catch <- dplyr::bind_rows(aged_catch_dat)

sort(unique(aged_catch$year))

ggplot(aged_catch, aes(length_group, fill = factor(age))) +
  geom_histogram(binwidth = 1, position = "stack") + 
  facet_wrap(~ Area, scales = "free") + 
  scale_fill_brewer(palette = "Set1")

# Save data frame for catch curve analysis
write.csv(aged_catch, "data/aged_catch_BT_FM_1987-2003.csv", row.names = FALSE)


# 1 row 1 ind in catch data? Check data cleaning script!
