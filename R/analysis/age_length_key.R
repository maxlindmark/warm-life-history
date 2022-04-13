#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 2021.08.04: Max Lindmark
#
# Create age-length keys and assign ages to catch data
#
# A. Load libraries
# 
# B. Read data
# 
# C. Create area-specific ALK
# 
# D. Age the catch data, save
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


# B. READ DATA =====================================================================
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

# Now we need to use the same length categories as in the catch data
catch_data <- read.csv("data/catch_BT_FM_1987-2003.csv")

length_at_age <- length_at_age %>% 
  filter(length > 5.1) %>% # This is the smallest category in the catch data
  mutate(len_cat = lencat(length, w = 0.5, startcat = 5.1))
  
len_cats <- sort(unique(length_at_age$len_cat))

ggplot(length_at_age, aes(len_cat, length)) + geom_point()


# C. AREA-SPECIFIC ALK =============================================================
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

  bt_temp_df <- length_at_age %>% filter(area == "BT")
  fm_temp_df <- length_at_age %>% filter(area == "FM")
  
  # Create age-length key
  bt_alk_raw <- xtabs(~len_cat + back_calc_age, data = bt_temp_df)
  fm_alk_raw <- xtabs(~len_cat + back_calc_age, data = fm_temp_df)
  
  bt_alk <- prop.table(bt_alk_raw, margin = 1)
  fm_alk <- prop.table(fm_alk_raw, margin = 1)
  
# Now age the catch data by each area using the area-specific age length key
  # area
  bt_catch <- catch_data %>% filter(Area == "BT")
  fm_catch <- catch_data %>% filter(Area == "FM")
  
  bt_aged_catch <- alkIndivAge(bt_alk, ~length_group, data = bt_catch)
  fm_aged_catch <- alkIndivAge(fm_alk, ~length_group, data = fm_catch)
  
# Combine into single data frame
aged_catch <- rbind(bt_aged_catch, fm_aged_catch)


# D. PLOT AND SAVE =================================================================
aged_catch %>% 
  mutate(area2 = ifelse(Area == "BT", "Warm", "Cold")) %>% 
  ggplot(., aes(length_group, fill = factor(age))) +
  geom_histogram(binwidth = 1, position = "stack") + 
  facet_wrap(~ area2, scales = "free") + 
  scale_fill_brewer(palette = "Set1") + 
  theme(text = element_text(size = 12),
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 8),
        aspect.ratio = 1) + 
  coord_cartesian(expand = 0) + 
  labs(fill = "Age",
       x = "Length group",
       y = "Count")
  
ggsave("figures/supp/age_by_length.png", width = 6.5, height = 6.5, dpi = 600)

# Save data frame for catch curve analysis
write.csv(aged_catch, "data/aged_catch_BT_FM_1987-2003.csv", row.names = FALSE)
