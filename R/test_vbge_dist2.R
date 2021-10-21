# Testing different VBGE distributions

# Load libraries (install first if needed)
library(tidyverse); theme_set(theme_classic(base_size = 12))
library(brms)
library(nlstools)
library(viridis)
library(bayesplot)
library(tidylog)
library(tidybayes)
library(RColorBrewer)
library(patchwork)
library(modelr)

# Print package versions for versions
sessionInfo() 

# For parallel processing
options(mc.cores = parallel::detectCores()) 


# Read data 
# df <- read.csv("data/Growth_data_BT_FM_1970-2004.csv", sep = ";") # This is the 
# original data from Huss et al (2019)
df <- read.csv("data/size_at_age_BT_FM_1970-2004.csv", sep = ";")

# Add an area-specific ID
df$ID2 <- paste(df$ID, df$area, sep = ".")

# Filter years to match the rest of the data sets and analyses
df <- df %>% filter(birth_year > 1980 & birth_year < 1998,
                    catch_year < 2004 & catch_age == age)

# Plot data
p1 <- ggplot(df, aes(age, length, color = area)) +
  geom_point() +
  facet_wrap(~ birth_year)
p1

# Now we need to reshape the data frame a little bit to add in the dummy variable
# for area, and we'll use this dataframe for the models
bt <- filter(df, area == "BT")
fm <- filter(df, area == "FM")

dfm <- data.frame(rbind(cbind(bt, areaW=1, areaC=0), cbind(fm, areaW=0, areaC=1)))

dfm$age <- as.numeric(dfm$age)

ggplot(filter(dfm, birth_year == 1988), aes(age, length, color = area)) +
  geom_point()

# They seem to come from a different gear
filter(dfm, area == "FM" & birth_year == 1988 & length > 310)
ggplot(filter(dfm, area == "FM"), aes(age, length, color = factor(gear))) +
  geom_point(size = 0.4, alpha = 0.5) +
  facet_wrap(~ birth_year)

# Remove the odd gear
sort(unique(dfm$gear))
dfm <- dfm %>% filter(!gear == 32)
sort(unique(dfm$gear))

# Now it looks better
ggplot(filter(dfm, area == "FM"), aes(age, length, color = factor(gear))) +
  geom_point(size = 0.4, alpha = 0.5) +
  facet_wrap(~ birth_year)

# And in BT more or less only gear 9 is used
ggplot(filter(dfm, area == "BT"), aes(age, length, color = factor(gear))) +
  geom_point(size = 0.4, alpha = 0.5) +
  facet_wrap(~ birth_year)

# Lastly, fit the model using cm not mm
dfm$length_cm <- dfm$length / 10
dfm$log_length_cm <- log(dfm$length_cm)

# Change age to integer
dfm$age <- as.integer(dfm$age)

min(dfm$age)
min(dfm$birth_year) 

colnames(dfm)
d <- dfm %>% dplyr::select(length_cm, areaW, age, areaC, birth_year)
write.csv(d, "/Users/maxlindmark/Desktop/R_STUDIO_PROJECTS/stan_data/d.csv")


### Fit models
## Start stan forum
prior <-
  prior(normal(-0.5, 1), nlpar = "t0C") +
  prior(normal(-0.5, 1), nlpar = "t0W") +
  prior(normal(0.2, 0.1), nlpar = "KC") +
  #prior("uniform(0, 0.6)", nlpar = "KC", lb = 0, ub = 0.6) +
  prior(normal(0.2, 0.1), nlpar = "KW") +
  #prior("uniform(0, 0.6)", nlpar = "KW", lb = 0, ub = 0.6) +
  prior(normal(45, 20), nlpar = "LinfC") +
  prior(normal(45, 20), nlpar = "LinfW")

inits <- list(
  t0C = -0.5,
  t0W = -0.5,
  KC = 0.3,
  KW = 0.3
)

m1 <- 
  brm(
    bf(log(length_cm) ~ areaW*log(LinfW*(1-exp(-KW*(age-t0W)))) + areaC*log(LinfC*(1-exp(-KC*(age-t0C)))),
       t0C ~ 1,
       t0W ~ 1,
       KC ~ 1 + (1|birth_year),
       KW ~ 1 + (1|birth_year),
       LinfC ~ 1 + (1|birth_year),
       LinfW ~ 1 + (1|birth_year),
       nl = TRUE),
    data = d,
    family = student(),
    prior = prior,
    seed = 9, 
    iter = 50, thin = 1, cores = 1, chains = 1,
    inits = list(inits)
  )
plot(m1)



# More chains
list_of_inits <- list(inits, inits)

m2 <- 
  brm(
    bf(log(length_cm) ~ areaW*log(LinfW*(1-exp(-KW*(age-t0W)))) + areaC*log(LinfC*(1-exp(-KC*(age-t0C)))),
       t0C ~ 1,
       t0W ~ 1,
       KC ~ 1 + (1|birth_year),
       KW ~ 1 + (1|birth_year),
       LinfC ~ 1 + (1|birth_year),
       LinfW ~ 1 + (1|birth_year),
       nl = TRUE),
    data = dfm,
    family = student(),
    prior = prior,
    seed = 9, 
    iter = 50, thin = 1, cores = 2, chains = 2,
    inits = list_of_inits
  )

## End stan forums
str(m1)

inits1 <- m1$fit@inits

str(inits1)
inits

inits1_main <- inits1[[1]][1:9]

m2 <- 
  brm(
    bf(log(length_cm) ~ areaW*log(LinfW*(1-exp(-KW*(age-t0W)))) + areaC*log(LinfC*(1-exp(-KC*(age-t0C)))),
       t0C ~ 1, 
       t0W ~ 1,
       KC ~ 1 + (1|birth_year),
       KW ~ 1 + (1|birth_year),
       LinfC ~ 1 + (1|birth_year),
       LinfW ~ 1 + (1|birth_year),
       nl = TRUE),
    data = dfm,
    family = student(),
    prior = prior,
    seed = 9, 
    iter = 50, thin = 1, cores = 2, chains = 2,
    inits = list(inits1_main, inits1_main)
  )


m2b <- 
  brm(
    bf(log(length_cm) ~ areaW*log(LinfW*(1-exp(-KW*(age-t0W)))) + areaC*log(LinfC*(1-exp(-KC*(age-t0C)))),
       t0C ~ 1,
       t0W ~ 1,
       KC ~ 1 + (1|birth_year),
       KW ~ 1 + (1|birth_year),
       LinfC ~ 1 + (1|birth_year),
       LinfW ~ 1 + (1|birth_year),
       nl = TRUE),
    data = dfm,
    family = student(),
    prior = prior,
    seed = 9, 
    iter = 50, thin = 1, cores = 3, chains = 3,
    inits = list(inits1_main, inits1_main, inits1_main)
  )


# rnorm(n = 1, mean = 0, sd = 0.01)
inits1 <- inits1_main
inits2 <- lapply(inits1_main, "*", 1 + rnorm(n = 1, mean = 0, sd = 0.01))
inits3 <- lapply(inits1_main, "*", 1 + rnorm(n = 1, mean = 0, sd = 0.01))

inits_3_chain <- list(inits1, inits2, inits3)

save(inits_3_chain, file = "R/analysis/vbge_3_chain_inits.RData")
?save

m2c <- 
  brm(
    bf(log(length_cm) ~ areaW*log(LinfW*(1-exp(-KW*(age-t0W)))) + areaC*log(LinfC*(1-exp(-KC*(age-t0C)))),
       t0C ~ 1,
       t0W ~ 1,
       KC ~ 1 + (1|birth_year),
       KW ~ 1 + (1|birth_year),
       LinfC ~ 1 + (1|birth_year),
       LinfW ~ 1 + (1|birth_year),
       nl = TRUE),
    data = dfm,
    family = student(),
    prior = prior,
    seed = 9, 
    iter = 50, thin = 1, cores = 3, chains = 3,
    inits = initis_3_chain
  )


# log normal standard model
# prior <-
#   prior(normal(-0.5, 1), nlpar = "t0C") +
#   prior(normal(-0.5, 1), nlpar = "t0W") +
#   prior(normal(0.2, 0.1), nlpar = "KC") +
#   prior(normal(0.2, 0.1), nlpar = "KW") +
#   prior(normal(45, 20), nlpar = "LinfC") +
#   prior(normal(45, 20), nlpar = "LinfW")
# 
# m3 <-
#   brm(
#     bf(length_cm ~ areaW*LinfW*(1-exp(-KW*(age-t0W))) + areaC*LinfC*(1-exp(-KC*(age-t0C))),
#        t0C ~ 1 + (1|birth_year),
#        t0W ~ 1 + (1|birth_year),
#        KC ~ 1 + (1|birth_year),
#        KW ~ 1 + (1|birth_year),
#        LinfC ~ 1 + (1|birth_year),
#        LinfW ~ 1 + (1|birth_year),
#        nl = TRUE),
#     data = dfm,
#     family = lognormal(),
#     prior = prior,
#     seed = 9,
#     iter = 4000, thin = 1, cores = 3, chains = 3
#   )


## Current VBGE but with random t_0
prior <-
  prior(normal(-0.5, 1), nlpar = "t0C") +
  prior(normal(-0.5, 1), nlpar = "t0W") +
  prior("uniform(0, 0.6)", nlpar = "KC", lb = 0, ub = 0.6) +
  prior("uniform(0, 0.6)", nlpar = "KW", lb = 0, ub = 0.6) +
  prior(normal(45, 20), nlpar = "LinfC") +
  prior(normal(45, 20), nlpar = "LinfW")

start_time <- Sys.time()
m4 <- 
  brm(
    bf(length_cm ~ areaW*LinfW*(1-exp(-KW*(age-t0W))) + areaC*LinfC*(1-exp(-KC*(age-t0C))),
       t0C ~ 1 + (1|birth_year),
       t0W ~ 1 + (1|birth_year),
       KC ~ 1 + (1|birth_year),
       KW ~ 1 + (1|birth_year),
       LinfC ~ 1 + (1|birth_year),
       LinfW ~ 1 + (1|birth_year),
       nl = TRUE),
    data = dfm,
    family = student(),
    prior = prior,
    seed = 9, 
    iter = 4000, thin = 1, cores = 3, chains = 3, inits = "0",
    control = list(max_treedepth = 13, adapt_delta = 0.9))
end_time <- Sys.time()
end_time - start_time

ggplot(df, aes(sample = y)) +
  stat_qq(distribution = qt, dparams = params["df"]) +
  stat_qq_line(distribution = qt, dparams = params["df"])

summary(m4)$spec_pars # Extract "fixed" effects from m1 for plotting the equation 
nu <- summary(m4)$spec_pars[2, 1]
nu

dfm %>%
  add_residual_draws(m4) %>%
  median_qi() %>%
  ggplot(aes(sample = .residual)) +
  geom_qq_line(distribution = qt, dparams = nu) +
  geom_qq(alpha = 0.8, distribution = qt, dparams = nu) +
  theme(text = element_text(size = 12))


### student likelihood on log data (not only one chain works! make it count by setting 4k iterations... )
prior <-
  prior(normal(-0.5, 1), nlpar = "t0C") +
  prior(normal(-0.5, 1), nlpar = "t0W") +
  prior(normal(0.2, 0.1), nlpar = "KC") +
  prior(normal(0.2, 0.1), nlpar = "KW") +
  prior(normal(45, 20), nlpar = "LinfC") +
  prior(normal(45, 20), nlpar = "LinfW")

m5 <- 
  brm(
    bf(log(length_cm) ~ areaW*log(LinfW*(1-exp(-KW*(age-t0W)))) + areaC*log(LinfC*(1-exp(-KC*(age-t0C)))),
       t0C ~ 1,
       t0W ~ 1,
       KC ~ 1 + (1|birth_year),
       KW ~ 1 + (1|birth_year),
       LinfC ~ 1 + (1|birth_year),
       LinfW ~ 1 + (1|birth_year),
       nl = TRUE),
    data = dfm,
    family = student(),
    prior = prior,
    seed = 9, 
    iter = 4000, thin = 1, cores = 1, chains = 1,
    inits = list(inits)
  )
plot(m5)

saveRDS(m4, "output/vbge/m4test.rds")
saveRDS(m5, "output/vbge/m5test.rds")

summary(m5)$spec_pars # Extract "fixed" effects from m1 for plotting the equation 
nu <- summary(m5)$spec_pars[2, 1]
nu

# https://ggplot2.tidyverse.org/reference/geom_qq.html
dfm %>%
  add_residual_draws(m5) %>%
  median_qi() %>%
  ggplot(aes(sample = .residual)) +
  geom_qq_line(distribution = qt, dparams = nu) +
  geom_qq(alpha = 0.8, distribution = qt, dparams = nu) +
  theme(text = element_text(size = 12))

