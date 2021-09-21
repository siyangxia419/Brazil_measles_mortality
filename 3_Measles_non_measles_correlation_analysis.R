### Analysis of the correlation between measles & non-measles disease mortality
###
### Input: Mortality data prepared to the desired format
### The data shoule locate in the directory "prepared data/" inside the location directory
### 
### Siyang Xia
### 2021.4.29



library(reshape2)
library(ggplot2)
library(cowplot)
library(dplyr)
library(mcp)
library(stargazer)
library(relaimpo)




### GLOBAL OPTIONS:
LOC <- "Brazil"
AGE <- "1-9"
DZZ <- "all"  # global variable



### Some miscellaneous preparation:
# Color-blind friendly palette:
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")



# Function to calculate proportional location in a range:
position_in_range <- function(x, p) min(x) + p * (max(x) - min(x))



# 1. Data input -----------------------------------------------------------

master.dir <- "C:/Users/Siyang Xia/Dropbox/My laptop/Mina lab/Projects/Measles/Empirical/Brazil/"
data.dir <- paste0(master.dir, LOC, "/prepared data/")
output.dir <- paste0(master.dir, LOC, "/measles vs non-measles/")
if(!dir.exists(output.dir)) dir.create(output.dir)


# Options of geographic locations and non-measles disease:
cities <- c("Brazil", "B horizonte", "Belem", "Fortaleza", "Goiania", "Manaus", "Porto Alegre", 
            "Recife", "Rio de Janeiro", "Salvador", "Sau Paulo")
diseases <- c("all", "TB", "URT", "bronchitis", "diarrheal", "meningitis", "pneumonia", 
              "pneumococcus", "rubella", "septicemia", "LRT", "asthma", "emphysema", "meningococcal") 

# Read the data:
infection <- read.csv(paste0(data.dir, LOC, "_measles_non-measles_mortality.csv"), 
                      header = TRUE, stringsAsFactors = FALSE)

# Some double-check:
all(round(infection$meV_mort / infection$meV_CFR - infection$meV_inc, digits = 10) == 0)
all(round(infection$meV_mort / infection$pop * 100000 - infection$meV_mort_rate, digits = 10) == 0)
all(round(infection$meV_inc / infection$pop * 100000 - infection$meV_inc_rate, digits = 10) == 0)
all(round(infection$dzz_mort / infection$pop * 100000 - infection$dzz_rate, digits = 10) == 0)





# 2. Preparation and data subsetting --------------------------------------

# Output directory for each age group
adt.dir <- paste0(output.dir, AGE, "/")
if(!dir.exists(adt.dir)) dir.create(adt.dir)

# Select the data of the entire Brazil:
d.bz <- infection %>% dplyr::filter(city == LOC, age == AGE, dzzName == DZZ)
nyear <- nrow(d.bz)

# Index the year
d.bz$year.idx <- d.bz$year - min(d.bz$year) + 1   # Use year as a time series





# 3. Change-point analysis -----------------------------------

# Load the function:
source(paste0(master.dir, "Function - change-point analysis.R"))

# Output directory:
cp.dir <- paste0(adt.dir, "change-point/")
if(!dir.exists(cp.dir)) dir.create(cp.dir)



# (a) measles time series: count -------------------------------------------

fit.bz.count <- change_point(d = d.bz, write.dir = cp.dir, 
                             ts.var = "meV_mort", year.var = "year.idx", 
                             LOC = LOC, DZZ = DZZ, output.ts = FALSE)

# Examine the results:
fit.bz.count.summary <- summary(fit.bz.count)

# Mean and 95% CI of the change-point
cp.bz.count <- fit.bz.count.summary[fit.bz.count.summary$name == "cp_1", c("mean", "lower", "upper")]
cp.bz.count.year <- as.numeric(cp.bz.count + min(d.bz$year) - 1)
names(cp.bz.count.year) <- c("mean", "lower", "upper")
print(cp.bz.count.year)

# The general time trend
cp.bz.count.trend <- stats::predict(fit.bz.count, 
                                    newdata = data.frame(year.idx = pretty(range(d.bz$year.idx), n = 100)))

# Visualization:
plot.fit.bz.count <- plot(fit.bz.count, q_fit = c(0.025, 0.975)) + 
  geom_line(data = cp.bz.count.trend, aes(x = year.idx, y = predict), 
            col = "red", size = 0.8, alpha = 0.6) +
  scale_x_continuous(name = "year", 
                     breaks = seq(1, 16, 5), 
                     labels = seq(1980, 1995, 5)) + 
  ylab("measles mortality") + 
  theme_classic() + 
  theme(panel.grid.major = element_line(color = "gray95"),
        axis.title = element_text(size = 10),
        axis.text = element_text(size = 10))
print(plot.fit.bz.count)
ggsave(filename = paste0(cp.dir, "Change-point-measles-count.png"), 
       plot = plot.fit.bz.count,
       width = 6, height = 4, units = "in", dpi = 300)

# Calculate the time-detrended measles mortality rate:
d.bz$meV_mort_dt <- d.bz$meV_mort - fitted(fit.bz.count)[, "fitted"]



# (b) measles time series: rate -------------------------------------------

fit.bz.rate <- change_point(d = d.bz, write.dir = cp.dir, 
                            ts.var = "meV_mort_rate", year.var = "year.idx", 
                            LOC = LOC, DZZ = DZZ, output.ts = FALSE)

# Look at the results:
fit.bz.rate.summary <- summary(fit.bz.rate)

# Mean and 95% CI of the change-point
cp.bz.rate <- fit.bz.rate.summary[fit.bz.rate.summary$name == "cp_1", c("mean", "lower", "upper")]
cp.bz.rate.year <- as.numeric(cp.bz.rate + min(d.bz$year) - 1)
names(cp.bz.rate.year) <- c("mean", "lower", "upper")
print(cp.bz.rate.year)

# The general time trend
cp.bz.rate.trend <- stats::predict(fit.bz.rate, 
                                   newdata = data.frame(year.idx = pretty(range(d.bz$year.idx), n = 100)))

# Visualization:
plot.fit.bz.rate <- plot(fit.bz.rate, q_fit = c(0.025, 0.975)) + 
  geom_line(data = cp.bz.rate.trend, aes(x = year.idx, y = predict), 
            col = "red", size = 0.8, alpha = 0.6) +
  scale_x_continuous(name = "year", 
                     breaks = seq(1, 16, 5), 
                     labels = seq(1980, 1995, 5)) + 
  ylab("measles mortality/100,000 ppl") + 
  theme_classic() + 
  theme(panel.grid.major = element_line(color = "gray95"),
        axis.title = element_text(size = 10),
        axis.text = element_text(size = 10))
print(plot.fit.bz.rate)
ggsave(filename = paste0(cp.dir, "Change-point-measles-rate.png"), 
       plot = plot.fit.bz.rate,
       width = 6, height = 4, units = "in", dpi = 300)

# Calculate the time-detrended measles mortality rate:
d.bz$meV_mort_rate_dt <- d.bz$meV_mort_rate - fitted(fit.bz.rate)[, "fitted"]



# (c) non-measles time series: count --------------------------------------

fit.bz.dzz.count <- change_point(d = d.bz, write.dir = cp.dir, 
                                 ts.var = "dzz_mort", year.var = "year.idx", 
                                 LOC = LOC, DZZ = DZZ, output.ts = FALSE)

# Look at the results:
fit.bz.dzz.count.summary <- summary(fit.bz.dzz.count)

# Mean and 95% CI of the change-point
cp.bz.dzz.count <- fit.bz.dzz.count.summary[fit.bz.dzz.count.summary$name == "cp_1", c("mean", "lower", "upper")]
cp.bz.dzz.count.year <- as.numeric(cp.bz.dzz.count + min(d.bz$year) - 1)
names(cp.bz.dzz.count.year) <- c("mean", "lower", "upper")
print(cp.bz.dzz.count.year)

# The general time trend
cp.bz.dzz.count.trend <- stats::predict(fit.bz.dzz.count, 
                                        newdata = data.frame(year.idx = pretty(range(d.bz$year.idx), n = 100)))

# Visualization:
plot.fit.bz.dzz.count <- plot(fit.bz.dzz.count, q_fit = c(0.025, 0.975)) + 
  geom_line(data = cp.bz.dzz.count.trend, aes(x = year.idx, y = predict), 
            col = "red", size = 0.8, alpha = 0.6) +
  scale_x_continuous(name = "year", 
                     breaks = seq(1, 16, 5), 
                     labels = seq(1980, 1995, 5)) + 
  ylab("non-measles mortality") + 
  theme_classic() + 
  theme(panel.grid.major = element_line(color = "gray95"),
        axis.title = element_text(size = 10),
        axis.text = element_text(size = 10))
print(plot.fit.bz.dzz.count)
ggsave(filename = paste0(cp.dir, "Change-point-non-measles-count.png"), 
       plot = plot.fit.bz.dzz.count,
       width = 6, height = 4, units = "in", dpi = 300)

# Calculate the time-detrended measles mortality rate:
d.bz$dzz_mort_dt <- d.bz$dzz_mort - fitted(fit.bz.dzz.count)[, "fitted"]



# (d) non-measles time series: rate ---------------------------------------

fit.bz.dzz.rate <- change_point(d = d.bz, write.dir = cp.dir, 
                                ts.var = "dzz_rate", year.var = "year.idx", 
                                LOC = LOC, DZZ = DZZ, output.ts = FALSE)

# Look at the results:
fit.bz.dzz.rate.summary <- summary(fit.bz.dzz.rate)

# Mean and 95% CI of the change-point
cp.bz.dzz.rate <- fit.bz.dzz.rate.summary[fit.bz.dzz.rate.summary$name == "cp_1", c("mean", "lower", "upper")]
cp.bz.dzz.rate.year <- as.numeric(cp.bz.dzz.rate + min(d.bz$year) - 1)
names(cp.bz.dzz.rate.year) <- c("mean", "lower", "upper")
print(cp.bz.dzz.rate.year)

# The general time trend
cp.bz.dzz.rate.trend <- stats::predict(fit.bz.dzz.rate, 
                                       newdata = data.frame(year.idx = pretty(range(d.bz$year.idx), n = 100)))

# Visualization:
plot.fit.bz.dzz.rate <- plot(fit.bz.dzz.rate, q_fit = c(0.025, 0.975)) + 
  geom_line(data = cp.bz.dzz.rate.trend, aes(x = year.idx, y = predict), 
            col = "red", size = 0.8, alpha = 0.6) +
  scale_x_continuous(name = "year", 
                     breaks = seq(1, 16, 5), 
                     labels = seq(1980, 1995, 5)) + 
  ylab("non-measles mortality/100,000 ppl") + 
  theme_classic() + 
  theme(panel.grid.major = element_line(color = "gray95"),
        axis.title = element_text(size = 10),
        axis.text = element_text(size = 10))
print(plot.fit.bz.dzz.rate)
ggsave(filename = paste0(cp.dir, "Change-point-non-measles-rate.png"), 
       plot = plot.fit.bz.dzz.rate,
       width = 6, height = 4, units = "in", dpi = 300)

# Calculate the time-detrended measles mortality rate:
d.bz$dzz_rate_dt <- d.bz$dzz_rate - fitted(fit.bz.dzz.rate)[, "fitted"]





# 4. Visualization: time seires ---------------------------------------------------------

source(paste0(master.dir, "Function - time series and correlation plots.R"))


# (a) Figure 1.1: Time series of the absolute case count ------------------

f1.1 <- time_series(d.ts = d.bz, 
                    meV_var = "meV_mort", nonmeV_var = "dzz_mort", 
                    cp = cp.bz.count.year, 
                    data_type = c("count", "count"), 
                    data_manipulation = "original")
print(f1.1)
ggsave(filename = paste0(adt.dir, "Time series-original-count.png"), 
       plot = f1.1,
       width = 7, height = 4, units = "in", dpi = 300)





# (b) Figure 1.2: Time series of the rate ---------------------------------

f1.2 <- time_series(d.ts = d.bz, 
                    meV_var = "meV_mort_rate", nonmeV_var = "dzz_rate", 
                    cp = cp.bz.rate.year, 
                    data_type = c("rate", "rate"), 
                    data_manipulation = "original")
print(f1.2)
ggsave(filename = paste0(adt.dir, "Time series-original-rate.png"), 
       plot = f1.2,
       width = 7, height = 4, units = "in", dpi = 300)





# (c) Figure 2.1 Correlation between measles and non-measles mortality count --------

f2.1 <- correlation_plot(d.ts = d.bz, 
                         meV_var = "meV_mort", nonmeV_var = "dzz_mort",  
                         data_type = c("count", "count"), 
                         data_manipulation = "original")
print(f2.1)
ggsave(filename = paste0(adt.dir, "Correlation-original-count.png"), 
       plot = f2.1, 
       width = 6, height = 4, units = "in", dpi = 300)





# (d) Figure 2.2 Correlation between measles and non-measles mortality rate --------

f2.2 <- correlation_plot(d.ts = d.bz, 
                         meV_var = "meV_mort_rate", nonmeV_var = "dzz_rate",  
                         data_type = c("rate", "rate"), 
                         data_manipulation = "original")
print(f2.2)
ggsave(filename = paste0(adt.dir, "Correlation-original-rate.png"), 
       plot = f2.2,
       width = 6, height = 4, units = "in", dpi = 300)




# (e) Figure 3.1: Time series of the detrended mortality count ------------

f3.1 <- time_series(d.ts = d.bz, 
                    meV_var = "meV_mort_dt", nonmeV_var = "dzz_mort_dt", 
                    cp = cp.bz.count.year, 
                    data_type = c("count", "count"), 
                    data_manipulation = "detrend")
print(f3.1)
ggsave(filename = paste0(adt.dir, "Time series-detrend-count.png"), 
       plot = f3.1,
       width = 7, height = 4, units = "in", dpi = 300)



# (f) Figure 3.2: Time series of the detrended mortality rate ------------

f3.2 <- time_series(d.ts = d.bz, 
                    meV_var = "meV_mort_rate_dt", nonmeV_var = "dzz_rate_dt", 
                    cp = cp.bz.rate.year, 
                    data_type = c("rate", "rate"), 
                    data_manipulation = "detrend")
print(f3.2)
ggsave(filename = paste0(adt.dir, "Time series-detrend-rate.png"), 
       plot = f3.2,
       width = 7, height = 4, units = "in", dpi = 300)




# (g) Figure 4.1 Correlation between detrended measles and non-measles mortality count --------

f4.1 <- correlation_plot(d.ts = d.bz, 
                         meV_var = "meV_mort_dt", nonmeV_var = "dzz_mort_dt",  
                         data_type = c("count", "count"), 
                         data_manipulation = "detrend")
print(f4.1)
ggsave(filename = paste0(adt.dir, "Correlation-detrend-count.png"), 
       plot = f4.1, 
       width = 6, height = 4, units = "in", dpi = 300)



# (h) Figure 4.2 Correlation between detrended measles and non-measles mortality rate --------

f4.2 <- correlation_plot(d.ts = d.bz, 
                         meV_var = "meV_mort_rate_dt", nonmeV_var = "dzz_rate_dt",  
                         data_type = c("rate", "rate"), 
                         data_manipulation = "detrend")
print(f4.2)
ggsave(filename = paste0(adt.dir, "Correlation-detrended-rate.png"), 
       plot = f4.2,
       width = 6, height = 4, units = "in", dpi = 300)






# 5. Linear regression ----------------------------------------------------

print(names(d.bz))

# (a) Full time series: original data ------------------------

# 1) Use only measles mortality as the predictor:
lm.bz.1.rate <- lm(dzz_rate ~ meV_mort_rate, data = d.bz)
print(lm.bz.1.rate)
summary(lm.bz.1.rate)

lm.bz.1.count <- lm(dzz_mort ~ meV_mort, data = d.bz)
print(lm.bz.1.count)
summary(lm.bz.1.count)


# 2) Use only year as the predictor:
lm.bz.2.rate <- lm(dzz_rate ~ year.idx, data = d.bz)
print(lm.bz.2.rate)
summary(lm.bz.2.rate)

lm.bz.2.count <- lm(dzz_mort ~ year.idx, data = d.bz)
print(lm.bz.2.count)
summary(lm.bz.2.count)


# 3) Add year and measles mortality as the predictor:
lm.bz.3.rate <- lm(dzz_rate ~ year.idx + meV_mort_rate, data = d.bz)
print(lm.bz.3.rate)
summary(lm.bz.3.rate)
anova(lm.bz.3.rate, lm.bz.1.rate)  # Effects of time
anova(lm.bz.3.rate, lm.bz.2.rate)  # Effects of measlse mortality

lm.bz.3.count <- lm(dzz_mort ~ year.idx + meV_mort, data = d.bz)
print(lm.bz.3.count)
summary(lm.bz.3.count)
anova(lm.bz.3.count, lm.bz.1.count)  # Effects of time
anova(lm.bz.3.count, lm.bz.2.count)  # Effects of measles mortality


# 4) Add interaction term between year and measles mortality:
lm.bz.4.rate <- lm(dzz_rate ~ year.idx * meV_mort_rate, data = d.bz)
print(lm.bz.4.rate)
summary(lm.bz.4.rate)
anova(lm.bz.4.rate, lm.bz.3.rate)  # Effects of the interaction

lm.bz.4.count <- lm(dzz_mort ~ year.idx * meV_mort, data = d.bz)
print(lm.bz.4.count)
summary(lm.bz.4.count)
anova(lm.bz.4.count, lm.bz.3.count)  # Effects of the interaction


# Put all regression results in a list
lm.bz.list.rate <- list(lm.bz.1.rate, lm.bz.2.rate, lm.bz.3.rate, lm.bz.4.rate)
lm.bz.list.count <- list(lm.bz.1.count, lm.bz.2.count, lm.bz.3.count, lm.bz.4.count)


# Model selection:
step.bz.rate <- step(lm.bz.4.rate, direction="both")
step.bz.count <- step(lm.bz.4.count, direction="both")

lm.bz.rate <- lm.bz.list.rate[[which.min(sapply(X = lm.bz.list.rate, FUN = AIC))]]
lm.bz.count <- lm.bz.list.count[[which.min(sapply(X = lm.bz.list.count, FUN = AIC))]]



# Examine the best model: rate
coefficients(lm.bz.rate)        # model coefficients
confint(lm.bz.rate, level=0.95) # CIs for model parameters
summary(lm.bz.rate)             # summary table
anova(lm.bz.rate)               # anova table


# Diagnostics of the best model:
plot(lm.bz.rate$residuals ~ lm.bz.rate$fitted.values)
qqnorm(lm.bz.rate$residuals); abline(a = 0, b = 1, col = "red")


# Relative importance of the predictors
relimp.bz.est <- calc.relimp(lm.bz.3.rate, type=c("lmg"), rela=TRUE)
print(relimp.bz.est)

relimp.bz.count.est <- calc.relimp(lm.bz.3.count, type=c("lmg"), rela=TRUE)
print(relimp.bz.count.est)


### write out the key results: mortality rate
sink(paste0(adt.dir, "Linear regression - full time series - rate - original.txt"))

cat("Location of interest: ")
cat(LOC)
cat("\nDisease of interest: ")
cat(DZZ)
cat("\nTime series: full\n\n\n")

cat("Linear regression models:\n")
stargazer(lm.bz.1.rate, lm.bz.2.rate, lm.bz.3.rate, lm.bz.4.rate,
          type = "text", style = "all", out = NULL,
          ci = TRUE, ci.level = 0.95,
          single.row = FALSE,
          column.labels = NULL,
          order = c(1, 3, 2),
          covariate.labels = c("Measles mortality rate", "Year", "Measles x Year"),
          dep.var.caption = "",
          dep.var.labels = "Dependent variable: Non-measles mortality rate",
          model.names = FALSE, model.numbers = TRUE, object.names = FALSE,
          omit = NULL, notes.align = "c")
cat("\n\n\n")

cat("Model comparision: effect of year\n\n")
print(anova(lm.bz.3.rate, lm.bz.1.rate))
cat("\n\n\n")

cat("Model comparision: effect of measlse mortality rate\n\n")
print(anova(lm.bz.3.rate, lm.bz.2.rate))
cat("\n\n\n")

cat("Model comparision: effect of the interactive term\n\n")
print(anova(lm.bz.4.rate, lm.bz.3.rate))
cat("\n\n\n")

cat("Step-wise model selection:\n\n")
step(lm.bz.4.rate, direction="both")
cat("\n\n\n")

cat("Relative importance of predictors:\n\n")
print(relimp.bz.est)
cat("\n")

sink()





# (b) Full time series: time-detrended data ------------------------

# 1) Use only measles mortality as the predictor:
lm.bz.dt.1.rate <- lm(dzz_rate_dt ~ meV_mort_rate_dt, data = d.bz)
print(lm.bz.dt.1.rate)
summary(lm.bz.dt.1.rate)

lm.bz.dt.1.count <- lm(dzz_mort_dt ~ meV_mort_dt, data = d.bz)
print(lm.bz.dt.1.count)
summary(lm.bz.dt.1.count)


# 2) Use only year as the predictor:
lm.bz.dt.2.rate <- lm(dzz_rate_dt ~ year.idx, data = d.bz)
print(lm.bz.dt.2.rate)
summary(lm.bz.dt.2.rate)

lm.bz.dt.2.count <- lm(dzz_mort_dt ~ year.idx, data = d.bz)
print(lm.bz.dt.2.count)
summary(lm.bz.dt.2.count)


# 3) Add year and measles mortality as the predictor:
lm.bz.dt.3.rate <- lm(dzz_rate_dt ~ year.idx + meV_mort_rate_dt, data = d.bz)
print(lm.bz.dt.3.rate)
summary(lm.bz.dt.3.rate)
anova(lm.bz.dt.3.rate, lm.bz.dt.1.rate)  # Effects of time
anova(lm.bz.dt.3.rate, lm.bz.dt.2.rate)  # Effects of measles mortality

lm.bz.dt.3.count <- lm(dzz_mort_dt ~ year.idx + meV_mort_dt, data = d.bz)
print(lm.bz.dt.3.count)
summary(lm.bz.dt.3.count)
anova(lm.bz.dt.3.count, lm.bz.dt.1.count)  # Effects of time
anova(lm.bz.dt.3.count, lm.bz.dt.2.count)  # Effects of measles mortality


# 4) Add interaction term between year and measles mortality:
lm.bz.dt.4.rate <- lm(dzz_rate_dt ~ year.idx * meV_mort_rate_dt, data = d.bz)
print(lm.bz.dt.4.rate)
summary(lm.bz.dt.4.rate)
anova(lm.bz.dt.4.rate, lm.bz.dt.3.rate)  # Effects of the interaction

lm.bz.dt.4.count <- lm(dzz_mort_dt ~ year.idx * meV_mort_dt, data = d.bz)
print(lm.bz.dt.4.count)
summary(lm.bz.dt.4.count)
anova(lm.bz.dt.4.count, lm.bz.dt.3.count)  # Effects of the interaction


# Put all regression results in a list
lm.bz.dt.list.rate <- list(lm.bz.dt.1.rate, lm.bz.dt.2.rate, lm.bz.dt.3.rate, lm.bz.dt.4.rate)
lm.bz.dt.list.count <- list(lm.bz.dt.1.count, lm.bz.dt.2.count, lm.bz.dt.3.count, lm.bz.dt.4.count)


# Model selection:
step.bz.dt.rate <- step(lm.bz.dt.4.rate, direction="both")
step.bz.dt.count <- step(lm.bz.dt.4.count, direction="both")

lm.bz.dt.rate <- lm.bz.dt.list.rate[[which.min(sapply(X = lm.bz.dt.list.rate, FUN = AIC))]]
lm.bz.dt.count <- lm.bz.dt.list.count[[which.min(sapply(X = lm.bz.dt.list.count, FUN = AIC))]]


# Examine the best model: rate
coefficients(lm.bz.dt.rate)        # model coefficients
confint(lm.bz.dt.rate, level=0.95) # CIs for model parameters
summary(lm.bz.dt.rate)             # summary table
anova(lm.bz.dt.rate)               # anova table


# Diagnostics of the best model:
plot(lm.bz.dt.rate$residuals ~ lm.bz.dt.rate$fitted.values)
qqnorm(lm.bz.dt.rate$residuals); abline(a = 0, b = 1, col = "red")


# Relative importance of the predictors
relimp.bz.dt.est <- calc.relimp(lm.bz.dt.3.rate, type=c("lmg"), rela=TRUE)
print(relimp.bz.dt.est)

relimp.bz.dt.count.est <- calc.relimp(lm.bz.dt.3.count, type=c("lmg"), rela=TRUE)
print(relimp.bz.dt.count.est)


# write out the key results:
sink(paste0(adt.dir, "Linear regression - full time series - rate - detrended.txt"))

cat("Location of interest: ")
cat(LOC)
cat("\nDisease of interest: ")
cat(DZZ)
cat("\nTime series: full\n\n\n")

cat("Linear regression models:\n")
stargazer(lm.bz.dt.1.rate, lm.bz.dt.2.rate, lm.bz.dt.3.rate, lm.bz.dt.4.rate,
          type = "text", style = "all", out = NULL,
          ci = TRUE, ci.level = 0.95,
          single.row = FALSE,
          column.labels = NULL,
          order = c(1, 3, 2),
          covariate.labels = c("Measles mortality rate", "Year", "Measles x Year"),
          dep.var.caption = "",
          dep.var.labels = "Dependent variable: Non-measles mortality rate",
          model.names = FALSE, model.numbers = TRUE, object.names = FALSE,
          omit = NULL, notes.align = "c")
cat("\n\n\n")

cat("Model comparision: effect of year\n\n")
print(anova(lm.bz.dt.3.rate, lm.bz.dt.1.rate))
cat("\n\n\n")

cat("Model comparision: effect of measlse mortality rate\n\n")
print(anova(lm.bz.dt.3.rate, lm.bz.dt.2.rate))
cat("\n\n\n")

cat("Model comparision: effect of the interactive term\n\n")
print(anova(lm.bz.dt.4.rate, lm.bz.dt.3.rate))
cat("\n\n\n")

cat("Step-wise model selection:\n\n")
step(lm.bz.dt.4.rate, direction="both")
cat("\n")

cat("Relative importance of predictors:\n\n")
print(relimp.bz.dt.est)
cat("\n")

sink()





# (c) Before the change-point: original data -----------------------
print(cp.bz.rate.year["mean"])
d.bz.t1 <- d.bz[d.bz$year < cp.bz.rate.year["mean"], ]
print(nrow(d.bz.t1))
range(d.bz.t1$year)


# 1) Use only measles mortality rate as the predictor
lm.bz.t1.1.rate <- lm(dzz_rate ~ meV_mort_rate, data = d.bz.t1)
print(lm.bz.t1.1.rate)
summary(lm.bz.t1.1.rate)

lm.bz.t1.1.count <- lm(dzz_mort ~ meV_mort, data = d.bz.t1)
print(lm.bz.t1.1.count)
summary(lm.bz.t1.1.count)


# 2) Use only year as the predictor:
lm.bz.t1.2.rate <- lm(dzz_rate ~ year.idx, data = d.bz.t1)
print(lm.bz.t1.2.rate)
summary(lm.bz.t1.2.rate)

lm.bz.t1.2.count <- lm(dzz_mort ~ year.idx, data = d.bz.t1)
print(lm.bz.t1.2.count)
summary(lm.bz.t1.2.count)


# 3) Add year and measles mortality as the predictor:
lm.bz.t1.3.rate <- lm(dzz_rate ~ year.idx + meV_mort_rate, data = d.bz.t1)
print(lm.bz.t1.3.rate)
summary(lm.bz.t1.3.rate)
anova(lm.bz.t1.3.rate, lm.bz.t1.1.rate)  # Effects of time
anova(lm.bz.t1.3.rate, lm.bz.t1.2.rate)  # Effects of measles mortality

lm.bz.t1.3.count <- lm(dzz_mort ~ year.idx + meV_mort, data = d.bz.t1)
print(lm.bz.t1.3.count)
summary(lm.bz.t1.3.count)
anova(lm.bz.t1.3.count, lm.bz.t1.1.count)  # Effects of time
anova(lm.bz.t1.3.count, lm.bz.t1.2.count)  # Effects of measles mortality


# 4) Add interaction term between year and measles mortality:
lm.bz.t1.4.rate <- lm(dzz_rate ~ year.idx * meV_mort_rate, data = d.bz.t1)
print(lm.bz.t1.4.rate)
summary(lm.bz.t1.4.rate)
anova(lm.bz.t1.4.rate, lm.bz.t1.3.rate)    # Effects of the interaction

lm.bz.t1.4.count <- lm(dzz_mort ~ year.idx * meV_mort, data = d.bz.t1)
print(lm.bz.t1.4.count)
summary(lm.bz.t1.4.count)
anova(lm.bz.t1.4.count, lm.bz.t1.3.count)  # Effects of the interaction


# Put all regression results in a list
lm.bz.t1.list.rate <- list(lm.bz.t1.1.rate, lm.bz.t1.2.rate, lm.bz.t1.3.rate, lm.bz.t1.4.rate)
lm.bz.t1.list.count <- list(lm.bz.t1.1.count, lm.bz.t1.2.count, lm.bz.t1.3.count, lm.bz.t1.4.count)

# Model selection:
if(nrow(d.bz.t1) >= 5){
  step.bz.t1.rate <- step(lm.bz.t1.4.rate, direction="both")
  step.bz.t1.count <- step(lm.bz.t1.4.count, direction="both")
}

lm.bz.t1.rate <- lm.bz.t1.list.rate[[which.min(sapply(X = lm.bz.t1.list.rate, FUN = AIC))]]
lm.bz.t1.count <- lm.bz.t1.list.count[[which.min(sapply(X = lm.bz.t1.list.count, FUN = AIC))]]


# Examine the best model: rate
coefficients(lm.bz.t1.rate)        # model coefficients
confint(lm.bz.t1.rate, level=0.95) # CIs for model parameters
summary(lm.bz.t1.rate)             # summary table
anova(lm.bz.t1.rate)               # anova table


# Diagnostics of the best model:
plot(lm.bz.t1.rate$residuals ~ lm.bz.t1.rate$fitted.values)
qqnorm(lm.bz.t1.rate$residuals); abline(a = 0, b = 1, col = "red")


# Relative importance of the predictors
relimp.bz.t1.est <- calc.relimp(lm.bz.t1.3.rate, type=c("lmg"), rela=TRUE)
print(relimp.bz.t1.est)

relimp.bz.count.t1.est <- calc.relimp(lm.bz.t1.3.count, type=c("lmg"), rela=TRUE)
print(relimp.bz.count.t1.est)


# write out the key results:
sink(paste0(adt.dir, "Linear regression - before change-point - rate - original.txt"))

cat("Location of interest: ")
cat(LOC)
cat("\nDisease of interest: ")
cat(DZZ)
cat("\nTime series: before change-point\n\n\n")

cat("Linear regression models:\n")
stargazer(lm.bz.t1.1.rate, lm.bz.t1.2.rate, lm.bz.t1.3.rate, lm.bz.t1.4.rate,
          type = "text", style = "all", out = NULL,
          ci = TRUE, ci.level = 0.95,
          single.row = FALSE,
          column.labels = NULL,
          order = c(1, 3, 2),
          covariate.labels = c("Measles mortality rate", "Year", "Measles x Year"),
          dep.var.caption = "",
          dep.var.labels = "Dependent variable: Non-measles mortality rate",
          model.names = FALSE, model.numbers = TRUE, object.names = FALSE,
          omit = NULL, notes.align = "c")
cat("\n\n\n")

cat("Model comparision: effect of year\n\n")
print(anova(lm.bz.t1.3.rate, lm.bz.t1.1.rate))
cat("\n\n\n")

cat("Model comparision: effect of measlse mortality rate\n\n")
print(anova(lm.bz.t1.3.rate, lm.bz.t1.2.rate))
cat("\n\n\n")

cat("Model comparision: effect of the interactive term\n\n")
print(anova(lm.bz.t1.4.rate, lm.bz.t1.3.rate))
cat("\n\n\n")

if(nrow(d.bz.t1) >= 5){
  cat("Step-wise model selection:\n\n")
  step(lm.bz.t1.4.rate, direction="both")
  cat("\n\n\n")
}

cat("Relative importance of predictors:\n\n")
print(relimp.bz.t1.est)
cat("\n")

sink()





# (d) Before the change-point: detrended data -----------------------

# 1) Use only measles mortality rate as the predictor
lm.bz.t1.dt.1.rate <- lm(dzz_rate_dt ~ meV_mort_rate_dt, data = d.bz.t1)
print(lm.bz.t1.dt.1.rate)
summary(lm.bz.t1.dt.1.rate)

lm.bz.t1.dt.1.count <- lm(dzz_mort_dt ~ meV_mort_dt, data = d.bz.t1)
print(lm.bz.t1.dt.1.count)
summary(lm.bz.t1.dt.1.count)


# 2) Use only year as the predictor:
lm.bz.t1.dt.2.rate <- lm(dzz_rate_dt ~ year.idx, data = d.bz.t1)
print(lm.bz.t1.dt.2.rate)
summary(lm.bz.t1.dt.2.rate)

lm.bz.t1.dt.2.count <- lm(dzz_mort_dt ~ year.idx, data = d.bz.t1)
print(lm.bz.t1.dt.2.count)
summary(lm.bz.t1.dt.2.count)


# 3) Add year and measles mortality as the predictor:
lm.bz.t1.dt.3.rate <- lm(dzz_rate_dt ~ year.idx + meV_mort_rate_dt, data = d.bz.t1)
print(lm.bz.t1.dt.3.rate)
summary(lm.bz.t1.dt.3.rate)
anova(lm.bz.t1.dt.3.rate, lm.bz.t1.dt.1.rate)  # Effects of time
anova(lm.bz.t1.dt.3.rate, lm.bz.t1.dt.2.rate)  # Effects of measles mortality

lm.bz.t1.dt.3.count <- lm(dzz_mort_dt ~ year.idx + meV_mort_dt, data = d.bz.t1)
print(lm.bz.t1.dt.3.count)
summary(lm.bz.t1.dt.3.count)
anova(lm.bz.t1.dt.3.count, lm.bz.t1.dt.1.count)  # Effects of time
anova(lm.bz.t1.dt.3.count, lm.bz.t1.dt.2.count)  # Effects of measles mortality


# 4) Add interaction term between year and measles mortality:
lm.bz.t1.dt.4.rate <- lm(dzz_rate_dt ~ year.idx * meV_mort_rate_dt, data = d.bz.t1)
print(lm.bz.t1.dt.4.rate)
summary(lm.bz.t1.dt.4.rate)
anova(lm.bz.t1.dt.4.rate, lm.bz.t1.dt.3.rate)  # Effects of the interaction

lm.bz.t1.dt.4.count <- lm(dzz_mort_dt ~ year.idx * meV_mort_dt, data = d.bz.t1)
print(lm.bz.t1.dt.4.count)
summary(lm.bz.t1.dt.4.count)
anova(lm.bz.t1.dt.4.count, lm.bz.t1.dt.3.count)  # Effects of the interaction


# Put all regression results in a list
lm.bz.t1.dt.list.rate <- list(lm.bz.t1.dt.1.rate, lm.bz.t1.dt.2.rate, lm.bz.t1.dt.3.rate, lm.bz.t1.dt.4.rate)
lm.bz.t1.dt.list.count <- list(lm.bz.t1.dt.1.count, lm.bz.t1.dt.2.count, lm.bz.t1.dt.3.count, lm.bz.t1.dt.4.count)


# Model selection:
if(nrow(d.bz.t1) >= 5){
  step.bz.t1.dt.rate <- step(lm.bz.t1.dt.4.rate, direction="both")
  step.bz.t1.dt.count <- step(lm.bz.t1.dt.4.count, direction="both")
}

lm.bz.t1.dt.rate <- lm.bz.t1.dt.list.rate[[which.min(sapply(X = lm.bz.t1.dt.list.rate, FUN = AIC))]]
lm.bz.t1.dt.count <- lm.bz.t1.dt.list.count[[which.min(sapply(X = lm.bz.t1.dt.list.count, FUN = AIC))]]


# Examine the best model: rate
coefficients(lm.bz.t1.dt.rate)        # model coefficients
confint(lm.bz.t1.dt.rate, level=0.95) # CIs for model parameters
summary(lm.bz.t1.dt.rate)             # summary table
anova(lm.bz.t1.dt.rate)               # anova table


# Diagnostics of the best model:
plot(lm.bz.t1.dt.rate$residuals ~ lm.bz.t1.dt.rate$fitted.values)
qqnorm(lm.bz.t1.dt.rate$residuals); abline(a = 0, b = 1, col = "red")


# Relative importance of the predictors
relimp.bz.t1.dt.est <- calc.relimp(lm.bz.t1.dt.3.rate, type=c("lmg"), rela=TRUE)
print(relimp.bz.t1.dt.est)

relimp.bz.count.t1.dt.est <- calc.relimp(lm.bz.t1.dt.3.count, type=c("lmg"), rela=TRUE)
print(relimp.bz.count.t1.dt.est)


# write out the key results:
sink(paste0(adt.dir, "Linear regression - before change-point - rate - detrended.txt"))

cat("Location of interest: ")
cat(LOC)
cat("\nDisease of interest: ")
cat(DZZ)
cat("\nTime series: before change-point\n\n\n")

cat("Linear regression models:\n")
stargazer(lm.bz.t1.dt.1.rate, lm.bz.t1.dt.2.rate, lm.bz.t1.dt.3.rate, lm.bz.t1.dt.4.rate,
          type = "text", style = "all", out = NULL,
          ci = TRUE, ci.level = 0.95,
          single.row = FALSE,
          column.labels = NULL,
          order = c(1, 3, 2),
          covariate.labels = c("Measles mortality rate", "Year", "Measles x Year"),
          dep.var.caption = "",
          dep.var.labels = "Dependent variable: Non-measles mortality rate",
          model.names = FALSE, model.numbers = TRUE, object.names = FALSE,
          omit = NULL, notes.align = "c")
cat("\n\n\n")

cat("Model comparision: effect of year\n\n")
print(anova(lm.bz.t1.dt.3.rate, lm.bz.t1.dt.1.rate))
cat("\n\n\n")

cat("Model comparision: effect of measlse mortality rate\n\n")
print(anova(lm.bz.t1.dt.3.rate, lm.bz.t1.dt.2.rate))
cat("\n\n\n")

cat("Model comparision: effect of the interactive term\n\n")
print(anova(lm.bz.t1.dt.4.rate, lm.bz.t1.dt.3.rate))
cat("\n\n\n")

if(nrow(d.bz.t1) >= 5){
  cat("Step-wise model selection:\n\n")
  step(lm.bz.t1.dt.4.rate, direction="both")
  cat("\n")
}

cat("Relative importance of predictors:\n\n")
print(relimp.bz.t1.dt.est)
cat("\n")

sink()





# (e) After the change-point: original data ----------------
d.bz.t2 <- d.bz[d.bz$year >= cp.bz.rate.year["mean"], ]
print(nrow(d.bz.t2))
range(d.bz.t2$year)


# 1) Use only measles mortality rate as the predictor
lm.bz.t2.1.rate <- lm(dzz_rate ~ meV_mort_rate, data = d.bz.t2)
print(lm.bz.t2.1.rate)
summary(lm.bz.t2.1.rate)

lm.bz.t2.1.count <- lm(dzz_mort ~ meV_mort, data = d.bz.t2)
print(lm.bz.t2.1.count)
summary(lm.bz.t2.1.count)


# 2) Use only year as the predictor:
lm.bz.t2.2.rate <- lm(dzz_rate ~ year.idx, data = d.bz.t2)
print(lm.bz.t2.2.rate)
summary(lm.bz.t2.2.rate)

lm.bz.t2.2.count <- lm(dzz_mort ~ year.idx, data = d.bz.t2)
print(lm.bz.t2.2.count)
summary(lm.bz.t2.2.count)


# 3) Add year and measles mortality as the predictor:
lm.bz.t2.3.rate <- lm(dzz_rate ~ year.idx + meV_mort_rate, data = d.bz.t2)
print(lm.bz.t2.3.rate)
summary(lm.bz.t2.3.rate)
anova(lm.bz.t2.3.rate, lm.bz.t2.1.rate)  # Effects of time
anova(lm.bz.t2.3.rate, lm.bz.t2.2.rate)  # Effects of measles mortality

lm.bz.t2.3.count <- lm(dzz_mort ~ year.idx + meV_mort, data = d.bz.t2)
print(lm.bz.t2.3.count)
summary(lm.bz.t2.3.count)
anova(lm.bz.t2.3.count, lm.bz.t2.1.count)  # Effects of time
anova(lm.bz.t2.3.count, lm.bz.t2.2.count)  # Effects of measles mortality


# 4) Add interaction term between year and measles mortality:
lm.bz.t2.4.rate <- lm(dzz_rate ~ year.idx * meV_mort_rate, data = d.bz.t2)
print(lm.bz.t2.4.rate)
summary(lm.bz.t2.4.rate)
anova(lm.bz.t2.4.rate, lm.bz.t2.3.rate)  # Effects of the interaction

lm.bz.t2.4.count <- lm(dzz_mort ~ year.idx * meV_mort, data = d.bz.t2)
print(lm.bz.t2.4.count)
summary(lm.bz.t2.4.count)
anova(lm.bz.t2.4.count, lm.bz.t2.3.count)  # Effects of the interaction


# Put all regression results in a list
lm.bz.t2.list.rate <- list(lm.bz.t2.1.rate, lm.bz.t2.2.rate, lm.bz.t2.3.rate, lm.bz.t2.4.rate)
lm.bz.t2.list.count <- list(lm.bz.t2.1.count, lm.bz.t2.2.count, lm.bz.t2.3.count, lm.bz.t2.4.count)


# Model selection:
if(nrow(d.bz.t2) >= 5){
  step.bz.t2.rate <- step(lm.bz.t2.4.rate, direction="both")
  step.bz.t2.count <- step(lm.bz.t2.4.count, direction="both")
}

lm.bz.t2.rate <- lm.bz.t2.list.rate[[which.min(sapply(X = lm.bz.t2.list.rate, FUN = AIC))]]
lm.bz.t2.count <- lm.bz.t2.list.count[[which.min(sapply(X = lm.bz.t2.list.count, FUN = AIC))]]


# Examine the best model: rate
coefficients(lm.bz.t2.rate)        # model coefficients
confint(lm.bz.t2.rate, level=0.95) # CIs for model parameters
summary(lm.bz.t2.rate)             # summary table
anova(lm.bz.t2.rate)               # anova table


# Diagnostics of the best model:
plot(lm.bz.t2.rate$residuals ~ lm.bz.t2.rate$fitted.values)
qqnorm(lm.bz.t2.rate$residuals); abline(a = 0, b = 1, col = "red")


# Relative importance of the predictors
relimp.bz.t2.est <- calc.relimp(lm.bz.t2.3.rate, type=c("lmg"), rela=TRUE)
print(relimp.bz.t2.est)

relimp.bz.count.t2.est <- calc.relimp(lm.bz.t2.3.count, type=c("lmg"), rela=TRUE)
print(relimp.bz.count.t2.est)


# write out the key results:
sink(paste0(adt.dir, "Linear regression - after change-point - rate - original.txt"))

cat("Location of interest: ")
cat(LOC)
cat("\nDisease of interest: ")
cat(DZZ)
cat("\nTime series: after change-point\n\n\n")

cat("Linear regression models:\n")
stargazer(lm.bz.t2.1.rate, lm.bz.t2.2.rate, lm.bz.t2.3.rate, lm.bz.t2.4.rate,
          type = "text", style = "all", out = NULL,
          ci = TRUE, ci.level = 0.95,
          single.row = FALSE,
          column.labels = NULL,
          order = c(1, 3, 2),
          covariate.labels = c("Measles mortality rate", "Year", "Measles x Year"),
          dep.var.caption = "",
          dep.var.labels = "Dependent variable: Non-measles mortality rate",
          model.names = FALSE, model.numbers = TRUE, object.names = FALSE,
          omit = NULL, notes.align = "c")
cat("\n\n\n")

cat("Model comparision: effect of year\n\n")
print(anova(lm.bz.t2.3.rate, lm.bz.t2.1.rate))
cat("\n\n\n")

cat("Model comparision: effect of measlse mortality rate\n\n")
print(anova(lm.bz.t2.3.rate, lm.bz.t2.2.rate))
cat("\n\n\n")

cat("Model comparision: effect of the interactive term\n\n")
print(anova(lm.bz.t2.4.rate, lm.bz.t2.3.rate))
cat("\n\n\n")

if(nrow(d.bz.t2) >= 5){
  cat("Step-wise model selection:\n\n")
  step(lm.bz.t2.4.rate, direction="both")
  cat("\n\n\n")
}

cat("Relative importance of predictors:\n\n")
print(relimp.bz.t2.est)
cat("\n")

sink()




# (f) After the change-point: detrended data ----------------

# 1) Use only measles mortality rate as the predictor
lm.bz.t2.dt.1.rate <- lm(dzz_rate_dt ~ meV_mort_rate_dt, data = d.bz.t2)
print(lm.bz.t2.dt.1.rate)
summary(lm.bz.t2.dt.1.rate)

lm.bz.t2.dt.1.count <- lm(dzz_mort_dt ~ meV_mort_dt, data = d.bz.t2)
print(lm.bz.t2.dt.1.count)
summary(lm.bz.t2.dt.1.count)


# 2) Use only year as the predictor:
lm.bz.t2.dt.2.rate <- lm(dzz_rate_dt ~ year.idx, data = d.bz.t2)
print(lm.bz.t2.dt.2.rate)
summary(lm.bz.t2.dt.2.rate)

lm.bz.t2.dt.2.count <- lm(dzz_mort_dt ~ year.idx, data = d.bz.t2)
print(lm.bz.t2.dt.2.count)
summary(lm.bz.t2.dt.2.count)


# 3) Add year and measles mortality as the predictor:
lm.bz.t2.dt.3.rate <- lm(dzz_rate_dt ~ year.idx + meV_mort_rate_dt, data = d.bz.t2)
print(lm.bz.t2.dt.3.rate)
summary(lm.bz.t2.dt.3.rate)
anova(lm.bz.t2.dt.3.rate, lm.bz.t2.dt.1.rate)  # Effects of time
anova(lm.bz.t2.dt.3.rate, lm.bz.t2.dt.2.rate)  # Effects of measles mortality

lm.bz.t2.dt.3.count <- lm(dzz_mort_dt ~ year.idx + meV_mort_dt, data = d.bz.t2)
print(lm.bz.t2.dt.3.count)
summary(lm.bz.t2.dt.3.count)
anova(lm.bz.t2.dt.3.count, lm.bz.t2.dt.1.count)  # Effects of time
anova(lm.bz.t2.dt.3.count, lm.bz.t2.dt.2.count)  # Effects of measles mortality


# 4) Add interaction term between year and measles mortality:
lm.bz.t2.dt.4.rate <- lm(dzz_rate_dt ~ year.idx * meV_mort_rate_dt, data = d.bz.t2)
print(lm.bz.t2.dt.4.rate)
summary(lm.bz.t2.dt.4.rate)
anova(lm.bz.t2.dt.4.rate, lm.bz.t2.dt.3.rate)  # Effects of the interaction

lm.bz.t2.dt.4.count <- lm(dzz_mort_dt ~ year.idx * meV_mort_dt, data = d.bz.t2)
print(lm.bz.t2.dt.4.count)
summary(lm.bz.t2.dt.4.count)
anova(lm.bz.t2.dt.4.count, lm.bz.t2.dt.3.count)  # Effects of the interaction


# Put all regression results in a list
lm.bz.t2.dt.list.rate <- list(lm.bz.t2.dt.1.rate, lm.bz.t2.dt.2.rate, lm.bz.t2.dt.3.rate, lm.bz.t2.dt.4.rate)
lm.bz.t2.dt.list.count <- list(lm.bz.t2.dt.1.count, lm.bz.t2.dt.2.count, lm.bz.t2.dt.3.count, lm.bz.t2.dt.4.count)


# Model selection:
if(nrow(d.bz.t2) >= 5){
  step.bz.t2.dt.rate <- step(lm.bz.t2.dt.4.rate, direction="both")
  step.bz.t2.dt.count <- step(lm.bz.t2.dt.4.count, direction="both")
}

lm.bz.t2.dt.rate <- lm.bz.t2.dt.list.rate[[which.min(sapply(X = lm.bz.t2.dt.list.rate, FUN = AIC))]]
lm.bz.t2.dt.count <- lm.bz.t2.dt.list.count[[which.min(sapply(X = lm.bz.t2.dt.list.count, FUN = AIC))]]


# Examine the best model: rate
coefficients(lm.bz.t2.dt.rate)        # model coefficients
confint(lm.bz.t2.dt.rate, level=0.95) # CIs for model parameters
summary(lm.bz.t2.dt.rate)             # summary table
anova(lm.bz.t2.dt.rate)               # anova table


# Diagnostics of the best model:
plot(lm.bz.t2.dt.rate$residuals ~ lm.bz.t2.dt.rate$fitted.values)
qqnorm(lm.bz.t2.dt.rate$residuals); abline(a = 0, b = 1, col = "red")


# Relative importance of the predictors
relimp.bz.t2.dt.est <- calc.relimp(lm.bz.t2.dt.3.rate, type=c("lmg"), rela=TRUE)
print(relimp.bz.t2.dt.est)

relimp.bz.count.t2.dt.est <- calc.relimp(lm.bz.t2.dt.3.count, type=c("lmg"), rela=TRUE)
print(relimp.bz.count.t2.dt.est)


# write out the key results:
sink(paste0(adt.dir, "Linear regression - after change-point - rate - detrended.txt"))

cat("Location of interest: ")
cat(LOC)
cat("\nDisease of interest: ")
cat(DZZ)
cat("\nTime series: after change-point\n\n\n")

cat("Linear regression models:\n")
stargazer(lm.bz.t2.dt.1.rate, lm.bz.t2.dt.2.rate, lm.bz.t2.dt.3.rate, lm.bz.t2.dt.4.rate,
          type = "text", style = "all", out = NULL,
          ci = TRUE, ci.level = 0.95,
          single.row = FALSE,
          column.labels = NULL,
          order = c(1, 3, 2),
          covariate.labels = c("Measles mortality rate", "Year", "Measles x Year"),
          dep.var.caption = "",
          dep.var.labels = "Dependent variable: Non-measles mortality rate",
          model.names = FALSE, model.numbers = TRUE, object.names = FALSE,
          omit = NULL, notes.align = "c")
cat("\n\n\n")

cat("Model comparision: effect of year\n\n")
print(anova(lm.bz.t2.dt.3.rate, lm.bz.t2.dt.1.rate))
cat("\n\n\n")

cat("Model comparision: effect of measlse mortality rate\n\n")
print(anova(lm.bz.t2.dt.3.rate, lm.bz.t2.dt.2.rate))
cat("\n\n\n")

cat("Model comparision: effect of the interactive term\n\n")
print(anova(lm.bz.t2.dt.4.rate, lm.bz.t2.dt.3.rate))
cat("\n\n\n")

if(nrow(d.bz.t2) >= 5){
  cat("Step-wise model selection:\n\n")
  step(lm.bz.t2.dt.4.rate, direction="both")
  cat("\n")
}

cat("Relative importance of predictors:\n\n")
print(relimp.bz.t2.dt.est)
cat("\n")

sink()





# 6. Detrend non-measles mortality by measles mortality ----------------------------------------

# Non-measles mortality count:
d.bz$dzz_mort_dmev <- d.bz$dzz_mort - lm.bz.1.count$fitted.values

# Non-measles mortality rate:
d.bz$dzz_rate_dmev <- d.bz$dzz_rate - lm.bz.1.rate$fitted.values



# (a) Time series of the mortality count after removing measles  --------------

lm.bz.count.dmev <- lm(dzz_mort_dmev ~ year, data = d.bz)

# R-squared and the coefficient:
R2.dmev <- round(summary(lm.bz.count.dmev)$r.squared, digits = 3)
coeff.dmev <- round(lm.bz.count.dmev$coefficients[2], digits = 3)

# Model fit and confidence interval:
ci.dmev <- predict(lm.bz.count.dmev, 
                   newdata = data.frame(year = d.bz$year), 
                   interval="confidence", level = 0.95)
ci.dmev <- as.data.frame(ci.dmev)
ci.dmev$year <- d.bz$year

# Make the figure
plot.dmev.count <- ggplot() + 
  geom_ribbon(data = ci.dmev, aes(x = year, ymin = lwr, ymax = upr), fill = "gray90") + 
  geom_line(data = ci.dmev, aes(x = year, y = fit), size = 1) + 
  geom_path(data = d.bz, aes(x = year, y = dzz_mort_dmev),
            size = 0.5, linetype = "dashed", color = "gray") +
  geom_point(data = d.bz, aes(x = year, y = dzz_mort_dmev, color = year),
             size = 2, shape = 19) +
  xlab(label = "Year") + 
  ylab(label = "Non-measles mortality") + 
  scale_color_gradient(name = "Year", 
                       breaks = seq(from = 1980, to = 1995, by = 5), 
                       limits = c(1980, 1995),
                       low = cbPalette[2], high = cbPalette[6]) + 
  annotate(geom = "text", 
           x = position_in_range(x = d.bz$year, p = 1),
           y = position_in_range(x = c(ci.dmev$lwr, ci.dmev$upr, d.bz$dzz_mort_dmev), p = 0.05), 
           label = bquote(paste(R^2 ,"= ", .(R2.dmev))),
           size = 4, hjust = 1) + 
  annotate(geom = "text", 
           x = position_in_range(x = d.bz$year, p = 1),
           y = position_in_range(x = c(ci.dmev$lwr, ci.dmev$upr, d.bz$dzz_mort_dmev), p = 0.12), 
           label = paste0("coeff = ", coeff.dmev),
           size = 4, hjust = 1) + 
  theme_classic() + 
  theme(axis.title = element_text(size = 12),
        axis.text = element_text(size = 12),
        legend.text = element_text(size = 10))
print(plot.dmev.count)
ggsave(filename = paste0(adt.dir, "Detrended by measles - count.png"), 
       plot = plot.dmev.count, 
       width = 6, height = 4, units = "in", dpi = 300)





# (b) Time series of the mortality rate after removing measles eff --------

lm.bz.rate.dmev <- lm(dzz_rate_dmev ~ year, data = d.bz)

# R-squared and the coefficient:
R2.dmev <- round(summary(lm.bz.rate.dmev)$r.squared, digits = 3)
coeff.dmev <- round(lm.bz.rate.dmev$coefficients[2], digits = 3)

# Model fit and confidence interval:
ci.dmev <- predict(lm.bz.rate.dmev, 
                   newdata = data.frame(year = d.bz$year), 
                   interval="confidence", level = 0.95)
ci.dmev <- as.data.frame(ci.dmev)
ci.dmev$year <- d.bz$year

# Make the figure
plot.dmev.rate <- ggplot() + 
  geom_ribbon(data = ci.dmev, aes(x = year, ymin = lwr, ymax = upr), fill = "gray90") + 
  geom_line(data = ci.dmev, aes(x = year, y = fit), size = 1) + 
  geom_path(data = d.bz, aes(x = year, y = dzz_rate_dmev),
            size = 0.5, linetype = "dashed", color = "gray") +
  geom_point(data = d.bz, aes(x = year, y = dzz_rate_dmev, color = year),
             size = 2, shape = 19) +
  xlab(label = "Year") + 
  ylab(label = "Non-measles mortality rate") + 
  scale_color_gradient(name = "Year", 
                       breaks = seq(from = 1980, to = 1995, by = 5), 
                       limits = c(1980, 1995),
                       low = cbPalette[2], high = cbPalette[6]) + 
  annotate(geom = "text", 
           x = position_in_range(x = d.bz$year, p = 1),
           y = position_in_range(x = c(ci.dmev$lwr, ci.dmev$upr, d.bz$dzz_rate_dmev), p = 0.05), 
           label = bquote(paste(R^2 ,"= ", .(R2.dmev))),
           size = 4, hjust = 1) + 
  annotate(geom = "text", 
           x = position_in_range(x = d.bz$year, p = 1),
           y = position_in_range(x = c(ci.dmev$lwr, ci.dmev$upr, d.bz$dzz_rate_dmev), p = 0.12), 
           label = paste0("coeff = ", coeff.dmev),
           size = 4, hjust = 1) + 
  theme_classic() + 
  theme(axis.title = element_text(size = 12),
        axis.text = element_text(size = 12),
        legend.text = element_text(size = 10))
print(plot.dmev.rate)
ggsave(filename = paste0(adt.dir, "Detrended by measles - rate.png"), 
       plot = plot.dmev.rate, 
       width = 6, height = 4, units = "in", dpi = 300)




# (c) Linear regression ---------------------------------------------------

# 1) Use only measles mortality as the predictor:
lm.bz.dmev.1.rate <- lm(dzz_rate_dmev ~ meV_mort_rate, data = d.bz)
print(lm.bz.dmev.1.rate)
summary(lm.bz.dmev.1.rate)

lm.bz.dmev.1.count <- lm(dzz_mort_dmev ~ meV_mort, data = d.bz)
print(lm.bz.dmev.1.count)
summary(lm.bz.dmev.1.count)


# 2) Use only year as the predictor:
lm.bz.dmev.2.rate <- lm(dzz_rate_dmev ~ year.idx, data = d.bz)
print(lm.bz.dmev.2.rate)
summary(lm.bz.dmev.2.rate)

lm.bz.dmev.2.count <- lm(dzz_mort_dmev ~ year.idx, data = d.bz)
print(lm.bz.dmev.2.count)
summary(lm.bz.dmev.2.count)


# 3) Add year and measles mortality as the predictor:
lm.bz.dmev.3.rate <- lm(dzz_rate_dmev ~ year.idx + meV_mort_rate, data = d.bz)
print(lm.bz.dmev.3.rate)
summary(lm.bz.dmev.3.rate)
anova(lm.bz.dmev.3.rate, lm.bz.dmev.1.rate)
anova(lm.bz.dmev.3.rate, lm.bz.dmev.2.rate)

lm.bz.dmev.3.count <- lm(dzz_mort_dmev ~ year.idx + meV_mort, data = d.bz)
print(lm.bz.dmev.3.count)
summary(lm.bz.dmev.3.count)
anova(lm.bz.dmev.3.count, lm.bz.dmev.1.count)
anova(lm.bz.dmev.3.count, lm.bz.dmev.2.count)


# 4) Add interaction term between year and measles mortality:
lm.bz.dmev.4.rate <- lm(dzz_rate_dmev ~ year.idx * meV_mort_rate, data = d.bz)
print(lm.bz.dmev.4.rate)
summary(lm.bz.dmev.4.rate)
anova(lm.bz.dmev.4.rate, lm.bz.dmev.3.rate)

lm.bz.dmev.4.count <- lm(dzz_mort_dmev ~ year.idx * meV_mort, data = d.bz)
print(lm.bz.dmev.4.count)
summary(lm.bz.dmev.4.count)
anova(lm.bz.dmev.4.count, lm.bz.dmev.3.count)


# Put all regression results in a list
lm.bz.dmev.list.rate <- list(lm.bz.dmev.1.rate, lm.bz.dmev.2.rate, lm.bz.dmev.3.rate, lm.bz.dmev.4.rate)
lm.bz.dmev.list.count <- list(lm.bz.dmev.1.count, lm.bz.dmev.2.count, lm.bz.dmev.3.count, lm.bz.dmev.4.count)


# Model selection:
step.bz.dmev.rate <- step(lm.bz.dmev.4.rate, direction="both")
step.bz.dmev.count <- step(lm.bz.dmev.4.count, direction="both")

lm.bz.dmev.rate <- lm.bz.dmev.list.rate[[which.min(sapply(X = lm.bz.dmev.list.rate, FUN = AIC))]]
lm.bz.dmev.count <- lm.bz.dmev.list.count[[which.min(sapply(X = lm.bz.dmev.list.count, FUN = AIC))]]


# Results of the best model
coefficients(lm.bz.dmev.rate)         # model coefficients
confint(lm.bz.dmev.rate, level=0.95)  # CIs for model parameters
summary(lm.bz.dmev.rate)              # summary table
anova(lm.bz.dmev.rate)                # anova table


# Diagnostics of the best model:
plot(lm.bz.dmev.rate$residuals ~ lm.bz.dmev.rate$fitted.values)
qqnorm(lm.bz.dmev.rate$residuals); abline(a = 0, b = 1, col = "red")


# Relative importance of the predictors
relimp.bz.dmev.est <- calc.relimp(lm.bz.dmev.3.rate, type=c("lmg"), rela=TRUE)
print(relimp.bz.dmev.est)

relimp.bz.count.dmev.est <- calc.relimp(lm.bz.dmev.3.count, type=c("lmg"), rela=TRUE)
print(relimp.bz.count.dmev.est)



# write out the key results:
sink(paste0(adt.dir, "Linear regression - detrended by measles - rate.txt"))

cat("Location of interest: ")
cat(LOC)
cat("\nDisease of interest: ")
cat(DZZ)
cat("\nTime series: full\n\n\n")

cat("Linear regression models:\n")
stargazer(lm.bz.dmev.1.rate, lm.bz.dmev.2.rate, lm.bz.dmev.3.rate, lm.bz.dmev.4.rate,
          type = "text", style = "all", out = NULL,
          ci = TRUE, ci.level = 0.95,
          single.row = FALSE,
          column.labels = NULL,
          order = c(1, 3, 2),
          covariate.labels = c("Measles mortality rate", "Year", "Measles x Year"),
          dep.var.caption = "",
          dep.var.labels = "Dependent variable: Non-measles mortality rate",
          model.names = FALSE, model.numbers = TRUE, object.names = FALSE,
          omit = NULL, notes.align = "c")
cat("\n\n\n")

cat("Model comparision: effect of year\n\n")
print(anova(lm.bz.dmev.3.rate, lm.bz.dmev.1.rate))
cat("\n\n\n")

cat("Model comparision: effect of measlse mortality rate\n\n")
print(anova(lm.bz.dmev.3.rate, lm.bz.dmev.2.rate))
cat("\n\n\n")

cat("Model comparision: effect of the interactive term\n\n")
print(anova(lm.bz.dmev.4.rate, lm.bz.dmev.3.rate))
cat("\n\n\n")

cat("Step-wise model selection:\n\n")
step(lm.bz.dmev.4.rate, direction="both")
cat("\n")

cat("Relative importance of predictors:\n\n")
print(relimp.bz.dmev.est)
cat("\n")

sink()





# 7. Difference between consective years -----------------------------------

d.bz.matrix <- data.matrix(d.bz[, !(names(d.bz) %in% c("age", "city", "dzzName"))])
diff.bz <- diff(d.bz.matrix)
diff.bz <- as.data.frame(diff.bz, stringsAsFactors = FALSE)

diff.bz$start_year <- d.bz$year[1:nrow(diff.bz)]
diff.bz$end_year <- d.bz$year[2:nrow(d.bz)]
diff.bz$year <- diff.bz$end_year   # Use the latter year in the two consecutive years



# (a) Time series: mortality count ----------------------------------------

plot.diff.bz.ts.count <- time_series(d.ts = diff.bz, 
                                     meV_var = "meV_mort", nonmeV_var = "dzz_mort", 
                                     cp = cp.bz.count.year, 
                                     data_type = c("count", "count"), 
                                     data_manipulation = "difference", 
                                     year_breaks = seq(1980, 1995, 5))
print(plot.diff.bz.ts.count)
ggsave(filename = paste0(adt.dir, "Time series-difference-count.png"), 
       plot = plot.diff.bz.ts.count,
       width = 7, height = 4, units = "in", dpi = 300)



# (b) Time series: mortality rate ----------------------------------------

plot.diff.bz.ts.rate <- time_series(d.ts = diff.bz, 
                                    meV_var = "meV_mort_rate", nonmeV_var = "dzz_rate", 
                                    cp = cp.bz.rate.year, 
                                    data_type = c("rate", "rate"), 
                                    data_manipulation = "difference", 
                                    year_breaks = seq(1980, 1995, 5))
print(plot.diff.bz.ts.rate)
ggsave(filename = paste0(adt.dir, "Time series-difference-rate.png"), 
       plot = plot.diff.bz.ts.rate,
       width = 7, height = 4, units = "in", dpi = 300)





# (c) Visualization of correlation: mortality count -----------------------

plot.diff.bz.cor.count <- correlation_plot(d.ts = diff.bz, 
                                           meV_var = "meV_mort", nonmeV_var = "dzz_mort",  
                                           data_type = c("count", "count"), 
                                           data_manipulation = "difference", 
                                           year_breaks = seq(1980, 1995, by = 5), 
                                           year_range = c(1980, 1995))
print(plot.diff.bz.cor.count)
ggsave(filename = paste0(adt.dir, "Difference-correlation-count.png"), 
       plot = plot.diff.bz.cor.count,
       width = 6, height = 4, units = "in", dpi = 300)




# (d) Visualization of correlation: mortality rate -----------------------

plot.diff.bz.cor.rate <- correlation_plot(d.ts = diff.bz, 
                                          meV_var = "meV_mort_rate", nonmeV_var = "dzz_rate",  
                                          data_type = c("rate", "rate"), 
                                          data_manipulation = "difference", 
                                          year_breaks = seq(1980, 1995, by = 5), 
                                          year_range = c(1980, 1995))
print(plot.diff.bz.cor.rate)
ggsave(filename = paste0(adt.dir, "Difference-correlation-rate.png"), 
       plot = plot.diff.bz.cor.rate,
       width = 6, height = 4, units = "in", dpi = 300)





# (e) Linear regression ---------------------------------------------------

diff.bz$year.idx <- diff.bz$year - min(diff.bz$year) + 1

# 1) Use only measles mortality as the predictor:
lm.bz.diff.1.rate <- lm(dzz_rate ~ meV_mort_rate, data = diff.bz)
print(lm.bz.diff.1.rate)
summary(lm.bz.diff.1.rate)

lm.bz.diff.1.count <- lm(dzz_mort ~ meV_mort, data = diff.bz)
print(lm.bz.diff.1.count)
summary(lm.bz.diff.1.count)


# 2) Use only year as the predictor:
lm.bz.diff.2.rate <- lm(dzz_rate ~ year.idx, data = diff.bz)
print(lm.bz.diff.2.rate)
summary(lm.bz.diff.2.rate)

lm.bz.diff.2.count <- lm(dzz_mort ~ year.idx, data = diff.bz)
print(lm.bz.diff.2.count)
summary(lm.bz.diff.2.count)


# 3) Add year and measles mortality as the predictor:
lm.bz.diff.3.rate <- lm(dzz_rate ~ year.idx + meV_mort_rate, data = diff.bz)
print(lm.bz.diff.3.rate)
summary(lm.bz.diff.3.rate)
anova(lm.bz.diff.3.rate, lm.bz.diff.1.rate)  # Effects of year
anova(lm.bz.diff.3.rate, lm.bz.diff.2.rate)  # Effects of measles mortality

lm.bz.diff.3.count <- lm(dzz_mort ~ year.idx + meV_mort, data = diff.bz)
print(lm.bz.diff.3.count)
summary(lm.bz.diff.3.count)
anova(lm.bz.diff.3.count, lm.bz.diff.1.count)  # Effects of year
anova(lm.bz.diff.3.count, lm.bz.diff.2.count)  # Effects of measles mortality


# 4) Add interaction term between year and measles mortality:
lm.bz.diff.4.rate <- lm(dzz_rate ~ year.idx * meV_mort_rate, data = diff.bz)
print(lm.bz.diff.4.rate)
summary(lm.bz.diff.4.rate)
anova(lm.bz.diff.4.rate, lm.bz.diff.3.rate)  # Effects of the interaction

lm.bz.diff.4.count <- lm(dzz_mort ~ year.idx * meV_mort, data = diff.bz)
print(lm.bz.diff.4.count)
summary(lm.bz.diff.4.count)
anova(lm.bz.diff.4.count, lm.bz.diff.3.count)  # Effects of the interaction


# Put all regression results in a list
lm.bz.diff.list.rate <- list(lm.bz.diff.1.rate, lm.bz.diff.2.rate, lm.bz.diff.3.rate, lm.bz.diff.4.rate)
lm.bz.diff.list.count <- list(lm.bz.diff.1.count, lm.bz.diff.2.count, lm.bz.diff.3.count, lm.bz.diff.4.count)


# Model selection:
step.bz.diff.rate <- step(lm.bz.diff.4.rate, direction="both")
step.bz.diff.count <- step(lm.bz.diff.4.count, direction="both")

lm.bz.diff.rate <- lm.bz.diff.list.rate[[which.min(sapply(X = lm.bz.diff.list.rate, FUN = AIC))]]
lm.bz.diff.count <- lm.bz.diff.list.count[[which.min(sapply(X = lm.bz.diff.list.count, FUN = AIC))]]


# Results of the best model
coefficients(lm.bz.diff.rate)         # model coefficients
confint(lm.bz.diff.rate, level=0.95)  # CIs for model parameters
summary(lm.bz.diff.rate)              # summary table
anova(lm.bz.diff.rate)                # anova table


# Diagnostics of the best model:
plot(lm.bz.diff.rate$residuals ~ lm.bz.diff.rate$fitted.values)
qqnorm(lm.bz.diff.rate$residuals); abline(a = 0, b = 1, col = "red")


# Relative importance of the predictors
relimp.bz.diff.est <- calc.relimp(lm.bz.diff.3.rate, type=c("lmg"), rela=TRUE)
print(relimp.bz.diff.est)

relimp.bz.diff.count.est <- calc.relimp(lm.bz.diff.3.count, type=c("lmg"), rela=TRUE)
print(relimp.bz.diff.count.est)


# write out the key results:
sink(paste0(adt.dir, "Linear regression - full time series - rate - difference.txt"))

cat("Location of interest: ")
cat(LOC)
cat("\nDisease of interest: ")
cat(DZZ)
cat("\nTime series: full")
cat("\nDifference between consective years\n\n\n")

cat("Linear regression models:\n")
stargazer(lm.bz.diff.1.rate, lm.bz.diff.2.rate, lm.bz.diff.3.rate, lm.bz.diff.4.rate,
          type = "text", style = "all", out = NULL,
          ci = TRUE, ci.level = 0.95,
          single.row = FALSE,
          column.labels = NULL,
          order = c(1, 3, 2),
          covariate.labels = c("Measles mortality rate", "Year", "Measles x Year"),
          dep.var.caption = "",
          dep.var.labels = "Dependent variable: Non-measles mortality rate",
          model.names = FALSE, model.numbers = TRUE, object.names = FALSE,
          omit = NULL, notes.align = "c")
cat("\n\n\n")

cat("Model comparision: effect of year\n\n")
print(anova(lm.bz.diff.3.rate, lm.bz.diff.1.rate))
cat("\n\n\n")

cat("Model comparision: effect of measlse mortality rate\n\n")
print(anova(lm.bz.diff.3.rate, lm.bz.diff.2.rate))
cat("\n\n\n")

cat("Model comparision: effect of the interactive term\n\n")
print(anova(lm.bz.diff.4.rate, lm.bz.diff.3.rate))
cat("\n\n\n")

cat("Step-wise model selection:\n\n")
step(lm.bz.diff.4.rate, direction="both")
cat("\n")

cat("Relative importance of predictors:\n\n")
print(relimp.bz.diff.est)
cat("\n")

sink()





# 8. Negative binomial regression ----------------------------------------------------

# NB model takes count data and use offset() to account for population size variations

# Examine if the data is over-dispersed
mean(d.bz$dzz_mort)
var(d.bz$dzz_mort)   # variance is much larger than the mean


# Model 1: Only the measles mortality count as the input
nbglm.bz.1 <- glm.nb(dzz_mort ~ meV_mort_rate, data = d.bz)
summary(nbglm.bz.1)


# Model 2: Add a offset to control for the change of population size
nbglm.bz.2 <- glm.nb(dzz_mort ~ meV_mort_rate + offset(log(pop / 100000)), data = d.bz)
summary(nbglm.bz.2)


# Model 3: Only the year as the input
nbglm.bz.3 <- glm.nb(dzz_mort ~ year.idx, data = d.bz)
summary(nbglm.bz.3)


# Model 4: Add a offset to control for the change of population size
nbglm.bz.4 <- glm.nb(dzz_mort ~ year.idx + offset(log(pop / 100000)), data = d.bz)
summary(nbglm.bz.4)


# Model 5: Both year and measles mortality as predictors
nbglm.bz.5 <- glm.nb(dzz_mort ~ year.idx + meV_mort_rate, data = d.bz)
summary(nbglm.bz.5)


# Model 6: Add a offset
nbglm.bz.6 <- glm.nb(dzz_mort ~ year.idx + meV_mort_rate + offset(log(pop / 100000)), data = d.bz)
summary(nbglm.bz.6)


# Model 7: Add the interaction between year and measles mortality
nbglm.bz.7 <- glm.nb(dzz_mort ~ year.idx * meV_mort_rate, data = d.bz)
summary(nbglm.bz.7)


# Model 8: Add a offset
nbglm.bz.8 <- glm.nb(dzz_mort ~ year.idx * meV_mort_rate + offset(log(pop / 100000)), data = d.bz)
summary(nbglm.bz.8)


# Examine the effects of each variable:
anova(nbglm.bz.1, nbglm.bz.5)  # Effects of year
anova(nbglm.bz.2, nbglm.bz.6)
anova(nbglm.bz.3, nbglm.bz.5)  # Effects of measles mortality
anova(nbglm.bz.4, nbglm.bz.6)
anova(nbglm.bz.5, nbglm.bz.7)  # Effects of the interaction
anova(nbglm.bz.6, nbglm.bz.8)


# Model comparison:
nbglm.bz.aic <- c(model.1 = nbglm.bz.1$aic, 
                  model.2 = nbglm.bz.2$aic, 
                  model.3 = nbglm.bz.3$aic,
                  model.4 = nbglm.bz.4$aic, 
                  model.5 = nbglm.bz.5$aic, 
                  model.6 = nbglm.bz.6$aic,
                  model.7 = nbglm.bz.7$aic, 
                  model.8 = nbglm.bz.8$aic)
print(nbglm.bz.aic)  # The best model has the lowest AIC


# Find the best model:
nbglm.bz.list <- list(nbglm.bz.1, nbglm.bz.2, nbglm.bz.3, nbglm.bz.4,
                      nbglm.bz.5, nbglm.bz.6, nbglm.bz.7, nbglm.bz.8)
nbglm.bz.best <- nbglm.bz.list[[which.min(nbglm.bz.aic)]]
summary(nbglm.bz.best)


# Stepwise model selection:
nbglm.bz.step <- step(nbglm.bz.8, direction = "both")


# Some further exploration of the best model:
summary(nbglm.bz.best) # summary
anova(nbglm.bz.best)   # anova table


### 95% confidence interval of the NB model coefficient of the full model:
nbglm.bz.best.cofficient <- cbind(estimate = coef(nbglm.bz.best), confint(nbglm.bz.best))
print(nbglm.bz.best.cofficient)

# Convert back to the empirical scale: multiplicable
nbglm.bz.best.cofficient.ir <- exp(nbglm.bz.best.cofficient)
print(nbglm.bz.best.cofficient.ir)


### Calculate the pseudo-r2 using different methods:
nbglm.bz.pseudo.R2 <- c(pscl::pR2(nbglm.bz.best), unlist(modEvA::RsqGLM(nbglm.bz.best)))
print(nbglm.bz.pseudo.R2)


### Checking model assumptions about over-dispersion: compare with the Poission model
pslm.bz.full <- glm(dzz_mort ~ year * meV_mort_rate, offset = log(pop / 100000), 
                    family = poisson(link = "log"), data = d.bz)
summary(pslm.bz.full)

# The chi-square and p-value comparing the two models
print(2 * (logLik(nbglm.bz.8) - logLik(pslm.bz.full)))  
pchisq(2 * (logLik(nbglm.bz.8) - logLik(pslm.bz.full)), 
       df = 1, lower.tail = FALSE)  # The p-value


# write out the key results:
sink(paste0(adt.dir, "Negative binomial regression - full time series - rate.txt"))

cat("Location of interest: ")
cat(LOC)
cat("\nDisease of interest: ")
cat(DZZ)
cat("\nTime series: full\n\n\n")

cat("Negative Binomial models:\n")
stargazer(nbglm.bz.2, nbglm.bz.4, nbglm.bz.6, nbglm.bz.8,
          type = "text", style = "all", out = NULL,
          ci = TRUE, ci.level = 0.95,
          single.row = FALSE,
          column.labels = NULL,
          order = c(1, 3, 2),
          covariate.labels = c("Measles mortality rate", "Year", "Measles x Year"),
          dep.var.caption = "",
          dep.var.labels = "Dependent variable: Non-measles mortality",
          model.names = FALSE, model.numbers = TRUE, object.names = FALSE,
          omit = NULL, notes.align = "c")
cat("\n\n\n")

cat("Model comparision: effect of year\n\n")
print(anova(nbglm.bz.2, nbglm.bz.6))
cat("\n\n\n")

cat("Model comparision: effect of measlse mortality rate\n\n")
print(anova(nbglm.bz.4, nbglm.bz.6))
cat("\n\n\n")

cat("Model comparision: effect of the interactive term\n\n")
print(anova(nbglm.bz.6, nbglm.bz.8))
cat("\n\n\n")

cat("Stepwise model selction:\n\n")
step(nbglm.bz.8, direction="both")
cat("\n\n\n")

cat("Model comparison:\n\n")
print(nbglm.bz.aic)
cat("\n")

sink()





# 9. Correlation between non-measles diseases -----------------------------

### Diseases included in each category:
nonmeV.list <- list(
  asthma = c("asthma"), 
  bronchitis = c("acute bronchitis", "bronchiectasis", "chronic bronchitis", 
                 "bronchitis unsp", "bronchopneumonia"),
  diarrheal = c("helminths", "intestinal"),
  emphysema = c("emphysema"),
  LRT = c("acute bronchitis", "bronchiectasis", "chronic bronchitis",
          "bronchitis unsp","bronchopneumonia", "pneumococcal pneumonia", 
          "pneumonia other", "pneumonia unsp","resp infection"),
  meningitis = c("mening unsp","bacterial mening"),
  meningococcal = c("meningococcal infec"),
  pneumococcus = c("pneumococcal pneumonia"),
  pneumonia = c("pneumococcal pneumonia","pneumonia unsp", "pneumonia other"),
  respiratory = c("URT","acute bronchitis", "bronchiectasis", "chronic bronchitis",
                  "bronchitis unsp","bronchopneumonia", "pneumococcal pneumonia", 
                  "pneumonia other", "pneumonia unsp","resp infection"),
  rubella = c("rubella"),
  septicemia = c("septicemia"),
  TB = c("TB"),
  URT = c("URT")
)

# pneumonia = pneumococcus + "pneumonia unsp" + "pneumonia other"
# LRT = bronchitis + pneumonia + "resp infection"
# respiratory = LRT + URT
# rubella and emphysema both have very few mortality per year (~10 or <10)
# URT also had < 100 mortality per year

# Remove overlapping categories and disease with very few mortality
dzz_to_exclude <- c("all", "pneumococcus", "LRT", "respiratory", 
                    "rubella", "emphysema", "URT")

# Load the function:
source(paste0(master.dir, "Function - pairwise disease correlation.R"))

# Output directory:
dzz.dir <- paste0(adt.dir, "pairwise/")
if(!dir.exists(dzz.dir)) dir.create(dzz.dir)





# (a) Extract different diseases ---------------------------------------

### 1) disease count: original data
dzz.count <- infection %>% 
  dplyr::filter(city == LOC, age == AGE, !(dzzName %in% dzz_to_exclude)) %>%
  dplyr::select(year, meV_mort, dzz_mort, city, dzzName) %>%
  reshape2::dcast(year + city + meV_mort ~ dzzName, value.var = "dzz_mort") %>%
  dplyr::select(-city) %>%
  dplyr::rename(measles = meV_mort) %>%
  dplyr::mutate(year.idx = year - min(year) + 1)


### 2) disease count: detrended data
dzz.count.dt <- dzz.count
disease.name.temp <- names(dzz.count)[!(names(dzz.count) %in% c("year", "city", "year.idx"))]
cp.dzz.count <- vector(mode = "list", length = length(disease.name.temp))
names(cp.dzz.count) <- disease.name.temp

for(i in disease.name.temp){
  cp.dzz.count[[i]] <- change_point(d = dzz.count, write.dir = dzz.dir, 
                                    ts.var = i, year.var = "year.idx", 
                                    LOC = LOC, DZZ = DZZ, 
                                    output.txt = FALSE, output.ts = FALSE, output.coeff = FALSE)
  dzz.count.dt[, i] <- dzz.count[, i] - fitted(cp.dzz.count[[i]])[, "fitted"]
}


### 3) disease count: difference between consecutive years
dzz.count.diff <- data.frame(diff(data.matrix(dzz.count)))
dzz.count.diff$year <- dzz.count$year[2:nrow(dzz.count)]  # Use the latter year in the two consecutive years


### 4) disease rate: original data
dzz.rate <- infection %>% 
  dplyr::filter(city == LOC, age == AGE, !(dzzName %in% dzz_to_exclude)) %>%
  dplyr::select(year, meV_mort_rate, dzz_rate, city, dzzName) %>%
  reshape2::dcast(year + city + meV_mort_rate ~ dzzName, value.var = "dzz_rate") %>%
  dplyr::select(-city) %>%
  dplyr::rename(measles = meV_mort_rate) %>%
  dplyr::mutate(year.idx = year - min(year) + 1)


### 5) disease rate: detrended data
dzz.rate.dt <- dzz.rate
disease.name.temp <- names(dzz.rate)[!(names(dzz.rate) %in% c("year", "city", "year.idx"))]
cp.dzz.rate <- vector(mode = "list", length = length(disease.name.temp))
names(cp.dzz.rate) <- disease.name.temp

for(i in disease.name.temp){
  cp.dzz.rate[[i]] <- change_point(d = dzz.rate, write.dir = dzz.dir, 
                                   ts.var = i, year.var = "year.idx", 
                                   LOC = LOC, DZZ = DZZ,
                                   output.txt = FALSE, output.ts = FALSE, output.coeff = FALSE)
  dzz.rate.dt[, i] <- dzz.rate[, i] - fitted(cp.dzz.rate[[i]])[, "fitted"]
}


### 6) disease rate: difference between consecutive years
dzz.rate.diff <- data.frame(diff(data.matrix(dzz.rate)))
dzz.rate.diff$year <- dzz.rate$year[2:nrow(dzz.rate)]  # Use the latter year in the two consecutive years


### Remove the year.idx variable
dzz.count$year.idx <- NULL
dzz.count.dt$year.idx <- NULL
dzz.count.diff$year.idx <- NULL
dzz.rate$year.idx <- NULL
dzz.rate.dt$year.idx <- NULL
dzz.rate.diff$year.idx <- NULL


### Order diseases by their rate in the year 1980
dzz.order.fixed <- c("year", "measles", 
                     names(sort(dzz.rate[dzz.rate$year == 1980,!(names(dzz.rate) %in% c("year", "measles"))], 
                                decreasing = TRUE)))





# (b) Full time series ------------------------------------------------------

### 1) Mortality count: original data
pairwise.count.ori <- disease_pairwise_correlation(dzz = dzz.count, 
                                                   dzz.order = dzz.order.fixed, 
                                                   data.type = "mortality count",
                                                   time.period = "full time series", 
                                                   data.manipulation = "original data", 
                                                   include.year = FALSE)
print(pairwise.count.ori$ts)           # time series
print(pairwise.count.ori$heatmap)      # heatmap of correlation coefficients
print(pairwise.count.ori$histogram)    # histogram of correlation coefficients
print(pairwise.count.ori$sample_size)  # correlation coefficients ~ sample size
summary(pairwise.count.ori$sample_size_lm)

ggsave(filename = paste0(dzz.dir, "Full time series-count-original - time series.png"), 
       plot = pairwise.count.ori$ts, 
       width = 10, height = 7, units = "in", dpi = 300)

ggsave(filename = paste0(dzz.dir, "Full time series-count-original - heatmap of correlation coefficient.png"), 
       plot = pairwise.count.ori$heatmap, 
       width = 6, height = 4, units = "in", dpi = 300)

ggsave(filename = paste0(dzz.dir, "Full time series-count-original - Sample size.png"), 
       plot = pairwise.count.ori$sample_size,
       width = 6, height = 4, units = "in", dpi = 300)

sink(paste0(dzz.dir, "Full time series-count-original - Sample size.txt"))
summary(pairwise.count.ori$sample_size_lm)
sink()



### 2) Mortality count: detrended data
pairwise.count.dt <- disease_pairwise_correlation(dzz = dzz.count.dt, 
                                                  dzz.order = dzz.order.fixed, 
                                                  data.type = "mortality count",
                                                  time.period = "full time series", 
                                                  data.manipulation = "detrended data", 
                                                  include.year = FALSE)
print(pairwise.count.dt$ts)           # time series
print(pairwise.count.dt$heatmap)      # heatmap of correlation coefficients
print(pairwise.count.dt$histogram)    # histogram of correlation coefficients
print(pairwise.count.dt$sample_size)  # correlation coefficients ~ sample size
summary(pairwise.count.dt$sample_size_lm)

ggsave(filename = paste0(dzz.dir, "Full time series-count-detrended - time series.png"), 
       plot = pairwise.count.dt$ts, 
       width = 10, height = 7, units = "in", dpi = 300)

ggsave(filename = paste0(dzz.dir, "Full time series-count-detrended - heatmap of correlation coefficient.png"), 
       plot = pairwise.count.dt$heatmap, 
       width = 6, height = 4, units = "in", dpi = 300)

ggsave(filename = paste0(dzz.dir, "Full time series-count-detrended - Sample size.png"), 
       plot = pairwise.count.dt$sample_size,
       width = 6, height = 4, units = "in", dpi = 300)

sink(paste0(dzz.dir, "Full time series-count-detrended - Sample size.txt"))
summary(pairwise.count.dt$sample_size_lm)
sink()



### 3) Mortality count: difference between consecutive years
pairwise.count.diff <- disease_pairwise_correlation(dzz = dzz.count.diff, 
                                                    dzz.order = dzz.order.fixed, 
                                                    data.type = "mortality count",
                                                    time.period = "full time series", 
                                                    data.manipulation = "difference", 
                                                    include.year = FALSE)
print(pairwise.count.diff$ts)           # time series
print(pairwise.count.diff$heatmap)      # heatmap of correlation coefficients
print(pairwise.count.diff$histogram)    # histogram of correlation coefficients
print(pairwise.count.diff$sample_size)  # correlation coefficients ~ sample size
summary(pairwise.count.diff$sample_size_lm)

ggsave(filename = paste0(dzz.dir, "Full time series-count-difference - time series.png"), 
       plot = pairwise.count.diff$ts, 
       width = 10, height = 7, units = "in", dpi = 300)

ggsave(filename = paste0(dzz.dir, "Full time series-count-difference - heatmap of correlation coefficient.png"), 
       plot = pairwise.count.diff$heatmap, 
       width = 6, height = 4, units = "in", dpi = 300)

ggsave(filename = paste0(dzz.dir, "Full time series-count-difference - Sample size.png"), 
       plot = pairwise.count.diff$sample_size,
       width = 6, height = 4, units = "in", dpi = 300)

sink(paste0(dzz.dir, "Full time series-count-difference - Sample size.txt"))
summary(pairwise.count.diff$sample_size_lm)
sink()



### 4) Mortality rate: original data
pairwise.rate.ori <- disease_pairwise_correlation(dzz = dzz.rate, 
                                                  dzz.order = dzz.order.fixed, 
                                                  data.type = "mortality/100,000 ppl",
                                                  time.period = "full time series", 
                                                  data.manipulation = "original data", 
                                                  include.year = FALSE)
print(pairwise.rate.ori$ts)           # time series
print(pairwise.rate.ori$heatmap)      # heatmap of correlation coefficients
print(pairwise.rate.ori$histogram)    # histogram of correlation coefficients
print(pairwise.rate.ori$sample_size)  # correlation coefficients ~ sample size
summary(pairwise.rate.ori$sample_size_lm)

ggsave(filename = paste0(dzz.dir, "Full time series-rate-original - time series.png"), 
       plot = pairwise.rate.ori$ts, 
       width = 10, height = 7, units = "in", dpi = 300)

ggsave(filename = paste0(dzz.dir, "Full time series-rate-original - heatmap of correlation coefficient.png"), 
       plot = pairwise.rate.ori$heatmap, 
       width = 6, height = 4, units = "in", dpi = 300)

ggsave(filename = paste0(dzz.dir, "Full time series-rate-original - Sample size.png"), 
       plot = pairwise.rate.ori$sample_size,
       width = 6, height = 4, units = "in", dpi = 300)

sink(paste0(dzz.dir, "Full time series-rate-original - Sample size.txt"))
summary(pairwise.rate.ori$sample_size_lm)
sink()



### 5) Mortality rate: detrended data
pairwise.rate.dt <- disease_pairwise_correlation(dzz = dzz.rate.dt, 
                                                 dzz.order = dzz.order.fixed, 
                                                 data.type = "mortality/100,000 ppl",
                                                 time.period = "full time series", 
                                                 data.manipulation = "detrended data", 
                                                 include.year = FALSE)
print(pairwise.rate.dt$ts)           # time series
print(pairwise.rate.dt$heatmap)      # heatmap of correlation coefficients
print(pairwise.rate.dt$histogram)    # histogram of correlation coefficients
print(pairwise.rate.dt$sample_size)  # correlation coefficients ~ sample size
summary(pairwise.rate.dt$sample_size_lm)

ggsave(filename = paste0(dzz.dir, "Full time series-rate-detrended - time series.png"), 
       plot = pairwise.rate.dt$ts, 
       width = 10, height = 7, units = "in", dpi = 300)

ggsave(filename = paste0(dzz.dir, "Full time series-rate-detrended - heatmap of correlation coefficient.png"), 
       plot = pairwise.rate.dt$heatmap, 
       width = 6, height = 4, units = "in", dpi = 300)

ggsave(filename = paste0(dzz.dir, "Full time series-rate-detrended - Sample size.png"), 
       plot = pairwise.rate.dt$sample_size,
       width = 6, height = 4, units = "in", dpi = 300)

sink(paste0(dzz.dir, "Full time series-rate-detrended - Sample size.txt"))
summary(pairwise.rate.dt$sample_size_lm)
sink()



### 6) Mortality rate: difference between consecutive years
pairwise.rate.diff <- disease_pairwise_correlation(dzz = dzz.rate.diff, 
                                                   dzz.order = dzz.order.fixed, 
                                                   data.type = "mortality/100,000 ppl",
                                                   time.period = "full time series", 
                                                   data.manipulation = "difference", 
                                                   include.year = FALSE)
print(pairwise.rate.diff$ts)           # time series
print(pairwise.rate.diff$heatmap)      # heatmap of correlation coefficients
print(pairwise.rate.diff$histogram)    # histogram of correlation coefficients
print(pairwise.rate.diff$sample_size)  # correlation coefficients ~ sample size
summary(pairwise.rate.diff$sample_size_lm)

ggsave(filename = paste0(dzz.dir, "Full time series-rate-difference - time series.png"), 
       plot = pairwise.rate.diff$ts, 
       width = 10, height = 7, units = "in", dpi = 300)

ggsave(filename = paste0(dzz.dir, "Full time series-rate-difference - heatmap of correlation coefficient.png"), 
       plot = pairwise.rate.diff$heatmap, 
       width = 6, height = 4, units = "in", dpi = 300)

ggsave(filename = paste0(dzz.dir, "Full time series-rate-difference - Sample size.png"), 
       plot = pairwise.rate.diff$sample_size,
       width = 6, height = 4, units = "in", dpi = 300)

sink(paste0(dzz.dir, "Full time series-rate-difference - Sample size.txt"))
summary(pairwise.rate.diff$sample_size_lm)
sink()






# (c) Before the change-point --------------------------------------------------

dzz.count.t1 <- dzz.count %>% filter(year < cp.bz.count.year["mean"])
dzz.count.t1.dt <- dzz.count.dt %>% filter(year < cp.bz.count.year["mean"])
dzz.count.t1.diff <- dzz.count.diff %>% filter(year < cp.bz.count.year["mean"])

dzz.rate.t1 <- dzz.rate %>% filter(year < cp.bz.rate.year["mean"])
dzz.rate.t1.dt <- dzz.rate.dt %>% filter(year < cp.bz.rate.year["mean"])
dzz.rate.t1.diff <- dzz.rate.diff %>% filter(year < cp.bz.rate.year["mean"])


### 1) Mortality count: original data
pairwise.count.t1.ori <- disease_pairwise_correlation(dzz = dzz.count.t1, 
                                                      dzz.order = dzz.order.fixed, 
                                                      data.type = "mortality count",
                                                      time.period = "before change-point", 
                                                      data.manipulation = "original data", 
                                                      include.year = FALSE)
print(pairwise.count.t1.ori$ts)           # time series
print(pairwise.count.t1.ori$heatmap)      # heatmap of correlation coefficients
print(pairwise.count.t1.ori$histogram)    # histogram of correlation coefficients
print(pairwise.count.t1.ori$sample_size)  # correlation coefficients ~ sample size
summary(pairwise.count.t1.ori$sample_size_lm)

ggsave(filename = paste0(dzz.dir, "Before change point-count-original - time series.png"), 
       plot = pairwise.count.t1.ori$ts, 
       width = 10, height = 7, units = "in", dpi = 300)

ggsave(filename = paste0(dzz.dir, "Before change point-count-original - heatmap of correlation coefficient.png"), 
       plot = pairwise.count.t1.ori$heatmap, 
       width = 6, height = 4, units = "in", dpi = 300)

ggsave(filename = paste0(dzz.dir, "Before change point-count-original - Sample size.png"), 
       plot = pairwise.count.t1.ori$sample_size,
       width = 6, height = 4, units = "in", dpi = 300)

sink(paste0(dzz.dir, "Before change point-count-original - Sample size.txt"))
summary(pairwise.count.t1.ori$sample_size_lm)
sink()



### 2) Mortality count: detrended data
pairwise.count.t1.dt <- disease_pairwise_correlation(dzz = dzz.count.t1.dt, 
                                                     dzz.order = dzz.order.fixed, 
                                                     data.type = "mortality count",
                                                     time.period = "before change-point", 
                                                     data.manipulation = "detrended data", 
                                                     include.year = FALSE)
print(pairwise.count.t1.dt$ts)           # time series
print(pairwise.count.t1.dt$heatmap)      # heatmap of correlation coefficients
print(pairwise.count.t1.dt$histogram)    # histogram of correlation coefficients
print(pairwise.count.t1.dt$sample_size)  # correlation coefficients ~ sample size
summary(pairwise.count.t1.dt$sample_size_lm)

ggsave(filename = paste0(dzz.dir, "Before change point-count-detrended - time series.png"), 
       plot = pairwise.count.t1.dt$ts, 
       width = 10, height = 7, units = "in", dpi = 300)

ggsave(filename = paste0(dzz.dir, "Before change point-count-detrended - heatmap of correlation coefficient.png"), 
       plot = pairwise.count.t1.dt$heatmap, 
       width = 6, height = 4, units = "in", dpi = 300)

ggsave(filename = paste0(dzz.dir, "Before change point-count-detrended - Sample size.png"), 
       plot = pairwise.count.t1.dt$sample_size,
       width = 6, height = 4, units = "in", dpi = 300)

sink(paste0(dzz.dir, "Before change point-count-detrended - Sample size.txt"))
summary(pairwise.count.t1.dt$sample_size_lm)
sink()



### 3) Mortality count: difference between consecutive years
pairwise.count.t1.diff <- disease_pairwise_correlation(dzz = dzz.count.t1.diff, 
                                                       dzz.order = dzz.order.fixed, 
                                                       data.type = "mortality count",
                                                       time.period = "before change-point", 
                                                       data.manipulation = "difference", 
                                                       include.year = FALSE)
print(pairwise.count.t1.diff$ts)           # time series
print(pairwise.count.t1.diff$heatmap)      # heatmap of correlation coefficients
print(pairwise.count.t1.diff$histogram)    # histogram of correlation coefficients
print(pairwise.count.t1.diff$sample_size)  # correlation coefficients ~ sample size
summary(pairwise.count.t1.diff$sample_size_lm)

ggsave(filename = paste0(dzz.dir, "Before change point-count-difference - time series.png"), 
       plot = pairwise.count.t1.diff$ts, 
       width = 10, height = 7, units = "in", dpi = 300)

ggsave(filename = paste0(dzz.dir, "Before change point-count-difference - heatmap of correlation coefficient.png"), 
       plot = pairwise.count.t1.diff$heatmap, 
       width = 6, height = 4, units = "in", dpi = 300)

ggsave(filename = paste0(dzz.dir, "Before change point-count-difference - Sample size.png"), 
       plot = pairwise.count.t1.diff$sample_size,
       width = 6, height = 4, units = "in", dpi = 300)

sink(paste0(dzz.dir, "Before change point-count-difference - Sample size.txt"))
summary(pairwise.count.t1.diff$sample_size_lm)
sink()



### 4) Mortality rate: original data
pairwise.rate.t1.ori <- disease_pairwise_correlation(dzz = dzz.rate.t1, 
                                                     dzz.order = dzz.order.fixed, 
                                                     data.type = "mortality/100,000 ppl",
                                                     time.period = "before change-point", 
                                                     data.manipulation = "original data", 
                                                     include.year = FALSE)
print(pairwise.rate.t1.ori$ts)           # time series
print(pairwise.rate.t1.ori$heatmap)      # heatmap of correlation coefficients
print(pairwise.rate.t1.ori$histogram)    # histogram of correlation coefficients
print(pairwise.rate.t1.ori$sample_size)  # correlation coefficients ~ sample size
summary(pairwise.rate.t1.ori$sample_size_lm)

ggsave(filename = paste0(dzz.dir, "Before change point-rate-original - time series.png"), 
       plot = pairwise.rate.t1.ori$ts, 
       width = 10, height = 7, units = "in", dpi = 300)

ggsave(filename = paste0(dzz.dir, "Before change point-rate-original - heatmap of correlation coefficient.png"), 
       plot = pairwise.rate.t1.ori$heatmap, 
       width = 6, height = 4, units = "in", dpi = 300)

ggsave(filename = paste0(dzz.dir, "Before change point-rate-original - Sample size.png"), 
       plot = pairwise.rate.t1.ori$sample_size,
       width = 6, height = 4, units = "in", dpi = 300)

sink(paste0(dzz.dir, "Before change point-rate-original - Sample size.txt"))
summary(pairwise.rate.t1.ori$sample_size_lm)
sink()



### 5) Mortality rate: detrended data
pairwise.rate.t1.dt <- disease_pairwise_correlation(dzz = dzz.rate.t1.dt, 
                                                    dzz.order = dzz.order.fixed, 
                                                    data.type = "mortality/100,000 ppl",
                                                    time.period = "before change-point", 
                                                    data.manipulation = "detrended data",
                                                    include.year = FALSE)
print(pairwise.rate.t1.dt$ts)           # time series
print(pairwise.rate.t1.dt$heatmap)      # heatmap of correlation coefficients
print(pairwise.rate.t1.dt$histogram)    # histogram of correlation coefficients
print(pairwise.rate.t1.dt$sample_size)  # correlation coefficients ~ sample size
summary(pairwise.rate.t1.dt$sample_size_lm)

ggsave(filename = paste0(dzz.dir, "Before change point-rate-detrended - time series.png"), 
       plot = pairwise.rate.t1.dt$ts, 
       width = 10, height = 7, units = "in", dpi = 300)

ggsave(filename = paste0(dzz.dir, "Before change point-rate-detrended - heatmap of correlation coefficient.png"), 
       plot = pairwise.rate.t1.dt$heatmap, 
       width = 6, height = 4, units = "in", dpi = 300)

ggsave(filename = paste0(dzz.dir, "Before change point-rate-detrended - Sample size.png"), 
       plot = pairwise.rate.t1.dt$sample_size,
       width = 6, height = 4, units = "in", dpi = 300)

sink(paste0(dzz.dir, "Before change point-rate-detrended - Sample size.txt"))
summary(pairwise.rate.t1.dt$sample_size_lm)
sink()



### 6) Mortality rate: difference between consecutive years
pairwise.rate.t1.diff <- disease_pairwise_correlation(dzz = dzz.rate.t1.diff, 
                                                      dzz.order = dzz.order.fixed, 
                                                      data.type = "mortality/100,000 ppl",
                                                      time.period = "before change-point", 
                                                      data.manipulation = "difference",
                                                      include.year = FALSE)
print(pairwise.rate.t1.diff$ts)           # time series
print(pairwise.rate.t1.diff$heatmap)      # heatmap of correlation coefficients
print(pairwise.rate.t1.diff$histogram)    # histogram of correlation coefficients
print(pairwise.rate.t1.diff$sample_size)  # correlation coefficients ~ sample size
summary(pairwise.rate.t1.diff$sample_size_lm)

ggsave(filename = paste0(dzz.dir, "Before change point-rate-difference - time series.png"), 
       plot = pairwise.rate.t1.diff$ts, 
       width = 10, height = 7, units = "in", dpi = 300)

ggsave(filename = paste0(dzz.dir, "Before change point-rate-difference - heatmap of correlation coefficient.png"), 
       plot = pairwise.rate.t1.diff$heatmap, 
       width = 6, height = 4, units = "in", dpi = 300)

ggsave(filename = paste0(dzz.dir, "Before change point-rate-difference - Sample size.png"), 
       plot = pairwise.rate.t1.diff$sample_size,
       width = 6, height = 4, units = "in", dpi = 300)

sink(paste0(dzz.dir, "Before change point-rate-difference - Sample size.txt"))
summary(pairwise.rate.t1.diff$sample_size_lm)
sink()





# (d) After the change-point --------------------------------------------------

dzz.count.t2 <- dzz.count %>% filter(year >= cp.bz.count.year["mean"])
dzz.count.t2.dt <- dzz.count.dt %>% filter(year >= cp.bz.count.year["mean"])
dzz.count.t2.diff <- dzz.count.diff %>% filter(year >= (cp.bz.count.year["mean"] + 1))
  # Remove the first year after the change-point as it is between the two period

dzz.rate.t2 <- dzz.rate %>% filter(year >= cp.bz.rate.year["mean"])
dzz.rate.t2.dt <- dzz.rate.dt %>% filter(year >= cp.bz.rate.year["mean"])
dzz.rate.t2.diff <- dzz.rate.diff %>% filter(year >= (cp.bz.rate.year["mean"] + 1))
  # Remove the first year after the change-point as it is between the two period


### 1) Mortality count: original data
pairwise.count.t2.ori <- disease_pairwise_correlation(dzz = dzz.count.t2, 
                                                      dzz.order = dzz.order.fixed, 
                                                      data.type = "mortality count",
                                                      time.period = "after change-point", 
                                                      data.manipulation = "original data",
                                                      include.year = FALSE)
print(pairwise.count.t2.ori$ts)           # time series
print(pairwise.count.t2.ori$heatmap)      # heatmap of correlation coefficients
print(pairwise.count.t2.ori$histogram)    # histogram of correlation coefficients
print(pairwise.count.t2.ori$sample_size)  # correlation coefficients ~ sample size
summary(pairwise.count.t2.ori$sample_size_lm)

ggsave(filename = paste0(dzz.dir, "After change point-count-original - time series.png"), 
       plot = pairwise.count.t2.ori$ts, 
       width = 10, height = 7, units = "in", dpi = 300)

ggsave(filename = paste0(dzz.dir, "After change point-count-original - heatmap of correlation coefficient.png"), 
       plot = pairwise.count.t2.ori$heatmap, 
       width = 6, height = 4, units = "in", dpi = 300)

ggsave(filename = paste0(dzz.dir, "After change point-count-original - Sample size.png"), 
       plot = pairwise.count.t2.ori$sample_size,
       width = 6, height = 4, units = "in", dpi = 300)

sink(paste0(dzz.dir, "After change point-count-original - Sample size.txt"))
summary(pairwise.count.t2.ori$sample_size_lm)
sink()



### 2) Mortality count: detrended data
pairwise.count.t2.dt <- disease_pairwise_correlation(dzz = dzz.count.t2.dt, 
                                                     dzz.order = dzz.order.fixed, 
                                                     data.type = "mortality count",
                                                     time.period = "after change-point", 
                                                     data.manipulation = "detrended data",
                                                     include.year = FALSE)
print(pairwise.count.t2.dt$ts)           # time series
print(pairwise.count.t2.dt$heatmap)      # heatmap of correlation coefficients
print(pairwise.count.t2.dt$histogram)    # histogram of correlation coefficients
print(pairwise.count.t2.dt$sample_size)  # correlation coefficients ~ sample size
summary(pairwise.count.t2.dt$sample_size_lm)

ggsave(filename = paste0(dzz.dir, "After change point-count-detrended - time series.png"), 
       plot = pairwise.count.t2.dt$ts, 
       width = 10, height = 7, units = "in", dpi = 300)

ggsave(filename = paste0(dzz.dir, "After change point-count-detrended - heatmap of correlation coefficient.png"), 
       plot = pairwise.count.t2.dt$heatmap, 
       width = 6, height = 4, units = "in", dpi = 300)

ggsave(filename = paste0(dzz.dir, "After change point-count-detrended - Sample size.png"), 
       plot = pairwise.count.t2.dt$sample_size,
       width = 6, height = 4, units = "in", dpi = 300)

sink(paste0(dzz.dir, "After change point-count-detrended - Sample size.txt"))
summary(pairwise.count.t2.dt$sample_size_lm)
sink()



### 3) Mortality count: difference between consecutive years
pairwise.count.t2.diff <- disease_pairwise_correlation(dzz = dzz.count.t2.diff, 
                                                       dzz.order = dzz.order.fixed, 
                                                       data.type = "mortality count",
                                                       time.period = "after change-point", 
                                                       data.manipulation = "difference",
                                                       include.year = FALSE)
print(pairwise.count.t2.diff$ts)           # time series
print(pairwise.count.t2.diff$heatmap)      # heatmap of correlation coefficients
print(pairwise.count.t2.diff$histogram)    # histogram of correlation coefficients
print(pairwise.count.t2.diff$sample_size)  # correlation coefficients ~ sample size
summary(pairwise.count.t2.diff$sample_size_lm)

ggsave(filename = paste0(dzz.dir, "After change point-count-difference - time series.png"), 
       plot = pairwise.count.t2.diff$ts, 
       width = 10, height = 7, units = "in", dpi = 300)

ggsave(filename = paste0(dzz.dir, "After change point-count-difference - heatmap of correlation coefficient.png"), 
       plot = pairwise.count.t2.diff$heatmap, 
       width = 6, height = 4, units = "in", dpi = 300)

ggsave(filename = paste0(dzz.dir, "After change point-count-difference - Sample size.png"), 
       plot = pairwise.count.t2.diff$sample_size,
       width = 6, height = 4, units = "in", dpi = 300)

sink(paste0(dzz.dir, "After change point-count-difference - Sample size.txt"))
summary(pairwise.count.t2.diff$sample_size_lm)
sink()



### 4) Mortality rate: original data
pairwise.rate.t2.ori <- disease_pairwise_correlation(dzz = dzz.rate.t2, 
                                                     dzz.order = dzz.order.fixed, 
                                                     data.type = "mortality/100,000 ppl",
                                                     time.period = "after change-point",
                                                     data.manipulation = "original data",
                                                     include.year = FALSE)
print(pairwise.rate.t2.ori$ts)           # time series
print(pairwise.rate.t2.ori$heatmap)      # heatmap of correlation coefficients
print(pairwise.rate.t2.ori$histogram)    # histogram of correlation coefficients
print(pairwise.rate.t2.ori$sample_size)  # correlation coefficients ~ sample size
summary(pairwise.rate.t2.ori$sample_size_lm)

ggsave(filename = paste0(dzz.dir, "After change point-rate-original - time series.png"), 
       plot = pairwise.rate.t2.ori$ts, 
       width = 10, height = 7, units = "in", dpi = 300)

ggsave(filename = paste0(dzz.dir, "After change point-rate-original - heatmap of correlation coefficient.png"), 
       plot = pairwise.rate.t2.ori$heatmap, 
       width = 6, height = 4, units = "in", dpi = 300)

ggsave(filename = paste0(dzz.dir, "After change point-rate-original - Sample size.png"), 
       plot = pairwise.rate.t2.ori$sample_size,
       width = 6, height = 4, units = "in", dpi = 300)

sink(paste0(dzz.dir, "After change point-rate-original - Sample size.txt"))
summary(pairwise.rate.t2.ori$sample_size_lm)
sink()



### 5) Mortality rate: detrended data
pairwise.rate.t2.dt <- disease_pairwise_correlation(dzz = dzz.rate.t2.dt, 
                                                    dzz.order = dzz.order.fixed, 
                                                    data.type = "mortality/100,000 ppl",
                                                    time.period = "after change-point", 
                                                    data.manipulation = "detrended data",
                                                    include.year = FALSE)
print(pairwise.rate.t2.dt$ts)           # time series
print(pairwise.rate.t2.dt$heatmap)      # heatmap of correlation coefficients
print(pairwise.rate.t2.dt$histogram)    # histogram of correlation coefficients
print(pairwise.rate.t2.dt$sample_size)  # correlation coefficients ~ sample size
summary(pairwise.rate.t2.dt$sample_size_lm)

ggsave(filename = paste0(dzz.dir, "After change point-rate-detrended - time series.png"), 
       plot = pairwise.rate.t2.dt$ts, 
       width = 10, height = 7, units = "in", dpi = 300)

ggsave(filename = paste0(dzz.dir, "After change point-rate-detrended - heatmap of correlation coefficient.png"), 
       plot = pairwise.rate.t2.dt$heatmap, 
       width = 6, height = 4, units = "in", dpi = 300)

ggsave(filename = paste0(dzz.dir, "After change point-rate-detrended - Sample size.png"), 
       plot = pairwise.rate.t2.dt$sample_size,
       width = 6, height = 4, units = "in", dpi = 300)

sink(paste0(dzz.dir, "After change point-rate-detrended - Sample size.txt"))
summary(pairwise.rate.t2.dt$sample_size_lm)
sink()



### 6) Mortality rate: difference between consecutive years
pairwise.rate.t2.diff <- disease_pairwise_correlation(dzz = dzz.rate.t2.diff, 
                                                      dzz.order = dzz.order.fixed, 
                                                      data.type = "mortality/100,000 ppl",
                                                      time.period = "after change-point", 
                                                      data.manipulation = "difference",
                                                      include.year = FALSE)
print(pairwise.rate.t2.diff$ts)           # time series
print(pairwise.rate.t2.diff$heatmap)      # heatmap of correlation coefficients
print(pairwise.rate.t2.diff$histogram)    # histogram of correlation coefficients
print(pairwise.rate.t2.diff$sample_size)  # correlation coefficients ~ sample size
summary(pairwise.rate.t2.diff$sample_size_lm)

ggsave(filename = paste0(dzz.dir, "After change point-rate-difference - time series.png"), 
       plot = pairwise.rate.t2.diff$ts, 
       width = 10, height = 7, units = "in", dpi = 300)

ggsave(filename = paste0(dzz.dir, "After change point-rate-difference - heatmap of correlation coefficient.png"), 
       plot = pairwise.rate.t2.diff$heatmap, 
       width = 6, height = 4, units = "in", dpi = 300)

ggsave(filename = paste0(dzz.dir, "After change point-rate-difference - Sample size.png"), 
       plot = pairwise.rate.t2.diff$sample_size,
       width = 6, height = 4, units = "in", dpi = 300)

sink(paste0(dzz.dir, "After change point-rate-difference - Sample size.txt"))
summary(pairwise.rate.t2.diff$sample_size_lm)
sink()





# 10. Generate report ------------------------------------------------------


# (a) Save key figure objects ----------------------------------------------

Brazil_1_9_meV_nonmeV <- list(ts_ori = f1.2,                            # time series
                              ts_dt = f3.2, 
                              ts_diff = plot.diff.bz.ts.rate,
                              cor_ori = f2.2,                           # correlation plot
                              cor_dt = f4.2, 
                              cor_diff = plot.diff.bz.cor.rate,
                              cp_meV = plot.fit.bz.rate,                # change-point time series
                              cp_dzz = plot.fit.bz.dzz.rate,
                              pairwise_full_ori = pairwise.rate.ori,    # disease pairwise correlation
                              pairwise_full_dt = pairwise.rate.dt,
                              pairwise_full_diff = pairwise.rate.diff,
                              pairwise_t1_ori = pairwise.rate.t1.ori, 
                              pairwise_t1_dt = pairwise.rate.t1.dt,
                              pairwise_t1_diff = pairwise.rate.t1.diff,
                              pairwise_t2_ori = pairwise.rate.t2.ori, 
                              pairwise_t2_dt = pairwise.rate.t2.dt,
                              pairwise_t2_diff = pairwise.rate.t2.diff,
                              lm_1 <- lm.bz.1.rate,                     # selected regression results  
                              lm_2 <- lm.bz.3.rate, 
                              lm_3 <- lm.bz.dt.1.rate,
                              lm_4 <- lm.bz.diff.1.rate, 
                              lm_5 <- nbglm.bz.6)
save(Brazil_1_9_meV_nonmeV, file = paste0(adt.dir, LOC, "_", AGE, "_meV_dzz_export_", Sys.Date(), ".Rdata"))





# (b) Combined figures: times series and correlation ---------------------------

# Adjust figure margins
ts_cor_adjust <- theme(legend.position = "none") + 
  theme(plot.margin = unit(c(15, 20, 5, 10), units = "pt")) + 
  theme(plot.title = element_text(hjust = 0.5, vjust = 1, size = 14))
ts_cor_legend_adjust <- theme(legend.box.margin = ggplot2::margin(0, 0, 0, 0))

# Extract the legends
ts_legend <- cowplot::get_legend(Brazil_1_9_meV_nonmeV$ts_ori + ts_cor_legend_adjust)
cor_legend <- cowplot::get_legend(Brazil_1_9_meV_nonmeV$cor_ori + ts_cor_legend_adjust)

# Generate the combined figures
ts_cor_combined_plot <- plot_grid(Brazil_1_9_meV_nonmeV$ts_ori + ts_cor_adjust + 
                                    ggtitle("Original mortality\n "),
                                  Brazil_1_9_meV_nonmeV$ts_dt + ts_cor_adjust + 
                                    ggtitle("Detrended mortality\n "),
                                  Brazil_1_9_meV_nonmeV$ts_diff + ts_cor_adjust + 
                                    ggtitle("Mortality difference\n "),
                                  ts_legend, 
                                  Brazil_1_9_meV_nonmeV$cor_ori + ts_cor_adjust,
                                  Brazil_1_9_meV_nonmeV$cor_dt + ts_cor_adjust,
                                  Brazil_1_9_meV_nonmeV$cor_diff + ts_cor_adjust,
                                  cor_legend, 
                                  nrow = 2, ncol = 4, 
                                  rel_widths = c(1, 1, 1, 0.4), 
                                  rel_heights = c(1, 1),
                                  labels = c("a", "c", "e", "", "b", "d", "f", ""), 
                                  hjust = -0.5, vjust = 1.5, label_size = 20, 
                                  align = "hv", axis = "tbl")
print(ts_cor_combined_plot)
save_plot(plot = ts_cor_combined_plot, 
          filename = paste0(adt.dir, LOC, "_", AGE, "_combined_time_series_correlation.pdf"),
          base_height = 8, base_width = 15)





# (c) Combined figure: heatmap of correlation coefficients ----------------

# Adjust figure margins
hm_adjust <- theme(legend.position = "none") + 
  theme(plot.margin = unit(c(5, 10, 5, 10), units = "pt")) + 
  theme(plot.title = element_text(hjust = 0.5, vjust = 1, size = 14, margin = margin(0, 0, 10, 0)))
hm_legend_adjust <- theme(legend.box.margin = ggplot2::margin(0, 0, 0, 0), 
                          legend.position="bottom", 
                          legend.key.height = unit(15, units = "pt"),
                          legend.key.width = unit(40, units = "pt"), 
                          legend.title = element_text(vjust = 0.9, hjust = 0)) 

# Extract the legends
hm_legend <- cowplot::get_legend(Brazil_1_9_meV_nonmeV$pairwise_t1_ori$heatmap + hm_legend_adjust)

# Add text to explain the top and bottom row
year_range_before <- paste0(range(d.bz$year)[1], "-", floor(cp.bz.rate.year["mean"]))
year_range_after <- paste0(ceiling(cp.bz.rate.year["mean"]), "-", range(d.bz$year)[2])

hm_lab_top <- ggdraw() + 
  draw_text(paste0("Before\nchange-point\n(", year_range_before, ")"), 
            x = 0, y = 0.5, size = 14) + 
  theme(plot.margin = unit(c(5, 10, 5, 10), units = "pt"))
hm_lab_bottom <- ggdraw() + 
  draw_text(paste0("After\nchange-point\n(", year_range_after, ")"),
            x = 0, y = 0.5, size = 14) +
  theme(plot.margin = unit(c(5, 10, 5, 10), units = "pt"))

# Generate the combined figures
hm_plot <- plot_grid(hm_lab_top, 
                     Brazil_1_9_meV_nonmeV$pairwise_t1_ori$heatmap + hm_adjust + 
                       ggtitle("Original mortality"),
                     Brazil_1_9_meV_nonmeV$pairwise_t1_dt$heatmap + hm_adjust + 
                       ggtitle("Detrended mortality"),
                     Brazil_1_9_meV_nonmeV$pairwise_t1_diff$heatmap + hm_adjust + 
                       ggtitle("Mortality difference"),
                     hm_lab_bottom,
                     Brazil_1_9_meV_nonmeV$pairwise_t2_ori$heatmap + hm_adjust + 
                       ggtitle(""),
                     Brazil_1_9_meV_nonmeV$pairwise_t2_dt$heatmap + hm_adjust + 
                       ggtitle(""),
                     Brazil_1_9_meV_nonmeV$pairwise_t2_diff$heatmap + hm_adjust + 
                       ggtitle(""),
                     nrow = 2, ncol = 4, 
                     rel_widths = c(0.4, 1, 1, 1), 
                     rel_heights = c(1, 1),
                     labels = c("", "a", "c", "e", "", "b", "d", "f"), 
                     hjust = -1, vjust = 1.5, label_size = 20,
                     align = "hv", axis = "tblr")

# Combine the plot with the legend
hm_combined_plot <- plot_grid(hm_plot, hm_legend, nrow = 2, ncol = 1, rel_heights = c(0.9, 0.1))
print(hm_combined_plot)

# Export the combined plot:
save_plot(hm_combined_plot, 
          filename = paste0(adt.dir, LOC, "_", AGE, "_combined_correlation_heatmap.pdf"),
          base_height = 9, base_width = 15)





# (d) Save the workspace --------------------------------------------------

# rm(cp.dzz.count, cp.dzz.rate)  # The two change point analysis list took too much space
save.image(paste0(adt.dir, LOC, "_", AGE, "_meV_dzz_", Sys.Date(), ".Rdata"))





# (e) R markdown to generate a summary report ----------------------------

library(rmarkdown)

render(paste0(master.dir, LOC, "/Summarize analysis results.Rmd"), 
       output_file = paste0(adt.dir, LOC, "_", AGE, "_summary_rate_", Sys.Date(), ".pdf"), 
       params = list(Rdata = paste0(adt.dir, LOC, "_", AGE, "_meV_dzz_", Sys.Date(), ".Rdata"),
                     LOC = LOC,
                     AGE = AGE,
                     DZZ = DZZ, 
                     rc = "rate"))

render(paste0(master.dir, LOC, "/Summarize analysis results.Rmd"), 
       output_file = paste0(adt.dir, LOC, "_", AGE, "_summary_count_", Sys.Date(), ".pdf"), 
       params = list(Rdata = paste0(adt.dir, LOC, "_", AGE, "_meV_dzz_", Sys.Date(), ".Rdata"),
                     LOC = LOC,
                     AGE = AGE,
                     DZZ = DZZ, 
                     rc = "count"))
