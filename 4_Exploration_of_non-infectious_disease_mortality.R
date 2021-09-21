### Analyze the correlation between measles and non-infectious disease mortality in
### the Brazilian mortality data
###
### Input: Mortality data download from Brazil Ministry of Health, DATASUS
### http://tabnet.datasus.gov.br/cgi/deftohtm.exe?sim/cnv/obt09br.def
### The data shoule locate in the directory "raw data/" inside the location directory
### 
### Siyang Xia
### 2021.4.25


library(dplyr)
library(reshape2)
library(ggplot2)
library(cowplot)
library(mcp)
library(relaimpo)
library(stargazer)



### GLOBAL OPTIONS:
LOC <- "Brazil"
TIME_PERIOD <- "1979_1995"
AGE <- "1-9"


# Color-blind friendly palette:
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")


# Function to calculate proportional location in a range:
position_in_range <- function(x, p) min(x) + p * (max(x) - min(x))





# 1. Preparation and data input -------------------------------------------

# Directories input and output
master.dir <- "C:/Users/Siyang Xia/Dropbox/My laptop/Mina lab/Projects/Measles/Empirical/Brazil/"
rawdata.dir <- paste0(master.dir, LOC, "/raw data/")        # directory of the raw data
prepared.dir <- paste0(master.dir, LOC, "/prepared data/")  # directory to store the prepared data


# Output directory for each age group
output.dir <- paste0(master.dir, LOC, "/other disease categories/")
if(!dir.exists(output.dir)) dir.create(output.dir)


# Data input: the age group of interest determined by "AGE"
d <- read.csv(paste0(rawdata.dir, "Mortality_1979_1995_both_Brazil_", 
                     gsub(pattern = "-", replacement = "_", x = AGE), ".csv"), 
              skip = 4, sep = ";", header = TRUE, stringsAsFactors = FALSE)

# Infectious disease of interest
DZZ.list <- read.csv(paste0(master.dir, "Disease of interest.csv"), 
                     header = TRUE, stringsAsFactors = FALSE)

# Population data
pop.long <- read.csv(paste0(prepared.dir, "Brazil_population_formatted.csv"), 
                     header = TRUE, stringsAsFactors = FALSE)



### Summarize the mortality count/rate of different disease categories
# Extract the ICD-9 number
d$ICD9.number <- sapply(X = d$Categoria.CID.9, 
                        FUN = function(x) as.numeric(strsplit(x, split = " - ")[[1]][1]))

# Remove the last two rows (total number and a format note)
d <- d[1:(nrow(d)-2), ]


# Disease categories and ICD-9 number
Neoplasm.ICD9 <- 140:239
Circulatory.ICD9 <- 390:459
Musculoskeletal.ICD <- 710:739
Injury.ICD <- 800:999
Infection.ICD <- DZZ.list$ICD.9_Number

d$category <- "other"
d$category[d$ICD9.number %in% Neoplasm.ICD9] <- "neoplasm"
d$category[d$ICD9.number %in% Circulatory.ICD9] <- "circulatory"
d$category[d$ICD9.number %in% Musculoskeletal.ICD] <- "musculoskeleta"
d$category[d$ICD9.number %in% Injury.ICD] <- "injury"
d$category[d$ICD9.number %in% Infection.ICD] <- "infection"


# Format the dataset
d.long <- d %>% 
  reshape2::melt(id.var = c("category", "ICD9.number", "Categoria.CID.9"), 
                 variable.name = "Year_ori", value.name = "Mortality_ori") %>%
  dplyr::filter(Year_ori != "Total") %>%
  dplyr::mutate(Year = as.numeric(gsub(pattern = "X", replacement = "", x = Year_ori)),
                Mortality = as.numeric(gsub(pattern = "-", replacement = "0", x = Mortality_ori))) %>%
  dplyr::select(-c(Year_ori, Mortality_ori))


# Add population information
pop.long <- pop.long %>% 
  dplyr::filter(Age == "1-9") %>% 
  dplyr::select(-Age)
d.long <- d.long %>% left_join(pop.long)


# Calculate mortality rate
d.long$Mortality.rate <- d.long$Mortality / d.long$Population * 100000


# Summary by disease category
d.long.sum <- d.long %>% 
  dplyr::group_by(category, Year) %>%
  dplyr::summarise(Mortality = sum(Mortality), 
                   Population = mean(Population),
                   Mortality.rate = sum(Mortality.rate)) %>%
  dplyr::ungroup()


# Remove year 1979 as it does not have population information
d.long.sum <- d.long.sum[d.long.sum$Year >= 1980, ]

# Calculate mortality count/rate as porportions of that in 1980 (for easy visualization)
d.long.sum$Mortality.prop <- 1
d.long.sum$Mortality.rate.prop <- 1
for(i in unique(d.long.sum$category)){
  d.long.sum$Mortality.prop[d.long.sum$category == i] <- 
    d.long.sum$Mortality[d.long.sum$category == i] / d.long.sum$Mortality[d.long.sum$category == i & d.long.sum$Year == 1980]
  d.long.sum$Mortality.rate.prop[d.long.sum$category == i] <- 
    d.long.sum$Mortality.rate[d.long.sum$category == i] / d.long.sum$Mortality.rate[d.long.sum$category == i & d.long.sum$Year == 1980]
}





# 2. Time series of different disease categories --------------------------

ts_mortality_count <- ggplot(data = d.long.sum,
                             aes(x = Year, y = Mortality, linetype = category, color = category)) +
  geom_line(size = 1) +
  scale_y_continuous(name = "Mortality") + 
  scale_linetype_manual(values = 1:5, name = "") +
  scale_color_manual(values = cbPalette[c(2, 3, 4, 6, 7)], name = "") +
  xlab("Year") +
  theme_bw() +
  theme(axis.title = element_text(size = 12),
        axis.text = element_text(size = 12),
        axis.line = element_line(colour = "black"),
        axis.title.y.right = element_text( angle = 90), 
        legend.text = element_text(size = 12),
        legend.key.height = unit(0.4, "in"),
        legend.key.width = unit(0.4, "in"),
        panel.border = element_blank())
print(ts_mortality_count)
ggsave(filename = paste0(output.dir, "Time series-mortality count absolute number.png"), 
       plot = ts_mortality_count,
       width = 7, height = 4, units = "in", dpi = 300)


ts_mortality_rate <- ggplot(data = d.long.sum,
                            aes(x = Year, y = Mortality.rate, linetype = category, color = category)) +
  geom_line(size = 1) +
  scale_y_continuous(name = "Mortality rate") + 
  scale_linetype_manual(values = 1:5, name = "") +
  scale_color_manual(values = cbPalette[c(2, 3, 4, 6, 7)], name = "") +
  xlab("Year") +
  theme_bw() +
  theme(axis.title = element_text(size = 12),
        axis.text = element_text(size = 12),
        axis.line = element_line(colour = "black"),
        axis.title.y.right = element_text( angle = 90), 
        legend.text = element_text(size = 12),
        legend.key.height = unit(0.4, "in"),
        legend.key.width = unit(0.4, "in"),
        panel.border = element_blank())
print(ts_mortality_rate)
ggsave(filename = paste0(output.dir, "Time series-mortality rate absolute number.png"), 
       plot = ts_mortality_rate,
       width = 7, height = 4, units = "in", dpi = 300)


ts_mortality_count_prop <- ggplot(data = d.long.sum,
                                  aes(x = Year, y = Mortality.prop, linetype = category, color = category)) +
  geom_line(size = 1) +
  scale_y_continuous(name = "Mortality (proportional to 1980)") + 
  scale_linetype_manual(values = 1:5, name = "") +
  scale_color_manual(values = cbPalette[c(2, 3, 4, 6, 7)], name = "") +
  xlab("Year") +
  theme_bw() +
  theme(axis.title = element_text(size = 12),
        axis.text = element_text(size = 12),
        axis.line = element_line(colour = "black"),
        axis.title.y.right = element_text( angle = 90), 
        legend.text = element_text(size = 12),
        legend.key.height = unit(0.4, "in"),
        legend.key.width = unit(0.4, "in"),
        panel.border = element_blank())
print(ts_mortality_count_prop)
ggsave(filename = paste0(output.dir, "Time series-mortality count proportion.png"), 
       plot = ts_mortality_count_prop,
       width = 7, height = 4, units = "in", dpi = 300)


ts_mortality_rate_prop <- ggplot(data = d.long.sum,
                                 aes(x = Year, y = Mortality.rate.prop, linetype = category, color = category)) +
  geom_line(size = 1) +
  scale_y_continuous(name = "Mortality rate (proportional to 1980)") + 
  scale_linetype_manual(values = 1:5, name = "") +
  scale_color_manual(values = cbPalette[c(2, 3, 4, 6, 7)], name = "") +
  xlab("Year") +
  theme_bw() +
  theme(axis.title = element_text(size = 12),
        axis.text = element_text(size = 12),
        axis.line = element_line(colour = "black"),
        axis.title.y.right = element_text( angle = 90), 
        legend.text = element_text(size = 12),
        legend.key.height = unit(0.4, "in"),
        legend.key.width = unit(0.4, "in"),
        panel.border = element_blank())
print(ts_mortality_rate_prop)
ggsave(filename = paste0(output.dir, "Time series-mortality rate proportion.png"), 
       plot = ts_mortality_rate_prop,
       width = 7, height = 4, units = "in", dpi = 300)

save.image(paste0(output.dir, "Comparing multiple disease categories_", Sys.Date(), ".RData"))





# 3. Correlation between measles and non-infectious disease mortality --------
# As a negative control, we perform a correlation analysis between measles mortality
# and circulatory or neoplasm diseases, which should be independent.

# Disease category of interest
DZZ <- "circulatory"

# A sub-directory
adt.dir <- paste0(output.dir, "measles_", DZZ, "/")
if(!dir.exists(adt.dir)) dir.create(adt.dir)



# Read in the measles data
measles_data_full <- read.csv(paste0(prepared.dir, LOC, "_measles_non-measles_mortality.csv"), 
                              header = TRUE, stringsAsFactors = FALSE)
measles_data <- measles_data_full %>%
  dplyr::filter(age == AGE, city == LOC, dzzName == "all") %>%
  dplyr::select(age, year, city, dzzName, meV_mort, meV_mort_rate, pop)

# Combine the measlse data with the disease of interest
d.bz <- d.long.sum %>% 
  dplyr::filter(category == all_of(DZZ)) %>%
  dplyr::select(year = Year, dzz_mort = Mortality, dzz_rate = Mortality.rate) %>%
  dplyr::right_join(measles_data) %>%
  dplyr::mutate(dzzName = DZZ) %>%
  as.data.frame()

d.bz$year.idx <- d.bz$year - min(d.bz$year) + 1   # Use year as a time series
nyear <- nrow(d.bz)





# 4. Change-point analysis -----------------------------------

# Load the function:
source(paste0(master.dir, "Function - change-point analysis.R"))



# (a) measles time series: rate -------------------------------------------

fit.bz.rate <- change_point(d = d.bz, write.dir = adt.dir, 
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
ggsave(filename = paste0(adt.dir, "Change-point-measles-rate.png"), 
       plot = plot.fit.bz.rate,
       width = 6, height = 4, units = "in", dpi = 300)

# Calculate the time-detrended measles mortality rate:
d.bz$meV_mort_rate_dt <- d.bz$meV_mort_rate - fitted(fit.bz.rate)[, "fitted"]



# (b) non-measles time series: rate ---------------------------------------

fit.bz.dzz.rate <- change_point(d = d.bz, write.dir = adt.dir, 
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
  ylab(paste0(DZZ, " mortality/100,000 ppl")) +
  theme_classic() + 
  theme(panel.grid.major = element_line(color = "gray95"),
        axis.title = element_text(size = 10),
        axis.text = element_text(size = 10))
print(plot.fit.bz.dzz.rate)
ggsave(filename = paste0(adt.dir, "Change-point-", DZZ, "-rate.png"), 
       plot = plot.fit.bz.dzz.rate,
       width = 6, height = 4, units = "in", dpi = 300)

# Calculate the time-detrended measles mortality rate:
d.bz$dzz_rate_dt <- d.bz$dzz_rate - fitted(fit.bz.dzz.rate)[, "fitted"]





# 5. Visualization: time seires ---------------------------------------------------------

source(paste0(master.dir, "Function - time series and correlation plots.R"))

# (a) Time series of the rate ---------------------------------

ts_ori_rate <- time_series(d.ts = d.bz, 
                           meV_var = "meV_mort_rate", nonmeV_var = "dzz_rate", 
                           cp = cp.bz.rate.year, 
                           data_type = c("rate", "rate"), 
                           data_manipulation = "original", 
                           dzz = DZZ)
print(ts_ori_rate)
ggsave(filename = paste0(adt.dir, "Time series-original-rate.png"), 
       plot = ts_ori_rate,
       width = 7, height = 4, units = "in", dpi = 300)




# (b) Correlation between measles and non-measles mortality rate --------

cor_ori_rate <- correlation_plot(d.ts = d.bz, 
                                 meV_var = "meV_mort_rate", nonmeV_var = "dzz_rate",  
                                 data_type = c("rate", "rate"), 
                                 data_manipulation = "original", 
                                 dzz = DZZ)
print(cor_ori_rate)
ggsave(filename = paste0(adt.dir, "Correlation-original-rate.png"), 
       plot = cor_ori_rate,
       width = 6, height = 4, units = "in", dpi = 300)




# (c) Time series of the detrended mortality rate ------------

ts_dt_rate <- time_series(d.ts = d.bz, 
                          meV_var = "meV_mort_rate_dt", nonmeV_var = "dzz_rate_dt", 
                          cp = cp.bz.rate.year, 
                          data_type = c("rate", "rate"), 
                          data_manipulation = "detrend", 
                          dzz = DZZ)
print(ts_dt_rate)
ggsave(filename = paste0(adt.dir, "Time series-detrend-rate.png"),
       plot = ts_dt_rate,
       width = 7, height = 4, units = "in", dpi = 300)





# (d) Correlation between detrended measles and non-measles mortality rate --------

cor_dt_rate <- correlation_plot(d.ts = d.bz, 
                                meV_var = "meV_mort_rate_dt", nonmeV_var = "dzz_rate_dt",  
                                data_type = c("rate", "rate"), 
                                data_manipulation = "detrend", 
                                dzz = DZZ)
print(cor_dt_rate)
ggsave(filename = paste0(adt.dir, "Correlation-detrended-rate.png"), 
       plot = cor_dt_rate,
       width = 6, height = 4, units = "in", dpi = 300)





# 6. Linear regression ----------------------------------------------------

print(names(d.bz))

# (a) Full time series: original data ------------------------

# 1) Use only measles mortality as the predictor:
lm.bz.1.rate <- lm(dzz_rate ~ meV_mort_rate, data = d.bz)
print(lm.bz.1.rate)
summary(lm.bz.1.rate)

# 2) Use only year as the predictor:
lm.bz.2.rate <- lm(dzz_rate ~ year.idx, data = d.bz)
print(lm.bz.2.rate)
summary(lm.bz.2.rate)

# 3) Add year and measles mortality as the predictor:
lm.bz.3.rate <- lm(dzz_rate ~ year.idx + meV_mort_rate, data = d.bz)
print(lm.bz.3.rate)
summary(lm.bz.3.rate)
anova(lm.bz.3.rate, lm.bz.1.rate)  # Effects of time
anova(lm.bz.3.rate, lm.bz.2.rate)  # Effects of measles mortality

# 4) Add interaction term between year and measles mortality:
lm.bz.4.rate <- lm(dzz_rate ~ year.idx * meV_mort_rate, data = d.bz)
print(lm.bz.4.rate)
summary(lm.bz.4.rate)
anova(lm.bz.4.rate, lm.bz.3.rate)  # Effects of the interaction

# Put all regression results in a list
lm.bz.list.rate <- list(lm.bz.1.rate, lm.bz.2.rate, lm.bz.3.rate, lm.bz.4.rate)

# Model selection:
step.bz.rate <- step(lm.bz.4.rate, direction = "both")
lm.bz.rate <- lm.bz.list.rate[[which.min(sapply(X = lm.bz.list.rate, FUN = AIC))]]

# Examine the best model: rate
coefficients(lm.bz.rate)
confint(lm.bz.rate, level=0.95)
summary(lm.bz.rate)
anova(lm.bz.rate)

# Diagnostics of the best model:
plot(lm.bz.rate$residuals ~ lm.bz.rate$fitted.values)
qqnorm(lm.bz.rate$residuals); abline(a = 0, b = 1, col = "red")

# Relative importance of the predictors
relimp.bz.est <- calc.relimp(lm.bz.3.rate, type=c("lmg"), rela=TRUE)
print(relimp.bz.est)


# write out the key results:
sink(paste0(adt.dir, "Linear regression - full time series - original.txt"))

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
          dep.var.labels = paste0("Dependent variable: ", DZZ, " mortality rate"),
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

# 2) Use only year as the predictor:
lm.bz.dt.2.rate <- lm(dzz_rate_dt ~ year.idx, data = d.bz)
print(lm.bz.dt.2.rate)
summary(lm.bz.dt.2.rate)

# 3) Add year and measles mortality as the predictor:
lm.bz.dt.3.rate <- lm(dzz_rate_dt ~ year.idx + meV_mort_rate_dt, data = d.bz)
print(lm.bz.dt.3.rate)
summary(lm.bz.dt.3.rate)
anova(lm.bz.dt.3.rate, lm.bz.dt.1.rate)  # Effects of time
anova(lm.bz.dt.3.rate, lm.bz.dt.2.rate)  # Effects of measles mortality

# 4) Add interaction term between year and measles mortality:
lm.bz.dt.4.rate <- lm(dzz_rate_dt ~ year.idx * meV_mort_rate_dt, data = d.bz)
print(lm.bz.dt.4.rate)
summary(lm.bz.dt.4.rate)
anova(lm.bz.dt.4.rate, lm.bz.dt.3.rate)  # Effects of the interaction

# Put all regression results in a list
lm.bz.dt.list.rate <- list(lm.bz.dt.1.rate, lm.bz.dt.2.rate, lm.bz.dt.3.rate, lm.bz.dt.4.rate)

# Model selection:
step.bz.dt.rate <- step(lm.bz.dt.4.rate, direction="both")
lm.bz.dt.rate <- lm.bz.dt.list.rate[[which.min(sapply(X = lm.bz.dt.list.rate, FUN = AIC))]]

# Examine the best model: rate
coefficients(lm.bz.dt.rate)
confint(lm.bz.dt.rate, level=0.95)
summary(lm.bz.dt.rate)
anova(lm.bz.dt.rate)

# Diagnostics of the best model:
plot(lm.bz.dt.rate$residuals ~ lm.bz.dt.rate$fitted.values)
qqnorm(lm.bz.dt.rate$residuals); abline(a = 0, b = 1, col = "red")

# Relative importance of the predictors
relimp.bz.dt.est <- calc.relimp(lm.bz.dt.3.rate, type=c("lmg"), rela=TRUE)
print(relimp.bz.dt.est)



# write out the key results:
sink(paste0(adt.dir, "Linear regression - full time series - detrended.txt"))

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
          dep.var.labels = paste0("Dependent variable: detrended ", DZZ, " mortality rate"),
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





# 7. Difference between consective years -----------------------------------

d.bz.matrix <- data.matrix(d.bz[, !(names(d.bz) %in% c("age", "city", "dzzName"))])
diff.bz <- diff(d.bz.matrix)
diff.bz <- as.data.frame(diff.bz, stringsAsFactors = FALSE)

diff.bz$start_year <- d.bz$year[1:nrow(diff.bz)]
diff.bz$end_year <- d.bz$year[2:nrow(d.bz)]
diff.bz$year <- diff.bz$end_year   # Use the latter year in the two consecutive years


# (a) Time series: mortality rate ----------------------------------------

plot.diff.bz.ts.rate <- time_series(d.ts = diff.bz, 
                                    meV_var = "meV_mort_rate", nonmeV_var = "dzz_rate", 
                                    cp = cp.bz.rate.year, 
                                    data_type = c("rate", "rate"), 
                                    data_manipulation = "difference", 
                                    year_breaks = seq(1980, 1995, 5),
                                    dzz = DZZ)
print(plot.diff.bz.ts.rate)
ggsave(filename = paste0(adt.dir, "Time series-difference-rate.png"), 
       plot = plot.diff.bz.ts.rate,
       width = 7, height = 4, units = "in", dpi = 300)




# (b) Visualization of correlation: mortality rate -----------------------

plot.diff.bz.cor.rate <- correlation_plot(d.ts = diff.bz, 
                                          meV_var = "meV_mort_rate", nonmeV_var = "dzz_rate",  
                                          data_type = c("rate", "rate"), 
                                          data_manipulation = "difference", 
                                          year_breaks = seq(1980, 1995, by = 5), 
                                          year_range = c(1980, 1995), 
                                          dzz = DZZ)
print(plot.diff.bz.cor.rate)
ggsave(filename = paste0(adt.dir, "Difference-correlation-rate.png"), 
       plot = plot.diff.bz.cor.rate,
       width = 6, height = 4, units = "in", dpi = 300)




# (c) Linear regression ---------------------------------------------------

diff.bz$year.idx <- diff.bz$year - min(diff.bz$year) + 1

# 1) Use only measles mortality as the predictor:
lm.bz.diff.1.rate <- lm(dzz_rate ~ meV_mort_rate, data = diff.bz)
print(lm.bz.diff.1.rate)
summary(lm.bz.diff.1.rate)

# 2) Use only year as the predictor:
lm.bz.diff.2.rate <- lm(dzz_rate ~ year.idx, data = diff.bz)
print(lm.bz.diff.2.rate)
summary(lm.bz.diff.2.rate)

# 3) Add year and measles mortality as the predictor:
lm.bz.diff.3.rate <- lm(dzz_rate ~ year.idx + meV_mort_rate, data = diff.bz)
print(lm.bz.diff.3.rate)
summary(lm.bz.diff.3.rate)
anova(lm.bz.diff.3.rate, lm.bz.diff.1.rate)  # Effects of time
anova(lm.bz.diff.3.rate, lm.bz.diff.2.rate)  # Effects of measles mortality

# 4) Add interaction term between year and measles mortality:
lm.bz.diff.4.rate <- lm(dzz_rate ~ year.idx * meV_mort_rate, data = diff.bz)
print(lm.bz.diff.4.rate)
summary(lm.bz.diff.4.rate)
anova(lm.bz.diff.4.rate, lm.bz.diff.3.rate)  # Effects of the interaction

# Put all regression results in a list
lm.bz.diff.list.rate <- list(lm.bz.diff.1.rate, lm.bz.diff.2.rate, lm.bz.diff.3.rate, lm.bz.diff.4.rate)

# Model selection:
step.bz.diff.rate <- step(lm.bz.diff.4.rate, direction="both")
lm.bz.diff.rate <- lm.bz.diff.list.rate[[which.min(sapply(X = lm.bz.diff.list.rate, FUN = AIC))]]

# Examine the best model: rate
coefficients(lm.bz.diff.rate)
confint(lm.bz.diff.rate, level=0.95)
summary(lm.bz.diff.rate)
anova(lm.bz.diff.rate)

# Diagnostics of the best model:
plot(lm.bz.diff.rate$residuals ~ lm.bz.diff.rate$fitted.values)
qqnorm(lm.bz.diff.rate$residuals); abline(a = 0, b = 1, col = "red")

# Relative importance of the predictors
relimp.bz.diff.est <- calc.relimp(lm.bz.diff.3.rate, type=c("lmg"), rela=TRUE)
print(relimp.bz.diff.est)



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
          dep.var.labels = paste0("Dependent variable: ", DZZ, " mortality rate"),
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





# 8. Negative binomial regression ------------------------------------------

# The data is over-dispersed: use NB model instead of Poisson model
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


# Examine the effects of each variable
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
print(nbglm.bz.aic)   # The best model has the lowest AIC


# Find the best model:
nbglm.bz.list <- list(nbglm.bz.1, nbglm.bz.2, nbglm.bz.3, nbglm.bz.4,
                      nbglm.bz.5, nbglm.bz.6, nbglm.bz.7, nbglm.bz.8)
nbglm.bz.best <- nbglm.bz.list[[which.min(nbglm.bz.aic)]]
summary(nbglm.bz.best)


# Stepwise model selection:
nbglm.bz.step <- step(nbglm.bz.8, direction="both")


# Some further exploration of the best model:
summary(nbglm.bz.best) # summary
anova(nbglm.bz.best)   # anova table


# 95% confidence interval of the NB model coefficient of the full model:
nbglm.bz.best.cofficient <- cbind(estimate = coef(nbglm.bz.best), confint(nbglm.bz.best))
print(nbglm.bz.best.cofficient)


# Convert back to the empirical scale: multiplicable
nbglm.bz.best.cofficient.ir <- exp(nbglm.bz.best.cofficient)
print(nbglm.bz.best.cofficient.ir)


# Calculate the pseudo-r2 using ifferent methods:
nbglm.bz.pseudo.R2 <- c(pscl::pR2(nbglm.bz.best), unlist(modEvA::RsqGLM(nbglm.bz.best)))
print(nbglm.bz.pseudo.R2)


# Checking model assumptions about over-dispersion: compare with the Poission model
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
          dep.var.labels = paste0("Dependent variable: ", DZZ, " mortality"),
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





# 9. Summarizing results and save the workspace --------------------------------


# (a) Save key figure objects ----------------------------------------------

Brazil_1_9_meV_noninf <- list(ts_ori = ts_ori_rate,                            # time series
                              ts_dt = ts_dt_rate, 
                              ts_diff = plot.diff.bz.ts.rate,
                              cor_ori = cor_ori_rate,                           # correlation plot
                              cor_dt = cor_dt_rate, 
                              cor_diff = plot.diff.bz.cor.rate,
                              cp_meV = plot.fit.bz.rate,                # change-point time series
                              cp_dzz = plot.fit.bz.dzz.rate,
                              lm_1 <- lm.bz.1.rate,                     # selected regression results  
                              lm_2 <- lm.bz.3.rate, 
                              lm_3 <- lm.bz.dt.1.rate,
                              lm_4 <- lm.bz.diff.1.rate, 
                              lm_5 <- nbglm.bz.6)
save(Brazil_1_9_meV_noninf, 
     file = paste0(adt.dir, LOC, "_", AGE, "_meV_", DZZ, "_export_", Sys.Date(), ".Rdata"))





# (b) Combined figures: times series and correlation ---------------------------

# Adjust figure margins
ts_cor_adjust <- theme(legend.position = "none") + 
  theme(plot.margin = unit(c(15, 20, 5, 10), units = "pt")) + 
  theme(plot.title = element_text(hjust = 0.5, vjust = 1, size = 14))
ts_cor_legend_adjust <- theme(legend.box.margin = ggplot2::margin(0, 0, 0, 0))

# Extract the legends
ts_legend <- cowplot::get_legend(Brazil_1_9_meV_noninf$ts_ori + ts_cor_legend_adjust)
cor_legend <- cowplot::get_legend(Brazil_1_9_meV_noninf$cor_ori + ts_cor_legend_adjust)

# Generate the combined figures
ts_cor_combined_plot <- plot_grid(Brazil_1_9_meV_noninf$ts_ori + ts_cor_adjust + 
                                    ggtitle("Original mortality\n "),
                                  Brazil_1_9_meV_noninf$ts_dt + ts_cor_adjust + 
                                    ggtitle("Detrended mortality\n "),
                                  Brazil_1_9_meV_noninf$ts_diff + ts_cor_adjust + 
                                    ggtitle("Mortality difference\n "),
                                  ts_legend, 
                                  Brazil_1_9_meV_noninf$cor_ori + ts_cor_adjust,
                                  Brazil_1_9_meV_noninf$cor_dt + ts_cor_adjust,
                                  Brazil_1_9_meV_noninf$cor_diff + ts_cor_adjust,
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





# (c) Save the workspace --------------------------------------------------
save.image(paste0(adt.dir, LOC, "_", AGE, "_meV_", DZZ, "_", Sys.Date(), ".Rdata"))
