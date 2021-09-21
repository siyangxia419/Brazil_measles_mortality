### Analysis of optimal duration of immunomodulation after measles infection
###
### Input: Mortality data prepared to the desired format
### The data shoule locate in the directory "prepared data/" inside the location directory
### 
### Siyang Xia
### 2021.4.21



library(reshape2)
library(ggplot2)
library(dplyr)



### GLOBAL OPTIONS:
LOC <- "Brazil"
AGE <- "1-9"
DZZ <- "all"  # For immunomodulation analysis, sum all non-measles diseases




### Some miscellaneous preparation:
# Color-blind friendly palette:
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")



# Function to calculate proportional location in a range:
position_in_range <- function(x, p) min(x) + p * (max(x) - min(x))



# 1. Data input -----------------------------------------------------------

master.dir <- "C:/Users/Siyang Xia/Dropbox/My laptop/Mina lab/Projects/Measles/Empirical/Brazil/"
data.dir <- paste0(master.dir, LOC, "/prepared data/")
output.dir <- paste0(master.dir, LOC, "/immunomodulation analysis/")
if(!dir.exists(output.dir)) dir.create(output.dir)


# Options of geographic locations and non-measles disease:
cities <- c("Brazil", "B horizonte", "Belem", "Fortaleza", "Goiania", "Manaus", "Porto Alegre", 
            "Recife", "Rio de Janeiro", "Salvador", "Sau Paulo")
diseases <- c("all", "TB", "URT", "bronchitis", "diarrheal", "meningitis", "pneumonia", 
              "pneumococcus", "rubella", "septicemia", "LRT", "asthma", "emphysema", "meningococcal") 

# Read the data:
infection <- read.csv(paste0(data.dir, LOC, "_measles_non-measles_mortality.csv"), 
                      header = TRUE, stringsAsFactors = FALSE)

# Some double-check: all should be TRUE
all(round(infection$meV_mort / infection$meV_CFR - infection$meV_inc, digits = 10) == 0)
all(round(infection$meV_mort / infection$pop * 100000 - infection$meV_mort_rate, digits = 10) == 0)
all(round(infection$meV_inc / infection$pop * 100000 - infection$meV_inc_rate, digits = 10) == 0)
all(round(infection$dzz_mort / infection$pop * 100000 - infection$dzz_rate, digits = 10) == 0)





# 2. Immunomodulation analysis --------------------------------------------- 

# Load the function for immunomodulation analysis:
source(paste0(master.dir, "Function - immunomodulation.R"))


# Select the data with the target location, age group and non-measles disease category:
d.im <- infection %>% dplyr::filter(city == LOC, age == AGE, dzzName == DZZ)
nyear <- nrow(d.im)


### Immunomodulation analysis with different options

op1 <- "mortality" 
op2 <- "rate"

# Focus on measles mortality instead of measles incidence as data of latter is not available

for(op3 in c(TRUE, FALSE)){
  for(op4 in c(TRUE, FALSE)){
    IM(d = d.im,
       target.dir = output.dir,
       target_disease = DZZ,
       target_city = LOC,
       incidence_or_mortality = op1,
       count_or_rate = op2,
       control_year = op3, 
       max_duration = 5, 
       difference = op4)
    print(paste(op1, op2, op3, op4))
  }
}
