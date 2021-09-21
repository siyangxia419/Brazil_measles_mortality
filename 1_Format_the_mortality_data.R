### Format the Brazilian mortality data
###
### Input: Mortality data download from Brazil Ministry of Health, DATASUS
### http://tabnet.datasus.gov.br/cgi/deftohtm.exe?sim/cnv/obt09br.def
### The data shoule locate in the directory "raw data/" inside the location directory
### 
### Siyang Xia
### 2021.4.28



library(dplyr)
library(reshape2)
library(tidyr)



### GLOBAL OPTIONS:
LOC <- "Brazil"
TIME_PERIOD <- "1979_1995"  # alternative: "1996_2019"



### directories for the specific location
rawdata.dir <- "./raw data/"   # directory of the downloaded data
out.dir <- "./prepared data/"  # directory to store the prepared data
if(!dir.exists(out.dir)) dir.create(out.dir)
global.dir <- "../"            # directory of some global input files



# 1. Demographic information ----------------------------------------------

# Population size was downloaded from Brazil Ministry of Health DATASUS
# http://tabnet.datasus.gov.br/cgi/tabcgi.exe?ibge/cnv/popbr.def

pop <- read.csv(paste0(rawdata.dir, "Population_1980_2012_both_", LOC, "_all.csv"), skip = 4,
                sep = ";", header = FALSE, stringsAsFactors = FALSE)
names(pop) <- c("Age", 1980:2012)
pop[, -1] <- apply(X = pop[, -1], MARGIN = 2, FUN = as.numeric)


# For some years, only age 0-4 was available instead of age 0-1 and 1-4
# We estimate the population size in the two age groups by fitting a linear 
# change of the decrease between 1980-1991 and 1991-1993
ratio_1980 <- pop$`1980`[2] / (pop$`1980`[2] + pop$`1980`[3])
ratio_1991 <- pop$`1991`[2] / (pop$`1991`[2] + pop$`1991`[3])
ratio_1993 <- pop$`1993`[2] / (pop$`1993`[2] + pop$`1993`[3])
ratio_1981_1990 <- ratio_1980 + (1:10) * (ratio_1991 - ratio_1980) / 11
ratio_1992 <- (ratio_1991 + ratio_1993) / 2
ratio_combined <- c(ratio_1981_1990, ratio_1992)
names(ratio_combined) <- c(1981:1990, 1992)
for(i in 2:ncol(pop)){
  if(is.na(pop[1, i])) pop[1, i] <- pop[2, i] + pop[3, i]
  if(is.na(pop[2, i])){
    pop[2, i] <- pop[1, i] * ratio_combined[names(pop)[i]]
    pop[3, i] <- pop[1, i] * (1 - ratio_combined[names(pop)[i]])
  }
}

# Translate the age categories:
pop$Age <- gsub(pattern = "Idade ignorada", replacement = "Unknown", x = 
                  gsub(pattern = "80 anos e mais", replacement = "80+", 
                       gsub(pattern = "Menor 1 ano", replacement = "0-1", 
                            gsub(pattern = " a ", replacement = "-", 
                                 gsub("([0-9] a [0-9]+).*$", "\\1", pop$Age)))))

# Add additional age categories: 1-9, 10-19, 10+, 20+
print(pop$Age)
pop.add <- pop[1:5, ]
pop.add$Age <- c("1-9", "10-19", "10+", "20+", "20-59")
pop.add[1, -1] <- colSums(pop[pop$Age %in% c("1-4", "5-9"), -1])
pop.add[2, -1] <- colSums(pop[pop$Age %in% c("10_14", "15-19"), -1])
pop.add[3, -1] <- colSums(pop[pop$Age %in% c("10-14", "15-19", "20-29", "30-39", "40-49", 
                                             "50-59", "60-69", "70-79", "80+"), -1])
pop.add[4, -1] <- colSums(pop[pop$Age %in% c("20-29", "30-39", "40-49", "50-59", "60-69", "70-79", "80+"), -1])
pop.add[5, -1] <- colSums(pop[pop$Age %in% c("20-29", "30-39", "40-49", "50-59"), -1])
pop <- rbind(pop, pop.add)
rm(pop.add)

# Convert to the long form
pop.long <- melt(pop, id.vars = "Age", variable.name = "Year", value.name = "Population") %>%
  relocate(Year, Age, Population) %>%
  mutate(Year = as.numeric(as.character(Year)))
pop.long$Age <- gsub(pattern = "Total", replacement = "all", x = pop.long$Age)

write.csv(x = pop.long, file = paste0(out.dir, LOC, "_population_formatted.csv"), 
          quote = FALSE, row.names = FALSE)





# 2. Disease mortality ----------------------------------------------------


# Read all data files (each file is for a different age group)
data.file.list <- list.files(path = rawdata.dir, pattern = "Mortality_1979_1995_*", full.names = TRUE)

# The first four lines in the data file for individual age group are titles, while for the data file
# that includes all age group, the top three lines are titles
data.file.allage <- data.file.list[grep(pattern = "_all", x = data.file.list)]
data.file.eachage <- data.file.list[-grep(pattern = "_all", x = data.file.list)]

data.list <- lapply(X = data.file.eachage, 
                    FUN = function(x) read.csv(x, header = TRUE, skip = 4,
                                               sep = ";", stringsAsFactors = FALSE))
data.list[[length(data.list) + 1]] <- read.csv(data.file.allage, header = TRUE, skip = 3,
                                               sep = ";", stringsAsFactors = FALSE)

names(data.list) <- gsub(pattern = "_", replacement = "-", 
                         x = gsub(pattern = paste0("./raw data/Mortality_1979_1995_both_", LOC, "_"), 
                                  replacement = "", 
                                  x = gsub(pattern = ".csv", replacement = "", 
                                           x = c(data.file.eachage, data.file.allage))))

# Combine into a single data frame
data.full <- bind_rows(data.list, .id = "Age")
data.full$ICD9.number <- sapply(X = data.full$Categoria.CID.9, 
                                FUN = function(x) as.numeric(strsplit(x, split = " - ")[[1]][1]))


# Disease of interest:
DZZ.list <- read.csv(paste0(global.dir, "Disease of interest.csv"), header = TRUE, stringsAsFactors = FALSE)


### Select the diseases of interest in the full data set:
d <- data.full %>% filter(ICD9.number %in% DZZ.list$ICD.9_Number)
d$disease.name <- DZZ.list$Infection_Label[match(d$ICD9.number, DZZ.list$ICD.9_Number)]


### Some cleaning and formatting
print(names(d))
d.long <- d %>% melt(id.var = c("Age", "disease.name", "ICD9.number", "Categoria.CID.9"), 
                     variable.name = "Year_ori", value.name = "Mortality_ori") %>%
  filter(Year_ori != "Total") %>%
  mutate(Year = as.numeric(gsub(pattern = "X", replacement = "", x = Year_ori)),
         Mortality = as.numeric(Mortality_ori)) %>%
  select(-c(Year_ori, Mortality_ori))

# Full combinations of age, year and diseases (i.e. ICD-9)
d.long2 <- d.long %>% expand(Age, nesting(ICD9.number, Categoria.CID.9, disease.name), Year)

# Convert NAs to zeros
d.long <- d.long %>% right_join(d.long2) %>%
  mutate(Mortality = ifelse(is.na(Mortality), 0, Mortality))


### Add population information
d.long <- d.long %>% left_join(pop.long)


### Calculate mortality rate
d.long$Mortality.rate <- d.long$Mortality / d.long$Population * 100000





# 3. Summarize by disease categories --------------------------------------

d.long.sum <- d.long %>% group_by(Age, disease.name, Year) %>%
  summarise(Mortality = sum(Mortality), 
            Population = mean(Population),
            Mortality.rate = sum(Mortality.rate)) %>%
  ungroup()

write.csv(x = d.long.sum, file = paste0(out.dir, LOC, "_individual_disease_mortality.csv"), 
          quote = FALSE, row.names = FALSE)





# 4. Convert to the desired format ----------------------------------------
# For the downstream analysis between measles and non-measles mortality

# a) Measles
measles.long <- d.long.sum %>% filter(disease.name == "measles")
measles.format <- measles.long %>% select(-disease.name) %>%
  rename(age = Age, year = Year, meV_mort = Mortality, 
         pop = Population, meV_mort_rate = Mortality.rate)


# b) Other diseases
dzz.long <- d.long.sum %>% filter(disease.name != "measles")

# All disease categories
dzz.category <- list(
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
  URT = c("URT"), 
  all = unique(dzz.long$disease.name)
)


# A list to store the data for each non-measles disease group
dzz.category.long.list <- vector(mode = "list", length = length(dzz.category))
names(dzz.category.long.list) <- names(dzz.category)

# Summarize the non-measles diseases for each disease group:
for(dzz in names(dzz.category)){
  # Select the diseases in the category
  dzz.temp <- dzz.long %>% filter(disease.name %in% dzz.category[[dzz]])
  
  # Sum the mortality
  dzz.category.long.list[[dzz]] <- dzz.temp %>% 
    group_by(Age, Year) %>%
    summarise(Mortality = sum(Mortality), 
              Population = mean(Population),
              Mortality.rate = sum(Mortality.rate)) %>%
    ungroup()
  
  # Bind with the measles information
  dzz.category.long.list[[dzz]] <- dzz.category.long.list[[dzz]] %>%
    rename(age = Age, year = Year, dzz_mort = Mortality, 
           pop = Population, dzz_rate = Mortality.rate) %>%
    mutate(dzzName = dzz, city = LOC) %>% 
    left_join(measles.format) %>%
    relocate("age", "year", "meV_mort", "dzz_mort", "pop", "meV_mort_rate", "dzz_rate", 
             "city", "dzzName")
}
rm(dzz.temp)


### Combine the data list:
measles_dzz <- bind_rows(dzz.category.long.list)
measles_dzz <- measles_dzz %>% filter(!is.na(pop))


### Export the final data
write.csv(x = measles_dzz, file = paste0(out.dir, LOC, "_measles_non-measles_mortality.csv"), 
          quote = FALSE, row.names = FALSE)






# 5. Alternative wide format of the data ----------------------------------

dzz_to_exclude <- c("pneumococcus", "LRT", "respiratory")

d.wide.count <- measles_dzz %>%
  dplyr::select(age, year, meV_mort, dzz_mort, city, pop, dzzName) %>%
  dplyr::filter(!(dzzName %in% dzz_to_exclude)) %>%
  reshape2::dcast(age + city + year + pop + meV_mort ~ dzzName, value.var = "dzz_mort") %>%
  dplyr::rename(measles = meV_mort, all_non_measles = all) %>%
  dplyr::mutate(age = paste("age", age, sep = "_"))

d.wide.rate <- measles_dzz %>%
  dplyr::select(age, year, meV_mort_rate, dzz_rate, city, pop, dzzName) %>%
  dplyr::filter(!(dzzName %in% dzz_to_exclude)) %>%
  reshape2::dcast(age + city + year + pop + meV_mort_rate ~ dzzName, value.var = "dzz_rate") %>%
  dplyr::rename(measles = meV_mort_rate, all_non_measles = all) %>%
  dplyr::mutate(age = paste("age", age, sep = "_"))

save(d.wide.count, d.wide.rate, 
     file = paste0(out.dir, "Brazil wide form prepared data.RData"))


write.csv(x = d.wide.count, 
          file = paste0(out.dir, "Brazil_disease_mortality_count.csv"),
          quote = FALSE, row.names = FALSE)
write.csv(x = d.wide.rate, 
          file = paste0(out.dir, "Brazil_disease_mortality_rate.csv"),
          quote = FALSE, row.names = FALSE)



##### Save the workspace
save.image(paste0(out.dir, LOC, "_format_", Sys.Date(), ".RData"))
