IM <- function(d,
               target.dir,
               target_disease = "all",
               target_city = "Brazil",
               incidence_or_mortality = "mortality",
               count_or_rate = "rate",
               control_year = TRUE,
               max_duration = 5,
               difference = FALSE){

  # d: the main data as a data frame (subset of the data.frame "infection")
  # target.dir: directory to export the results
  # target_disease: the category of non-measles disaese of interest (choose from the vector diseases)
  # target_city: the city of interest or the entire Brazil (choose from the vector "cities")
  # incidence_or_mortality: use incidence or mortality of the measles data
  # count_or_rate: use rate or absolute case count
  # control_year: whether to include year as a covariant: TRUE = control year
  # max_duration: maximum number of years for immunomodulation
  # difference: if TRUE, use the difference between consecutive years

  
  
  # a) Preparation -------------------------------------------------------------
  
  nyear <- nrow(d)

  # Labels for output
  year_label <- ifelse(control_year, "control_year", "ignore_year")
  output_label <- paste(incidence_or_mortality, count_or_rate, ifelse(difference, "diff", "ori"),  year_label, sep = "_")

  # Calculate immune anmesia: depending on "incidence_or_mortality"
  if(incidence_or_mortality == "incidence"){
    for(i in 1:max_duration) d[, paste0("ia", i)] <- c(rep(0, i-1), d$meV_inc[1:(nyear-i+1)])
  }else{
    for(i in 1:max_duration) d[, paste0("ia", i)] <- c(rep(0, i-1), d$meV_mort[1:(nyear-i+1)])
  }

  # Calculation preparation:
  d.ia.matrix <- as.matrix(d[, paste0("ia", 1:max_duration)])

  # gamma weight:
  gamma_adjustment <- 1 - pgamma(seq(from = 6, by = 12, length.out = max_duration),
                                 shape = 35, rate = 1.236)



  # b) Calculate immunomodulation --------------------------------------------

  # A data frame to record the R2 with different duration of immune amnesia:
  R2_table <- data.frame(duration = 1:max_duration,
                         start_year = rep(-1, max_duration),
                         end_year = rep(-1, max_duration),
                         n = rep(-1, max_duration),
                         beta_additive = rep(-1, max_duration),
                         p_additive = rep(-1, max_duration),
                         R2_additive = rep(-1, max_duration),
                         beta_gamma = rep(-1, max_duration),
                         p_gamma = rep(-1, max_duration),
                         R2_gamma = rep(-1, max_duration),
                         beta_annual = rep(-1, max_duration),
                         p_annual = rep(-1, max_duration),
                         R2_annual = rep(-1, max_duration),
                         optimal = rep(-1, max_duration))
  
  # A list to store the full analysis results of each duration 
  ia.list <- vector(mode = "list", length = max_duration)


  for(i in 1:max_duration){

    ### Calculate the accumulated immune amnesia:
    d.ia <- d
    inclusion <- c(rep(1, i), rep(0, max_duration-i))
    d.ia$additive <- as.numeric(d.ia.matrix %*% inclusion)
    d.ia$gamma <- as.numeric(d.ia.matrix %*% (inclusion * gamma_adjustment))

    # Remove the first "max_duration" year where the accumulation of IA is not completed for all possible durations of immunomodulation
    d.ia <- d.ia[max_duration:nyear, ]

    # Convert to rate per 100,000 people
    if(count_or_rate == "rate"){
      d.ia$additive <- d.ia$additive / d.ia$pop * 100000
      d.ia$gamma <- d.ia$gamma / d.ia$pop * 100000

      if(incidence_or_mortality == "incidence"){
        d.ia$annual <- d.ia$meV_inc_rate
      }else{
        d.ia$annual <- d.ia$meV_mort_rate
      }
      d.ia$dzz <- d.ia$dzz_rate
    }else{
      if(incidence_or_mortality == "incidence"){
        d.ia$annual <- d.ia$meV_inc
      }else{
        d.ia$annual <- d.ia$meV_mort
      }
      d.ia$dzz <- d.ia$dzz_mort
    }
    
    
    
    ### Using the original quantities or difference between consecutive years
    if(difference){
      d.ia.diff <- d.ia %>% 
        dplyr::select(-c(age, city, dzzName)) %>% 
        data.matrix() %>%
        diff() %>%
        as.data.frame()
      d.ia.diff$year <- d.ia$year[2:nrow(d.ia)]
      d.ia <- d.ia.diff
    }


    ### Linear models to estimate R-squared:
    # 1. Annual:
    if(control_year){
      lm.ia.annual <- lm(dzz ~ year + annual, data = d.ia)
    }else{
      lm.ia.annual <- lm(dzz ~ annual, data = d.ia)
    }
    lm.temp <- summary(lm.ia.annual)
    if("annual" %in% row.names(lm.temp$coefficients)){
      R2_table$beta_annual[i] <- lm.temp$coefficients["annual", "Estimate"]
      R2_table$p_annual[i] <- lm.temp$coefficients["annual", "Pr(>|t|)"]
      R2_table$R2_annual[i] <- summary(lm.ia.annual)$r.squared
    }else{
      R2_table$beta_annual[i] <- NA
      R2_table$p_annual[i] <- NA
      R2_table$R2_annual[i] <- summary(lm.ia.annual)$r.squared
    }

    # 2. Additive IM:
    if(control_year){
      lm.ia.additive <- lm(dzz ~ year + additive, data = d.ia)
    }else{
      lm.ia.additive <- lm(dzz ~ additive, data = d.ia)
    }
    lm.temp <- summary(lm.ia.additive)
    if("additive" %in% row.names(lm.temp$coefficients)){
      R2_table$beta_additive[i] <- lm.temp$coefficients["additive", "Estimate"]
      R2_table$p_additive[i] <- lm.temp$coefficients["additive", "Pr(>|t|)"]
      R2_table$R2_additive[i] <- summary(lm.ia.additive)$r.squared
    }else{
      R2_table$beta_additive[i] <- NA
      R2_table$p_additive[i] <- NA
      R2_table$R2_additive[i] <- summary(lm.ia.additive)$r.squared
    }

    # 3. Gamma IM:
    if(control_year){
      lm.ia.gamma <- lm(dzz ~ year + gamma, data = d.ia)
    }else{
      lm.ia.gamma <- lm(dzz ~ gamma, data = d.ia)
    }
    lm.temp <- summary(lm.ia.gamma)
    if("gamma" %in% row.names(lm.temp$coefficients)){
      R2_table$beta_gamma[i] <- lm.temp$coefficients["gamma", "Estimate"]
      R2_table$p_gamma[i] <- lm.temp$coefficients["gamma", "Pr(>|t|)"]
      R2_table$R2_gamma[i] <- summary(lm.ia.gamma)$r.squared
    }else{
      R2_table$beta_gamma[i] <- NA
      R2_table$p_gamma[i] <- NA
      R2_table$R2_gamma[i] <- summary(lm.ia.gamma)$r.squared
    }

    
    # Determine which method results in the highest R2
    R2_table$optimal[i] <- c("additive", "gamma", "annual")[which.max(R2_table[i, c("R2_additive", "R2_gamma", "R2_annual")])]
    
    # Add some basic information
    R2_table$n[i] <- nrow(d.ia)
    R2_table$start_year[i] <- min(d.ia$year)
    R2_table$end_year[i] <- max(d.ia$year)

    # Record the analysis results:
    ia.list[[i]] <- d.ia



    ### Preparation for the figures
    xlabels <- c("Measles", "Additive", "Gamma")
    if(incidence_or_mortality == "incidence"){
      xlabels <- paste(xlabels, "incidence", sep = " ")
    }else{
      xlabels <- paste(xlabels, "mortality", sep = " ")
    }
    ylabels <- "Non-measlse mortality"
    if(difference){
      xlabels <- paste0(xlabels, " difference")
      ylabels <- paste0(ylabels, " difference")
    }
    if(count_or_rate == "rate"){
      xlabels <- paste0(xlabels, "/100,000 ppl")
      ylabels <- paste0(ylabels, "/100,000 ppl")
    }
    names(xlabels) <- c("annual", "additive", "gamma")
    if(control_year) ylabels <- paste0(ylabels, " (residual)")
  }


  ### Figure: Change of R2 with different duration of immunomodulation
  R2_table_long <- R2_table %>%
    dplyr::filter(duration <= max_duration) %>%
    dplyr::select(duration, R2_additive, R2_gamma, R2_annual) %>%
    reshape2::melt(id.var = "duration") %>%
    dplyr::mutate(variable = ordered(variable,
                              levels = c("R2_annual", "R2_additive", "R2_gamma"),
                              labels = c("annual measles mortality", 
                                         "additive immunomodulation",
                                         "gamma immunomodulation")))

  # Make the figure:
  f <- ggplot(data = R2_table_long, aes(x = duration, y = value, color = variable)) +
    geom_line(size = 1.2) +
    scale_color_manual(name = "non-measles mortality vs.", values = cbPalette[1:4]) +
    scale_x_continuous(name = "Duration of immunomodulation (year)") +
    ylab("R-squared") +
    theme_bw() +
    theme(panel.border = element_blank(),
          panel.grid.minor = element_blank(),
          axis.line = element_line(),
          axis.text = element_text(size = 10),
          axis.title = element_text(size = 12))
  print(f)
  ggsave(filename = paste0(target.dir, "Change of R2-", output_label, ".png"), plot = f,
         width = 5, height = 4, units = "in", dpi = 320)


  save.image(paste0(target.dir, output_label, "_", Sys.Date(), ".RData"))
  write.csv(x = R2_table, file = paste0(target.dir, "Regression results-", output_label,".csv"), row.names = FALSE)
  
  return(list(R2_table, f))
}