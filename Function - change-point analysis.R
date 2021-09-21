change_point <- function(d, write.dir, ts.var, year.var, LOC, DZZ, 
                         output.txt = TRUE, output.ts = TRUE, output.coeff = TRUE){
  # d: input data that contains at least two columns: incidence/mortality & time/year
  # write.dir: directory to write the analysis results
  # ts.var: name of the column that contains the incidence or mortality data
  # year.var: name of the column that contains time or year
  # LOC: location of interest
  # DZZ: disease of interest
  # output.txt: whether to write the key results into a text file
  # output.ts: whether to generate a time series plot
  # output.coeff: whether to generate a plot with posterior distribution and trace plots
  
  if(!require("mcp")){
    install.packages("mcp"); library(mcp)
  }
  
  d$year.idx <- d[, year.var]
  d$ts.var <- d[, ts.var]
  
  # Model 1: decrease + jointed platau, no temporal auto-correlation
  model.bz1 <- list(
    ts.var ~ 1 + year.idx,  # linear change with year
    ts.var ~ 1 ~ 0          # Jointed platau
  )
  fit.bz1 <- mcp(model = model.bz1, data = d)
  plot(fit.bz1)   # The three posterior density plot represents three independent chains
  summary(fit.bz1)                                                                                                                                                                                                                                                                                                                                                                                                                         
  fit.bz1$loo <- loo(fit.bz1)
  
  
  # Model 2: decrease + jointed platau, with temporal auto-correlation
  model.bz2 <- list(
    ts.var ~ 1 + year.idx + ar(1),  # linear change with year and temporal auto-correlation
    ts.var ~ 1 ~ 0                  # Jointed platau          
  )
  fit.bz2 <- mcp(model = model.bz2, data = d)
  plot(fit.bz2)
  summary(fit.bz2)
  fit.bz2$loo <- loo(fit.bz2)
  
  
  # Model 3: decrease + disjointed platau, no temporal auto-correlation
  model.bz3 <- list(
    ts.var ~ 1 + year.idx,  # linear change with year
    ts.var ~ 1 ~ 1          # disjointed platau
  )
  fit.bz3 <- mcp(model = model.bz3, data = d)
  plot(fit.bz3)   # The three posterior density plot represents three independent chains
  summary(fit.bz3)                                                                                                                                                                                                                                                                                                                                                                                                                         
  fit.bz3$loo <- loo(fit.bz3)
  
  
  # Model 4: decrease + disjointed platau, with temporal auto-correlation
  model.bz4 <- list(
    ts.var ~ 1 + year.idx + ar(1),  # linear change with year and temporal auto-correlation
    ts.var ~ 1 ~ 1                  # Disjointed platau          
  )
  fit.bz4 <- mcp(model = model.bz4, data = d)
  plot(fit.bz4)
  summary(fit.bz4)
  fit.bz4$loo <- loo(fit.bz4)
  
  
  # Model 5: two segments with jointed intercept
  model.bz5 <- list(
    ts.var ~ 1 + year.idx,
    ts.var ~ 1 ~ 0 + year.idx
  )
  fit.bz5 <- mcp(model = model.bz5, data = d)
  plot(fit.bz5)   # The three posterior density plot represents three independent chains
  summary(fit.bz5)                                                                                                                                                                                                                                                                                                                                                                                                                         
  fit.bz5$loo <- loo(fit.bz5)
  
  
  # Model 6: two segments with jointed intercept and temporal autocorrelation
  model.bz6 <- list(
    ts.var ~ 1 + year.idx + ar(1), 
    ts.var ~ 1 ~ 0 + year.idx
  )
  fit.bz6 <- mcp(model = model.bz6, data = d)
  plot(fit.bz6)   # The three posterior density plot represents three independent chains
  summary(fit.bz6)                                                                                                                                                                                                                                                                                                                                                                                                                         
  fit.bz6$loo <- loo(fit.bz6)
  
  
  # Model 7: two segments with disjointed intercept
  model.bz7 <- list(
    ts.var ~ 1 + year.idx,
    ts.var ~ 1 ~ 1 + year.idx
  )
  fit.bz7 <- mcp(model = model.bz7, data = d)
  plot(fit.bz7)   # The three posterior density plot represents three independent chains
  summary(fit.bz7)                                                                                                                                                                                                                                                                                                                                                                                                                         
  fit.bz7$loo <- loo(fit.bz7)
  
  
  # Model 8: two segments with disjointed intercept and temporal autocorrelation
  model.bz8 <- list(
    ts.var ~ 1 + year.idx + ar(1), 
    ts.var ~ 1 ~ 1 + year.idx
  )
  fit.bz8 <- mcp(model = model.bz8, data = d)
  plot(fit.bz8)   # The three posterior density plot represents three independent chains
  summary(fit.bz8)                                                                                                                                                                                                                                                                                                                                                                                                                         
  fit.bz8$loo <- loo(fit.bz8)
  
  
  # Model 9: two segments with jointed intercept
  model.bz9 <- list(
    ts.var ~ 1 + year.idx,
    ts.var ~ 1 ~ 0 + year.idx + ar(1)
  )
  fit.bz9 <- mcp(model = model.bz9, data = d)
  plot(fit.bz9)   # The three posterior density plot represents three independent chains
  summary(fit.bz9)                                                                                                                                                                                                                                                                                                                                                                                                                         
  fit.bz9$loo <- loo(fit.bz9)
  
  
  # Model 10: two segments with jointed intercept and temporal autocorrelation
  model.bz10 <- list(
    ts.var ~ 1 + year.idx + ar(1), 
    ts.var ~ 1 ~ 0 + year.idx + ar(1)
  )
  fit.bz10 <- mcp(model = model.bz10, data = d)
  plot(fit.bz10)   # The three posterior density plot represents three independent chains
  summary(fit.bz10)                                                                                                                                                                                                                                                                                                                                                                                                                         
  fit.bz10$loo <- loo(fit.bz10)
  
  
  # Model 11: two segments with disjointed intercept
  model.bz11 <- list(
    ts.var ~ 1 + year.idx,
    ts.var ~ 1 ~ 1 + year.idx + ar(1)
  )
  fit.bz11 <- mcp(model = model.bz11, data = d)
  plot(fit.bz11)   # The three posterior density plot represents three independent chains
  summary(fit.bz11)                                                                                                                                                                                                                                                                                                                                                                                                                         
  fit.bz11$loo <- loo(fit.bz11)
  
  
  # Model 12: two segments with disjointed intercept and temporal autocorrelation
  model.bz12 <- list(
    ts.var ~ 1 + year.idx + ar(1), 
    ts.var ~ 1 ~ 1 + year.idx + ar(1)
  )
  fit.bz12 <- mcp(model = model.bz12, data = d)
  plot(fit.bz12)   # The three posterior density plot represents three independent chains
  summary(fit.bz12)                                                                                                                                                                                                                                                                                                                                                                                                                         
  fit.bz12$loo <- loo(fit.bz12)
  
  
  # Model 13: No change point but consider temporal autocorrelation
  model.bz13 <- list(
    ts.var ~ 1 + year.idx
  )
  fit.bz13 <- mcp(model = model.bz13, data = d)
  plot(fit.bz13)   # The three posterior density plot represents three independent chains
  summary(fit.bz13)                                                                                                                                                                                                                                                                                                                                                                                                                         
  fit.bz13$loo <- loo(fit.bz13)
  
  
  # Model 14: No change point but consider temporal autocorrelation
  model.bz14 <- list(
    ts.var ~ 1 + year.idx + ar(1)
  )
  fit.bz14 <- mcp(model = model.bz14, data = d)
  plot(fit.bz14)   # The three posterior density plot represents three independent chains
  summary(fit.bz14)                                                                                                                                                                                                                                                                                                                                                                                                                         
  fit.bz14$loo <- loo(fit.bz14)
  
  
  # Put results in a list
  fit.bz.list <- list(fit.bz1, fit.bz2, fit.bz3, fit.bz4,
                      fit.bz5, fit.bz6, fit.bz7, fit.bz8,
                      fit.bz9, fit.bz10, fit.bz11, fit.bz12, 
                      fit.bz13, fit.bz14)
  
  
  # Compare multiple models
  fit.bz.compare <- loo::loo_compare(fit.bz1$loo, fit.bz2$loo, 
                                     fit.bz3$loo, fit.bz4$loo,
                                     fit.bz5$loo, fit.bz6$loo, 
                                     fit.bz7$loo, fit.bz8$loo,
                                     fit.bz9$loo, fit.bz10$loo, 
                                     fit.bz11$loo, fit.bz12$loo,
                                     fit.bz13$loo, fit.bz14$loo)
  print(fit.bz.compare)
  fit.bz.index <- as.numeric(gsub("model", "", row.names(fit.bz.compare)[1]))
  fit.bz <- fit.bz.list[[fit.bz.index]]
  fit.bz.summary <- summary(fit.bz)
  
  # Convert year back to the real scale:
  cp.bz <- fit.bz.summary[fit.bz.summary$name == "cp_1", c("mean", "lower", "upper")]
  cp.bz.year <- as.numeric(cp.bz + min(d$year) - 1)
  names(cp.bz.year) <- c("mean", "lower", "upper")

  
  
  ### Write the key results out:
  if(output.txt){
    sink(paste0(write.dir, "Change-point analysis-", ts.var, ".txt"))
    
    cat("Location of interest: ")
    cat(LOC)
    cat("\nDisease of interest: ")
    cat(DZZ)
    cat("\nVariable of interest: ")
    cat(ts.var)
    cat("\nTime series: full\n\n\n")
    
    summary(fit.bz)
    cat("\n\n")
    print(fit.bz.compare)
    cat("\n\n")
    print(cp.bz.year)
    
    sink()
  }
    
    
  ### Time series:
  if(output.ts){
    fit.bz.ts <- plot(fit.bz, q_fit = c(0.025, 0.975)) + 
      scale_x_continuous(name = "year") + 
      ylab(ts.var) + 
      theme_classic() + 
      theme(panel.grid.major = element_line(color = "gray95"),
            axis.title = element_text(size = 10),
            axis.text = element_text(size = 10))
    ggsave(filename = paste0(write.dir, "Change-point analysis-time series-", ts.var, ".png"), 
           plot = fit.bz.ts, width = 6, height = 4, units = "in", dpi = 300)
  }
  
  ### trace plot and change-point posterior
  if(output.coeff){
    if(fit.bz.index <= 12){
      fit.bz.cp.plot <- plot_pars(fit.bz, regex_pars = "cp_", type = "combo") + 
        theme_classic() + 
        theme(panel.grid.major = element_line(color = "gray95"),
              axis.title = element_text(size = 10),
              axis.text = element_text(size = 10))
      ggsave(filename = paste0(write.dir, "Change-point analysis-coefficient-", ts.var, ".png"), 
             plot = fit.bz.cp.plot, width = 8, height = 4, units = "in", dpi = 300)
    }
  }
  
  return(fit.bz)
}
