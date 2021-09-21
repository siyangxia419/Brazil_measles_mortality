### Function to plot time series of measles and non-measles disease mortality
time_series <- function(d.ts, 
                        meV_var = "meV_mort_rate", nonmeV_var = "dzz_rate", 
                        cp = NULL, 
                        data_type = c("rate", "rate"), data_manipulation = "original", 
                        year_breaks = NULL, year_range = NULL, 
                        dzz = "non-measles", 
                        col_choice = c("#D55E00", "#0072B2"), 
                        line_choice = c("dashed", "solid")){
  # d.ts: time series data frame with each row contains mortality data from one year
  # meV_var: name of the column that container measles mortality data
  # nonmeV_var: name of the column that container measles mortality data
  # cp: a vector of change-point analysis results: mean, 95% CI lower, 95% CI higher
  # data_type: whether the measles and nonmeasles data are mortality rate or count, respectively. 
  #            Note this choice has to be consistent with meV_var and nonmeV_var
  # data_manipulation: "original", "detrend", "difference"
  # year_breaks and year_range: breaks and ranges for the year legend
  # dzz: name of the non-measles disease category
  # col_choice: colors used for non-measles and measles data, respectively
  # line_choice: line types used for non-measles and measles data, respectively
  
  ### scale the non-measles mortality:
  raw.scale <- max(d.ts[, nonmeV_var]) / max(d.ts[, meV_var])
  if(raw.scale >= 1){
    sec.axis.scale <- floor(raw.scale)
  }else{
    sec.axis.scale <- 1/ceiling(1/raw.scale)
  }
  d.ts$nonmeV_scaled <- d.ts[, nonmeV_var] / sec.axis.scale
  
  ### Some plotting preparation
  # x axis breaks
  if(is.null(year_breaks)) year_breaks <- seq(min(d.ts$year), max(d.ts$year), 5)
  
  # y axis breaks
  if(raw.scale >= 1){
    meV_breaks <- pretty(c(d.ts[, meV_var], d.ts[, "nonmeV_scaled"]))
    nonmeV_breaks <- meV_breaks * sec.axis.scale
  }else{
    nonmeV_breaks <- pretty(c(d.ts[, nonmeV_var], d.ts[, meV_var] * sec.axis.scale))
    meV_breaks <- nonmeV_breaks / sec.axis.scale
  }
  
  # axis titles
  axis_title <- c("measles mortality", paste0(dzz, " mortality"))
  if(data_manipulation == "detrend") axis_title <- paste("detrended", axis_title)
  if(data_manipulation == "difference") axis_title <- paste(axis_title, "difference")
  if(data_type[1] == "rate") axis_title[1] <- paste0(axis_title[1], "/100,000 ppl")
  if(data_type[2] == "rate") axis_title[2] <- paste0(axis_title[2], "/100,000 ppl")
  
  # Convert the data to the long form to feed ggplot
  d.ts.plot <- d.ts %>%
    dplyr::select(year, all_of(meV_var), nonmeV_scaled) %>%
    reshape2::melt(id.var = "year") %>%
    dplyr::mutate(variable = ordered(variable,
                                     levels = c("nonmeV_scaled", meV_var)))
  
  
  ### Generating the plot:
  ts_plot <- ggplot(data = d.ts.plot,
                    aes(x = year, y = value, linetype = variable, color = variable))
  
  # Add the change-point analysis results:
  if(!is.null(cp)){
    ts_plot <- ts_plot + 
      geom_vline(xintercept = cp["mean"], size = 1, color = "darkgray") + 
      annotate(geom = "rect",
               xmin = cp["lower"], xmax = cp["upper"],
               ymin = -Inf, ymax = Inf, 
               alpha = 0.3, fill = "gray")
  }
  
  # Plot the time series
  ts_plot <- ts_plot +
    geom_line(size = 0.8) +
    geom_point(size = 2, shape = 19) + 
    scale_x_continuous(name = "year", breaks = year_breaks) + 
    scale_y_continuous(name = axis_title[1], breaks = meV_breaks,
                       sec.axis = sec_axis(~ . * sec.axis.scale,
                                           name = axis_title[2], 
                                           breaks = nonmeV_breaks, 
                                           labels = round(nonmeV_breaks, digits = 2))) +
    scale_linetype_manual(values = line_choice[1:2],
                          labels = c(dzz, "measles"),
                          name = "") +
    scale_color_manual(values = col_choice[1:2],
                       labels = c(dzz, "measles"),
                       name = "") +
    theme_bw() +
    theme(axis.title = element_text(size = 12),
          axis.text = element_text(size = 12),
          axis.line = element_line(colour = "black"),
          axis.title.y.right = element_text( angle = 90), 
          legend.text = element_text(size = 12),
          legend.key.height = unit(0.8, "in"),
          legend.key.width = unit(0.4, "in"),
          panel.border = element_blank(),
          panel.grid.minor = element_blank())
  
  # print(ts_plot)
  return(ts_plot)
}



### Function to plot correlation between measles and non-measles disease mortality
correlation_plot <- function(d.ts, 
                             meV_var = "meV_mort_rate", nonmeV_var = "dzz_rate", 
                             data_type = c("rate", "rate"), data_manipulation = "original",
                             year_breaks = NULL, year_range = NULL, 
                             dzz = "non-measles", 
                             col_choice = c("#E69F00", "#0072B2")){
  
  # d.ts: time series data frame with each row contains mortality data from one year
  # meV_var: name of the column that container measles mortality data
  # nonmeV_var: name of the column that container measles mortality data
  # data_type: whether the measles and nonmeasles data are mortality rate or count, respectively. 
  #            Note this choice has to be consistent with meV_var and nonmeV_var
  # data_manipulation: "original", "detrend", "difference"
  # year_breaks and year_range: breaks and ranges for the year legend
  # dzz: name of the non-measles disease category
  # col_choice: the two end of the color spectral that reflects year
  
  ### Linear regression
  d.lm <- d.ts %>% dplyr::select(meV = all_of(meV_var), nonmeV = all_of(nonmeV_var), year)
  lm.cor <- lm(nonmeV ~ meV, data = d.lm)
  
  # R-squared and the coefficient:
  R2 <- round(summary(lm.cor)$r.squared, digits = 3)
  coeff <- round(lm.cor$coefficients[2], digits = 3)
  
  # Model fit and confidence interval:
  temp <- pretty(range(d.lm$meV), n = 100)
  ci <- predict(lm.cor, newdata = data.frame(meV = temp), interval = "confidence", level = 0.95)
  ci <- as.data.frame(ci)
  ci$meV <- temp
  
  
  ### Plotting preparation
  axis_title <- c("measles mortality", paste0(dzz, " mortality"))
  if(data_manipulation == "detrend") axis_title <- paste("detrended", axis_title)
  if(data_manipulation == "difference") axis_title <- paste(axis_title, "difference")
  if(data_type[1] == "rate") axis_title[1] <- paste0(axis_title[1], "/100,000 ppl")
  if(data_type[2] == "rate") axis_title[2] <- paste0(axis_title[2], "/100,000 ppl")
  
  if(is.null(year_breaks)) year_breaks <- pretty(d.lm$year)
  if(is.null(year_range)) year_range <- range(d.lm$year)
  
  ### Make the figure
  cr_plot <- ggplot() + 
    geom_ribbon(data = ci, aes(x = meV, ymin = lwr, ymax = upr), fill = "gray90") + 
    geom_line(data = ci, aes(x = meV, y = fit), size = 1) + 
    geom_path(data = d.lm, aes(x = meV, y = nonmeV),
              size = 0.5, linetype = "dashed", color = "gray") +
    geom_point(data = d.lm, aes(x = meV, y = nonmeV, color = year),
               size = 2, shape = 19) +
    xlab(label = axis_title[1]) + 
    ylab(label = axis_title[2]) + 
    scale_color_gradient(name = "Year", 
                         breaks = year_breaks, limits = year_range,
                         low = col_choice[1], high = col_choice[2]) + 
    annotate(geom = "text", 
             x = position_in_range(x = d.lm$meV, p = 1),
             y = position_in_range(x = c(d.lm$nonmeV, ci$lwr, ci$upr), p = 0.05), 
             label = bquote(paste(R^2 ,"= ", .(R2))),
             size = 4, hjust = 1) + 
    annotate(geom = "text", 
             x = position_in_range(x = d.lm$meV, p = 1),
             y = position_in_range(x = c(d.lm$nonmeV, ci$lwr, ci$upr), p = 0.12), 
             label = paste0("coeff = ", coeff),
             size = 4, hjust = 1) + 
    theme_classic() + 
    theme(axis.title = element_text(size = 12),
          axis.text = element_text(size = 12),
          legend.text = element_text(size = 10))
  
  # print(cr_plot)
  return(cr_plot)
}
