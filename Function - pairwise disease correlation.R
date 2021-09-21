disease_pairwise_correlation <- function(dzz, dzz.order, 
                                         data.type = "mortality/100,000 ppl",
                                         time.period = "full time series",
                                         data.manipulation = "original data",
                                         include.year = TRUE,
                                         color.choice = c("#0072B2", "#FFFFFF", "#D55E00")){
  # dzz: a matrix of disease count or rate with each row represents a year and each column is a disease.
  #      The first column is year
  # dzz.order: order of all disease categories appearing on the axes
  
  
  # 1. Time series ----------------------------------------------------------
  
  # Convert dzz into a long form:
  dzz.long <- dzz %>%
    reshape2::melt(id.var = "year") %>%
    dplyr::mutate(color = ifelse(variable == "measles", "meV", "dzz"))
  
  
  # Visualization of the time series
  ts.dzz <- ggplot(data = dzz.long, aes(x = year, y = value, color = color)) + 
    geom_line(size = 1) +
    facet_wrap(vars(variable), scales = "free_y") + 
    scale_y_continuous(name = data.type) +
    scale_color_manual(values = c("black", "red")) + 
    xlab("year") +
    theme_bw() +
    theme(axis.title = element_text(size = 12),
          axis.text = element_text(size = 12),
          axis.line = element_line(colour = "black"),
          axis.title.y.right = element_text( angle = 90),
          panel.border = element_blank(), 
          strip.background = element_blank(),
          strip.text = element_text(size = 12),
          legend.position = "none")  
  # print(ts.dzz)
  
  
  # 2. Heatmap of pairwise correlation coefficients -------------------------
  
  # If include.year is false
  if(!include.year){
    dzz$year <- NULL
    dzz.order <- dzz.order[dzz.order != "year"]
  }
  
  # Calculate the correlation between disease categories:
  dzz <- data.matrix(dzz)
  dzzcor <- cor(dzz)
  
  # Sort the diseases by their correlation with measles (minimum -> maximum)
  dzzcor.sorted <- dzzcor[dzz.order, dzz.order]
  
  # Convert into the long form:
  dzzcor.sorted[lower.tri(dzzcor.sorted, diag = TRUE)] <- NA
  dzzcor.long <- melt(dzzcor.sorted)
  dzzcor.long <- dzzcor.long[!is.na(dzzcor.long$value), ]
  dzzcor.long <- dzzcor.long[order(dzzcor.long$value), ]
  dzzcor.long$value2 <- round(dzzcor.long$value, digits = 2)
  
  # Order of the diseases
  dzzcor.long$Var1 <- ordered(dzzcor.long$Var1, levels = dzz.order)
  dzzcor.long$Var2 <- ordered(dzzcor.long$Var2, levels = dzz.order)
  
  
  # Some visualization:
  heatmap.dzzcor <- ggplot(data = dzzcor.long) + 
    geom_tile(aes(x = Var1, y = Var2, fill = value)) + 
    geom_text(aes(x = Var1, y = Var2, label = value2), size = 2) + 
    scale_fill_gradient2(name = "Correlation", 
                         low = color.choice[1], mid = color.choice[2], high = color.choice[3],
                         breaks = seq(from = -1, to = 1, by = 0.5),
                         limits = c(-1, 1)) + 
    ggtitle(paste(data.type, time.period, data.manipulation, sep = " - ")) + 
    theme_bw() + 
    theme(panel.border = element_blank(),
          panel.grid = element_blank(),
          axis.line = element_line(), 
          axis.text.x = element_text(size = 10, angle = 60, hjust = 1),
          axis.text.y = element_text(size = 10),
          axis.title = element_blank(),
          plot.title = element_text(size = 10, hjust = 0.5),
          legend.text = element_text(size = 10))
  # print(heatmap.dzzcor)
  
  
  
  # 3. Histogram of all correlation coefficients ----------------------------
  
  dzzcor.long2 <- dzzcor.long[dzzcor.long$Var1 != "year" & dzzcor.long$Var2 != "year", ]
  hist.dzzcor <- ggplot(data = dzzcor.long2) + 
    geom_histogram(aes(x = value)) + 
    xlab("Correlation coefficient") + 
    ylab("Frequency") + 
    theme_bw() + 
    theme(panel.border = element_blank(),
          axis.line = element_line(), 
          axis.text = element_text(size = 10),
          axis.title = element_text(size = 12))
  # print(hist.dzzcor)
  
  
  
  # 4. Sample size and correlation coefficients -----------------------------
  
  # Calculate the mean mortality rate or count across all years as a proxy for sample sizes
  dzz.mean <- sort(colMeans(dzz), decreasing = TRUE)
  samplesize.dzz.order <- names(dzz.mean)
  dzzcor.n.sorted <- dzzcor[samplesize.dzz.order, samplesize.dzz.order]
  
  # Convert into the long form:
  dzzcor.n.sorted[lower.tri(dzzcor.n.sorted, diag = TRUE)] <- NA
  dzzcor.n.long <- melt(dzzcor.n.sorted)
  dzzcor.n.long <- dzzcor.n.long[!is.na(dzzcor.n.long$value), ]
  
  # Sample size and correlation coefficients:
  dzzcor.n.long$n1 <- dzz.mean[match(as.character(dzzcor.n.long$Var1), 
                                     names(dzz.mean))]
  dzzcor.n.long$n2 <- dzz.mean[match(as.character(dzzcor.n.long$Var2), 
                                     names(dzz.mean))]
  dzzcor.n.long <- dzzcor.n.long %>% filter(Var1 != "year" & Var2 != "year")
  
  dzzcor.n.range <- range(c(dzzcor.n.long$n1, dzzcor.n.long$n2))
  dzzcor.n.jitter <- diff(dzzcor.n.range) * 0.01
  
  # Visualization
  samplesize.dzzcor <- ggplot(data = dzzcor.n.long) + 
    geom_jitter(aes(x = n1, y = n2, color = value), 
                size = 2.5, shape = 16, stroke = 2, 
                width = dzzcor.n.jitter, height = dzzcor.n.jitter) + 
    scale_color_gradient2(name = "Correlation", 
                          low = color.choice[1], mid = "white", high = color.choice[3],
                          breaks = seq(from = -1, to = 1, by = 0.5),
                          limits = c(-1, 1)) + 
    ggtitle("Sample size and correlation") + 
    xlab(paste0(data.type, " of disease 1")) + 
    ylab(paste0(data.type, " of disease 2")) + 
    theme_bw() + 
    theme(panel.border = element_blank(),
          axis.line = element_line(), 
          axis.text.x = element_text(size = 10, angle = 0, hjust = 1),
          axis.text.y = element_text(size = 10),
          axis.title = element_text(size = 10),
          plot.title = element_text(size = 10, hjust = 0.5),
          legend.text = element_text(size = 10))
  # print(samplesize.dzzcor)
  
  
  # Some statistical analysis
  samplesize.lm.dzzcor <- lm(value ~ n1 + n2, data = dzzcor.n.long)
  summary(samplesize.lm.dzzcor)
  
  return(list(ts = ts.dzz, 
              heatmap = heatmap.dzzcor, 
              histogram = hist.dzzcor, 
              sample_size = samplesize.dzzcor, 
              sample_size_lm = samplesize.lm.dzzcor))
}