#############################################
#                                           #
#      COPD Plot functions                  #
#      -Mean Cumulative Count (MCC)         #
#      -Cumulative Incidence (CI)           #
#      -KM curve for death                  #
#                                           #
#############################################       
# Function to run MCC functions from Dong et. al (MCC, SCI)
MCC.run <- function(data){
  id    <- as.numeric(data[,"id"])
  time  <- as.numeric(data[,"day"])
  cause <- as.numeric(data[,"event"])
  MCC   <- MCC(id = id, time = time, cause = cause, type = "MCC")
  return(MCC)
}

# Function to plot MCC  
MCC.plot <- function(event){
  for(i in 1:2){
    Xi <- i-1
    fit.df <- matrix(NA, n.sim, cycles)
    for(k in 1:n.sim) { 
      data.list <- list.files(path = output.path, pattern = paste0(event,".", Xi,".",k,".*\\.csv"), recursive = T)
      data.sim  <- do.call("rbind",lapply(data.list, function(x) read.csv(paste(output.path, x, sep=''), stringsAsFactors = FALSE)))
      data.sim  <- subset(data.sim, sim.no == k)
      
      Xi <- as.numeric(data.sim$X[1])
      k <- data.sim$sim.no[1]
      col <- Xi + 1
      
      
      MCC.data <- data.sim %>%
        group_by(sim.no) %>%
        do(data.frame(MCC=MCC.run(data.sim)))
      MCC.data$sim.no <- k
      MCC.data$X <- Xi
      
      
      # fit spline and extract datapoints 
      m   <- gam(MCC.MCC ~ s(MCC.time), data = MCC.data)
      pd  <- data.frame(MCC.time=seq(1,cycles))
      fit <- predict.gam(m, newdata=pd)
      
      fit.df[k,] <-  fit
      print(k)
    }
    if(Xi==0){
      fit.df.0 <- fit.df
      fit.df.sort.0 <- apply(fit.df.0, 2, sort)
      fit.df.ci.0 <- apply(fit.df.sort.0, 2, quantile, probs = c(0.025, 0.975))
      fit.df.med.0 <- apply(fit.df.sort.0, 2, quantile, probs = 0.5)
      fit.df.mean.0 <- apply(fit.df.sort.0, 2, mean)
      fit.df.ci.l.0 <- fit.df.ci.0[1,]
      fit.df.ci.u.0 <- fit.df.ci.0[2,]
      print(Xi)
    }
    else{
      fit.df.1 <- fit.df
      fit.df.sort.1 <- apply(fit.df.1, 2, sort)
      fit.df.ci.1 <- apply(fit.df.sort.1, 2, quantile, probs = c(0.025, 0.975))
      fit.df.med.1 <- apply(fit.df.sort.1, 2, quantile, probs = 0.5)
      fit.df.mean.1 <- apply(fit.df.sort.1, 2, mean)
      fit.df.ci.l.1 <- fit.df.ci.1[1,]
      fit.df.ci.u.1 <- fit.df.ci.1[2,]
      print(Xi)
    }
  }
  
  # plots 
  x <- seq_along(1:ncol(fit.df))
  
  ymin <- 0
  ymax <- max(max(fit.df.ci.u.1), max(fit.df.ci.u.0)) + 0.1
  
  plot(x, fit.df.ci.l.0, type ="l", bty ="L", col = rgb(0, 1, 0, 0.1), ylim = c(ymin, ymax), xlab = "Days", ylab = "MCC")
  points(x, fit.df.ci.u.0, type ="l" , col = rgb(0, 1, 0, 0.1))
  polygon(c(x,rev(x)), c(fit.df.ci.u.0, rev(fit.df.ci.l.0)), col = rgb(0, 1, 0, 0.3), border = "white")
  lines(x, fit.df.mean.0, col = "green")
  par(new = T)
  plot(x, fit.df.ci.l.1, type ="l", bty ="L", col = rgb(0, 0, 1, 0.1), ylim = c(ymin, ymax), xlab = "Days", ylab = "MCC", main = paste0(event," MCC")) 
  points(x,fit.df.ci.u.1, type ="l", col = rgb(0, 0, 1, 0.1))
  polygon(c(x,rev(x)),c(fit.df.ci.u.1,rev(fit.df.ci.l.1)), col = rgb(0, 0, 1, 0.3), border = NA)  
  lines(x, fit.df.mean.1, col = "purple")
  legend('topleft', c("0", "1"), lty=1, col=c(rgb(0, 1, 0, 0.3), rgb(0, 0, 1, 0.3)), bty='n')
  
  par(new = F)
  fit.df.diff <- fit.df.1 - fit.df.0
  fit.df.sort <- apply(fit.df.diff, 2, sort)
  fit.df.ci   <- apply(fit.df.sort, 2, quantile, probs = c(0.025, 0.975))
  fit.df.med  <- apply(fit.df.sort, 2, quantile, probs = 0.5)
  fit.df.mean <- apply(fit.df.sort, 2, mean)
  fit.df.diff.ci.l <- fit.df.ci[1,]
  fit.df.diff.ci.u <- fit.df.ci[2,]
  
  ymin <- min(fit.df.diff.ci.l) - 0.1
  ymax <- max(fit.df.diff.ci.u) + 0.1
  
  plot(x, fit.df.diff.ci.l, type="l", bty="L", col = rgb(0, 1, 1, 0.1), ylim = c(ymin, ymax), xlab = "Days", ylab = "MCC", main = paste0(event," MCC Abs"))
  points(x, fit.df.diff.ci.u, type="l" ,  col = rgb(0, 1, 1, 0.1))
  polygon(c(x,rev(x)), c(fit.df.diff.ci.u, rev(fit.df.diff.ci.l)), col = rgb(0, 1, 1, 0.3), border = "white")
  lines(x, fit.df.mean, col = "turquoise3")
  legend('topleft', "MCC diff", lty=1, col="turquoise3", bty='n')      
  
  par(new = F)
  fit.df.rel <- fit.df.1 / fit.df.0
  fit.df.sort <- apply(fit.df.rel,  2, sort)
  fit.df.ci <- apply(fit.df.sort, 2, quantile, probs = c(0.025, 0.975))
  fit.df.med <- apply(fit.df.sort, 2, quantile, probs = 0.5)
  fit.df.mean <- apply(fit.df.sort, 2, mean)
  fit.df.rel.ci.l <- fit.df.ci[1,]
  fit.df.rel.ci.u <- fit.df.ci[2,]
  
  ymin <- min(fit.df.rel.ci.l) - 0.1
  ymax <- max(fit.df.rel.ci.u) + 0.1
  
  plot(x, fit.df.rel.ci.l, type="l", bty="L", col = rgb(1, 0, 0, 0.1), ylim = c(ymin, ymax), xlab = "Days", ylab = "MCC", main = paste0(event," MCC Rel"))
  points(x, fit.df.rel.ci.u, type="l" ,  col = rgb(1, 0, 0, 0.1))
  polygon(c(x,rev(x)), c(fit.df.rel.ci.u, rev(fit.df.rel.ci.l)), col = rgb(1, 0, 0, 0.3), border = "white")
  lines(x, fit.df.mean, col = "pink3")
  legend('topleft', "MCC Ratio", lty=1, col="pink3", bty='n')      
  
  
}

########################################### 
#                                         #
#   START OF:                             #
#   -Cumulative Incidence and Survival    #
#   -Plots                                #
#                                         #
###########################################        
KM.plot <- function(event){
  for(i in 1:2){
    Xi <- i-1
    surv.df <- matrix(NA, n.sim, cycles)
    for(k in 1:n.sim) { 
      data.list <- list.files(path = output.path, pattern = paste0(event,".", Xi,".",k,".*\\.csv"), recursive = T)
      data.sim <- do.call("rbind",lapply(data.list, function(x) read.csv(paste(output.path, x, sep=''), stringsAsFactors = FALSE)))
      data.sim <- subset(data.sim, sim.no == k)
      
      Xi <- data.sim$X[1]
      k <- data.sim$sim.no[1]
      col <- Xi + 1
      
      event.data <- subset(data.sim, event != 2) 
      km <- survfit(Surv(day, event == 1) ~ 1, data = event.data)
      tidy_km <- tidy(km)
      
      # fit spline and extract datapoints 
      m   <- gam(estimate ~ s(time), data = tidy_km)
      pd  <- data.frame(time = seq(1,cycles))
      fit <- predict.gam(m, newdata = pd)
      
      surv.df[k,] <- fit
    }
    
    if(Xi==0){
      surv.df.0 <- surv.df
      fit.df.sort.0 <- apply(surv.df.0, 2, sort)
      fit.df.ci.0 <- apply(fit.df.sort.0, 2, quantile, probs = c(0.025, 0.975))
      fit.df.med.0 <- apply(fit.df.sort.0, 2, quantile, probs = 0.5)
      fit.df.mean.0 <- apply(fit.df.sort.0, 2, mean)
      fit.df.ci.l.0 <- fit.df.ci.0[1,]
      fit.df.ci.u.0 <- fit.df.ci.0[2,]
    }
    else{
      surv.df.1 <- surv.df
      fit.df.sort.1 <- apply(surv.df.1, 2, sort)
      fit.df.ci.1 <- apply(fit.df.sort.1, 2, quantile, probs = c(0.025, 0.975))
      fit.df.med.1 <- apply(fit.df.sort.1, 2, quantile, probs = 0.5)
      fit.df.mean.1 <- apply(fit.df.sort.1, 2, mean)
      fit.df.ci.l.1 <- fit.df.ci.1[1,]
      fit.df.ci.u.1 <- fit.df.ci.1[2,]
    }
  }
  # plots 
  x <- seq_along(1:ncol(surv.df))
  
  ymin <- 0
  ymax <- max(max(fit.df.ci.u.1), max(fit.df.ci.u.0)) + 0.1
  
  plot(x, fit.df.ci.l.0, type="l", bty="L", col = rgb(0, 1, 0, 0.1), ylim = c(ymin, ymax), xlab = "Days", ylab = "Survival Probability")
  points(x, fit.df.ci.u.0, type="l" , col = rgb(0, 1, 0, 0.1))
  polygon(c(x,rev(x)), c(fit.df.ci.u.0, rev(fit.df.ci.l.0)), col = rgb(0, 1, 0, 0.3), border = "white")
  lines(x, fit.df.mean.0, col = "green")
  par(new = T)
  plot(x, fit.df.ci.l.1, type="l", bty="L", col = rgb(0, 0, 1, 0.1), ylim = c(ymin, ymax), xlab = "Days", ylab = "Survival Probability", main = paste0(event," KM"))
  points(x,fit.df.ci.u.1, type="l", col = rgb(0, 0, 1, 0.1))
  polygon(c(x,rev(x)),c(fit.df.ci.u.1,rev(fit.df.ci.l.1)), col = rgb(0, 0, 1, 0.3), border = NA)  
  lines(x, fit.df.mean.1, col = "purple")
  legend('bottomleft', c("0", "1"), lty=1, col=c(rgb(0, 1, 0, 0.3), rgb(0, 0, 1, 0.3)), bty='n')  
  
  
  par(new = F)
  surv.df.rel <- surv.df.1/surv.df.0
  fit.df.rel <- apply(surv.df.rel, 2, sort)
  fit.df.rel.ci <- apply(fit.df.rel,  2, quantile, probs = c(0.025, 0.975))
  fit.df.rel.med <- apply(fit.df.rel,  2, quantile, probs = 0.5)
  fit.df.rel.mean <- apply(fit.df.rel,  2, mean)
  fit.df.rel.ci.l <- fit.df.rel.ci[1,]
  fit.df.rel.ci.u <- fit.df.rel.ci[2,]
  
  ymin <- 0
  ymax <- max(fit.df.rel.ci.u) + 0.1
  
  plot(x, fit.df.rel.ci.l, type="l", bty="L", col = rgb(0, 1, 0, 0.1), ylim = c(ymin, ymax), xlab = "Days", ylab = "Survival Probability", main = paste0(event," KM Rel"))
  points(x, fit.df.rel.ci.u, type="l" , col = rgb(0, 1, 0, 0.1))
  polygon(c(x,rev(x)), c(fit.df.rel.ci.u, rev(fit.df.rel.ci.l)), col = rgb(0, 1, 0, 0.3), border = "white")
  lines(x, fit.df.rel.mean, col = "green")
  legend('bottomleft', "ratio", lty=1, col=c(rgb(0, 1, 0, 0.3), rgb(0, 1, 0, 0.3)), bty='n')   
  
  
  par(new = F)
  surv.df.abs <- surv.df.1 - surv.df.0
  fit.df.abs <- apply(surv.df.abs, 2, sort)
  fit.df.abs.ci <- apply(fit.df.abs,  2, quantile, probs = c(0.025, 0.975))
  fit.df.abs.med <- apply(fit.df.abs,  2, quantile, probs = 0.5)
  fit.df.abs.mean <- apply(fit.df.abs,  2, mean)
  fit.df.abs.ci.l <- fit.df.abs.ci[1,]
  fit.df.abs.ci.u <- fit.df.abs.ci[2,]
  
  ymin <- min(fit.df.abs.ci.l) - 0.1
  ymax <- max(fit.df.abs.ci.u) + 0.1
  
  plot(x, fit.df.abs.ci.l, type="l", bty="L", col = rgb(0, 0, 1, 0.1), ylim = c(ymin, ymax), xlab = "Days", ylab = "Survival Probability", main = paste0(event," KM Abs"))
  points(x, fit.df.abs.ci.u, type="l" , col = rgb(0, 0, 1, 0.1))
  polygon(c(x,rev(x)), c(fit.df.abs.ci.u, rev(fit.df.abs.ci.l)), col = rgb(0, 0, 1, 0.3), border = "white")
  lines(x, fit.df.abs.mean, col = "pink3")
  legend('bottomleft', "diff", lty=1, col=c(rgb(0, 0, 1, 0.3), rgb(0, 1, 0, 0.3)), bty='n')  
  
} 


CI.plot <- function(event){
  for(i in 1:2){
    Xi <- i-1
    surv.df <- matrix(NA, n.sim, cycles)
    for(k in 1:n.sim) { 
      data.list <- list.files(path = output.path, pattern = paste0(event,".", Xi,".",k,".*\\.csv"), recursive = T)
      data.sim <- do.call("rbind",lapply(data.list, function(x) read.csv(paste(output.path, x, sep=''), stringsAsFactors = FALSE)))
      data.sim <- subset(data.sim, sim.no == k)
      
      Xi <- data.sim$X[1]
      k <- data.sim$sim.no[1]
      col <- Xi + 1
      
      event.data <- subset(data.sim, event != 2) 
      km <- survfit(Surv(day, event==1) ~ 1, data = event.data)
      tidy_km <- tidy(km)
      tidy_km$ci <- 1-tidy_km$estimate
      
      # fit spline and extract datapoints 
      m   <- gam(ci ~ s(time), data = tidy_km)
      pd  <- data.frame(time = seq(1, cycles))
      fit <- predict.gam(m, newdata = pd)
      
      surv.df[k,] <- fit
    }     
    if(Xi==0){
      surv.df.0 <- surv.df
      fit.df.sort.0 <- apply(surv.df.0, 2, sort)
      fit.df.ci.0 <- apply(fit.df.sort.0, 2, quantile, probs = c(0.025, 0.975))
      fit.df.med.0 <- apply(fit.df.sort.0, 2, quantile, probs = 0.5)
      fit.df.mean.0 <- apply(fit.df.sort.0, 2, mean)
      fit.df.ci.l.0 <- fit.df.ci.0[1,]
      fit.df.ci.u.0 <- fit.df.ci.0[2,]
    }
    else{
      surv.df.1 <- surv.df
      fit.df.sort.1 <- apply(surv.df.1, 2, sort)
      fit.df.ci.1 <- apply(fit.df.sort.1, 2, quantile, probs = c(0.025, 0.975))
      fit.df.med.1 <- apply(fit.df.sort.1, 2, quantile, probs = 0.5)
      fit.df.mean.1 <- apply(fit.df.sort.1, 2, mean)
      fit.df.ci.l.1 <- fit.df.ci.1[1,]
      fit.df.ci.u.1 <- fit.df.ci.1[2,]
      
    }
  }
  # plots
  x <- seq_along(1:ncol(surv.df))
  
  ymin <- 0
  ymax <- max(max(fit.df.ci.u.1), max(fit.df.ci.u.0)) + 0.1
  
  plot(x, fit.df.ci.l.0, type="l", bty="L", col = rgb(0, 1, 0, 0.1), ylim = c(ymin, ymax), xlab = "Days", ylab = "Cumulative Incidence")
  points(x, fit.df.ci.u.0, type="l" , col = rgb(0, 1, 0, 0.1))
  polygon(c(x,rev(x)), c(fit.df.ci.u.0, rev(fit.df.ci.l.0)), col = rgb(0, 1, 0, 0.3), border = "white")
  lines(x, fit.df.mean.0, col = "green")
  par(new = T)
  plot(x, fit.df.ci.l.1, type="l", bty="L", col = rgb(0, 0, 1, 0.1), ylim = c(ymin, ymax), xlab = "Days", ylab = "Cumulative Incidence", main = paste0(event," CI"))
  points(x,fit.df.ci.u.1, type="l", col = rgb(0, 0, 1, 0.1))
  polygon(c(x,rev(x)),c(fit.df.ci.u.1,rev(fit.df.ci.l.1)), col = rgb(0, 0, 1, 0.3), border = NA)  
  lines(x, fit.df.mean.1, col = "purple")
  legend('topleft', c("0", "1"), lty=1, col=c(rgb(0, 1, 0, 0.3), rgb(0, 0, 1, 0.3)), bty='n')  
  
  
  par(new = F)
  surv.df.rel <- surv.df.1/surv.df.0
  fit.df.rel  <- apply(surv.df.rel, 2, sort)
  fit.df.rel.ci   <- apply(fit.df.rel,  2, quantile, probs = c(0.025, 0.975))
  fit.df.rel.med  <- apply(fit.df.rel,  2, quantile, probs = 0.5)
  fit.df.rel.mean <- apply(fit.df.rel,  2, mean)
  fit.df.rel.ci.l <- fit.df.rel.ci[1,]
  fit.df.rel.ci.u <- fit.df.rel.ci[2,]
  
  ymin <- min(fit.df.rel.ci.l) - 0.1
  ymax <- max(fit.df.rel.ci.u) + 0.1
  
  plot(x, fit.df.rel.ci.l, type="l", bty="L", col = rgb(0, 1, 0, 0.1), ylim = c(0, ymax), xlab = "Days", ylab = "Cumulative Incidence", main = paste0(event," CI Rel"))
  points(x, fit.df.rel.ci.u, type="l" , col = rgb(0, 1, 0, 0.1))
  polygon(c(x,rev(x)), c(fit.df.rel.ci.u, rev(fit.df.rel.ci.l)), col = rgb(0, 1, 0, 0.3), border = "white")
  lines(x, fit.df.rel.mean, col = "green")
  legend('bottomleft', "ratio", lty=1, col=c(rgb(0, 1, 0, 0.3), rgb(0, 1, 0, 0.3)), bty='n')   
  
  
  par(new = F)
  surv.df.abs <- surv.df.1 - surv.df.0
  fit.df.abs  <- apply(surv.df.abs, 2, sort)
  fit.df.abs.ci   <- apply(fit.df.abs,  2, quantile, probs = c(0.025, 0.975))
  fit.df.abs.med  <- apply(fit.df.abs,  2, quantile, probs = 0.5)
  fit.df.abs.mean <- apply(fit.df.abs,  2, mean)
  fit.df.abs.ci.l <- fit.df.abs.ci[1,]
  fit.df.abs.ci.u <- fit.df.abs.ci[2,]
  
  ymin <- min(fit.df.abs.ci.l) - 0.1
  ymax <- max(fit.df.abs.ci.u) + 0.1
  
  plot(x, fit.df.abs.ci.l, type ="l", bty="L", col = rgb(0, 0, 1, 0.1), ylim = c(ymin, ymax), xlab = "Days", ylab = "Cumulative Incidence", main = paste0(event," CI Abs"))
  points(x, fit.df.abs.ci.u, type ="l" , col = rgb(0, 0, 1, 0.1))
  polygon(c(x,rev(x)), c(fit.df.abs.ci.u, rev(fit.df.abs.ci.l)), col = rgb(0, 0, 1, 0.3), border = "white")
  lines(x, fit.df.abs.mean, col = "pink3")
  legend('bottomleft', "diff", lty = 1, col = c(rgb(0, 0, 1, 0.3), rgb(0, 1, 0, 0.3)), bty = 'n')  
  
}