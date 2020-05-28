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
MCC.plot <- function(event, ymax_all,  ymax_diff, ymin_diff, ymax_rel, ymin_rel){
  for(i in 1:1){
    
    Xi <- i-1
    fit.df.00 <- fit.df.01 <- fit.df.10 <- matrix(NA, n.sim, cycles)

    
    for(k in 1:n.sim) { 
      data.list <- list.files(path = output.path, pattern = paste0(event,".", Xi,".", k,".*\\.csv"), recursive = T)
      data.sim  <- do.call("rbind",lapply(data.list, function(x) read.csv(paste(output.path, x, sep=''), stringsAsFactors = FALSE)))
      data.sim  <- subset(data.sim, sim.no == k)

      k <- data.sim$sim.no[1]
      
      data.sim.00 <- subset(data.sim, X == 0)
      data.sim.01 <- subset(data.sim, X == 1)
      data.sim.10 <- subset(data.sim, X == 10)
      
      MCC.data.00 <- data.sim.00 %>%
                     group_by(sim.no, X) %>%
                     do(data.frame(MCC = MCC.run(data.sim.00)))
      
      MCC.data.01 <- data.sim.01 %>%
                     group_by(sim.no, X) %>%
                     do(data.frame(MCC = MCC.run(data.sim.01)))
      
      MCC.data.10 <- data.sim.10 %>%
                     group_by(sim.no, X) %>%
                     do(data.frame(MCC = MCC.run(data.sim.10)))

      # fit spline and extract datapoints 
        m.00   <- gam(MCC.MCC ~ s(MCC.time), data = MCC.data.00)
        pd.00  <- data.frame(MCC.time = seq(1, cycles))
        fit.00 <- predict.gam(m.00, newdata = pd.00)
        fit.df.00[k,] <-  fit.00
        
        m.01   <- gam(MCC.MCC ~ s(MCC.time), data = MCC.data.01)
        pd.01  <- data.frame(MCC.time = seq(1, cycles))
        fit.01 <- predict.gam(m.01, newdata = pd.01)
        fit.df.01[k, ] <-  fit.01
        
        m.10   <- gam(MCC.MCC ~ s(MCC.time), data = MCC.data.10)
        pd.10  <- data.frame(MCC.time = seq(1, cycles))
        fit.10 <- predict.gam(m.10, newdata = pd.10)
        fit.df.10[k, ] <-  fit.10

    }
        fit.df.sort.00 <- apply(fit.df.00, 2, sort)
        fit.df.ci.00   <- apply(fit.df.sort.00, 2, quantile, probs = c(0.025, 0.975))
        fit.df.med.00  <- apply(fit.df.sort.00, 2, quantile, probs = 0.5)
        fit.df.mean.00 <- apply(fit.df.sort.00, 2, mean)
        fit.df.ci.l.00 <- fit.df.ci.00[1,]
        fit.df.ci.u.00 <- fit.df.ci.00[2,]
        
        fit.df.sort.01 <- apply(fit.df.01, 2, sort)
        fit.df.ci.01   <- apply(fit.df.sort.01, 2, quantile, probs = c(0.025, 0.975))
        fit.df.med.01  <- apply(fit.df.sort.01, 2, quantile, probs = 0.5)
        fit.df.mean.01 <- apply(fit.df.sort.01, 2, mean)
        fit.df.ci.l.01 <- fit.df.ci.01[1,]
        fit.df.ci.u.01 <- fit.df.ci.01[2,]
        
        fit.df.sort.10 <- apply(fit.df.10, 2, sort)
        fit.df.ci.10   <- apply(fit.df.sort.10, 2, quantile, probs = c(0.025, 0.975))
        fit.df.med.10  <- apply(fit.df.sort.10, 2, quantile, probs = 0.5)
        fit.df.mean.10 <- apply(fit.df.sort.10, 2, mean)
        fit.df.ci.l.10 <- fit.df.ci.10[1,]
        fit.df.ci.u.10 <- fit.df.ci.10[2,]
        

  
  # plots 
    x <- seq_along(1:ncol(fit.df.00))
    
    ymin <- 0
    ymax <- ymax_all
    
    plot  (x, fit.df.ci.l.00, type ="l", bty ="L", col = rgb(0, 1, 0, 0.1), ylim = c(ymin, ymax), xlab = "Days", ylab = "MCC")
    points(x, fit.df.ci.u.00, type ="l" , col = rgb(0, 1, 0, 0.1))
    polygon(c(x,rev(x)), c(fit.df.ci.u.00, rev(fit.df.ci.l.00)), col = rgb(0, 1, 0, 0.3), border = NA)
    lines(x, fit.df.mean.00, col = "green")
    par(new = T)
    plot  (x, fit.df.ci.l.01, type ="l", bty ="L", col = rgb(0, 0, 1, 0.1), ylim = c(ymin, ymax), xlab = "Days", ylab = "MCC", main = paste0(event," MCC")) 
    points(x, fit.df.ci.u.01, type ="l", col = rgb(0, 0, 1, 0.1))
    polygon(c(x,rev(x)),c(fit.df.ci.u.01, rev(fit.df.ci.l.01)), col = rgb(0, 0, 1, 0.3), border = NA)  
    lines(x, fit.df.mean.01, col = "purple")
    par(new = T)
    plot  (x, fit.df.ci.l.10, type ="l", bty ="L", col = rgb(0, 1, 1, 0.1), ylim = c(ymin, ymax), xlab = "Days", ylab = "MCC", main = paste0(event," MCC")) 
    points(x, fit.df.ci.u.10, type ="l", col = rgb(0, 1, 1, 0.1))
    polygon(c(x,rev(x)),c(fit.df.ci.u.10, rev(fit.df.ci.l.10)), col = rgb(0, 1, 1, 0.3), border = NA)  
    lines(x, fit.df.mean.10, col = "turquoise")
    legend('topleft', c("never", "former", "current"), lty=1, col=c("green", "purple", "turquoise"), bty='n')
    
    par(new = F)
    fit.df.diff.10.00 <- fit.df.10 - fit.df.00
    fit.df.sort.10.00 <- apply(fit.df.diff.10.00, 2, sort)
    fit.df.ci.10.00   <- apply(fit.df.sort.10.00, 2, quantile, probs = c(0.025, 0.975))
    fit.df.med.10.00  <- apply(fit.df.sort.10.00, 2, quantile, probs = 0.5)
    fit.df.mean.10.00 <- apply(fit.df.sort.10.00, 2, mean)
    fit.df.diff.ci.l.10.00 <- fit.df.ci.10.00[1,]
    fit.df.diff.ci.u.10.00 <- fit.df.ci.10.00[2,]
    
    fit.df.diff.01.00 <- fit.df.01 - fit.df.00
    fit.df.sort.01.00 <- apply(fit.df.diff.01.00, 2, sort)
    fit.df.ci.01.00   <- apply(fit.df.sort.01.00, 2, quantile, probs = c(0.025, 0.975))
    fit.df.med.01.00  <- apply(fit.df.sort.01.00, 2, quantile, probs = 0.5)
    fit.df.mean.01.00 <- apply(fit.df.sort.01.00, 2, mean)
    fit.df.diff.ci.l.01.00 <- fit.df.ci.01.00[1,]
    fit.df.diff.ci.u.01.00 <- fit.df.ci.01.00[2,]
    
    ymin <- ymin_diff
    ymax <- ymax_diff 
    
    plot  (x, fit.df.diff.ci.l.10.00, type="l", bty="L", col = rgb(0, 1, 1, 0.1), ylim = c(ymin, ymax), xlab = "Days", ylab = "MCC", main = paste0(event," MCC Abs"))
    points(x, fit.df.diff.ci.u.10.00, type="l" ,  col = rgb(0, 1, 1, 0.1))
    polygon(c(x,rev(x)), c(fit.df.diff.ci.u.10.00, rev(fit.df.diff.ci.l.10.00)), col = rgb(0, 1, 1, 0.3), border = NA)
    lines(x, fit.df.mean.10.00, col = "turquoise1")
    par(new = T)
    plot  (x, fit.df.diff.ci.l.01.00, type="l", bty="L", col = rgb(1, 0, 1, 0.1), ylim = c(ymin, ymax), xlab = "Days", ylab = "MCC", main = paste0(event," MCC Abs"))
    points(x, fit.df.diff.ci.u.01.00, type="l" ,  col = rgb(1, 0, 1, 0.1))
    polygon(c(x,rev(x)), c(fit.df.diff.ci.u.01.00, rev(fit.df.diff.ci.l.01.00)), col = rgb(1, 0, 1, 0.3), border = NA)
    lines(x, fit.df.mean.01.00, col = "purple")
    legend('topleft', c("current - never", "former - never"), lty=1, col= c("turquoise1", "purple"), bty='n') 
    
    
    
    par(new = F)
         fit.df.rel.10.00 <- fit.df.10 / fit.df.00
        fit.df.sort.10.00 <- apply(fit.df.rel.10.00,  2, sort)
          fit.df.ci.10.00 <- apply(fit.df.sort.10.00, 2, quantile, probs = c(0.025, 0.975))
         fit.df.med.10.00 <- apply(fit.df.sort.10.00, 2, quantile, probs = 0.5)
        fit.df.mean.10.00 <- apply(fit.df.sort.10.00, 2, mean)
    fit.df.rel.ci.l.10.00 <- fit.df.ci.10.00[1,]
    fit.df.rel.ci.u.10.00 <- fit.df.ci.10.00[2,]
    
         fit.df.rel.01.00 <- fit.df.01 / fit.df.00
        fit.df.sort.01.00 <- apply(fit.df.rel.01.00,  2, sort)
          fit.df.ci.01.00 <- apply(fit.df.sort.01.00, 2, quantile, probs = c(0.025, 0.975))
         fit.df.med.01.00 <- apply(fit.df.sort.01.00, 2, quantile, probs = 0.5)
        fit.df.mean.01.00 <- apply(fit.df.sort.01.00, 2, mean)
    fit.df.rel.ci.l.01.00 <- fit.df.ci.01.00[1,]
    fit.df.rel.ci.u.01.00 <- fit.df.ci.01.00[2,]
    
    
    ymin <- ymin_rel
    ymax <- ymax_rel 
    
    plot  (x, fit.df.rel.ci.l.10.00,   type="l", bty="L", col = rgb(1, 0, 0, 0.1), ylim = c(ymin, ymax), xlab = "Days", ylab = "MCC", main = paste0(event," MCC Rel"))
    points(x, fit.df.rel.ci.u.10.00, type="l" ,  col = rgb(1, 0, 0, 0.1))
    polygon(c(x,rev(x)), c(fit.df.rel.ci.u.10.00, rev(fit.df.rel.ci.l.10.00)), col = rgb(1, 0, 0, 0.3), border = NA)
    lines(x, fit.df.mean.10.00, col = "red")
    par(new = T)
    plot  (x, fit.df.rel.ci.l.01.00, type="l", bty="L", col = rgb(1, 1, 0, 0.1), ylim = c(ymin, ymax), xlab = "Days", ylab = "MCC", main = paste0(event," MCC Rel"))
    points(x, fit.df.rel.ci.u.01.00, type="l" ,  col = rgb(1, 1, 0, 0.1))
    polygon(c(x,rev(x)), c(fit.df.rel.ci.u.01.00, rev(fit.df.rel.ci.l.01.00)), col = rgb(1, 1, 0, 0.3), border = NA)
    lines(x, fit.df.mean.01.00, col = "orange")
    legend('topleft', c("current / never", "former / never"), lty=1, col=c("red", "orange"), bty='n')      
    
    print(paste0(event, 
      "5  year mean MCC 10~~ ", fit.df.mean.10[5  * 365], 
      "5  year CI L MCC 10~~ ", fit.df.ci.l.10[5  * 365],
      "5  year CI U MCC 10~~ ", fit.df.ci.u.10[5  * 365],
      "14 year mean MCC 10~~ ", fit.df.mean.10[14 * 365],
      "14 year CI L MCC 10~~ ", fit.df.ci.l.10[14 * 365],
      "14 year CI U MCC 10~~ ", fit.df.ci.u.10[14 * 365],
      "5  year mean MCC 01~~ ", fit.df.mean.01[5  * 365], 
      "5  year CI L MCC 01~~ ", fit.df.ci.l.01[5  * 365],
      "5  year CI U MCC 01~~ ", fit.df.ci.u.01[5  * 365],
      "14 year mean MCC 01~~ ", fit.df.mean.01[14 * 365],
      "14 year CI L MCC 01~~ ", fit.df.ci.l.01[14 * 365],
      "14 year CI U MCC 01~~ ", fit.df.ci.u.01[14 * 365],
      "5  year mean MCC 00~~ ", fit.df.mean.00[5  * 365], 
      "5  year CI L MCC 00~~ ", fit.df.ci.l.00[5  * 365],
      "5  year CI U MCC 00~~ ", fit.df.ci.u.00[5  * 365],
      "14 year mean MCC 00~~ ", fit.df.mean.00[14 * 365],
      "14 year CI L MCC 00~~ ", fit.df.ci.l.00[14 * 365],
      "14 year CI U MCC 00~~ ", fit.df.ci.u.00[14 * 365]))
    
  
  } 
}

########################################### 
#                                         #
#   START OF:                             #
#   -Cumulative Incidence and Survival    #
#   -Plots                                #
#                                         #
###########################################        
KM.plot <- function(event){
  for(i in 1:1){
    Xi <- i-1
    surv.df.00 <- surv.df.10 <- surv.df.01  <- matrix(NA, n.sim, cycles)

    for(k in 1:n.sim) { 
      data.list <- list.files(path = output.path, pattern = paste0(event,".", Xi,".",k,".*\\.csv"), recursive = T)
      data.sim <- do.call("rbind",lapply(data.list, function(x) read.csv(paste(output.path, x, sep=''), stringsAsFactors = FALSE)))
      data.sim <- subset(data.sim, sim.no == k)
      
      event.data <- subset(data.sim, event != 2) 
      event.data.00 <- subset(event.data, X == 0)
      event.data.01 <- subset(event.data, X == 1)
      event.data.10 <- subset(event.data, X == 10)
      
      km.00 <- survfit(Surv(day, event == 1) ~ 1, data = event.data.00)
      tidy_km.00 <- tidy(km.00)

      km.10 <- survfit(Surv(day, event == 1) ~ 1, data = event.data.10)
      tidy_km.10 <- tidy(km.10)

      km.01 <- survfit(Surv(day, event == 1) ~ 1, data = event.data.01)
      tidy_km.01 <- tidy(km.01)

      # fit spline and extract datapoints 
        m.00     <- gam(estimate ~ s(time), data = tidy_km.00)
        pd.00    <- data.frame(time = seq(1,cycles))
        fit.00   <- predict.gam(m.00, newdata = pd.00)
        
        m.10     <- gam(estimate ~ s(time), data = tidy_km.10)
        pd.10    <- data.frame(time = seq(1,cycles))
        fit.10   <- predict.gam(m.10, newdata = pd.10)
        
        m.01     <- gam(estimate ~ s(time), data = tidy_km.01)
        pd.01    <- data.frame(time = seq(1,cycles))
        fit.01   <- predict.gam(m.01, newdata = pd.01)
        
        surv.df.00[k,] <- fit.00
        surv.df.10[k,] <- fit.10
        surv.df.01[k,] <- fit.01

    }
    
      fit.df.sort.00 <- apply(surv.df.00, 2, sort)
        fit.df.ci.00 <- apply(fit.df.sort.00, 2, quantile, probs = c(0.025, 0.975))
       fit.df.med.00 <- apply(fit.df.sort.00, 2, quantile, probs = 0.5)
      fit.df.mean.00 <- apply(fit.df.sort.00, 2, mean)
      fit.df.ci.l.00 <- fit.df.ci.00[1,]
      fit.df.ci.u.00 <- fit.df.ci.00[2,]

      fit.df.sort.10 <- apply(surv.df.10, 2, sort)
        fit.df.ci.10 <- apply(fit.df.sort.10, 2, quantile, probs = c(0.025, 0.975))
       fit.df.med.10 <- apply(fit.df.sort.10, 2, quantile, probs = 0.5)
      fit.df.mean.10 <- apply(fit.df.sort.10, 2, mean)
      fit.df.ci.l.10 <- fit.df.ci.10[1,]
      fit.df.ci.u.10 <- fit.df.ci.10[2,]
      
      fit.df.sort.01 <- apply(surv.df.01, 2, sort)
        fit.df.ci.01 <- apply(fit.df.sort.01, 2, quantile, probs = c(0.025, 0.975))
       fit.df.med.01 <- apply(fit.df.sort.01, 2, quantile, probs = 0.5)
      fit.df.mean.01 <- apply(fit.df.sort.01, 2, mean)
      fit.df.ci.l.01 <- fit.df.ci.01[1,]
      fit.df.ci.u.01 <- fit.df.ci.01[2,]
      
  # KM plots 
    x <- seq_along(1:ncol(surv.df.00))
    
    ymin <- 0
    ymax <- 1
    
    plot  (x, fit.df.ci.l.00, type="l", bty="L", col = rgb(0, 1, 0, 0.1), ylim = c(ymin, ymax), xlab = "Days", ylab = "Survival Probability")
    points(x, fit.df.ci.u.00, type="l" , col = rgb(0, 1, 0, 0.1))
    polygon(c(x,rev(x)), c(fit.df.ci.u.00, rev(fit.df.ci.l.00)), col = rgb(0, 1, 0, 0.3), border = NA)
    lines(x, fit.df.mean.00, col = "green")
    par(new = T)
    plot  (x, fit.df.ci.l.01, type="l", bty="L", col = rgb(0, 0, 1, 0.1), ylim = c(ymin, ymax), xlab = "Days", ylab = "Survival Probability")
    points(x, fit.df.ci.u.01, type="l", col = rgb(0, 0, 1, 0.1))
    polygon(c(x,rev(x)),c(fit.df.ci.u.01, rev(fit.df.ci.l.01)), col = rgb(0, 0, 1, 0.3), border = NA)  
    lines(x, fit.df.mean.01, col = "purple")
    par(new = T)
    plot  (x, fit.df.ci.l.10, type="l", bty="L", col = rgb(1, 0, 0, 0.1), ylim = c(ymin, ymax), xlab = "Days", ylab = "Survival Probability")
    points(x, fit.df.ci.u.10, type="l", col = rgb(1, 0, 0, 0.1))
    polygon(c(x,rev(x)),c(fit.df.ci.u.10, rev(fit.df.ci.l.10)), col = rgb(1, 0, 0, 0.3), border = NA)  
    lines(x, fit.df.mean.10, col = "red")
    legend('bottomleft', c("never", "former", "current"), lty=1, col=c("green", "purple", "red"), bty='n')  
    
    par(new = F)
        surv.df.rel.10.00 <- surv.df.10/surv.df.00
         fit.df.rel.10.00 <- apply(surv.df.rel.10.00, 2, sort)
      fit.df.rel.ci.10.00 <- apply(fit.df.rel.10.00,  2, quantile, probs = c(0.025, 0.975))
     fit.df.rel.med.10.00 <- apply(fit.df.rel.10.00,  2, quantile, probs = 0.5)
    fit.df.rel.mean.10.00 <- apply(fit.df.rel.10.00,  2, mean)
    fit.df.rel.ci.l.10.00 <- fit.df.rel.ci.10.00[1,]
    fit.df.rel.ci.u.10.00 <- fit.df.rel.ci.10.00[2,]
    
        surv.df.rel.01.00 <- surv.df.01/surv.df.00
         fit.df.rel.01.00 <- apply(surv.df.rel.01.00, 2, sort)
      fit.df.rel.ci.01.00 <- apply(fit.df.rel.01.00,  2, quantile, probs = c(0.025, 0.975))
     fit.df.rel.med.01.00 <- apply(fit.df.rel.01.00,  2, quantile, probs = 0.5)
    fit.df.rel.mean.01.00 <- apply(fit.df.rel.01.00,  2, mean)
    fit.df.rel.ci.l.01.00 <- fit.df.rel.ci.01.00[1,]
    fit.df.rel.ci.u.01.00 <- fit.df.rel.ci.01.00[2,]
    
    ymin <- 0.5
    ymax <- 2
    
    plot  (x, fit.df.rel.ci.l.10.00 , type="l", bty="L", col = rgb(0, 1, 0, 0.1), ylim = c(ymin, ymax), xlab = "Days", ylab = "Survival Probability", main = paste0(event," KM Relative"))
    points(x, fit.df.rel.ci.u.10.00 , type="l" , col = rgb(0, 1, 0, 0.1))
    polygon(c(x,rev(x)), c(fit.df.rel.ci.u.10.00 , rev(fit.df.rel.ci.l.10.00 )), col = rgb(0, 1, 0, 0.3), border = NA)
    lines(x, fit.df.rel.mean.10.00 , col = "darkgreen")
    par(new = T)
    plot  (x, fit.df.rel.ci.l.01.00, type="l", bty="L", col = rgb(0, 1, 1, 0.1), ylim = c(ymin, ymax), xlab = "Days", ylab = "Survival Probability", main = paste0(event," KM Relative"))
    points(x, fit.df.rel.ci.u.01.00, type="l" , col = rgb(0, 1, 1, 0.1))
    polygon(c(x,rev(x)), c(fit.df.rel.ci.u.01.00, rev(fit.df.rel.ci.l.01.00)), col = rgb(0, 1, 1, 0.3), border = NA)
    lines(x, fit.df.rel.mean.01.00, col = "blue")
    legend('bottomleft', c("current/never", "former/never"), lty=1, col=c("darkgreen", "blue"), bty='n')   
    
    par(new = F)
        surv.df.abs.10.00 <- surv.df.10 - surv.df.00
         fit.df.abs.10.00 <- apply(surv.df.abs.10.00, 2, sort)
      fit.df.abs.ci.10.00 <- apply(fit.df.abs.10.00,  2, quantile, probs = c(0.025, 0.975))
     fit.df.abs.med.10.00 <- apply(fit.df.abs.10.00,  2, quantile, probs = 0.5)
    fit.df.abs.mean.10.00 <- apply(fit.df.abs.10.00,  2, mean)
    fit.df.abs.ci.l.10.00 <- fit.df.abs.ci.10.00[1,]
    fit.df.abs.ci.u.10.00 <- fit.df.abs.ci.10.00[2,]
    
        surv.df.abs.01.00 <- surv.df.01 - surv.df.00
         fit.df.abs.01.00 <- apply(surv.df.abs.01.00, 2, sort)
      fit.df.abs.ci.01.00 <- apply(fit.df.abs.01.00,  2, quantile, probs = c(0.025, 0.975))
     fit.df.abs.med.01.00 <- apply(fit.df.abs.01.00,  2, quantile, probs = 0.5)
    fit.df.abs.mean.01.00 <- apply(fit.df.abs.01.00,  2, mean)
    fit.df.abs.ci.l.01.00 <- fit.df.abs.ci.01.00[1,]
    fit.df.abs.ci.u.01.00 <- fit.df.abs.ci.01.00[2,]
    
    ymin <- -2
    ymax <- 2
    
    plot  (x, fit.df.abs.ci.l.10.00 , type="l", bty="L", col = rgb(0, 0, 1, 0.1), ylim = c(ymin, ymax), xlab = "Days", ylab = "Survival Probability", main = paste0(event," KM Absolute"))
    points(x, fit.df.abs.ci.u.10.00 , type="l" , col = rgb(0, 0, 1, 0.1))
    polygon(c(x,rev(x)), c(fit.df.abs.ci.u.10.00 , rev(fit.df.abs.ci.l.10.00)), col = rgb(0, 0, 1, 0.3), border = NA)
    lines(x, fit.df.abs.mean.10.00, col = "purple")
    par(new = T)
    plot  (x, fit.df.abs.ci.l.01.00, type="l", bty="L", col = rgb(0, 1, 1, 0.1), ylim = c(ymin, ymax), xlab = "Days", ylab = "Survival Probability", main = paste0(event," KM Absolute"))
    points(x, fit.df.abs.ci.u.01.00, type="l" , col = rgb(0, 1, 1, 0.1))
    polygon(c(x,rev(x)), c(fit.df.abs.ci.u.01.00, rev(fit.df.abs.ci.l.01.00)), col = rgb(0, 1, 1, 0.3), border = NA)
    lines(x, fit.df.abs.mean.01.00, col = "blue")
    legend('bottomleft', c("current-never", "former-never"), lty=1, col=c("purple", "blue"), bty='n')  
    
    print(paste0(
      "5  year mean KM 10~~ ", fit.df.mean.10[5  * 365],
      "5  year CI L KM 10~~ ", fit.df.ci.l.10[5  * 365],
      "5  year CI U KM 10~~ ", fit.df.ci.u.10[5  * 365],
      "14 year mean KM 10~~ ", fit.df.mean.10[14 * 365],
      "14 year CI L KM 10~~ ", fit.df.ci.l.10[14 * 365],
      "14 year CI U KM 10~~ ", fit.df.ci.u.10[14 * 365],
      "5  year mean KM 01~~ ", fit.df.mean.01[5  * 365],
      "5  year CI L KM 01~~ ", fit.df.ci.l.01[5  * 365],
      "5  year CI U KM 01~~ ", fit.df.ci.u.01[5  * 365],
      "14 year mean KM 01~~ ", fit.df.mean.01[14 * 365],
      "14 year CI L KM 01~~ ", fit.df.ci.l.01[14 * 365],
      "14 year CI U KM 01~~ ", fit.df.ci.u.01[14 * 365],
      "5  year mean KM 00~~ ", fit.df.mean.00[5  * 365],
      "5  year CI L KM 00~~ ", fit.df.ci.l.00[5  * 365],
      "5  year CI U KM 00~~ ", fit.df.ci.u.00[5  * 365],
      "14 year mean KM 00~~ ", fit.df.mean.00[14 * 365],
      "14 year CI L KM 00~~ ", fit.df.ci.l.00[14 * 365],
      "14 year CI U KM 00~~ ", fit.df.ci.u.00[14 * 365]))

  }  
} 


