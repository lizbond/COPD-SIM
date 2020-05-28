###################################################################
###                                                             ###
###    Make the baseline cohort based on real world cohort      ###
###                                                             ###
###################################################################
#library(dplyr) 
#library(openxlsx)

#setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Number of individuals in the simulated cohort
  #n.i <- 2000
  
  set.seed(9)

# Dataframe to store the results
  bl.dat <- data.frame(INT         = rep(1, n.i),                 # intercept
                       AGE         = rnorm(n.i, 63.98, 14),       # based on real world cohort
                       SEX         = rbinom(n.i, 1, 0.449),       # based on real world cohort
                       SMOKINGF    = rep(0, n.i), 
                       SMOKINGN    = rep(0, n.i),
                       CHF         = rep(0, n.i), 
                       IHD         = rep(0, n.i),       
                       CANCER      = rep(0, n.i),      
                       DM          = rep(0, n.i),        
                       ASTHMA      = rep(0, n.i),
                       DEMENTIA    = rep(0, n.i),  
                       DEPRESSION  = rep(0, n.i),  
                       ANXIETY     = rep(0, n.i),   
                       HTN         = rep(0, n.i),
                       RIO2        = rep(0, n.i), 
                       RIO3        = rep(0, n.i),             
                       DEP2        = rep(0, n.i), 
                       DEP3        = rep(0, n.i), 
                       DEP4        = rep(0, n.i), 
                       DEP5        = rep(0, n.i),
                       SMK_TIME    = rep(0, n.i))
  
# Esimates from real world data
  est.m <- read.csv(paste0(fp,"Data/ICES baseline 20200110.csv"))
  est.m[is.na(est.m)] <- 0
  
  
# Loop throuh the columns
  # SMOKING CATEGORICAL
    for(i in 5){
        odds <- data.frame(SMK1 = 1, 
                           SMK2 = exp(as.matrix(bl.dat) %*% as.matrix(est.m[i - 4])),
                           SMK3 = exp(as.matrix(bl.dat) %*% as.matrix(est.m[i - 3])))
        
        probs.m <- odds/rowSums(odds)
        
        smk.prb <- matrix(NA, nrow = n.i, ncol = 1)
        
        for(j in 1:n.i){
          smk.prb[j,] <- sample(ncol(probs.m), 1, prob = probs.m[j,])
        }
        smk.dat <- as.data.frame(smk.prb)
        smk.dat$SMOKINGF <- ifelse(smk.prb == 2, 1, 0)
        smk.dat$SMOKINGN <- ifelse(smk.prb == 3, 1, 0)
        bl.dat[,4:5] <- smk.dat[,2:3]
    }
  
  # BINARY VARAIBLES
    for(i in 6:14){
      
        odds1 <- exp(as.matrix(bl.dat) %*% as.matrix(est.m[i - 3])) 
        prob1 <- odds1 / (1 + odds1)
        bl.dat[,i] <- rbinom(n.i, 1, prob = prob1)
    
    }
  # RIO CATEGORICAL
    for(i in 16){
        odds2 <- data.frame(Q1 = 1, 
                            Q2 = exp(as.matrix(bl.dat) %*% as.matrix(est.m[i - 3])),
                            Q3 = exp(as.matrix(bl.dat) %*% as.matrix(est.m[i - 2])))
      
        probs.m2 <- odds2/rowSums(odds2)
        
        rio.dat <- matrix(NA, nrow = n.i, ncol = 1)
        
        for(j in 1:n.i){
            rio.dat[j,] <- sample(ncol(probs.m2), 1, prob = probs.m2[j,])
        }
        rio.dat <- as.data.frame(rio.dat)
        rio.dat$RIO2 <- ifelse(rio.dat[,1] == 2, 1, 0)
        rio.dat$RIO3 <- ifelse(rio.dat[,1] == 3, 1, 0)
        bl.dat[,15:16] <- rio.dat[,2:3]
    }

  # DEPRIVATION CATEGORICAL
    for(i in 17){
        odds3 <- data.frame(Q1 = 1, 
                            Q2 = exp(as.matrix(bl.dat) %*% as.matrix(est.m[i - 3])),
                            Q3 = exp(as.matrix(bl.dat) %*% as.matrix(est.m[i - 2])),
                            Q4 = exp(as.matrix(bl.dat) %*% as.matrix(est.m[i - 1])),
                            Q5 = exp(as.matrix(bl.dat) %*% as.matrix(est.m[i    ])))
        
        probs.m3 <- odds3/rowSums(odds3)
        
        dep.dat <- matrix(NA, nrow = n.i, ncol = 1)
        
        for(j in 1:n.i){
            dep.dat[j,] <- sample(ncol(probs.m3), 1, prob = probs.m3[j,])
        }
        dep.dat <- as.data.frame(dep.dat)
        dep.dat$DEP2 <- ifelse(dep.dat[,1] == 2, 1, 0)
        dep.dat$DEP3 <- ifelse(dep.dat[,1] == 3, 1, 0)
        dep.dat$DEP4 <- ifelse(dep.dat[,1] == 4, 1, 0)
        dep.dat$DEP5 <- ifelse(dep.dat[,1] == 5, 1, 0)
        
        bl.dat[,17:20] <- dep.dat[,2:5]
    }
  
# SMK_TIME CONTINUOUS (0 if SMOKINGN == 1)  
  for(i in 18){
      smk_time <- data.frame(ST = as.matrix(bl.dat) %*% as.matrix(est.m[i]))

      st.dat <- matrix(NA, nrow = n.i, ncol = 1)
      
      for(j in 1:n.i){
        st.dat[j,] <- ifelse(bl.dat$SMOKINGN[j] != 1, rnorm(1, smk_time[j,], 1), 0)
      }
      
       bl.dat[,21] <- st.dat
  }

  prop.table(table(smk.prb))
  #hist(bl.dat$SMK_TIME)
  
  # Add final variables for simulation dataset
    bl.dat$EDcount_C   <- 0
    bl.dat$HOSPcount_C <- 0
    bl.dat$EDcount_O   <- 0
    bl.dat$HOSPcount_O <- 0
    bl.dat$SMOKINGC <- ifelse(bl.dat$SMOKINGF == 0 & bl.dat$SMOKINGN == 0, 1, 0)
    bl.dat$SMOKINGC_TIME <- bl.dat$SMOKINGC * bl.dat$SMK_TIME 
    bl.dat$SMOKINGF_TIME <- bl.dat$SMOKINGF * bl.dat$SMK_TIME 

    bl.dat <- subset(bl.dat, select = -c(SMOKINGN, SMK_TIME))
    bl.dat <- bl.dat[,c(1:3, 5:24, 4, 25:26)]

    # Intercept	age_dx	sexM	c_chf1	c_ihd1	c_cancer1	c_dm1	c_asthma1	c_dementia1	c_depression1	
    # c_anxiety1	c_htn1	rio_cat2.10-44	rio_cat3.45+	deprivationQ2	deprivationQ3	deprivationQ4	deprivationQ5	
    # ned_c	nhosp_c	ned_o	nhosp_o	smoking_current	smoking_former	smoking_current:smk_time_dx	smoking_former:smk_time_dx
    


