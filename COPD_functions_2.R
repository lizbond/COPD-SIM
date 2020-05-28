#####################################################################################################################
###                                                                                                               ###
###   Functions for COPD Simulation                                                                               ### 
###                                                                                                               ###
###   Citations: Jalal H, et al. An Overview of R in Health Decision Sciences.                                    ### 
###                         Med. Decis. Making. 2017; 37(3): 735-746.                                             ###      
###              Krijkamp EM, et al. Microsimulation modeling for health decision sciences using R: a tutorial.   ### 
###                         Med. Decis. Making. 2018; (in press).                                                 ###
###                                                                                                               ###
#####################################################################################################################

options(scipen=999)

############################### The Probs() function is used to update transition probabilities based on the health state occupied at cycle j
Probs <- function(M_it, j) {
  # M_it: current health state
  # j:    current cycle
  
  p.it <- matrix(0, nrow = n.s, ncol = n.ind)   # create vector of state transition probabilities
  rownames(p.it) <- v.n                       # name the rows
  
  pos.DX  <- M_it == "DX"       # the time points for the individuals in DX
  pos.EDC <- M_it == "EDC"      # the time points for the individuals in EDC
  pos.EDO <- M_it == "EDO"      # the time points for the individuals in EDO
  pos.HC  <- M_it == "HC"       # the time points for the individuals in HC
  pos.HO  <- M_it == "HO"       # the time points for the individuals in HO
  
  # Update p.it with appropriate transition probabilities
  p.it[, M_it == "DX"]  <- rbind(1 - p.DXEDC[pos.DX]  - p.DXEDO[pos.DX]  - p.DXHC[pos.DX] - p.DXHO[pos.DX] - p.DXMO[pos.DX], p.DXEDC[pos.DX], p.DXEDO[pos.DX], p.DXHC[pos.DX], p.DXHO[pos.DX], p.DXMO[pos.DX]) # trans pr from DX
  p.it[, M_it == "EDC"] <- rbind(1 - p.EDCHC[pos.EDC] - p.EDCHO[pos.EDC] - p.EDCMO[pos.EDC], 0, 0, p.EDCHC[pos.EDC], p.EDCHO[pos.EDC], p.EDCMO[pos.EDC])                                                       # trans pr from EDC
  p.it[, M_it == "EDO"] <- rbind(1 - p.EDOHC[pos.EDO] - p.EDOHO[pos.EDO] - p.EDOMO[pos.EDO], 0, 0, p.EDOHC[pos.EDO], p.EDOHO[pos.EDO], p.EDOMO[pos.EDO])                                                       # trans pr from EDO
  p.it[, M_it == "HC"]  <- rbind(p.HCDX[pos.HC], 0, 0, 1 - p.HCDX[pos.HC] - p.HCMO[pos.HC], 0, p.HCMO[pos.HC])                                                                                                 # trans pr from HC
  p.it[, M_it == "HO"]  <- rbind(p.HODX[pos.HO], 0, 0, 0, 1 - p.HODX[pos.HO] - p.HOMO[pos.HO], p.HOMO[pos.HO])                                                                                                 # trans pr from HO
  p.it[, M_it == "MO"]  <-     c(0, 0, 0, 0, 0, 1)                                                                                                                                                             # trans pr from MO
  
  return(t(p.it))# return the transition probabilities
}

################################# samplev solution ############################################################################################
# samplev : efficient implementation of the rMultinom() function of the Hmisc package 
# probs: the matrix of probabilities for each individual for each state at each time point
# m: the number of values that need to be sampled at a time per individual

samplev <- function (probs, m) {
  d <- dim(probs)
  n <- d[1]
  k <- d[2]
  lev <- dimnames(probs)[[2]]
  if (!length(lev)) 
    lev <- 1:k
  ran <- matrix(lev[1], ncol = m, nrow = n)
  U <- t(probs)
  for(i in 2:k) {
    U[i, ] <- U[i, ] + U[i - 1, ]
  }
  if (any((U[k, ] - 1) > 1e-05))
    stop("error in multinom: probabilities do not sum to 1")
  for (j in 1:m) {
    un <- rep(runif(n), rep(k, n))
    ran[, j] <- lev[1 + colSums(un > U)]
  }
  ran
}



#########################################
###### Flexsurv Output Functions ########
######     Transform Output      ########
#########################################
# WEIBULL 
  # prob
    WeibP <- function (q = 3600, shape, scale, lower.tail= FALSE, log.p=FALSE) {
      WeibP <- pweibull(q = q, shape , scale, lower.tail = lower.tail, log.p = log.p)
      return(WeibP)
    }
  
  # transfrom shape and scale and output transprob
    WeibTP <- function(times, shapePar, scalePar, B, X) 
    {
      logSC <- log(scalePar)
      SH <- as.numeric(shapePar)
      SC <- as.numeric(exp(logSC + X%*%B))
      
      WeibTP <- 1 -  WeibP(times, SH, SC) /
                     WeibP(times - 1, SH, SC) 
      return(WeibTP)
    }

  
# LOGLOGISTIC
  # prob
    LLP <- function (q = 3600, shape, scale, lower.tail= FALSE, log.p=FALSE) {
           LLP <- pllogis(q = q, shape , scale, lower.tail = lower.tail, log.p = log.p)
                  return(LLP)
    }
    
  # transfrom shape and scale and output transprob
    LLTP <- function(times, shapePar, scalePar, B, X) 
    {
      logSC <- log(scalePar)
      SH    <- as.numeric(shapePar)
      SC    <- as.numeric(exp(logSC + X%*%B))
      
      LLTP  <- 1 -  LLP(times, SH, SC) /
                    LLP(times - 1, SH, SC) 
      return(LLTP)
    }
    
# LOGNORMAL
  # prob
    LNP <- function (q = 3600, mean, sd, lower.tail= FALSE, log.p=FALSE) {
      LNP <- plnorm(q = q, mean, sd, lower.tail = lower.tail, log.p = log.p)
      return(LNP)
    }
  
  # transfrom mean and SD and output transprob
    LNTP <- function(times, meanlog, sdlog, B, X) 
    {
      M  <- as.numeric(meanlog + X%*%B)
      SD <- as.numeric(sdlog)
      
      LNTP <- 1 - LNP(times, M, SD) / 
                  LNP(times - 1, M, SD)
      return(LNTP)
    }


########################################
########  Import Model Output   ######## 
########################################
# MSM
  MSM <- read.csv(paste0(fp, "Data/MSM.csv"), header = TRUE)
  as.data.frame(MSM)
 
# From ED_C 
  MNEDC <- t(read.csv(paste0(fp, "Data/EDC.csv"), header = TRUE))
  MNEDC 
 
# From ED_O 
  MNEDO <- t(read.csv(paste0(fp, "Data/EDO.csv"), header = TRUE))
  MNEDO 

###################################
### MNM ED probability function ### 
###################################
  
  EDcycle  <- function(betas1, betas2, betas3, X){
                        exp.resp <- data.frame(Q1  = 1, 
                                               Q2 = exp(X %*% betas1),
                                               Q3 = exp(X %*% betas2),
                                               Q4 = exp(X %*% betas3))
                        
                        res <- exp.resp/rowSums(exp.resp)
    return(res)
  }

################################################################################## 
############################# COVAR Functions ####################################  
##################################################################################   

# read in model output CSV files (MSM and MNM)  
  DXEDC_covar <- as.matrix(read.csv(paste0(fp, "Data/DXEDC_covar.csv")))
  DXEDO_covar <- as.matrix(read.csv(paste0(fp, "Data/DXEDO_covar.csv")))
  DXHC_covar  <- as.matrix(read.csv(paste0(fp, "Data/DXHC_covar.csv")))
  DXHO_covar  <- as.matrix(read.csv(paste0(fp, "Data/DXHO_covar.csv")))
  DXMO_covar  <- as.matrix(read.csv(paste0(fp, "Data/DXMO_covar.csv")))
  HCDX_covar  <- as.matrix(read.csv(paste0(fp, "Data/HCDX_covar.csv")))
  HCMO_covar  <- as.matrix(read.csv(paste0(fp, "Data/HCMO_covar.csv")))
  HODX_covar  <- as.matrix(read.csv(paste0(fp, "Data/HODX_covar.csv")))
  HOMO_covar  <- as.matrix(read.csv(paste0(fp, "Data/HOMO_covar.csv")))
  EDC_covar   <- as.matrix(read.csv(paste0(fp, "Data/EDC_covar.csv")))
  EDO_covar   <- as.matrix(read.csv(paste0(fp, "Data/EDO_covar.csv")))
  
  MSM_csv <- read.csv(paste0(fp, "Data/MSM.csv"))
      DXEDC_mean  <- t(MSM_csv[,2])
      DXEDO_mean  <- t(MSM_csv[,3])
      DXHC_mean   <- t(MSM_csv[,4])
      DXHO_mean   <- t(MSM_csv[,5])
      DXMO_mean   <- t(MSM_csv[,6])
      HCDX_mean   <- t(MSM_csv[,7])
      HCMO_mean   <- t(MSM_csv[,8])
      HODX_mean   <- t(MSM_csv[,9])
      HOMO_mean   <- t(MSM_csv[,10])
  MNM_csv <- read.csv(paste0(fp, "Data/MNM.csv"))
      EDC_mean  <- t(MNM_csv[,1])
      EDO_mean  <- t(MNM_csv[,2])
      
#####################################################
######  Run rmvnorm                            ######
######  1) Transform MSM model output          ######
######  2) Run MSM output in rmvnorm           ######
######  3) Run MNM output in rmvnorm           ######
######  4) Transform rmvnorm MSM model output  ###### 
#####################################################
      
# Transform MSM model output and concatenate with Betas
  # WEIBULL and LOGLOGISTIC      
    trans.NB <- function(par1, par2, B) 
                      {
                        logSH <- log(par1)
                        logSC <- log(par2)
                        B <- B
                        results <- c(logSH, logSC, B)
                        return(results)
                      }
  # LOGNORMAL   
    lnorm.NB <- function(par1, par2, B) 
                      {
                        M  <- par1
                        SD <- log(par2)
                        B  <- B
                        results <- c(M, SD, B)
                        return(results)
                      }
                       
    # Run mvrnorm for multistate models
      # Transform MSM output
        DXEDC.op.par <- as.numeric(trans.NB(par1 = MSM_csv[1,2],  par2 = MSM_csv[2,2],  B = as.matrix(MSM_csv[3:nrow(MSM_csv),2])))    # DX to EDC
        DXEDO.op.par <- as.numeric(trans.NB(par1 = MSM_csv[1,3],  par2 = MSM_csv[2,3],  B = as.matrix(MSM_csv[3:nrow(MSM_csv),3])))    # DX to EDO
        DXHC.op.par  <- as.numeric(lnorm.NB(par1 = MSM_csv[1,4],  par2 = MSM_csv[2,4],  B = as.matrix(MSM_csv[3:nrow(MSM_csv),4])))    # DX to HC
        DXHO.op.par  <- as.numeric(trans.NB(par1 = MSM_csv[1,5],  par2 = MSM_csv[2,5],  B = as.matrix(MSM_csv[3:nrow(MSM_csv),5])))    # DX to HO
        DXMO.op.par  <- as.numeric(lnorm.NB(par1 = MSM_csv[1,6],  par2 = MSM_csv[2,6],  B = as.matrix(MSM_csv[3:nrow(MSM_csv),6])))    # DX to MO
        HCDX.op.par  <- as.numeric(trans.NB(par1 = MSM_csv[1,7],  par2 = MSM_csv[2,7],  B = as.matrix(MSM_csv[3:nrow(MSM_csv),7])))    # HC to DX
        HCMO.op.par  <- as.numeric(trans.NB(par1 = MSM_csv[1,8],  par2 = MSM_csv[2,8],  B = as.matrix(MSM_csv[3:nrow(MSM_csv),8])))    # HC to MO
        HODX.op.par  <- as.numeric(trans.NB(par1 = MSM_csv[1,9],  par2 = MSM_csv[2,9],  B = as.matrix(MSM_csv[3:nrow(MSM_csv),9])))    # HO to DX
        HOMO.op.par  <- as.numeric(trans.NB(par1 = MSM_csv[1,10], par2 = MSM_csv[2,10], B = as.matrix(MSM_csv[3:nrow(MSM_csv),10])))   # HO to MO
      # Run MSM in rmvnorm
        DXEDC.rmv.op <- mvtnorm::rmvnorm(n.sim, DXEDC.op.par, DXEDC_covar)
        DXEDO.rmv.op <- mvtnorm::rmvnorm(n.sim, DXEDO.op.par, DXEDO_covar)
        DXHC.rmv.op  <- mvtnorm::rmvnorm(n.sim, DXHC.op.par,  DXHC_covar)  # Lognormal
        DXHO.rmv.op  <- mvtnorm::rmvnorm(n.sim, DXHO.op.par,  DXHO_covar)
        DXMO.rmv.op  <- mvtnorm::rmvnorm(n.sim, DXMO.op.par,  DXMO_covar)  # Lognormal
        HCDX.rmv.op  <- mvtnorm::rmvnorm(n.sim, HCDX.op.par,  HCDX_covar)
        HCMO.rmv.op  <- mvtnorm::rmvnorm(n.sim, HCMO.op.par,  HCMO_covar)
        HODX.rmv.op  <- mvtnorm::rmvnorm(n.sim, HODX.op.par,  HODX_covar)
        HOMO.rmv.op  <- mvtnorm::rmvnorm(n.sim, HOMO.op.par,  HOMO_covar)
      # Run MNM in rmvnorm
        EDC.rmv.op <- mvtnorm::rmvnorm(n.sim, EDC_mean, EDC_covar)
        EDO.rmv.op <- mvtnorm::rmvnorm(n.sim, EDO_mean, EDO_covar)

# Transform rmvnorm output  
  # Function for Weibull and Loglogistic   
    nb.op <- function(sim.op, nb, n.sim = n.sim){
              for(k in 1:n.sim){
                    par1 <- exp(sim.op[k,1])
                    par2 <- exp(sim.op[k,2])
                    B <-  sim.op[k,3:ncol(sim.op)]
                    nb[k,] <- c(par1, par2, B)
              }
              return(nb)
     }
     
  # Function for Lognormal  
    nb.op.ln <- function(sim.op, nb, n.sim = n.sim){ 
                 for(k in 1:n.sim){
                     par1 <- sim.op[k,1] 
                     par2 <- exp(sim.op[k,2]) 
                     B <- sim.op[k,3:ncol(sim.op)]
                     nb[k,] <- c(par1, par2, B)
                   }
                   return(nb)
    }
    

 # RUN rmvnorm for n.sim
 # Matrices for MSM output
   DXEDC.norm.mat <- DXEDO.norm.mat <- DXHC.norm.mat <- DXHO.norm.mat <- DXMO.norm.mat <-
                     HCDX.norm.mat  <- HCMO.norm.mat <- HODX.norm.mat <- HOMO.norm.mat <- matrix(NA, ncol = ncol(DXEDC.rmv.op), nrow = n.sim)
 # MSM rmvnorm output 
   DXEDC.norm.mat <- nb.op   (sim.op = DXEDC.rmv.op, nb = DXEDC.norm.mat, n.sim)
   DXEDO.norm.mat <- nb.op   (sim.op = DXEDO.rmv.op, nb = DXEDO.norm.mat, n.sim)
   DXHC.norm.mat  <- nb.op.ln(sim.op = DXHC.rmv.op,  nb = DXHC.norm.mat,  n.sim) # Lognormal
   DXHO.norm.mat  <- nb.op   (sim.op = DXHO.rmv.op,  nb = DXHO.norm.mat,  n.sim)
   DXMO.norm.mat  <- nb.op.ln(sim.op = DXMO.rmv.op,  nb = DXMO.norm.mat,  n.sim) # Lognormal
   HCDX.norm.mat  <- nb.op   (sim.op = HCDX.rmv.op,  nb = HCDX.norm.mat,  n.sim)
   HCMO.norm.mat  <- nb.op   (sim.op = HCMO.rmv.op,  nb = HCMO.norm.mat,  n.sim)
   HODX.norm.mat  <- nb.op   (sim.op = HODX.rmv.op,  nb = HODX.norm.mat,  n.sim)
   HOMO.norm.mat  <- nb.op   (sim.op = HOMO.rmv.op,  nb = HOMO.norm.mat,  n.sim)
 # MNM rmvnorm output
   EDC.norm.mat   <- EDC.rmv.op
   EDC.dim <- ncol(EDC.norm.mat) / 3
       EDCHC.norm.mat <- EDC.rmv.op[,1                :  EDC.dim]
       EDCHO.norm.mat <- EDC.rmv.op[,(EDC.dim + 1)    : (EDC.dim * 2)]
       EDCMO.norm.mat <- EDC.rmv.op[,(EDC.dim * 2 + 1): ncol(EDC.norm.mat)]
   EDO.norm.mat   <- EDO.rmv.op
   EDO.dim <- ncol(EDO.norm.mat) / 3
       EDOHC.norm.mat <- EDO.rmv.op[,1                :  EDO.dim]
       EDOHO.norm.mat <- EDO.rmv.op[,(EDO.dim + 1)    : (EDO.dim * 2)]
       EDOMO.norm.mat <- EDO.rmv.op[,(EDO.dim * 2 + 1): ncol(EDO.norm.mat)]
       

       



  

