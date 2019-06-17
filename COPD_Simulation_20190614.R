                                                                                      # Authors: Elizabeth Bond, Petros Pechlivanoglou                                                  #
  #########        #########      #########    #########                              # Funding Provided By: Canadian Respiratory Research Network                                      #
 ###     ###      ###     ###     ###    ###   ###    ###                             #                                                                                                 #
###       ###    ###       ###    ###     ###  ###     ###                            # Citations:                                                                                      #
###             ###         ###   ###     ###  ###      ###                           # -Jalal H, et al. An Overview of R in Health Decision Sciences.                                  # 
###             ###         ###   ###    ###   ###      ###   #####  ##  ###   ###    #                  Med. Decis. Making. 2017; 37(3): 735-746.                                      #
###             ###         ###   ########     ###      ###   ##     ##  ## # # ##    # -Krijkamp EM, et al. Microsimulation modeling for health decision sciences using R: a tutorial. #    
###       ###    ###       ###    ###          ###     ###    #####  ##  ##  #  ##    #                      Med. Decis. Making. 2018; (in press).                                      #
 ###     ###      ###     ###     ###          ###    ###        ##  ##  ##     ##    #                                                                                                 #
  #########        #########      ###          #########      #####  ##  ##     ##    #                                                                                                 #            
                                                                                      #                                                                                                 #
#########################################################################################################################################################################################
####                                                                               ####
####  COPD SIMULATION  (Parallel)                                                  ####
####                                                                               ####
####        Six States:  Post-COPD Diagnosis                                       ####
####                     COPD-Related ED Visit                                     ####
####                     Non-COPD Related ED Visit                                 ####
####                     COPD-Related Hospitalization                              ####
####                     Non-COPD Related Hospitalization                          #### 
####                     All-Cause Mortality                                       ####
####                                                                               ####
#######################################################################################
####                                                                               ####
#### FINAL OUTPUT FILES:                                                           ####
####    -Mean Cumulative Count (MCC) Plots:                                        ####
####                     COPD-Related ED Visit                                     ####
####                     Non-COPD Related ED Visit                                 ####
####                     COPD-Related Hospitalization                              ####
####                     Non-COPD Related Hospitalization                          #### 
####    -Surival and Cumulative Incidence Plots:                                   ####
####                     All-Cause Mortality                                       ####
####    -Absolute and Relative difference plots for above                          ####
####                                                                               ####
#######################################################################################
####                                                                               ####
#### USER INPUTS:                                                                  ####
####    -Update where to save .pdf and .csv                                        ####
####            see line                                                           ####
####    -Update simulation inputs (default 14 years, 10,000 people)                ####
####            see line                                                           ####
####    -Update your covariate of interest (default smokers vs. non-smokers)       ####
####            see line                                                           ####
####    -Update maximum age alive (default 110 years)                              ####
####            see line                                                           ####
####                                                                               ####
#######################################################################################
rm(list=ls())
###############################################
##                                           ##  
##        Install and load packages          ##
##                                           ##
###############################################
install.packages(c("survival", "demography", "dplyr", "flexsurv", "DAAG", "msm", "gems", "cmprsk", "MASS",
                   "tidyverse", "readxl", "reshape", "etm", "rlang", "mgcv", "mvtnorm", "corpcor", "gridExtra",
                   "data.table", "plyr", "ggplot2", "purrr", "survminer", "foreach", "doParallel", "splines", "broom"))
 library(plyr) 
 library(survival) 
 library(demography)
 library(dplyr) 
 library(flexsurv) 
 library(DAAG) 
 library(msm)
 library(gems) 
 library(cmprsk)
 library(MASS)
 library(tidyverse) 
 library(readxl) 
 library(reshape)
 library(etm) 
 library(rlang) 
 library(mgcv) 
 library(mvtnorm) 
 library(corpcor)
 library(gridExtra)
 library(data.table) 
 library(ggplot2) 
 library(purrr)
 library(survminer) 
 library(foreach) 
 library(doParallel) 
 library(splines)
 library(broom)

###############################################
##                                           ##  
##          Set working directory            ##
##                                           ##
###############################################
  #setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) # set working directory to the folder where the file is saved 

###############################################
##                                           ##  
##  Set where output files should be saved   ##
##                                           ##
###############################################
  #output.path <- "XXXXXXXXXXX"
  

                                                 
######################################                    ###################################### 
###################################### Simulation Inputs  ######################################
###################################### Update Accordingly ######################################
######################################                    ######################################
   
v.n     <- c("DX","EDC","EDO","HC", "HO","MO")  # model states
n.s     <- length(v.n)                          # number of states
n.t     <- 14                                   # number of years to run
c.l     <- 1/365                                # cycle length (a day)
n.i     <- 10000                                # number of individuals (for microsimulation)
times   <- seq(from = 0, to = n.t*365, by = 1)  # sequence of times to be considered in the model
cycles  <- length(times)                        # number of cycles
n.sim   <- 200                                  # number of simulations

source(file = "COPD_functions.R")                # source simulation functions
source(file = "MCC_functions_dong_20190131.R")   # source MCC functions


######################################                    ###################################### 
######################################   MICROSIMULATION  ######################################
######################################        START       ######################################
######################################                    ######################################

p = Sys.time()

## OPEN LOOP (1) 
## FOR PATIENT SUBGROUPS 
   for(i in 1:2) {    

     Xi <- i-1   # Xi value (1 = YES, 0 = NO)
     set.seed(2) # Set Seed
    
   # Baseline Individual Charateristics                                                                                                           ### UPDATING YOUR COVARIATE OF INTEREST                   ###
     X_init <- data.frame(1, rnorm(n.i, 63.98, 14), rbinom(n.i, 1, 0.449), rbinom(n.i, 1, Xi),     rbinom(n.i, 1, 0.117),                         ###                                                       ###           
                             rbinom(n.i, 1, 0.181), rbinom(n.i, 1, 0.137), rbinom(n.i, 1, 0.197),  rbinom(n.i, 1, 0.190),                         ### Covariates appear in the same order as colnames below ###           
                             rbinom(n.i, 1, 0.044), rbinom(n.i, 1, 0.048), rbinom(n.i, 1, 0.180),  rbinom(n.i, 1, 0.528),                         ### Xi is your covariate(s) or subgroup of interest       ###
                             t(rmultinom(n.i, 1, c( 0.46,0.398, 0.142))[-1,]), t(rmultinom(n.i, 1, c(0.136, 0.176, 0.207, 0.232, 0.249))[-1,]),   ### Xi (1 = YES, 0 = NO)                                  ### 
                             0 , 0 , 0 , 0)                                                                                                       ### Current setting: smokers vs. non-smokers              ###             
     
         colnames(X_init) <- c("INT", "AGE",       "SEX",         "SMOKING",   "CHF",
                                      "IHD",       "CANCER",      "DM",        "ASTHMA",
                                      "DEMENTIA",  "DEPRESSION",  "ANXIETY",   "HTN",
                                      "RIO2", "RIO3",             "DEP2", "DEP3", "DEP4", "DEP5",
                                      "EDcount_C", "HOSPcount_C", "EDcount_O", "HOSPcount_O")

 # Initialize storage for simuation output
   MCC.data.edc <- MCC.data.edo <- MCC.data.hc <- MCC.data.ho  <- MCC.data.mo <- data.frame()

 # Set up parallel
    Cores <- detectCores()
       cl <- makeCluster(Cores)
       registerDoParallel(cl)

## OPEN LOOP (2)      
## FOR NUMBER OF SIMULATIONS
   for(k in 1:n.sim){    

       results.sim.par <- 0                                                            # Initialize storage for parallel loop
                  inds <- split(seq_len(n.i), sort(rep_len(seq_len(Cores), n.i)))      # Slipt number of indivuals into groups for parallel loop
                  
## OPEN LOOP (3)                 
## FOR NUMBER OF INDIVIDUALS (PARALLEL)
   results.sim.par <- foreach(l = seq_along(inds), .packages=c("mvtnorm",   "corpcor",  "rlang",  "survival",
                                                               "flexsurv",  "MASS",     "msm",    "gems", 
                                                               "dplyr",     "DAAG",     "cmprsk", "reshape", 
                                                               "splines",   "survival", "etm",    "data.table", 
                                                               "plyr",      "ggplot2",  "purrr",  "survminer", 
                                                               "gridExtra", "foreach",  "doParallel")) %dopar% {   

     n.ind <- length(inds[[l]])  # number of groups 

     v.M_Init   <- rep("DX", n.ind)                                                  # initial state for all individuals
     v.EDC_Init <- v.EDO_Init <- v.HC_Init <- v.HO_Init <- rep(0, times = n.ind)     # no illness onset at start of model
     v.DX_Init  <- rep(1, times = n.ind)                                             # no illness onset at start of model

   # Total ED and Hosp count matrices
     EDCcountMat <- HCcountMat <- EDOcountMat <- HOcountMat <- matrix(0, nrow = n.ind, ncol = length(times),
                                                                      dimnames = list(paste("ind", 1:n.ind, sep = " "),  # name the rows ind1, ind2, ind3, etc.
                                                                                      paste(times, sep = " ")))          # name the columns cycle1, cycle2, cycle3, etc.
   # ED, Hosp, and Mortality indicatior matrices
     EDCcountMat.ind  <- HCcountMat.ind <- EDOcountMat.ind <- HOcountMat.ind <- MOcountMat.ind <- matrix(0, nrow = n.ind, ncol = length(times),
                                                                                                         dimnames = list(paste("ind", 1:n.ind, sep = " "), # name the rows ind1, ind2, ind3, etc.
                                                                                                                         paste(times, sep = " ")))         # name the columns cycle1, cycle2, cycle3, etc.

   # Store results from every cycle loop
     event.mat.edc <- event.mat.edo <- event.mat.hc <- event.mat.ho <- event.mat.mo <- matrix(NA , ncol = 5, nrow = n.ind*length(times))

   # Assign X values
     X = X_init[inds[[l]],]

     # CREATE RESULT MATRIX
     # m.M: health state for each patient at each cycle
       m.M <-  matrix(nrow = n.ind, ncol = length(times),
                      dimnames = list(paste("ind", 1:n.ind, sep = " "),     # name the rows ind1, ind2, ind3, etc.
                                      paste("cycle", times, sep = " ")))    # name the columns cycle1, cycle2, cycle3, etc.

     # CREATE ATTRIBUTE MATRICES
     # m.DX, m.EDC, m.EDO, m.HC, m.HO: track time in state
       m.DX <- m.EDC <- m.EDO <- m.HC <- m.HO <- matrix(0, nrow = n.ind, ncol = length(times),
                                                        dimnames = list(paste("ind", 1:n.ind, sep = " "),     # name the rows ind1, ind2, ind3, etc.
                                                                        paste("cycle", times, sep = " ")))    # name the columns cycle1, cycle2, cycle3, etc.

     # Initialize health-state matrices
       m.M  [,1] <- v.M_Init    # initial health state for individual i
       m.DX [,1] <- v.DX_Init   # initialize time since DX  for individual i
       m.EDC[,1] <- v.EDC_Init  # initialize time since EDC for individual i
       m.EDO[,1] <- v.EDO_Init  # initialize time since EDO for individual i
       m.HC [,1] <- v.HC_Init   # initialize time since HC  for individual i
       m.HO [,1] <- v.HO_Init   # initialize time since HO  for individual i

     # Import parameters and define new vars in function and simulation
       j = 0
       
       
## OPEN LOOP (4)                 
## FOR TIME (number of cycles)
   for (t in  times[-length(times)]) { 

     j = j + 1

     # TRANSITION PROBABILITIES for n.sim
       # MNM
         p.EDCHC <- t(EDcycle(betas=EDCHC.norm.mat[k,], as.matrix(X)))     # EDC to HC
         p.EDCHO <- t(EDcycle(betas=EDCHO.norm.mat[k,], as.matrix(X)))     # EDC to HO
         p.EDCMO <- t(EDcycle(betas=EDCMO.norm.mat[k,], as.matrix(X)))     # EDC to MO
         p.EDOHC <- t(EDcycle(betas=EDOHC.norm.mat[k,], as.matrix(X)))     # EDO to HC
         p.EDOHO <- t(EDcycle(betas=EDOHO.norm.mat[k,], as.matrix(X)))     # EDO to HO
         p.EDOMO <- t(EDcycle(betas=EDOMO.norm.mat[k,], as.matrix(X)))     # EDO to MO
       # MSM
         p.DXEDC <- t(LLTP  (m.DX[,j], shapePar= DXEDC.norm.mat[k,1], scalePar= DXEDC.norm.mat[k,2], B= as.matrix(DXEDC.norm.mat[k,3:24]), X= as.matrix(X[-1])))    # DX to EDC
         p.DXEDO <- t(WeibTP(m.DX[,j], shapePar= DXEDO.norm.mat[k,1], scalePar= DXEDO.norm.mat[k,2], B= as.matrix(DXEDO.norm.mat[k,3:24]), X= as.matrix(X[-1])))    # DX to EDO
         p.DXHC  <- t(LNTP  (m.DX[,j], meanlog=  DXHC.norm.mat [k,1], sdlog=    DXHC.norm.mat [k,2], B= as.matrix(DXHC.norm.mat [k,3:24]), X= as.matrix(X[-1])))    # DX to HC
         p.DXHO  <- t(WeibTP(m.DX[,j], shapePar= DXHO.norm.mat [k,1], scalePar= DXHO.norm.mat [k,2], B= as.matrix(DXHO.norm.mat [k,3:24]), X= as.matrix(X[-1])))    # DX to HO
         p.DXMO  <- t(LNTP  (m.DX[,j], meanlog=  DXMO.norm.mat [k,1], sdlog=    DXMO.norm.mat [k,2], B= as.matrix(DXMO.norm.mat [k,3:24]), X= as.matrix(X[-1])))    # DX to MO
         p.HCDX  <- t(LLTP  (m.HC[,j], shapePar= HCDX.norm.mat [k,1], scalePar= HCDX.norm.mat [k,2], B= as.matrix(HCDX.norm.mat [k,3:24]), X= as.matrix(X[-1])))    # HC to DX
         p.HCMO  <- t(LLTP  (m.HC[,j], shapePar= HCMO.norm.mat [k,1], scalePar= HCMO.norm.mat [k,2], B= as.matrix(HCMO.norm.mat [k,3:24]), X= as.matrix(X[-1])))    # HC to MO
         p.HODX  <- t(LLTP  (m.HO[,j], shapePar= HODX.norm.mat [k,1], scalePar= HODX.norm.mat [k,2], B= as.matrix(HODX.norm.mat [k,3:24]), X= as.matrix(X[-1])))    # HO to DX
         p.HOMO  <- t(LLTP  (m.HO[,j], shapePar= HOMO.norm.mat [k,1], scalePar= HOMO.norm.mat [k,2], B= as.matrix(HOMO.norm.mat [k,3:24]), X= as.matrix(X[-1])))    # HO to MO

       # Get transition probabilities based on current health state and time since illness
         m.p <- Probs(m.M[,j], j)

       # Sample the  health state at t + 1 based on transition probabilities p
         #m.M[,j + 1] <- samplev(probs = m.p, m = 1)                              # no maximum age alive 
         m.M[,j + 1] <- ifelse(X$AGE >= 110, "MO", samplev(probs = m.p, m = 1))   # maximum age alive is 110

       # Update time since illness onset for t + 1
         m.DX [,j + 1] <- ifelse(m.M[,j + 1] == "DX" , m.DX [,j] + 1, 0)
         m.EDC[,j + 1] <- ifelse(m.M[,j + 1] == "EDC", m.EDC[,j] + 1, 0)
         m.EDO[,j + 1] <- ifelse(m.M[,j + 1] == "EDO", m.EDO[,j] + 1, 0)
         m.HC [,j + 1] <- ifelse(m.M[,j + 1] == "HC" , m.HC [,j] + 1, 0)
         m.HO [,j + 1] <- ifelse(m.M[,j + 1] == "HO" , m.HO [,j] + 1, 0)

       # Update ED and Hosp counts in X
         X$EDcount_C  [(m.M[,j + 1] == "EDC") & (m.M[,j] != "EDC")] <- X$EDcount_C  [(m.M[,j + 1] == "EDC") & (m.M[,j] != "EDC")] + 1
         X$EDcount_O  [(m.M[,j + 1] == "EDO") & (m.M[,j] != "EDO")] <- X$EDcount_C  [(m.M[,j + 1] == "EDO") & (m.M[,j] != "EDO")] + 1
         X$HOSPcount_C[(m.M[,j + 1] == "HC")  & (m.M[,j] != "HC")]  <- X$HOSPcount_C[(m.M[,j + 1] == "HC")  & (m.M[,j] != "HC")]  + 1
         X$HOSPcount_O[(m.M[,j + 1] == "HO")  & (m.M[,j] != "HO")]  <- X$HOSPcount_O[(m.M[,j + 1] == "HO")  & (m.M[,j] != "HO")]  + 1

       # Update ED and HOSP count
         EDCcountMat[,j + 1]  <-  ifelse(m.M[,j + 1] == "EDC" & m.M[,j] != "EDC", EDCcountMat[,j] + 1, EDCcountMat [,j])
         EDOcountMat[,j + 1]  <-  ifelse(m.M[,j + 1] == "EDO" & m.M[,j] != "EDO", EDOcountMat[,j] + 1, EDOcountMat [,j])
         HCcountMat [,j + 1]  <-  ifelse(m.M[,j + 1] == "HC"  & m.M[,j] != "HC",  HCcountMat [,j] + 1,  HCcountMat [,j])
         HOcountMat [,j + 1]  <-  ifelse(m.M[,j + 1] == "HO"  & m.M[,j] != "HO",  HOcountMat [,j] + 1,  HOcountMat [,j])

       # Update X.event.list for MCC (event of interest 1, competing event 2, death occurrence 3, death continuance 4)
         event.edc <- data.frame(out = ifelse(m.M[,j + 1] == "EDC" & m.M[,j] != "EDC", 1,
                                       ifelse(m.M[,j + 1] == "EDO" & m.M[,j] != "EDO", 2,
                                       ifelse(m.M[,j + 1] == "HC"  & m.M[,j] != "HC",  2,
                                       ifelse(m.M[,j + 1] == "HO"  & m.M[,j] != "HO",  2,
                                       ifelse(m.M[,j + 1] == "MO"  & m.M[,j] != "MO",  3,
                                       ifelse(m.M[,j + 1] == "MO"  & m.M[,j] == "MO",  4, 0)))))))

         event.edo <- data.frame(out = ifelse(m.M[,j + 1] == "EDO" & m.M[,j] != "EDO", 1,
                                       ifelse(m.M[,j + 1] == "EDC" & m.M[,j] != "EDC", 2,
                                       ifelse(m.M[,j + 1] == "HC"  & m.M[,j] != "HC",  2,
                                       ifelse(m.M[,j + 1] == "HO"  & m.M[,j] != "HO",  2,
                                       ifelse(m.M[,j + 1] == "MO"  & m.M[,j] != "MO",  3,
                                       ifelse(m.M[,j + 1] == "MO"  & m.M[,j] == "MO",  4, 0)))))))

         event.hc  <- data.frame(out = ifelse(m.M[,j + 1] == "HC"  & m.M[,j] != "HC",  1,
                                       ifelse(m.M[,j + 1] == "HC"  & m.M[,j] == "HC",  2,   # added so that subjects cannot go to the hospital if they are already in the hospital
                                       ifelse(m.M[,j + 1] == "EDO" & m.M[,j] != "EDO", 2,
                                       ifelse(m.M[,j + 1] == "EDC" & m.M[,j] != "EDC", 2,
                                       ifelse(m.M[,j + 1] == "HO"  & m.M[,j] != "HO",  2,
                                       ifelse(m.M[,j + 1] == "MO"  & m.M[,j] != "MO",  3,
                                       ifelse(m.M[,j + 1] == "MO"  & m.M[,j] == "MO",  4, 0))))))))

         event.ho  <- data.frame(out = ifelse(m.M[,j + 1] == "HO"  & m.M[,j] != "HO",  1,
                                       ifelse(m.M[,j + 1] == "HO"  & m.M[,j] == "HO",  2,   # added so that subjects cannot go to the hospital if they are already in the hospital
                                       ifelse(m.M[,j + 1] == "EDO" & m.M[,j] != "EDO", 2,
                                       ifelse(m.M[,j + 1] == "HC"  & m.M[,j] != "HC",  2,
                                       ifelse(m.M[,j + 1] == "EDC" & m.M[,j] != "EDC", 2,
                                       ifelse(m.M[,j + 1] == "MO"  & m.M[,j] != "MO",  3,
                                       ifelse(m.M[,j + 1] == "MO"  & m.M[,j] == "MO",  4, 0))))))))

         event.mo  <- data.frame(out = ifelse(m.M[,j + 1] == "MO"  & m.M[,j] != "MO",  1,
                                       ifelse(m.M[,j + 1] == "MO"  & m.M[,j] == "MO",  4, 0)))

         event.edc$j  <- event.edo$j  <- event.hc$j  <- event.ho$j  <- event.mo$j  <- j
         event.edc$Xi <- event.edo$Xi <- event.hc$Xi <- event.ho$Xi <- event.mo$Xi <- Xi
         event.edc$k  <- event.edo$k  <- event.hc$k  <- event.ho$k  <- event.mo$k  <- k
         event.edc$id <- event.edo$id <- event.hc$id <- event.ho$id <- event.mo$id <- inds[[l]][1:n.ind]

         event.mat.edc[(1 + (j - 1) * n.ind):(j * n.ind),] <- as.matrix(event.edc)
         event.mat.edo[(1 + (j - 1) * n.ind):(j * n.ind),] <- as.matrix(event.edo)
         event.mat.hc [(1 + (j - 1) * n.ind):(j * n.ind),] <- as.matrix(event.hc)
         event.mat.ho [(1 + (j - 1) * n.ind):(j * n.ind),] <- as.matrix(event.ho)
         event.mat.mo [(1 + (j - 1) * n.ind):(j * n.ind),] <- as.matrix(event.mo)

       # Update X Age
         X$AGE <- X$AGE + 1/365.25

       # Organize events and export CSV for last cycle of j
         if (j == length(times)-1)  {

           if(any(grepl("package:plyr", search()))) detach("package:plyr") else message("plyr not loaded")
           MCC.data.sim <- function(data){
                           colnames(data) <- c("event", "day", "X", "sim.no", "id")
             # add indicator variables
               one <- rep(1, nrow(data))
               data <- as.data.frame(cbind(data, one))
               data.n <- data %>%
                         group_by(sim.no, id) %>%
                         mutate(obs.num = rev(cumsum(rev(one)))) %>%
                         ungroup()
             # drop 0 values that arent the last value
               data.n$event <- ifelse(data.n$obs.num == 1 & data.n$event == 0, 99, data.n$event)
               data.n <- subset(data.n, event != 0)
               data.n$event[data.n$event == 99] <- 0
             # drop death event if not first recorded death event
               data.n <- subset(data.n, event != 4)
             # drop added indicator columns
               data.n <- data.n %>%
                         dplyr::select(-c(one, obs.num))
             # run event.data
               data.n <- as.data.frame(data.n)
               data.n$ind.X.sim <- paste(data.n$X, data.n$sim.no, sep = "")
               data.n <- data.n[order(data.n[,"id"]), ]

               return(data.n)
           }

         # Save output to CSV files
           MCC.data.edc <- MCC.data.sim(data = event.mat.edc)
           MCC.data.edo <- MCC.data.sim(data = event.mat.edo)
           MCC.data.hc  <- MCC.data.sim(data = event.mat.hc)
           MCC.data.ho  <- MCC.data.sim(data = event.mat.ho)
           MCC.data.mo  <- MCC.data.sim(data = event.mat.mo)

           write.csv(MCC.data.edc, paste0(output.path, "edc.",Xi,".", k, ".grp.", l, ".csv"), row.names = FALSE)
           write.csv(MCC.data.edo, paste0(output.path, "edo.",Xi,".", k, ".grp.", l, ".csv"), row.names = FALSE)
           write.csv(MCC.data.hc,  paste0(output.path, "hc.", Xi,".", k, ".grp.", l, ".csv"), row.names = FALSE)
           write.csv(MCC.data.ho,  paste0(output.path, "ho.", Xi,".", k, ".grp.", l, ".csv"), row.names = FALSE)
           write.csv(MCC.data.mo,  paste0(output.path, "mo.", Xi,".", k, ".grp.", l, ".csv"), row.names = FALSE)

         }

     } # CLOSE LOOP (4) 
         rbind (X$AGE)

    } # CLOSE LOOP (3)
   } # CLOSE LOOP (2)
  } # CLOSE LOOP (1)
       comp.time = Sys.time() - p
       print(n.i)
       print(n.sim)
       print(comp.time)
       print(Cores)
       stopCluster(cl)
       
#############################################
#                                           #
#      PLOTS:                               #
#      -Mean Cumulative Count (MCC)         #
#      -Cumulative Incidence (CI)           #
#      -KM curve for death                  #
#                                           #
#############################################
source(file = "COPD_plot_functions.R")                # source simulation plot functions
# Plot MCC  
  par(mfrow=c(3,1))
  pdf(file = paste0("MCC.edc.", n.i,".pdf"))
  MCC.plot(event="edc")
  dev.off()
  pdf(file = paste0("MCC.edo.", n.i,".pdf"))
  MCC.plot(event="edo")
  dev.off()
  pdf(file = paste0("MCC.hc.", n.i,".pdf"))
  MCC.plot(event="hc")
  dev.off()
  pdf(file = paste0("MCC.ho.", n.i,".pdf"))
  MCC.plot(event="ho")
  dev.off()

# Plot CI and KM   
  m <- rbind(c(1, 1), c(2, 3))
  layout(m)
  par(mar = c(4, 4, 1, 1), xpd = NA)
  pdf(file = paste0("KM.mo.", n.i,".pdf"))
  KM.plot(event = "mo")
  dev.off()
  pdf(file = paste0("CI.mo.", n.i,".pdf"))
  CI.plot(event = "mo")
  dev.off() 



  
  
  
  
print("finish")
  
  

  





