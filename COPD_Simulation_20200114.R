#                                                                                                   #
#  AUTHORS: Elizabeth Bond, Petros Pechlivanoglou                                                   #
#                                                                                                   #
#  FUNDING Provided By: Canadian Respiratory Research Network                                       #
#                                                                                                   #
#  CITATIONS:                                                                                       #
#  -Jalal H, et al. An Overview of R in Health Decision Sciences.                                   # 
#   Med. Decis. Making. 2017; 37(3): 735-746.                                                       #
#  -Krijkamp EM, et al. Microsimulation modeling for health decision sciences using R: a tutorial.  #    
#   Med. Decis. Making. 2018; (in press).                                                           #
#                                                                                                   #
#####################################################################################################
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
####    -Surival Plots:                                                            ####
####                     All-Cause Mortality                                       ####
####    -Absolute and Relative difference plots for above                          ####
####    -Length of Stay table                                                      ####  
####                                                                               ####
#######################################################################################
####                                                                               ####
#### USER INPUTS:                                                                  ####
####    -Update fp and output.path                                                 ####
####    -Update simulation inputs (default 14 years, 10,000 people)                ####
####    -Update your covariate of interest                                         ####
####    -Update maximum age alive (default 110 years)                              ####
####                                                                               ####
#######################################################################################
#rm(list=ls())
###############################################
##                                           ##  
##        Install and load packages          ##
##                                           ##
###############################################
# install.packages(c("survival", "demography", "dplyr", "flexsurv", "DAAG", "msm", "gems", "cmprsk", "MASS",
#                    "tidyverse", "readxl", "reshape", "etm", "rlang", "mgcv", "mvtnorm", "corpcor", "gridExtra",
#                    "data.table", "plyr", "ggplot2", "purrr", "survminer", "foreach", "doParallel", "splines", "broom", "gmodels"))
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
 library(gmodels)

###############################################
##                                           ##  
##          Set working directory            ##
##                                           ##
###############################################
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) # set working directory to the folder where the file is saved 
  fp <- "XXXXXXXXXXX"
  
###############################################
##                                           ##  
##  Set where output files should be saved   ##
##                                           ##
###############################################
  output.path <- "XXXXXXXXXXX"
                                      
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

source(file = paste0(fp, "COPD_functions_2.R"))              # source simulation functions
source(file = paste0(fp, "MCC_functions_dong_20190131.R"))   # source MCC functions

######################################                    ###################################### 
######################################   MICROSIMULATION  ######################################
######################################        START       ######################################
######################################                    ######################################

# Baseline covariate generation
# covaraites: AGE SEX CHF IHD CANCER DM ASTHMA DEMENTIA DEPRESSION ANXIETY HTN RIO2 RIO3 DEP2 DEP3 DEP4 DEP5
#             EDcount_C HOSPcount_C EDcount_O HOSPcount_O SMOKINGC SMOKINGF SMOKINGC_TIME SMOKINGF_TIME

  cov_int <- c("SMOKINGC", "SMOKINGF")  # 10 = Current, 01 = Former, 00 = Never

  source(file = paste0(fp, "Cov gen code_2.R"))

p = Sys.time()

## OPEN LOOP (1)
## FOR PATIENT SUBGROUPS
   for(i in 1:1) {

     Xi <- i - 1   # Xi value (1 = YES, 0 = NO)
     set.seed(2) # Set Seed

   # Baseline Individual Charateristics
     X_init <- bl.dat[-c(ncol(bl.dat)-1:-ncol(bl.dat))]

 # Initialize storage for simuation output
   MCC.data.edc <- MCC.data.edo <- MCC.data.hc <- MCC.data.ho  <- MCC.data.mo <- data.frame()

 # Set up parallel
    Cores <- detectCores() - 1
       cl <- makeCluster(Cores)
       registerDoParallel(cl)

## OPEN LOOP (2)
## FOR NUMBER OF SIMULATIONS
   for(k in 1:n.sim){

       results.sim.par <- 0                                                            # Initialize storage for parallel loop
                  inds <- split(seq_len(n.i), sort(rep_len(seq_len(Cores), n.i)))      # Split number of indivuals into groups for parallel loop

## OPEN LOOP (3)
## FOR NUMBER OF INDIVIDUALS (PARALLEL)
   results.sim.par <- foreach(l = seq_along(inds), .combine = rbind, .packages=c("mvtnorm",   "corpcor",  "rlang",  "survival",
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
     event.mat.edc <- event.mat.edo <- event.mat.hc <- event.mat.ho <- event.mat.mo <- matrix(NA , ncol = 5, nrow = n.ind * length(times))

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
         p.EDC <- t(EDcycle(betas1 = EDCHC.norm.mat[k,],                    # EDC to HC
                            betas2 = EDCHO.norm.mat[k,],                    # EDC to HO
                            betas3 = EDCMO.norm.mat[k,], X = as.matrix(X))) # EDC to MO

         p.EDO <- t(EDcycle(betas1 = EDOHC.norm.mat[k,],                    # EDO to HC
                            betas2 = EDOHO.norm.mat[k,],                    # EDO to HO
                            betas3 = EDOMO.norm.mat[k,], X = as.matrix(X))) # EDO to MO
         p.EDCHC <- p.EDC[2, ]
         p.EDCHO <- p.EDC[3, ]
         p.EDCMO <- p.EDC[4, ]
         p.EDOHC <- p.EDO[2, ]
         p.EDOHO <- p.EDO[3, ]
         p.EDOMO <- p.EDO[4, ]
       # MSM
         p.DXEDC <- t(LLTP  (m.DX[,j], shapePar= DXEDC.norm.mat[k,1], scalePar= DXEDC.norm.mat[k,2], B= as.matrix(DXEDC.norm.mat[k,3:ncol(DXEDC.norm.mat)]), X= as.matrix(X[-1])))    # DX to EDC
         p.DXEDO <- t(WeibTP(m.DX[,j], shapePar= DXEDO.norm.mat[k,1], scalePar= DXEDO.norm.mat[k,2], B= as.matrix(DXEDO.norm.mat[k,3:ncol(DXEDO.norm.mat)]), X= as.matrix(X[-1])))    # DX to EDO
         p.DXHC  <- t(LNTP  (m.DX[,j], meanlog=  DXHC.norm.mat [k,1], sdlog=    DXHC.norm.mat [k,2], B= as.matrix(DXHC.norm.mat [k,3:ncol(DXHC.norm.mat)]),  X= as.matrix(X[-1])))    # DX to HC
         p.DXHO  <- t(WeibTP(m.DX[,j], shapePar= DXHO.norm.mat [k,1], scalePar= DXHO.norm.mat [k,2], B= as.matrix(DXHO.norm.mat [k,3:ncol(DXHO.norm.mat)]),  X= as.matrix(X[-1])))    # DX to HO
         p.DXMO  <- t(LNTP  (m.DX[,j], meanlog=  DXMO.norm.mat [k,1], sdlog=    DXMO.norm.mat [k,2], B= as.matrix(DXMO.norm.mat [k,3:ncol(DXMO.norm.mat)]),  X= as.matrix(X[-1])))    # DX to MO
         p.HCDX  <- t(LLTP  (m.HC[,j], shapePar= HCDX.norm.mat [k,1], scalePar= HCDX.norm.mat [k,2], B= as.matrix(HCDX.norm.mat [k,3:ncol(HCDX.norm.mat)]),  X= as.matrix(X[-1])))    # HC to DX
         p.HCMO  <- t(LLTP  (m.HC[,j], shapePar= HCMO.norm.mat [k,1], scalePar= HCMO.norm.mat [k,2], B= as.matrix(HCMO.norm.mat [k,3:ncol(HCMO.norm.mat)]),  X= as.matrix(X[-1])))    # HC to MO
         p.HODX  <- t(LLTP  (m.HO[,j], shapePar= HODX.norm.mat [k,1], scalePar= HODX.norm.mat [k,2], B= as.matrix(HODX.norm.mat [k,3:ncol(HODX.norm.mat)]),  X= as.matrix(X[-1])))    # HO to DX
         p.HOMO  <- t(LLTP  (m.HO[,j], shapePar= HOMO.norm.mat [k,1], scalePar= HOMO.norm.mat [k,2], B= as.matrix(HOMO.norm.mat [k,3:ncol(HOMO.norm.mat)]),  X= as.matrix(X[-1])))    # HO to MO

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

         cov_int_col <- as.numeric(apply(X[ ,cov_int], 1, paste, collapse = ""))

         event.edc$j  <- event.edo$j  <- event.hc$j  <- event.ho$j  <- event.mo$j  <- j
         event.edc$Xi <- event.edo$Xi <- event.hc$Xi <- event.ho$Xi <- event.mo$Xi <- cov_int_col
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

           # write.csv(MCC.data.edc, paste0(output.path, "edc.", Xi,".", k, ".grp.", l, ".csv"), row.names = FALSE)
           # write.csv(MCC.data.edo, paste0(output.path, "edo.", Xi,".", k, ".grp.", l, ".csv"), row.names = FALSE)
           # write.csv(MCC.data.hc,  paste0(output.path, "hc.",  Xi,".", k, ".grp.", l, ".csv"), row.names = FALSE)
           # write.csv(MCC.data.ho,  paste0(output.path, "ho.",  Xi,".", k, ".grp.", l, ".csv"), row.names = FALSE)
           # write.csv(MCC.data.mo,  paste0(output.path, "mo.",  Xi,".", k, ".grp.", l, ".csv"), row.names = FALSE)

        # Organize data for LENGTH OF STAY
        #   -Only keep hosp that are discharged alive
        # HO
          data.M  <- t(m.M)
          data.M  <- melt(data.M)

          data.HO <- t(m.HO)
          data.HO <- melt(data.HO)
          data.HO <- cbind(data.HO, data.M[,3])

                   st.HO  <- data.HO[,2:4]
          colnames(st.HO) <- c("id", "time", "state")
          st.HO$id <- colsplit(st.HO$id, " ", names = c("ind", "id"))[,2]
          st.HO <- as.data.frame(st.HO)
          st.HO <- st.HO %>%
                   group_by(id) %>%
                   mutate(TO = dplyr::lead(state)) %>%
                   ungroup()

          st.HO <- left_join(st.HO, event.mo[,c("id", "Xi")], by = "id")

          time.HO <- subset(st.HO, time != 0)
          time.HO <- time.HO %>%
                     mutate(keep = ifelse(time > lead(time), 1, 0))
          time.HO <- subset(time.HO, is.na(keep) | keep == 1)
          time.HO <- subset(time.HO, TO != "MO")  # Only keep hosp that are discharged alive
          time.HO <- time.HO %>%
                     group_by(id, Xi)%>%
                     summarise_at(vars(time), funs(c.time = sum(time),
                                                       no = n()))
          hosp.HO <- sum(time.HO$c.time) / sum(time.HO$no)

        # HC
          data.HC <- t(m.HC)
          data.HC <- melt(data.HC)
          data.HC <- cbind(data.HC, data.M[,3])

          st.HC  <- data.HC[,2:4]
          colnames(st.HC) <- c("id", "time", "state")
          st.HC$id <- colsplit(st.HC$id, " ", names = c("ind", "id"))[,2]
          st.HC <- as.data.frame(st.HC)
          st.HC <- st.HC %>%
                   group_by(id) %>%
                   mutate(TO = dplyr::lead(state)) %>%
                   ungroup()

          st.HC <- left_join(st.HC, event.mo[,c("id", "Xi")], by = "id")

          time.HC <- subset(st.HC, time != 0)
          time.HC <- time.HC %>%
                     mutate(keep = ifelse(time > lead(time), 1, 0))
          time.HC <- subset(time.HC, is.na(keep) | keep == 1)
          time.HC <- subset(time.HC, TO != "MO")  # Only keep hosp that are discharged alive
          time.HC <- time.HC %>%
                     group_by(id, Xi)%>%
                     summarise_at(vars(time), funs(c.time = sum(time),
                                                       no = n()))
          hosp.HC <- sum(time.HC$c.time) / sum(time.HC$no)

        # Merge and Export
          hosp.data <- as.data.frame(rbind(hosp.HO, hosp.HC))
          hosp.data$l    <- l
          hosp.data$sim  <- k
          return(hosp.data)
         }
     } # CLOSE LOOP (4)
    } # CLOSE LOOP (3)
        write.csv(results.sim.par,  paste0(output.path, "hosp.", Xi, ".", k, ".csv"), row.names = T)
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
#   Length of Stay for:                     #
#         -COPD-Related and                 #
#         -Non-COPD Related Hosp            #
#                                           #
#############################################

Xi = 0

      data.los.list <- list.files(path = output.path, pattern = paste0("hosp.", Xi,".*\\.csv"), recursive = T)
      data.los.sim  <- do.call("rbind",lapply(data.los.list, function(x) read.csv(paste(output.path, x, sep=''), stringsAsFactors = FALSE)))
      
      colnames(data.los.sim) <- c("hosp", "time", "l", "sim")
      data.los.sim[,"hosp"] <- rep(c("HO","HC"), (nrow(data.los.sim)/2))
      
      data.los.sim      <- as.data.frame(data.los.sim)
      
      data.los <- data.los.sim %>%
                  group_by(hosp) %>%
                  summarise(mean  = ci(time, na.rm = TRUE)[1],
                            lowCI = ci(time, na.rm = TRUE)[2],
                            hiCI  = ci(time, na.rm = TRUE)[3],
                            sd    = ci(time, na.rm = TRUE)[4],
                            med   = quantile(time,0.50, na.rm = TRUE),
                            IQRl  = quantile(time,0.25, na.rm = TRUE),
                            IQRu  = quantile(time,0.75, na.rm = TRUE))
      
      write.csv(data.los,  paste0(output.path, "HospSum.csv"), row.names = T)

#############################################
#                                           #
#      PLOTS:                               #
#      -Mean Cumulative Count (MCC)         #
#      -KM curve for death                  #
#                                           #
#############################################
# source(file = paste0(fp, "COPD_plot_functions_4.R"))                # source simulation plot functions
# # Plot MCC
#   pdf(file = paste0(output.path, "MCC.edc.", n.i,".pdf"))
#   MCC.plot(event="edc", ymax_all = 2,  ymax_diff = 2, ymin_diff = -1, ymax_rel = 2.5, ymin_rel = 0.5)
#   dev.off()
#   pdf(file = paste0(output.path, "MCC.edo.", n.i,".pdf"))
#   MCC.plot(event="edo", ymax_all = 10, ymax_diff = 2, ymin_diff = -1, ymax_rel = 2.5, ymin_rel = 0.5)
#   dev.off()
#   pdf(file = paste0(output.path, "MCC.hc.", n.i,".pdf"))
#   MCC.plot(event="hc",  ymax_all = 2,  ymax_diff = 2, ymin_diff = -1, ymax_rel = 2.5, ymin_rel = 0.5)
#   dev.off()
#   pdf(file = paste0(output.path, "MCC.ho.", n.i,".pdf"))
#   MCC.plot(event="ho",  ymax_all = 4,  ymax_diff = 2, ymin_diff = -1, ymax_rel = 2.5, ymin_rel = 0.5)
#   dev.off()
# 
# # Plot KM
#   pdf(file = paste0(output.path, "KM.mo.", n.i,".pdf"))
#   KM.plot(event = "mo")
#   dev.off()

print("finish")
  
  

  





