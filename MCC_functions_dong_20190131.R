#################################################################################
#                                                                               
#                        MEAN CUMULATIVE COUNT                                  
#                                                                               
# Last updated February 26, 2015.                                                  
#################################################################################
# R code to accompany:
#
#   Dong H, Robison LL, Leisenring WM, Martin LJ, Armstrong GT, Yasui Y. (in press).
#   Estimating the burden of recurrent events in the presence of competing risks: 
#   The method of mean cumulative count. American Journal of Epidemiology. 
#
# NOTE:
#   Our purpose of providing a simple hypothetical example and the computation code is 
#   that it would serve as a useful tutorial for researchers who want to learn how to 
#   apply the method of Mean Cumulative Count. Should you have any question, please conctact 
#   Huiru Dong by email: huiru@ualberta.ca.
#
# R version 3.1.1 (2014-07-10)
# Platform: x86_64-w64-mingw32/x64 (64-bit)


############ -------------  PARAMETERS -------------- #########
############----- id:    Participants's ID
############----- time:  Follow up time of participant to interest/competing-risk event/censoring
############----- cause: Event of interest is coded as 1; competing-risk event as 2; censoring as 0 
############----- type:  Character string specifying the type of calcultion. 
#				 Possible values are "MCC" (MCC calculation by equation) and "SCI" (sum of cumulative incidence). 


#################################################################################

######### Run this only once to install the package #########
#install.packages("reshape")   
#install.packages("etm")

# library(reshape)
# library(splines)
# library(survival)
# library(etm)

###################   This is the main fucntion to run  ###################

MCC=function(id,time,cause,type)
{
  if(type=="MCC"){
    outdata=MCC_eq(id,time,cause)
  }
  if(type=="SCI"){
    outdata=SCI(id,time,cause)
  }
  return(outdata)
}



###################  Calculate the Sum of Cumulative Incidence  ########## 


SCI=function(id,time,cause)
{
  indata=data.frame(id,time,cause)
  ii=order(indata$id,indata$time,-indata$cause)  ##if time is the same for 0/1 or 0/2 of the same person, make 0 last
  sort_data=indata[ii,]
  sort_id=sort_data$id
  
  #The time interval we are interested#
  freq=1
  time.interval=seq(min(time),max(time),freq)
  
  #Calculate the total number of participants: N#
  Nodup_id=unique(id)
  N=length(unique(id)) 
  
  
  ######## Huiru's code has loop for each person and loop for each id, which was very slow (I actually cannot 
  ### have it work when there in my CCSS --too much time) I will revise it to avoid looping pfr person 1: no. of person
  
  sort_data$first <- ave(sort_data$time, sort_data$id, FUN = seq_along)
  M=max(sort_data$first) 		#maximum numer of events
  
  ### Take the last row so we know the maximum number of events per person
  event.number=do.call(rbind, lapply(split(sort_data, sort_data$id), tail, 1))[,c("id","first")]
  colnames(event.number)=c("id","maxE")
  
  alldata=merge(sort_data,event.number,by.x="id",by.y="id")
  
  #### make data for cumulative incidence, with dimension M*N
  data.MCC=NULL
  for(i in 1:M){
    if(i==1){
      data_temp=alldata[alldata$first==1,]
      data_temp$m_event=i
      data.MCC=rbind(data.MCC,data_temp)
    }
    if(i>1){
      ## for those with i or more records, take the ith record
      the_ith=alldata[alldata$first==i,]
      ##for those with <i records, take the last record. If the last record is an event 1, then change it to 0;
      ## because it contributed in data (i-1) already.
      the_last=alldata[alldata$maxE<i & alldata$first==alldata$maxE,]
      the_last$cause[the_last$cause==1]=0
      CIdata_ith=rbind(the_ith,the_last)
      CIdata_ith$m_event=i
      data.MCC=rbind(data.MCC,CIdata_ith)
    }
  }
  #Calculate cumulative incidence for M times event#
  #Store the results in matrix MCC.base#
  MCC.base=NULL
  for (j in 1:M)
  {
    dataj=data.MCC[data.MCC$m_event==j,]
    if(max(dataj$cause!=0)){  ##in example simple2, the 2nd time all cause=0
      #Calculate cumulative incidence
      P=summary(etmCIF(Surv(time, cause!=0)~ 1,data=dataj,etype=cause,failcode=1))[[1]]$'CIF 1'$P
      #Calculate the change/increase
      Deta_P=c(P[1],diff(P))  ## at time 1, the difference is P-0
      
      #Keep the time points
      Time=summary(etmCIF(Surv(time, cause!=0)~ 1,data=dataj,etype=cause,failcode=1))[[1]]$'CIF 1'$time
      #Combine all the information
      if (is.null(P)==FALSE & is.null(Time)==FALSE)
      { MCC.base=rbind(MCC.base,cbind(Time,P, Deta_P, CumI=j))}
    }
  }
  
  #Only keep the time points that affect MCC#
  Nodup_MCC.base=MCC.base[!duplicated(MCC.base[,c(2,4)]),]
  
  #Sort by event dates#
  jj=order(Nodup_MCC.base[,1])
  sort_MCC.base=Nodup_MCC.base[jj,]
  
  #MCC is calculated by using the summation of cumulative incidences for all event of interest (first and recurrent)#
  MCC=cumsum(sort_MCC.base[,3])
  
  combine_MCC=cbind(sort_MCC.base,MCC)
  
  #Only show the time points that have MCC change, and remove irrelevant information#
  MCC.final=aggregate(combine_MCC[,5],list(combine_MCC[,1]),max)###
  colnames(MCC.final)=c("Time","SumCIs")###
  
  return(MCC.final)
}
#-------------------------------------------------------------------------------------------------#



######### -----------------------  MCC function by equation ---------------##################
MCC_eq=function(id,time,cause)
{
  
  ntotal=length(unique(id)) 
  count=rep(1,length(id))
  indata=data.frame(id,time,cause,count)
  
  freq_cause=aggregate(count~time+cause, data=indata,sum)
  
  lifetable_1 <- cast(freq_cause, time~cause,value="count",fill=0)  
  colnames(lifetable_1)[colnames(lifetable_1)=="1"]="event"
  colnames(lifetable_1)[colnames(lifetable_1)=="0"]="censor"
  colnames(lifetable_1)[colnames(lifetable_1)=="2"]="cmprk" 
  colnames(lifetable_1)[colnames(lifetable_1)=="3"]="death"            ############################################## CHANGED                                                       
  lifetable_1
  
  ### need to consider the situation that there is no censor 0, or event 1, competing risk 2 in the data at all.
  cause_in=unique(freq_cause$cause)  
  if(0 %in% cause_in==FALSE){
    censor=rep(0,dim(lifetable_1)[1])
    lifetable_1=data.frame(lifetable_1,censor)
  }
  if(1 %in% cause_in==FALSE){
    event=rep(0,dim(lifetable_1)[1])
    lifetable_1=data.frame(lifetable_1,event)
  }
  if(2 %in% cause_in==FALSE){
    cmprk=rep(0,dim(lifetable_1)[1])
    lifetable_1=data.frame(lifetable_1,cmprk)
  }
  if(3 %in% cause_in==FALSE){
    death=rep(0,dim(lifetable_1)[1])
    lifetable_1=data.frame(lifetable_1,death)                          ############################################## CHANGED
  }
  
  
  
  
  
  ## n at risk at j = n at risk at j-1 -C-R, so get the running sum of C and R over time
  sum_censor=cumsum(lifetable_1[,"censor"])
  cmprk=lifetable_1[,"cmprk"]                                          ############################################## CHANGED
  sum_death=cumsum(lifetable_1[,"death"])                              ############################################## CHANGED
  lifetable_2=cbind(lifetable_1,sum_censor,cmprk,sum_death)            ############################################## CHANGED
  nrisk=ntotal-(sum_censor+sum_death+cmprk)                            ############################################## CHANGED
  nrisk_previous=c(ntotal,nrisk[1:(length(nrisk)-1)]) ## at the first time point, n at risk is the original number
  
  lifetable=data.frame(time=lifetable_1$time,nrisk=nrisk_previous,lifetable_1[,c("censor","event","cmprk","death")])    ############################################## CHANGED
  lifetable
  
  surv_prob=1-lifetable$death/lifetable$nrisk                          ############################################## CHANGED               
  overall_surv=cumprod(surv_prob)
  overall_surv_previous=c(1,overall_surv[1:(length(overall_surv)-1)])   ###KM(Tj-1) is used in the MCC equation.
  Ave_events=overall_surv_previous*lifetable$event/lifetable$nrisk
  MCC=cumsum(Ave_events)
  
  MCCtable=data.frame(lifetable,MCC)
  MCC.final=do.call(rbind, lapply(split(MCCtable, MCCtable$MCC), head, 1))[,c("time","MCC")]
  rownames(MCC.final)=NULL
  return(MCC.final)
}
#-------------------------------------------------------------------------------------------------#


# # ------------- Run the function with examples ------------#
# 
# ### ----- Simpel example with 3 ids, and everyone had an event ----one row per person, this equals cumulative incidence
# Participant_ID=c(1,2,3)
# Exit_Time=c(8,9,10)
# Cause=c(1,1,1)
# simple1=data.frame(Participant_ID, Exit_Time, Cause)
# simple1
# MCC(id=Participant_ID, time=Exit_Time,cause=Cause,type="MCC")
# MCC(id=Participant_ID, time=Exit_Time,cause=Cause,type="SCI")
# 
# ### -----Simpel example with 3 ids, and everyone had an event and then immedicatly censored
# Participant_ID=c(1,1,2,2,3,3)
# Exit_Time=c(8,8,9,9,10,10)
# Cause=c(1,0,1,0,1,0)
# simple2=data.frame(Participant_ID, Exit_Time, Cause)
# simple2
# 
# MCC(id=Participant_ID, time=Exit_Time,cause=Cause,type="MCC")
# MCC(id=Participant_ID, time=Exit_Time,cause=Cause,type="SCI")
# 
# ### ---- Example data used in the Paper ---- ###
# Participant_ID=c(1,2,3,4,4,4,4,5,5)
# Exit_Time=c(15.2,1.3,7.4,3.7,9.6,11.5,15.2,4.5,5.3)
# Cause=c(0,0,2,1,1,1,0,1,2)
# MCC(id=Participant_ID, time=Exit_Time,cause=Cause,type="MCC")
# MCC(id=Participant_ID, time=Exit_Time,cause=Cause,type="SCI")
# 
# 
# ### ---- Example data: multiple people had events at 4.2 ---- ###
# id=c(1,2,2,2,2,3,4,4,4,5,5,6,6,6,6,7,7,7,7)
# time=c(0.4,1.3,3,4.8,5.3,1.5,1,3,4.8,4.2,5.3,2.6,4.2,5.4,5.8,1.3,2,4.2,5.7)
# cause=c(0,1,1,1,0,2,1,1,2,1,0,1,1,1,1,1,1,1,2)
# data.frame(id, time, cause)
# MCC(id=id, time=time,cause=cause,type="MCC")
# MCC(id=id, time=time,cause=cause,type="SCI")
# 
# 
# 
# 
# 
# 
# 
# 
