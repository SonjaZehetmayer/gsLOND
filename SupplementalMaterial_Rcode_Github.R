##############################################################################################################################
# Last modificiation: September 2022                                                                                         #
#                                                                                                                            #
# This is a source code in the R programming language (http://cran.r-project.org/) for the manuscript                        #
# "Online control of the False Discovery Rate in group-sequential platform trials"                                           #
# by Sonja Zehetmayer, Martin Posch, Franz Koenig (Statistical Methods in Medical Research, To Appear).                      #   
#                                                                                                                            #
# Implements group-sequential LOND procedure (gsLOND) for online FDR control.                                                #
#                                                                                                                            #
# An input vector of p-values is given together with the specification of hypothesis, stage and time (for details see below).#
#                                                                                                                            #
#                                                                                                                            #
#                                                                                                                            #
##############################################################################################################################
#                                                                                                                            #
# Mainfunctions:                                                                                                             #
# gsLOND                                                                                                                     #
# gsLOND.II                                                                                                                  #
# gsLOND.III                                                                                                                 #
#                                                                                                                            #
# Parameters of gsLOND functions:                                                                                            #
# ----------------------------                                                                                               #
# pvals: vector of p-values of interim and final analyses                                                                    #
# hyp: vector defining which p-values belong to which hypothesis                                                             #
# stage: vector defining if a p-values was calculated for interim or for final data.                                         #
# time: timing of test decisions                                                                                             #
# alpha: efficacy boundary                                                                                                   #
# alpha1: futility boundary                                                                                                  #
# N: upper bound for total number of hypotheses/treatments. Either an integer or NA (no upper bound).                        #
#                     for beta="equal", N must be defined as integer                                                         #
# iuse.val: type of alpha spending: 1...O'Brien-Fleming type spending function                                               #
#                                   2...Pocock type spending function                                                        #    
# beta: type of allocation: "desc": descending beta values according to Javanmard and Montanari (2018), equation 31          # 
#                           "equal": equal allocation                                                                        #
#                                                                                                                            #
# Output of gsLOND functions:                                                                                                #
# ------------------------                                                                                                   #  
# list of input elements and results:                                                                                        #
#                                                                                                                            #
# inputs: matrix of input parameters pvals, stage, hyp, time                                                                 #
#                                                                                                                            #
# results: matrix of output parameters R,stopstage,betai,gsalpha.1,gsalpha.2 for each hypothesis:                            #
#    R: 0 (if hypothesis is not rejected) or 1 (if hypothesis is rejected)                                                   #
#    stopstage: 1 if hypothesis stops in the interim analysis either for futility or efficacy or                             #
#               2 if final test decision is performed in the final analysis.                                                 #
#    betai: allocation of significance level alpha                                                                           #
#    alphai: individual significance thresholds for each hypothesis                                                          #
#    gsalpha.1: group-sequential boundaries for test decision in interim analysis                                            #
#    gsalpha.2: group-sequential boundaries for test decision in final analysis; 0 if hypothesis stopped interim analysis    #
#                                                                                                                            #
#                                                                                                                            #
#                                                                                                                            #
##############################################################################################################################
# No warranty is given about the validity of this program!                                                                   #
##############################################################################################################################



library(ldbounds)
library(rpact)

gsLOND<-function(pvals,hyp,stage,time,alpha=0.025,alpha1=0.5,N=1000,iuse.val=1,beta="desc")
{
  K<-max(hyp)

  mat<-cbind(pvals,stage,hyp,time)              
  ordering<-mat[order(mat[,4]),]                
  colnames(ordering)<-c("pvals","stage","hyp","time")
  
  if (beta=="desc")
  {
    jj<-1:max(N,K,na.rm=TRUE)
    betai <-  0.07720838 * alpha * log(pmax(jj, 2))/(jj * exp(sqrt(log(jj))))
    if (!is.na(N))  betai<-betai/sum(betai) *alpha    #values of beta, if N is defined
  }
  
  if (beta=="equal")        
    betai<-rep(alpha/N,N)
  
  
  betai<-betai[1:K]
  R<-rep(0,K)
  stopstage<-rep(0,K)
  gsalpha.1<-rep(0,K)
  gsalpha.2<-rep(0,K)
  alphai<-rep(0,K)
  
  if (iuse.val==1) meth="asOF"
  if (iuse.val==2) meth="asP"
  
  for (i in seq_len(nrow(ordering) )) 
  {
    
    if (ordering[i,2]==1)   #stage 1
    {
      alphai[ordering[i,3]]<-betai[ordering[i,3]]*(1+sum(R[1:(ordering[i,3]-1)]))
      #gsb<-ldBounds(c(0.5,1),iuse=iuse.val,alpha=alphai[ordering[i,3]],sides=1)
      #gsalpha.1[ordering[i,3]]<-1-pnorm(gsb$upper.bounds)[1] #interim level
      
      gsb<-getDesignGroupSequential(typeOfDesign = meth, informationRates = c(0.5, 1),alpha = alphai[ordering[i,3]], sided = 1, tolerance = 1e-08)[["stageLevels"]][1]
      
      gsalpha.1[ordering[i,3]]<-gsb #interim level
      R[ordering[i,3]]<-ifelse(ordering[i,1]<=gsalpha.1[ordering[i,3]],1,0)
      stopstage[ordering[i,3]]<-ifelse((ordering[i,1]<=gsalpha.1[ordering[i,3]])|(ordering[i,1]>=alpha1),1,0)
    }
    
    if (ordering[i,2]==2)   #stage 2
    {
      if (stopstage[ordering[i,3]]==0) #the hypothesis has not stopped in the interim analysis
      {
        alphai[ordering[i,3]]<-betai[ordering[i,3]]*(1+sum(R[1:(ordering[i,3]-1)]))
        gsb<-getDesignGroupSequential(typeOfDesign = meth, informationRates = c(0.5, 1),alpha = alphai[ordering[i,3]], sided = 1, tolerance = 1e-08)[["stageLevels"]][2]
        
        gsalpha.2[ordering[i,3]]<-gsb
        R[ordering[i,3]]<-ifelse(ordering[i,1]<=gsalpha.2[ordering[i,3]],1,0)
        stopstage[ordering[i,3]]<-2
      }
    }
  }
  
  result<-list(ordering,cbind(1:K,R,stopstage,betai,alphai,gsalpha.1,gsalpha.2))
  names(result)<-c("input","result")
  result
}



gsLOND.II<-function(pvals,hyp,stage,time,alpha=0.025,alpha1=0.5,N=1000,iuse.val=1,beta="desc")
{
  K<-max(hyp)
  
  mat<-cbind(pvals,stage,hyp,time)  
  ordering<-mat[order(mat[,4]),]    
  colnames(ordering)<-c("pvals","stage","hyp","time")
  
  if (iuse.val==1) meth="asOF"
  if (iuse.val==2) meth="asP"
  
  if (beta=="desc")
  {
    jj<-1:max(N,K,na.rm=TRUE)
    betai <-  0.07720838 * alpha * log(pmax(jj, 2))/(jj * exp(sqrt(log(jj))))
    if (!is.na(N))  betai<-betai/sum(betai) *alpha    #values of beta, if N is defined
  }
  
  if (beta=="equal")        
    betai<-rep(alpha/N,N)

  betai<-betai[1:K]
  R<-rep(0,K)
  stopstage<-rep(0,K)
  gsalpha.1<-rep(0,K)
  gsalpha.2<-rep(0,K)
  alphai<-rep(0,K)
  
  for (i in seq_len(nrow(ordering) )) 
  {
    
    if (ordering[i,2]==1)   #stage 1
    {
      alphai[ordering[i,3]]<-betai[ordering[i,3]]*(1+sum(R[1:(ordering[i,3]-1)]))
      gsb<-getDesignGroupSequential(typeOfDesign = meth, informationRates = c(0.5, 1),alpha = alphai[ordering[i,3]], sided = 1, tolerance = 1e-08)[["stageLevels"]][1]
      
      gsalpha.1[ordering[i,3]]<-gsb  #interim level
      R[ordering[i,3]]<-ifelse(ordering[i,1]<=gsalpha.1[ordering[i,3]],1,0)
      stopstage[ordering[i,3]]<-ifelse((ordering[i,1]<=gsalpha.1[ordering[i,3]])|(ordering[i,1]>=alpha1),1,0)
    }
    
    if (ordering[i,2]==2)   #stage 2
    {
      if (stopstage[ordering[i,3]]==0) #the hypothesis has not stopped in the interim analysis
      {
        alphai[ordering[i,3]]<-betai[ordering[i,3]]*(1+sum(R[1:(ordering[i,3]-1)]))

        gsalpha.2[ordering[i,3]]<-getDesignGroupSequential(typeOfDesign = "asUser", informationRates = c(0.5, 1),alpha = alphai[ordering[i,3]], 
                                 sided = 1, tolerance = 1e-08,userAlphaSpending = c(gsalpha.1[ordering[i,3]],alphai[ordering[i,3]]),
                                 typeBetaSpending = "none")[["stageLevels"]][[2]]
        
        R[ordering[i,3]]<-ifelse(ordering[i,1]<=gsalpha.2[ordering[i,3]],1,0)  
        stopstage[ordering[i,3]]<-2
      }
    }
  }
  
  result<-list(ordering,cbind(1:K,R,stopstage,betai,alphai,gsalpha.1,gsalpha.2))
  names(result)<-c("input","result")
  result
}




gsLOND.III<-function(pvals,hyp,stage,time,alpha,alpha1,N,iuse.val,beta="desc")
{
  K<-max(hyp)
  
  mat<-cbind(pvals,stage,hyp,time)                                      
  ordering<-mat[order(mat[,4]),]                
  colnames(ordering)<-c("pvals","stage","hyp","time")
  
  if (iuse.val==1) meth="asOF"
  if (iuse.val==2) meth="asP"
  
  if (beta=="desc")
  {
    jj<-1:max(N,K,na.rm=TRUE)
    betai <-  0.07720838 * alpha * log(pmax(jj, 2))/(jj * exp(sqrt(log(jj))))
    if (!is.na(N))  betai<-betai/sum(betai) *alpha    #values of beta, if N is defined
  }
  
  if (beta=="equal")        
    betai<-rep(alpha/N,N)

  betai<-betai[1:K]
  R<-rep(0,K)
  stopstage<-rep(0,K)
  gsalpha.1<-rep(0,K)
  gsalpha.2<-rep(0,K)
  alphai<-rep(0,K)
  
  for (i in seq_len(nrow(ordering) )) 
  {
    if (ordering[i,2]==1)   #stage 1
    {
      alphai[ordering[i,3]]<-betai[ordering[i,3]]*(1+sum(R))
      gsb<-getDesignGroupSequential(typeOfDesign = meth, informationRates = c(0.5, 1),alpha = alphai[ordering[i,3]], sided = 1, tolerance = 1e-08)[["stageLevels"]][1]
      
      gsalpha.1[ordering[i,3]]<-gsb#1-pnorm(gsb$upper.bounds)[1] #interim niveau
      R[ordering[i,3]]<-ifelse(ordering[i,1]<=gsalpha.1[ordering[i,3]],1,0)
      stopstage[ordering[i,3]]<-ifelse((ordering[i,1]<=gsalpha.1[ordering[i,3]])|(ordering[i,1]>=alpha1),1,0)
    }
    
    if (ordering[i,2]==2)   #stage2
    {
      if (stopstage[ordering[i,3]]==0) #the hypothesis has not stopped in the interim analysis
      {
        alphai[ordering[i,3]]<-betai[ordering[i,3]]*(1+sum(R))
        gsb<-getDesignGroupSequential(typeOfDesign = meth, informationRates = c(0.5, 1),alpha = alphai[ordering[i,3]], sided = 1, tolerance = 1e-08)[["stageLevels"]][2]
        
        gsalpha.2[ordering[i,3]]<-gsb 
        R[ordering[i,3]]<-ifelse(ordering[i,1]<=gsalpha.2[ordering[i,3]],1,0)  
        stopstage[ordering[i,3]]<-2
      }
    }
  }
  
  result<-list(ordering,cbind(1:K,R,stopstage,betai,alphai,gsalpha.1,gsalpha.2))
  names(result)<-c("input","result")
  result
}



#############################################################################################################
#                                                                                                           #
#                                                Example                                                    #
#                                                                                                           #
#############################################################################################################

#platform trial with currently 2 hypotheses, H1 with interim and final analysis, H2 only with interim analysis

pvals<-c(0.01,0.0001,0.00002)
hyp<-c(1,1,2)
stage<-c(1,2,1)
time<-c(1,3,2)

#O'Brien-Fleming design (rejects H1 at final analysis, but not H2 in interim analysis), no upper bound N defined
gsLOND(pvals,hyp,stage,time,alpha=0.025,alpha1=0.5,N=NA,iuse.val=1)
gsLOND.II(pvals,hyp,stage,time,alpha=0.025,alpha1=0.5,N=NA,iuse.val=1)
gsLOND.III(pvals,hyp,stage,time,alpha=0.025,alpha1=0.5,N=NA,iuse.val=1)

#Pocock design with upper bound N=100 on the total number of hypotheses
#(rejects H1 at final analysis and H2 in interim analysis)
gsLOND(pvals,hyp,stage,time,alpha=0.025,alpha1=0.5,N=100,iuse.val=2)
gsLOND.II(pvals,hyp,stage,time,alpha=0.025,alpha1=0.5,N=100,iuse.val=2)
gsLOND.III(pvals,hyp,stage,time,alpha=0.025,alpha1=0.5,N=100,iuse.val=2)






#values from toy example, table 1, Pockock design
pvals<-c(0.1,0.000001,0.2,0.000001,0.1,0.1)
hyp<-c(1,1,2,2,3,3)
stage<-c(1,2,1,2,1,2)
time<-c(1,2,3,4,5,6)

gsLOND(pvals,hyp,stage,time,alpha=0.05,alpha1=0.5,N=3,iuse.val=2,beta="equal")
gsLOND.II(pvals,hyp,stage,time,alpha=0.05,alpha1=0.5,N=3,iuse.val=2,beta="equal")
gsLOND.III(pvals,hyp,stage,time,alpha=0.05,alpha1=0.5,N=3,iuse.val=2,beta="equal")


#values from toy example, table 1, O'Brien-Fleming design
pvals<-c(0.1,0.000001,0.2,0.000001,0.1,0.1)
hyp<-c(1,1,2,2,3,3)
stage<-c(1,2,1,2,1,2)
time<-c(1,2,3,4,5,6)

gsLOND(pvals,hyp,stage,time,alpha=0.05,alpha1=0.5,N=3,iuse.val=1,beta="equal")
gsLOND.II(pvals,hyp,stage,time,alpha=0.05,alpha1=0.5,N=3,iuse.val=1,beta="equal")
gsLOND.III(pvals,hyp,stage,time,alpha=0.05,alpha1=0.5,N=3,iuse.val=1,beta="equal")


#values from toy example, where H1 is rejected in the final analysis table 2, Pockock design 
pvals<-c(0.1,0.000002,0.1,0.000001,0.1,0.1)
hyp<-c(1,1,2,2,3,3)
stage<-c(1,2,1,2,1,2)
time<-c(1,3,2,5,4,6)

gsLOND(pvals,hyp,stage,time,alpha=0.05,alpha1=0.5,N=3,iuse.val=2,beta="equal")
gsLOND.II(pvals,hyp,stage,time,alpha=0.05,alpha1=0.5,N=3,iuse.val=2,beta="equal")
gsLOND.III(pvals,hyp,stage,time,alpha=0.05,alpha1=0.5,N=3,iuse.val=2,beta="equal")

