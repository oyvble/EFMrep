rm(list=ls())
#EXAMPLE FILE TO HELP DEVELOPMENT
library(euroformix)
library(EFMrep)

#pkgpath = path.package("euroformix")
pkgpath = path.package("EFMrep")

#Other Settings:
AT = 50 #define 
pC = 0.05
lambda = 0.01
fst = 0

#Hypothesis settings:
nC = 2 #number of contributors
condOrder_hp=c(0,1)
condOrder_hd=c(0,0)
knownRef_hp = NULL
knownRef_hd = 2

#Import data:
setwd(paste0(pkgpath,"//examples"))

popfn = paste0("task1_freqs.csv")
popFreq = freqImport(popfn)[[1]] #must contain ALL loci to analyse (over all replicates)

#Replicates
samples =  sample_tableToList(tableReader("task1_evidESX17.csv"))

#getKit()
kit = c("ESX17")#,"Fusion 6C","GlobalFiler") #selected kit for each replicate (order must match!)

refData =  sample_tableToList(tableReader("task1_ref1.csv"))
refData =  append(refData, sample_tableToList(tableReader("task1_ref2.csv")) )
#SPECIFY MODEL OPTIONS:

#UNCOMMON PARAMS:
nReps = length(samples)
  
#PH models:
PHexp <- 1:nReps
PHvar <- 1:nReps

#Degradation setting:
DEG  = 1:nReps  #unequal per replicate

#Mixture proportions:
mixProp = 1:nReps 

#Stutter model settings
stuttBW = 1:nReps #unequal across replicates
stuttFW = 1:nReps #unequal across replicates

#condOrder=condOrder_hp;knownRef=knownRef_hp
#minF=NULL;normalize=FALSE;maxIter=1000; steptol=1e-6;nDone=3;maxThreads=128;delta=1;seed=NULL;verbose=TRUE;knownRel=NULL;ibd=NULL
steptol0 = 1e-6
mleHp = contLikMLE2(nC,samples,popFreq,refData, condOrder_hp, knownRef_hp, kit, AT,pC,lambda,fst, mixProp, PHexp, PHvar, DEG, stuttBW, stuttFW,steptol = steptol0)
mleHp$fit$loglik

#Compare with EFM:
dat = prepareData(samples,refData,popFreq[names(samples[[1]])])
mleHpEFM = contLikMLE(nC,dat$samples,dat$popFreq,dat$refData , condOrder_hp,knownRef_hp, xi = NULL,pC,nDone=3, AT,  fst, lambda, xiFW = NULL, kit=kit,steptol = steptol0)
mleHpEFM$fit$loglik

#Obtain validation plots
validhpEFM = validMLEmodel(mleHpEFM)
validhp = validMLEmodel2(mleHp)
range(validhpEFM$pvalue-validhp$pvalue)

