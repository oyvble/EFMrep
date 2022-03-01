#Testing that the numerical calculation of loglik is correct (also comparing against EFM)
#rm(list=ls());library(euroformix);library(testthat);library(EFMrep)
pkgpath = path.package("EFMrep")

#Settings:
#Other Settings:
AT <<- 40 #define 
pC <<- 0.11
lambda <<- 0.054
fst <<- 0.01
NOC <<- 3

#getKit()

#Import data:
setwd(paste0(pkgpath,"//examples"))
popfn <- paste0("NIST.csv")
popFreq <<- freqImport(popfn)[[1]] #must contain ALL loci to analyse (over all replicates)
samples <<- sample_tableToList(tableReader("PROVEDIt_evids.csv"))[1]
refData <<-  sample_tableToList(tableReader("PROVEDIt_refs.csv"))
nReps <<- length(samples)

kit <<- c(rep("Fusion 6C",1))
fstv <<- setNames( rep(fst,length(popFreq)), names(popFreq) )
ATv <<- setNames(c(AT),names(samples)) #specify AT for each replicate
lambdav <<- setNames( rep(lambda,nReps), names(samples))
pCv <<- setNames( rep(pC,nReps), names(samples))

#PARAM: Unequal params#
PHexp <<- 1:nReps
PHvar <<- 1:nReps
DEG  <<- 1:nReps  #unequal per replicate
stuttBW <<- 1:nReps #unequal across replicates
stuttFW <<- 1:nReps #unequal across replicates
mixProp <<- 1:nReps 

#Hypothesis 1: Cond on all refs
test_that("check hyp1", {
  testFuns(condOrder=1:3,compareEFM = TRUE,steptol0=1e-6,nDone0=3,seed0=1)
})



