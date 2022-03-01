#Testing that the numerical calculation of loglik is correct
#rm(list=ls());library(euroformix);library(testthat)
pkgpath = path.package("EFMrep")

#Settings:
pC = 0.05
lambda = 0.01
fst = 0.01 #marker specific
NOC <<- 2 #number of contributors

#Import data:
setwd(paste0(pkgpath,"//examples"))
popfn <- paste0("test2_freqs.csv")
popFreq <<- freqImport(popfn)[[1]] #must contain ALL loci to analyse (over all replicates)
samples <<- sample_tableToList(tableReader("test2_evids.csv"))[1]
refData <<-  sample_tableToList(tableReader("test2_refs.csv"))
nReps <<- length(samples)

kit <<- c(rep("SGMPlus",1))
fstv <<- setNames( rep(fst,length(popFreq)), names(popFreq) )
ATv <<- setNames(c(200),names(samples)) #specify AT for each replicate
lambdav <<- setNames( rep(lambda,nReps), names(samples))
pCv <<- setNames( rep(pC,nReps), names(samples))

#PARAM: Unequal params#
PHexp <<- 1:nReps
PHvar <<- 1:nReps
DEG  <<- 1:nReps  #unequal per replicate
stuttBW <<- 1:nReps #unequal across replicates
stuttFW <<- 1:nReps #unequal across replicates
mixProp <<- 1:nReps 

#Hypothesis 1: Cond on both ref1 and ref2
test_that("check hyp1", {
  testFuns(condOrder=1:2)
})

#Hypothesis 2: Cond on ref1
test_that("check hyp2", {
  testFuns(condOrder=c(1,0),knownRef=2)
})

#Hypothesis 3: Cond on ref2 (missing marker at TH01)
test_that("check hyp3", {
  testFuns(condOrder=c(0,1),knownRef=1)
})

#Hypothesis 4: 2 Unknowns
test_that("check hyp4", {
  testFuns(condOrder=c(0,0),knownRef=1:2)
})


