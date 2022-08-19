#Testing that the numerical calculation of loglik is correct
#Special test: Second replicate include an extra marker that is not used.
#CAUSES ERROR
#rm(list=ls());library(euroformix);library(testthat);library(EFMrep)

pkgpath = path.package("EFMrep")

#Other Settings:
pC = 0.05
lambda = 0.01
fst = 0 #marker specific
NOC <<- 1 #number of contributors

#Import data:
setwd(paste0(pkgpath,"//examples"))
popfn <- paste0("task1_freqs.csv")
popFreq <<- freqImport(popfn)[[1]] #must contain ALL loci to analyse (over all replicates)
samples <<- sample_tableToList(tableReader("test1_evids_extended.csv")) #notice that we evaluate both reps here
refData <<-  sample_tableToList(tableReader("test1_ref1_missing.csv"))
nReps <<- length(samples)

kit <<- c(rep("testkit",nReps))
fstv <<- setNames( rep(fst,length(popFreq)), names(popFreq) )
ATv <<- setNames(rep(50,nReps),names(samples)) #specify AT for each replicate
lambdav <<- setNames( rep(lambda,nReps), names(samples))
pCv <<- setNames( rep(pC,nReps), names(samples))

#PARAM: Unequal params#
PHexp <<- 1:nReps
PHvar <<- 1:nReps
DEG  <<- 1:nReps  #unequal per replicate
stuttBW <<- 1:nReps #unequal across replicates
stuttFW <<- 1:nReps #unequal across replicates
mixProp <<- 1:nReps 

#Hypothesis 1: Cond on ref1 
test_that("check hyp1", {
  expect_error( testFuns(condOrder=1) ) #should throw an error
})
