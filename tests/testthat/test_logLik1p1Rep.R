#Testing that the numerical calculation of loglik is correct
#rm(list=ls());library(euroformix);library(testthat)

pkgpath = path.package("EFMrep")

#Other Settings:
pC = 0.05
lambda = 0.01
fst = 0 #marker specific
NOC <<- 1 #number of contributors

#Import data:
setwd(paste0(pkgpath,"//examples"))
popfn <- paste0("test2_freqs.csv")
popFreq <<- freqImport(popfn)[[1]] #must contain ALL loci to analyse (over all replicates)
samples <<- sample_tableToList(tableReader("test1_evids.csv"))[1] #only 1 rep
refData <<-  sample_tableToList(tableReader("test1_ref1_missing.csv"))
nReps <<- length(samples)

kit <<- rep("testkit",nReps)
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
  testFuns(condOrder=1)
})
#sapply(ret,function(x) paste0("c(",paste0(round(x,4),collapse = ","),")"))

#Hypothesis 2: Unknown unrelated
test_that("check hyp2", {
  testFuns(condOrder=0,knownRef=1)
})

#Hypothesis 3: Related sibling
test_that("check hyp3", {
  ibd = list(c(1/4,1/2,1/4))
  testFuns(condOrder=0,knownRel = 1,ibd=ibd)
})

#Hypothesis 4: Related child
test_that("check hyp4", {
  ibd = list(c(0,1,0))
  testFuns(condOrder=0,knownRel = 1,ibd=ibd)
})


