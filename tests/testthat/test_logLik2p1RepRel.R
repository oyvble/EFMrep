#Testing that the numerical calculation of loglik is correct
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

kit <<- "SGMPlus" #c("ESX17","Fusion 6C","GlobalFiler") #selected kit for each replicate (order must match!)
fstv <<- setNames( rep(fst,length(popFreq)), names(popFreq) )
ATv <<- setNames(c(200),names(samples)) #specify AT for each replicate
lambdav <<- setNames( rep(lambda,nReps), names(samples))
pCv <<- setNames( rep(pC,nReps), names(samples))

#PARAM: Unequal params#
PHexp <<- 1:nReps
PHvar <<- 1:nReps
DEG  <<- 1:nReps  #unequal per replicate
mixProp <<- 1:nReps 
stuttBW <<- 1:nReps #unequal across replicates
stuttFW <<- 1:nReps #unequal across replicates

#Hypothesis 1: Cond on ref2 (missing marker on TH01), and 1st Unknown is related (sibling) to ref1
ibd0 = c(1/4,1/2,1/4) #sibling koefficients
test_that("check hyp1", {
# condOrder=c(0,1);knownRel=1;ibd=list(ibd0);knownRef=NULL;
  testFuns(condOrder=c(0,1), knownRel=1, ibd=list(ibd0))
})

#Hypothesis 2: 1st Unknown is related (sibling) to ref1, 2nd unknown is related (sibling) to ref2
test_that("check hyp2", {
# condOrder=c(0,0);knownRel=1:2;ibd=list(ibd0,ibd0);knownRef=NULL;
  testFuns(condOrder=c(0,0),knownRel=1:2,ibd=list(ibd0,ibd0))
})


