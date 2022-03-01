#CHECK THAT CALS ARE CORRECT FOR 2 unknown related (siblings)
#rm(list=ls())
condRelIndex = 1:2
refRlist = list(ref1=c("26","19"),ref2=c("22","24"))
ibd0=c(0.25,0.50,0.25)
ibdList = list(ibd0,ibd0)
freq=c(0.108794,0.0308214,0.0282059,0.10632,0.104411,0.103421,0.0200057,0.498021)#paste0(freq,collapse=",")
names(freq) = c("19","21","22","23","24","25","26","99")

ret = calcGjoint2(freq,nU=2,refK=NULL,refRlist=refRlist,ibdList=ibdList,condRelIndex=condRelIndex)
P = ret$Gprob

ret1a = calcGjoint2(freq,nU=1,refRlist=refRlist[1],ibdList=ibdList[1],condRelIndex=1)
ret1b = calcGjoint2(freq,nU=1,refRlist=refRlist[2],ibdList=ibdList[2],condRelIndex=1)
if(0) { #check that EFM gives same probs
  ret2a = euroformix::calcGjoint(freq,nU=1,refR=refRlist[[1]],ibd=ibd[[1]])
  ret2b = euroformix::calcGjoint(freq,nU=1,refR=refRlist[[2]],ibd=ibd[[2]])
  all(ret1a$Gprob==ret2a$Gprob)
  all(ret1b$Gprob==ret2b$Gprob)
}

test_that("check genProb of 2 unknown siblings", {
  nG = nrow(ret1a$G)
  for(i in 1:nG) {
    for(j in 1:nG) {
      expect_equal(ret1a$Gprob[i]*ret1b$Gprob[j],P[i,j])
    }
  }
})
