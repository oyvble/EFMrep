#rm(list=ls())
fst=0.1 #use high value
freq = c(0.1,0.2,0.7 )
names(freq) = seq_along(freq)
refK =  c(1,2,3,3) #typed
refR = c(1,2) #ref alleles
dec=4 #round off error

nU = 2 #2 unknowns
refRlist = list(refR)

############Siblings

ibd0=c(0.25,0.5,0.25)
ibdList = list(ibd0)

#Unknown first then related:
condRelIndex = 2
ret = calcGjoint2(freq,nU,fst,refK,refRlist,ibdList,condRelIndex)
P1 = ret$Gprob

#Related first then unknown:
condRelIndex = 1
ret = calcGjoint2(freq,nU,fst,refK,refRlist,ibdList,condRelIndex)
P2 = ret$Gprob
P2t = t(P1)

P1-P2t

############ParentCHild

ibd0=c(0,1,0)
ibdList = list(ibd0)

#Unknown first then related:
condRelIndex = 2
ret = calcGjoint2(freq,nU,fst,refK,refRlist,ibdList,condRelIndex)
P1 = ret$Gprob

#Related first then unknown:
condRelIndex = 1
ret = calcGjoint2(freq,nU,fst,refK,refRlist,ibdList,condRelIndex)
P2 = ret$Gprob
P2t = t(P1)

P1-P2t



