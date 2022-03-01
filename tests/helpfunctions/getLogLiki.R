#MANUAL DERIVATION FOR EACH GENOTYPE OUTCOME
getLogLiki = function(samples,refData,popFreq,condOrder,knownRef,par,NOC,ATv,pCv,lambdav,fstv, knownRel=NULL, ibd=NULL) {
  #Prepare data:
  Qallele = "99"
  dat = prepareData2(samples,refData,popFreq,AT=ATv) #obtain data to use for analysis
  repNames=names(samples) #obtain replicate names (this is the index order for model params)
  
  #COMPARE PER MARKER RESULTS WITH MANUAL DERIVED:
  #expect_equal(sum(logLikv),mle$fit$loglik)
  kitinfoList = lapply(unique(kit),euroformix::getKit)
  names(kitinfoList) = unique(kit)
  
  locs = names(dat) #obtain loci to traverse
  #logLiki = setNames( rep(NA,length(locs)), locs) #overall likelihood
  LikList = list() #store vals for each genotype outcome
  for(loc in locs) {
#  loc=locs[5]
    datloc = dat[[loc]]
    freq = datloc$freq 
    Aset = names(freq) #get allele outcome 
    Aset0 = setdiff(Aset,Qallele) #keep non-dropouts
    Astutt = c(as.numeric(Aset0)-1,as.numeric(Aset0)+1)
    Aset = c(Aset,setdiff(Astutt,Aset0)) #insert potential stutters last
    
    #Obtain number of conditionals and unkowns
    condRef =  datloc$refs[which(condOrder>0)] #get alleles of conds
    if(length(condRef)==0) {
      nCond = 0
    } else {
      nCond = sum(sapply(condRef,function(x) !is.null(x) && length(x)==2)) #number of contributors
    }
    nU = NOC - nCond
    
    #Replicates
    repNamesMarker = names(dat[[loc]]$samples)
    repIDs = which(repNames%in%repNamesMarker)
    
    #Prepare PH data (same order as Aset):
    Ylist = list() #list per samples
    for(sample in repNamesMarker) {
      Ylist[[sample]] = rep(0,length(Aset))  #preparing PHs after extending allele vector
      Ylist[[sample]][match( names(dat[[loc]]$samples[[sample]]) ,Aset)] = dat[[loc]]$samples[[sample]]  #insert PHs
    }
    
    #Obtain contribution matrix
    nAG = matrix(0,nrow=length(Aset),ncol=NOC,dimnames = list(Aset)) #create allele counting matrix for unknown contributors
    for(ind in seq_len(nCond)) nAG[,ind] = table(factor(condRef[[ind]],level=Aset)) #insert contribution
    
    nAG0 = nAG #temporary store
    if(nU>0) {
      refK = unlist(datloc$refs[c(which(condOrder>0),knownRef)]) #obtain references
      refRlist = NULL
      if(!is.null(knownRel)) refRlist = datloc$refs[knownRel] #obtain references to consider
      
      #obtain contributution index of related (must be given if non-regular references): 
      #IMPORTANT: RELATEDS ARE ALWAYS THE FIRST UNKNOWNS UNLESS THERE ARE KNOWN CONTRIBUTORS MISSING MARKER(S) 
      condRelIndex =  sum(condOrder>0) + seq_len(length(knownRel))  - nCond #subtract with number of found conditional
      
      Glist = calcGjoint2(freq=freq,nU=nU,fst=fstv[loc],refK=refK,refRlist=refRlist,ibdList=ibd,condRelIndex=condRelIndex  )
#      Glist2 = euroformix::calcGjoint(freq,nU,fstv[loc],refK,refRlist[[1]],ibd[[1]])
#      max(abs(c(Glist$Gprob)-c(Glist2$Gprob)))
      Gset = Glist$G #get allele out come of unknowns
      Gprob = c(Glist$Gprob ) #vectorize 
      numGenos1p = nrow(Gset)
      nJointComb = numGenos1p^nU
    } else {
      nJointComb = 1
      Gprob = 1
    }
    
    #TRAVERSING ALL COMBINATION OF unknown genotypes
    #Lik = 0 #obtain probabiliity of evidence (sum over all genotypes)
    LikList[[loc]] = rep(NA,nJointComb) #init vector to store vals
    for(gind in seq_len(nJointComb)) { #traverse all genorypes
#      gind=1
      nAG = nAG0 #get copy
      
      gind2 = gind - 1 #get C index 
      modrest = gind2
      inserted = FALSE
      jointGind = rep(0,nU)
      
      for(kk in seq_len(nU)) {
        if (!inserted) { 
          if ( kk>1 ) {
            modrest = as.integer( (modrest - jointGind[kk-1]) / numGenos1p); 
          }
          jointGind[kk] = modrest %% numGenos1p; 
          if (modrest < numGenos1p) inserted = TRUE  #then all digits are inserted
        }
        nAG[,nCond + kk] = table(factor(Gset[ jointGind[kk]+1,],level=Aset)) #insert contribution
      }

      logEvid = 0 #to be summed up (for marker)
      for(rr in repIDs) { #traverse each replicataes (outer loop)
#     rr=1
        sample = repNames[rr]
        kitinfo = kitinfoList[[kit[rr]]]
        if(!repNames[rr]%in%names(Ylist)) next #skip if no data  
        
        #Prepare base pair info for kits:
        bpv = rep(NA,length(Aset)) #obtain base pairs for kit
        subkit <- subset(kitinfo,toupper(kitinfo$Marker)==loc) #control on dyer
        for(aa in Aset) { #for each allele
          ind <- which(subkit$Allele==aa)
          if(length(ind)==0) ind <- which.min(abs(as.numeric(subkit$Allele) - as.numeric(aa))) #nearest neighbour (in allele name)-> maximum of the alleles
          bpv[aa==Aset] <- subkit$Size[ind] #insert base pair
        }
        
        #PARAM:
        mx = par$mixProp[NOC*(rr-1) + 1:NOC]
        shape0 = 1/(par$PHvar[rr]^2) #get shape param
        scale0 = par$PHexp[rr]/shape0 #get scale param
        beta =  par$DEG[rr] #degrad slope param
        xiB = par$stutt[2*rr - 1] #backwards stutter prop
        xiF = par$stutt[2*rr ] #forward stutter prop
        
        mui <-  c(nAG[,drop=FALSE]%*%mx)*shape0 #expected contribution
        mui  <- mui*exp( (bpv-125)/100*log(beta) )  #in case of degradation
        
        
        if( xiB>0 || xiF>0) {
          #Delegate stutterprod
          FWstutt = match(Aset,as.character(as.numeric(Aset)+1)) #alleleinds to stutter from
          BWstutt = match(Aset,as.character(as.numeric(Aset)-1)) #alleleinds to stutter from
          stuttB <- mui[BWstutt]*xiB #backward-stutter parts
          stuttF <- mui[FWstutt]*xiF #forward-stutter parts
          indBW = !is.na(stuttB)
          indFW = !is.na(stuttF)
          indLooseStutt = indBW & indFW #index of alleles which are assumed to not loose stutter product
          mui[indLooseStutt] = mui[indLooseStutt]*(1- (xiB+xiF)) #loose stutter products
          mui[indBW] = mui[indBW] + stuttB[indBW]
          mui[indFW] = mui[indFW] + stuttF[indFW]
          #      names(mui) = Aset
        }
        
        vali = 0
        #Divide set into dropin, contr and dropout
        psiDI <- which( Ylist[[sample]]>0 & mui==0 ) #drop-in (not a stutter). ANY DROPIN?
        psiYmu <- which( Ylist[[sample]]>0 & mui>0 )  #contributing to model and observed PH
        psiDO <- which( Ylist[[sample]]==0 & mui>0 )  #dropout elem
        
        if(length(psiYmu)>0) vali =  vali + sum(dgamma(Ylist[[sample]][psiYmu],shape=mui[psiYmu],scale=scale0,log=TRUE))   # CONTRIBUTION TO PH
        if(length(psiDO)>0) vali = vali + sum(pgamma(ATv[sample],shape=mui[psiDO],scale=scale0,log=TRUE)) #CONTR TO DROPOUT
        if(length(psiDI)>0) vali = vali + sum( dexp( Ylist[[sample]][psiDI]-ATv[sample], rate=lambdav[sample], log=TRUE) + log(pCv[sample]*freq[psiDI]) ) #CONTR TO DROPIN: Scaled for each dropin. THIS IS OK
        if(length(psiDI)==0) vali = vali + log(1-pCv[sample]) #in case of no dropin
        logEvid = logEvid + vali #add to likelihood of marker
      } #end for each replicate
     # print(logEvid)
      #Lik = Lik + exp(logEvid)*Gprob[gind] #SUM UP OVER ALL GENOTYPE
      LikList[[loc]][gind] = exp(logEvid)*Gprob[gind] #insert for genotype comb
    } #end for each genotype outcome
    #logLiki[loc] = log(Lik) #insert
  } #end for each marker
  return(LikList)
} #end function
