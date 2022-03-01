#' @title deconvolve2
#' @author Oyvind Bleka
#' @description deconvolve ranks the set of the most conditional posterior probability of genotypes the STR DNA mixture given a fitted model under a hypothesis.
#' @details The procedure calculates the likelihood for each single locus. Then it combines the most probable genotypes from each loci to produce a ranked list of deconvolved profiles.
#' 
#' @param mlefit Fitted object using contLikMLE function.
#' @param alpha Required sum of the listed posterior probabilities.
#' @param maxlist The ranked deconvolved profile list will not exceed this number (used to avoid endless search).
#' @return ret A list(table1,table2,table3,table4,rankGi,rankG,pG) where rankG is the ranked genotypes with corresponding probabilities in pG. rankgGi is the same, but per marker. table1 is rankG and pG combined (joint results). table2 uses rankGi to find marginal results for top-genotypes. table3 and table4 shows this marginalized on genotypes and alleles per contributor 
#' @export

# alpha=0.95;maxlist=1000
deconvolve2 = function(mlefit,alpha=0.95,maxlist=1000){ 
#  mlefit =  get("setEVID",envir=mmTK)$mlefit_hd
  
  par = mlefit$fit$par
  c <- mlefit$prepareC #returned from prepareC
  nC = c$NOC #number of contrs
  nM = c$nLocs #number of markers
  locs = c$locNames #loci to evaluate
  
  #Step 1) Calculate L(E|g,thetahat) for each marker
  dlist <- list()  #joint genotype probs for each combinations per marker
  GClist2 <- list() #get index of ranked combined genotypes for all genotypes
  Glist2 <- list() #used to store names of genotypes 
  
  #nrow(combGind)==c$nJointCombs[m] #must be the same
  loglikVEC = rep(0,sum(c$nJointCombs)) #init big calculation vector (likelihood for all outcome)
  loglikVEC = .C("loglikGamma_allcomb2", as.numeric(loglikVEC), c$nJointCombs, c$NOC, c$NOK,
            as.numeric(par$mixProp),  as.numeric(par$PHexp), as.numeric(par$PHvar), as.numeric(par$DEG), as.numeric(par$stutt), 
            as.numeric(c$AT),as.numeric(c$fst),as.numeric(c$dropinProb),as.numeric(c$dropinWeight),
            c$nReps, c$nLocs, c$nRepMarkers, c$nAlleles, c$nAlleles2, c$startIndMarker_nAlleles, c$startIndMarker_nAllelesReps,
            c$peakHeights, c$freq, c$nTyped, c$maTyped, c$basepairs, c$repID, c$startIndMarker_nRepMarkers, 
            c$nGenos, c$outG1allele, c$outG1contr, c$startIndMarker_outG1allele, c$startIndMarker_outG1contr, c$startIndMarker_nJointCombs,
            c$nStutters, c$stuttFromInd, c$stuttToInd, c$stuttParamInd , c$startIndMarker_nStutters,
            c$knownGind, as.integer(mlefit$maxThreads), c$relGind, c$ibdLong )[[1]] #obtain likelihood of all outcome

  #CHECK THAT COMPONENTS OBTAIN SAME:
  #for(m in 1:nM)  logLik = logLik + log( sum(exp( loglikVEC[ c$startIndMarker_nJointCombs[m] + 1:c$nJointCombs[m] ])) )

  #Step 1) Structuing  L(E|thetahat) for each marker
  for(m in 1:nM) { #extract info in c relevant for each markers:
#    print(m)
    loc =locs[m] #obtain locus name
    locDat = mlefit$data[[loc]] #obtain data for locus

    indKnownGind = (nC*(m-1)+1) : (nC*m) #get index where the genotype indices are
    knownGind = c$knownGind[indKnownGind] #extract the vector with genotype indices (-1 means unknown)
    uind <- which(knownGind==-1) #unknown genotype indices
    kind <- which(knownGind!=-1) #known genotype indices
    nU <- length(uind) #number of unknowns
    if(nU==0) { #if no unknowns
      dlist[[loc]] <- 0
      GClist2[[loc]] <- t(as.matrix(knownGind + 1)) #insert the conditional genotypes
      Glist2[[loc]] <- c$Gset[[loc]] #calcGjoint(freq=mlefit$popFreq[[loc]],nU=1)$G #Get genotypes
      next #skip to next marker
    }
    Gset <- Glist2[[loc]] <-  c$Gset[[loc]] #genotype possibilities
    Gset1 = paste0(Gset[,1],"/",Gset[,2])
    
    Glist <- list() #get genotype index list of all outcome
    for(k in 1:nU) Glist[[k]] <- 1:nrow(Gset)
    combGind <- expand.grid(Glist) #get all combinations
    combGind <- as.matrix(combGind,nrow=nrow(combGind))
    
    #obtain calculated likelihood direclty
    Gcombind <- c$startIndMarker_nJointCombs[m] + 1:c$nJointCombs[m] #get index of genotype outcome
    dvec <- loglikVEC[Gcombind] #obtain calculated values
    #Gset1[combGind]
    
    combGind2 <- numeric() #including the genotype of all contributors:
    for(k in 1:nC) { #for each contributors
      if(k%in%kind) combGind2 <- cbind(combGind2, rep(knownGind[k]+1,nrow(combGind)) ) #add genotype index of the reference
      if(k%in%uind) combGind2 <- cbind(combGind2, combGind[,which(k==uind)]) #add genotype index of the reference
    }
    isOK <- !is.infinite(dvec) #remove genotypes giving zero likelihood
    combGind2 <- combGind2[isOK,,drop=F] #removing genotypes which gives zero likelihood 
    dvec <- dvec[isOK]
    rank <- order(dvec,decreasing=TRUE) #order the genotypes wrt post prob values
    dlist[[loc]] <- dvec[rank]   #log-likelihood per per marker per
    GClist2[[loc]] <- combGind2[rank,,drop=F] #get index of ranked combined genotypes 
  } #end for each loci
 
  #POST PROCESSING:
  kvec <- 1:nC #index of contributors
  colN <- paste0("C",kvec) #column name of contributors
  
  #Step 3) Convert rank-list to list with allele-names
  deconvlisti <- list() #list per locus for all contributors
  for(loc in locs) {
   ii <- which(locs==loc)
   genv <-  paste0(Glist2[[loc]][,1],"/",Glist2[[loc]][,2]) #get vector of genotypes
   pGi <- exp(dlist[[loc]]) #convert to normal scale
   deconvlisti[[loc]] <- matrix(genv[ GClist2[[loc]] ],nrow=nrow(GClist2[[loc]])) #translate indices to genotype names
   deconvlisti[[loc]] <- cbind(deconvlisti[[loc]], pGi/sum(pGi) ) #add probabilities per makers
   colnames(deconvlisti[[loc]]) <- c(colN,"Probability")
  }
  #Step4) Create table layouts:
  
  #Helpfunctions to obtain marginal probabilities
  getMarg <- function(x,y) { #get marginal of genotypes
    agg <- aggregate(y,by=list(x),sum) #get probabilities
    ord <- order(agg[,2],decreasing=TRUE)
    agg2 <- agg[ord,,drop=F]
    colnames(agg2) <- c("Genotype","Probability")
    return(agg2)
  }
  getMarg2 <- function(x,y) { #get marginal of alleles
    tmp <- unlist(strsplit(x,"/"))
    unA <- unique(tmp) #unique alleles
    x2 <- t(matrix(tmp,nrow=2))
    prob <- rep(NA,length(unA))  
    for(aa in unA) prob[which(unA==aa)] <- sum(y[rowSums(x2==aa)>0]) #sum probabilities
    ord <- order(prob,decreasing=TRUE)
    agg <- data.frame(Allele=unA[ord],Probability=prob[ord])
    return(agg)
  }
  maxI <- function(p) min(min(which(cumsum(p)>=alpha)),maxlist,length(p))  #helpfunction to obtain a maximum size of a vector (bounded in both length and probability)
  
  #A) Calculate marginal probabilities for all contributors (genotypes and alleles):
  deconvlistic <- list() #genotype list per contributor
  deconvlistica <- list() #allele list per contributor
  cn <-  c("TopGenotype","probability","ratioToNextGenotype") #names for each contributor
  toplist <- list()
  for(loc in locs) {
    deconvlistica[[loc]] <- deconvlistic[[loc]] <- list()
    X <- deconvlisti[[loc]]
    nc <- ncol(X) #number of column
    nr <- nc - 1 #number of contributors
    tab <- matrix(,nrow=3,ncol=nr)
    rownames(tab) <- cn 
    colnames(tab) <- colN
    for(rr in 1:nr) {
      deconvlistic[[loc]][[colN[rr]]] <- tmp <- getMarg(x=X[,rr],y=as.numeric(X[,nc]))
      deconvlistica[[loc]][[colN[rr]]] <- getMarg2(x=X[,rr],y=as.numeric(X[,nc]))
      rat <- ifelse(nrow(tmp)>1, tmp[1,2]/tmp[2,2],NA) #get ratio from first to second genotype
      tab[,rr] <- c(tmp[1,1],signif(tmp[1,2],4),signif(rat,4)) 
    }
    toplist[[loc]] <- tab
  }
  
  
  #B) Create tables
  table1 <- table2 <- table3 <- table4 <- numeric()
  for(loc in locs) {
    combs <- deconvlisti[[loc]]
    prob <- as.numeric(combs[,ncol(combs)])
    combs[,ncol(combs)] <- signif(prob,4)
    maxind <- maxI(prob)
    combs <- combs[1:maxind,,drop=F]
    table1 <- rbind(table1,cbind(loc,1:nrow(combs),combs))
    table1 <- rbind(table1, rep("",ncol(table1)) )
  }
  colnames(table1)[1:2] <- c("Locus","Rank")
  
  for(loc in locs) table2 <- rbind(table2, c(toplist[[loc]]) )
  colnames(table2) <- paste0(cn,"_",c(t(replicate(length(cn),colN))))
  rownames(table2) <- locs
  
  maxI2 <- function(p) min(max(which(p>(1-alpha))),maxlist,length(p))
  space <- cbind("","","","")
  for(cc in colN) {
    for(loc in locs) {
     tmp <- deconvlistic[[loc]][[cc]]
     maxind <- maxI(p=tmp$Probability)
     newrow <- tmp[1:maxind,,drop=F]
     newrow[,2] <- signif(newrow[,2],4)
     newrows <- as.matrix(cbind(cc,loc,newrow))
     table3 <- rbind(table3,newrows,space)
    
     tmp <- deconvlistica[[loc]][[cc]]
     maxind <- maxI2(p=tmp$Probability)
     newrow <- tmp[1:maxind,,drop=F]
     newrow[,2] <- signif(newrow[,2],4)
     newrows <- as.matrix(cbind(cc,loc,newrow))
     table4 <- rbind(table4,newrows,space)
    }
  }
  colnames(table3)[1:2] <- colnames(table4)[1:2] <- c("Contr.","Locus")
  
  #new version adds rankGic and rankGica (genotype ranks per contributor in addition to per allele)
  return(list(table1=table1,table2=table2,table3=table3,table4=table4,toprankGi=toplist,rankGi=deconvlisti))
} #end function

