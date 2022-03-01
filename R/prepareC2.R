#' @title prepareC2
#' @author Oyvind Bleka
#' @description prepareC is used in the functions contLikMLE to prepare input to C-call
#' @details Assumes that all PH below threshold are removed. The function builds the data input to the C-code
#' @param dat A list which contains the output from function prepareData2 (samples,refs,freq)
#' @param repNames Names of replicate samples
#' @param nC Number of contributors in model (must be a consitant, same for all replicates)
#' @param condOrder Specify conditioning references from refData (must be consistent order). For instance condOrder=(0,2,1,0) means that we restrict the model such that Ref2 and Ref3 are respectively conditioned as 2. contributor and 1. contributor in the model. 
#' @param knownRef Specify known non-contributing references from refData (indices). For instance knownRef=(1,2) means that reference 1 and 2 is known non-contributor in the hypothesis. This affectes coancestry correction.
#' @param kit shortname of kit: Obtained from getKit(). 
#' @param AT Analytical threshold. A vector per sample, or a list with marker elements per sample
#' @param pC Dropin prob parameter. A vector per sample, or a list with marker elements per sample
#' @param lambda DropinPH parameter. A vector per sample, or a list with marker elements per sample
#' @param fst The co-ancestry coefficient. Default is 0. A vector per marker.
#' @param incBS A boolean whether potential BW stutters are included
#' @param incFS A boolean whether potential FW stutters are included
#' @param knownRel Specify the index of the related contributing reference from refData (one index). For instance knownRel=2 means that unknown1 is related to reference 2 with ibd specified relationship.
#' @param ibd the identical by decent coefficients list of the relationship for each unknowns (specifies the type of relationship). Default is NULL, meaning no related inds. Can be a list for each defined related pairs.
#' @return ret A list of data input to call the C-code with
#' @export 


prepareC2 = function(dat,repNames,nC,condOrder=NULL,knownRef=NULL,kit=NULL,AT=50, pC=0.05,lambda=0.01,fst=0,incBS=FALSE,incFS=FALSE, knownRel=NULL, ibd=NULL){
 Qallele="99" #Name of allele given if missing in evidence (drop-out allele). Defualt is 99. This is important when considering the degradation model since 99 is closest to maximum allelein a locus. 

 #An overview of which markers has the replicate:
 #sampleMarkerNames = sapply(dat, function(x) names(x$samples)) 

 modelDegradation = rep(FALSE,length(repNames))
 if(!is.null(kit) && any(!is.na(kit)) ) {
   kits = na.omit(unique(kit))
   kitinfoList = list()
   for(k in 1:length(kits)) {
     kitinfo = euroformix::getKit(kits[k])
     if(is.na(kitinfo)[1]) stop("Specified kit not found. Needs to be specified in order to model degradation.")
     kitinfoList[[ kits[k] ]] = kitinfo
   }
   modelDegradation = kit%in%names(kitinfoList) #indicate whether to model degradation for different reps
 }
 
 #obtain markers to evaluate
 locs = names(dat) #can also be obtained using this object
 nLocs = length(locs) #number of loci
 
 #Get list of genotypes for each markers (MUST BE SAME ORDER AS THOSE CREATED IN C++ code!!)
 Gset <- list() 
 for(loc in locs) {
#   loc=locs[1]
    freq=dat[[loc]]$freq
    alleles <- names(freq)   
    nn = length(alleles)
    
    #G matrix is the vectorized upper triangular (1,1),(1,2),...,(1,n),(2,2),(2,3),...,(2,n),....,(n,n)
    G = numeric()
    for(i in 1:nn) {
       G = rbind(G, cbind( alleles[rep(i,nn - i + 1)], alleles[i:nn] ))
    }
    Gset[[loc]] <- G #store genotype
 }
 nGenos = sapply(Gset,nrow) #number of genotypes per marker
 
 #Fix references as known contributors: Assign genotypes of known references to knownGind-matrix
 NOK = rep(0,nLocs) #number of known contributors per loci
 knownGind <- matrix(-1,ncol=nLocs,nrow=nC) #default is no references (=-1)
 #assign references to knownGind-matrix by values of Glist
 if(!is.null(condOrder) && any(condOrder>0)) {
   for(loc in locs) {
     locind = which(loc==locs) #get loc-index
     subRef <- dat[[loc]]$refs #take out reference profiles
     if(length(subRef)==0)  stop(paste('Missing locus (',loc,') in refData.',sep=''))
     for(k in which(condOrder>0) ) { #consider each reference to condition on
         if(length(subRef[[k]])==0) next #Allowing unknown contributors for missing markers (provided with empty element character())
         if(length(subRef[[k]])!=2) stop("References need to have exactly two alleles.") 
         
         NOK[locind] = NOK[locind] + 1 #add known contributor
         Gind1 <- subRef[[k]][1]==Gset[[loc]][,1] & subRef[[k]][2]==Gset[[loc]][,2]
         Gind2 <- subRef[[k]][2]==Gset[[loc]][,1] & subRef[[k]][1]==Gset[[loc]][,2]
         knownGind[condOrder[k],locind] = which(Gind1 | Gind2) - 1 #subtract with one since we work from 0-indice
     }
   }
 }
 #NOU = nC - NOK #number of unknowns per markers

 #RELATEDNESS:
 #Assign genotypes of related references to relGind-matrix
 relGind <- matrix(-1,ncol=nLocs,nrow=nC) #default is no references (=-1)
 #assign references to knownGind-matrix by values of Glist
 if(!is.null(knownRel) ) {
   for(loc in locs) {
     locind = which(loc==locs)
     subRef <- dat[[loc]]$refs #take out reference profiles
     if(length(subRef)==0)  stop(paste0("Missing locus ",loc," in refData."))
     for(rel in knownRel) { #for each related reference
       subRefRel = subRef[[rel]] #obtain alleles
       relUind =  max(NOK) + which(rel==knownRel) #index of unknown contributors (Note to check structure in C code)
       if(relUind > nC) stop("More related than unknowns were specified!")
       if(length(subRefRel)==0) { #allow for missing markers (exchanging with an unknown)
         #print(paste0("Note: Missing alleles at locus ",loc," for related reference ",names(subRef)[k]))
         next
       } else if(length(subRefRel)!=2) {
         stop("References must have exactly two alleles.")
       }
       
       #Obtain genotype of related
       Gind1 <- subRefRel[1]==Gset[[loc]][,1] & subRefRel[2]==Gset[[loc]][,2]
       Gind2 <- subRefRel[2]==Gset[[loc]][,1] & subRefRel[1]==Gset[[loc]][,2]
       relGind[relUind ,locind] = which(Gind1 | Gind2) - 1 #subtract with one since we work from 0-indice
     }
   }
 }
 ibd0 =  c(1,0,0) #unrelatedness is default for each contributors (also known)
 ibdList = list()
 for(cc in seq_len(nC) ) ibdList[[cc]] = ibd0 #insert default values for each contributors  
 if(!is.null(knownRel) && !is.null(ibd)){ #Relatedness was defined
    if(!is.list(ibd)) ibd = list(ibd) #force to be list
    if(length(knownRel)!=length(ibd)) stop("ibd must have same number of list elements as knownRel")
    for(rel in knownRel) { #for each related reference we extract the corresponding ibd-vector
      indRel =  which(rel==knownRel) #index of relatedness
      ibdList[[ max(NOK) + indRel]] = ibd[[indRel]] #index of unknown contributors]] 
    }
 }
 ibdLong = unlist(ibdList) #long vector over all contributors
 
 
 #TRAVERSE EACH MARKER AND STRUCTURE DATA FOR c++ input:
 #PREPARE FREQUENCY AND COVERAGE INFO
 datList = list() 
 for(m in 1:nLocs) { #m=5
   datList[[m]] = list()
   
   loc = locs[m] #for selected loc (already upper)
   freq = dat[[loc]]$freq #get allele freqs
   allelesALL = names(freq) #allele also including Q-alele
   alleles = setdiff(allelesALL,Qallele) #assumed order as in freqs (assume alleles already added!)
      
   datList[[m]]$freq = freq #STORE FREQS (Q-allele already included)
   datList[[m]]$alleles = allelesALL #alleles #store allele names (Include possible Q-alleles)
   
   #PREPARE OBSERVATIONS AND DROPIN WEIGHTS
   datList[[m]]$nAlleles <- nAlleles <- length(allelesALL) #number of allees to travers through (Including allele 99)
   locDat = dat[[loc]]$samples#obtain replicate info
   
   datList[[m]]$nReps <- nReps <- length(locDat) #obtain number of samples
   repNames_marker = names(dat[[loc]]$samples) #obtain replicate names for specific marker
   
   repID <-  which(repNames%in%repNames_marker) #obtain replicate ID of reps considered 
   datList[[m]]$repID = repID - 1 #MUST START WITH INDEX 0 for C code
   
   #obtain other params (different formats are allowed):
   getVal = function(x) { #helpfunction
     if(length(x)==1) {
       x0 = x #is constant
     } else { #x is a vector or a list
       if( is.vector(x) ) {
         x0 = x[sample] #assume a constant per replicate
       } else if( is.list(AT) ) {
         x0 = x[[sample]][[loc]] #assume specified per marker per replicate
       } else {
         stop("Non recognized format for argument")
       }
     }
   }
 
   ATs <- lambdas <- pCs <- fsts <- rep(NA,nReps)
   for(r in 1:nReps) { #traverse each replicate
    sample = repNames_marker[r] #obtain sample name
    ATs[r] = getVal(AT)
    lambdas[r] = getVal(lambda)
    pCs[r] = getVal(pC)
   }
   datList[[m]]$AT <- ATs #vector per marker
   datList[[m]]$pC <- pCs #vector per marker
   datList[[m]]$lambda <- lambdas #vector per marker
   
   fst0 = fst[loc] #obtain marker specific (if given as vector)
   if(is.na(fst0)) fst0 = fst[1] #not found so setting first (default)
   datList[[m]]$fst <- fst0 #this is constant per marker
   
   #Inesrt peak height vector   
   yv <- fv <- bv <- matrix(0,ncol=nAlleles,nrow=nReps) #create PH-matrix and drop-in weight matrix
   for(r in 1:nReps) { #for each replicates (following handle unordered loci under each sample)
     yv[r, match(names(locDat[[r]]),alleles)] = locDat[[r]] #insert PHs
     #fv[r,] = log(lambda0) - lambda0*(yv[r,] - AT0 ) + log(prC0) + log(freq)
     
     #drop-in weights vector (per allele depending on PH and freqs) 
     weight  = log(lambdas[r]) - lambdas[r]*(yv[r,] - ATs[r]) #obtain weight for each observation (on log-scale)
     #dexp( yv[r,] - ATs[r], rate=lambdas[r], log=TRUE)  #
     fv[r,] = weight + log(freq) #Important to include allele freq as part of drop-in weights
     
     #Obtain base pair info (used for degradation)
     if( !is.null(kit) ) {
       kit0 = kit[ repID[r] ] #kit for particular replicate
       if(!is.na(kit0)) { #must exist
         kitinfo = kitinfoList[[ kit0 ]]
         subkitinfo = kitinfo[toupper(kitinfo$Marker)==loc,,drop=FALSE]
         if(nrow(subkitinfo)==0) stop(paste0("No kit info found at marker ",loc," for ",kit0))
    
         kitinfo_size =  subkitinfo$Size 
         kitinfo_allele = subkitinfo$Allele #assume name as allele names
         bp <- numeric() #fragment length/size in bp
         for(allel in alleles) { #for each alleles (not Q-allele
           bind <- which( kitinfo_allele==allel ) #obtain index of position
           if(length(bind)==0) { #if missing alleles in kitinfo
               bind <- which.min(abs(as.numeric(allel) - as.numeric(kitinfo_allele))) 
           } 
           bp <- c(bp,kitinfo_size[bind]) #get fragment length
         } #end for each alleles
         if(Qallele%in%allelesALL) bp = c(bp, max(kitinfo_size) ) #set max size for Q-allele if considered
         bv[r,] = (bp-125)/100 #rescale before
       }
     }
   }
   datList[[m]]$PHobs = unlist(yv) #insert PHs for obs alleles
   datList[[m]]$PHdropin = unlist(fv) #insert drop-in weights (for each observations)
   datList[[m]]$basepairs = unlist(bv) #insert fragment length info (basepairs)
  
   #Check if alleles are numeric:
   suppressWarnings({ alleles2 = as.numeric(alleles)}) #convert to numeric
   isNum = ifelse(!any(is.na(alleles2)),TRUE,FALSE ) #check if all alleles are numeric

   #PREPARE STUTTER INFO	  
   #obtain potential stutters from observed alleles:
   stuttTypes = NULL #obtain stutter types
   if(incBS & isNum) stuttTypes = append(stuttTypes,"-1")
   if(incFS & isNum) stuttTypes = append(stuttTypes,"+1")
   
	 getStutteredAllele = function(allele,type) {
	   paste0( as.numeric(allele) + as.numeric(type) )
	 }
	    
	 potStutters = NULL
	 stuttList = list() #store stuttered info
	 for(type in stuttTypes) { #for each stutter type
		 stutt <- stuttList[[type]] <- getStutteredAllele(alleles,type)
		 potStutters = c(potStutters,stutt[!stutt%in%alleles]) #INsert potential stutters
	 }
	  
	 #obtain Alleles to include as data:
	 alleles2 = c(allelesALL,unique(potStutters))  #insert potential stutters to allele list     
	 nAlleles2 = length(alleles2) #total number of alleles (including potential stutters)
	  
	 #Obtain indices of each stutters types
	 stuttFromInd <- stuttToInd <- stuttParamInd <- NULL #init vector
	 for(type in stuttTypes) { #for each stutter type
#	   type=stuttTypes[1]
	   if(length(alleles)==0) next #skip if no alleles
		 fromInd = 1:length(alleles) - 1 #from index. NB; subtract index to start from 0 (c++)
		 toInd = match(stuttList[[type]],alleles2) - 1 #obtain index which each allele in alleles stutter to. NB: subtract index
		 paramInd = rep(which(type==stuttTypes) ,length(toInd)) - 1 #obtain param index NB; subtract index to start from 0 (c++)
		 
		 stuttFromInd = c(stuttFromInd,fromInd) 
		 stuttToInd = c(stuttToInd,toInd)
		 stuttParamInd = c(stuttParamInd,paramInd) 
	 }
   datList[[m]]$nStutters = length(stuttParamInd) #total number of stutters to traverse
   datList[[m]]$stuttFromInd = stuttFromInd
   datList[[m]]$stuttToInd = stuttToInd 
   datList[[m]]$stuttParamInd = stuttParamInd 
   
   nPS = length(alleles2) - length(allelesALL) #number of potential stutters (not observed)
   datList[[m]]$nAlleles2 = nAlleles + nPS #number of alleles including the potential ones
   
   #Obtain contribution matrices (per locus)
   numGenos1p = nGenos[loc] #obtain number of genotypes (1 contr)
   outG1contr = matrix(0,nrow=numGenos1p, ncol = nAlleles, 0.0); #init nG1xnA matrix (contribution matrix). Indicating what alleles that are contributoed
   outG1allele = matrix(0,nrow=numGenos1p,ncol=2) #init nG1x2 matrix (allele names as indices 0,...,nA-1)
   cc = 1 #0; #counter oveer all
   for (i  in 0:(nAlleles-1)) {
      for (j in i:(nAlleles-1)) { 
         outG1allele[cc,1] = i; #include index
         outG1allele[cc,2] = j; #include index
         outG1contr[cc,i+1] = outG1contr[cc,i+1] + 1.0; #insert contr at allele i
         outG1contr[cc,j+1] = outG1contr[cc,j+1] + 1.0; #insert contr at allele j
         cc = cc + 1; #iterate to next genotype outcome
      }
   }
   datList[[m]]$outG1allele = outG1allele
   datList[[m]]$outG1contr = outG1contr
   
   #Get number of typed allelse (must count alleles)
   tmp <- rep(0, length(allelesALL))
   if(!is.null(dat[[loc]]$refs)) {
      typedRefs = unique( c(which(condOrder>0),knownRef) )#,knownRel) ) #get unique referneces
      for(k in typedRefs) { #for each typed refs
         ind <- which( allelesALL%in%dat[[loc]]$refs[[k]] )
         tmp[ind] = tmp[ind] + (length(ind)==1) + 1 #add twice sampled if homozygote 
      }
   }
   datList[[m]]$nTyped <- sum(tmp) #number of total sampled (for each loci)
   datList[[m]]$maTyped <- tmp  #add vector of typed
 } #end for each marker
      
 
 #INITS:
 #Variables to use in the C code:
#  names(datList[[1]])
 #nLocs #number of loci
 
 nReps = as.integer(sapply(datList,function(x) x$nReps)) #number of samples per locus
 repID = as.integer( unlist( sapply(datList,function(x) x$repID) )) #get replicate ID per locus
 nAlleles =  as.integer(sapply(datList,function(x) x$nAlleles)) #number of observed alleles (including Q) per locus
 nAlleles2 =  as.integer(sapply(datList,function(x) x$nAlleles2)) #number of observed alleles (including Q) per locus
 freq =  as.numeric(unlist(sapply(datList,function(x) x$freq)))

 peakHeights = as.numeric( unlist(sapply(datList,function(x) x$PHobs)) ) #NOte the vectorization makes each replicates after each other for each allele
 dropinWeight = as.numeric( unlist(sapply(datList,function(x) x$PHdropin)) ) #NOte the vectorization makes each replicates after each other for each allele
 basepairs =  as.numeric( unlist(sapply(datList,function(x) x$basepairs)) ) 
 
 if(length(peakHeights)!=length(dropinWeight)) stop("Missmatch length")
 if(length(peakHeights)!=length(basepairs)) stop("Missmatch length")

 nStutters = as.integer(sapply(datList,function(x) x$nStutters)) #number of stutters per locus
 stuttFromInd = as.integer( unlist(sapply(datList,function(x) x$stuttFromInd)) ) #- 1 #Note the -1 subtraction
 stuttToInd = as.integer( unlist(sapply(datList,function(x) x$stuttToInd)) ) #- 1 #Note the -1 subtraction
 stuttParamInd = as.integer( unlist(sapply(datList,function(x) x$stuttParamInd)) ) #- 1 #Note the -1 subtraction
 
 #Genotype traversing:
 nGenos =  as.integer(nGenos) #number of genotypes to traverse
 outG1allele =  as.integer( unlist(sapply(datList,function(x) t(x$outG1allele))) ) #obtain genotype outcome
 outG1contr = as.integer( unlist(sapply(datList,function(x) t(x$outG1contr))) ) #obtain genotype outcome
 
 #Other info:
 maTyped = as.numeric( unlist(sapply(datList,function(x) x$maTyped)) )
 nTyped = as.numeric( sapply(datList,function(x) x$nTyped) )
 AT = as.numeric( unlist( sapply(datList,function(x) x$AT) ) )
 fst = as.numeric(sapply(datList,function(x) x$fst) )
 lambda = as.numeric( unlist( sapply(datList,function(x) x$lambda) ) )
 dropinProb = as.numeric( unlist( sapply(datList,function(x) x$pC) )) #dropin probs
 
 #CREATE CUMULATIVE vectors to make simple lookup in C++ code (per marker)
 startIndMarker_nAlleles = as.integer(c(0,cumsum(nAlleles))) #used to loop peak heights
 #startIndMarker_nAlleles2 = c(0,cumsum(nAlleles2)) #used to init shape vector
 startIndMarker_nStutters = as.integer(c(0,cumsum(nStutters))) #used for looping stutter indices
 startIndMarker_nGenos = as.integer(c(0,cumsum(nGenos)))
 startIndMarker_nAllelesReps = as.integer(c(0,cumsum(nAlleles*nReps))) #used to loop peak heights, fragment size etc (across Replicates)
 startIndMarker_nRepMarkers = as.integer(c(0,cumsum(nReps))) #used to indicate indicate which repID to use ()
 
 #Vectorize outG1 matrices and calculate start index (per locus) for these
 startIndMarker_outG1allele = as.integer(2*startIndMarker_nGenos) #start index for 1p contributor matrix
 startIndMarker_outG1contr = as.integer(c(0,cumsum(nAlleles*nGenos))) #start index for per-allele contribution matrix

 #REMAINING variables:
 NOU= as.integer(nC - NOK) #number of unknowns
 nJointCombs = as.integer(nGenos^NOU) #joint combinations
 startIndMarker_nJointCombs = as.integer(c(0,cumsum(nJointCombs)))
 
 #contribution, relatedness etc
 NOC=as.integer(nC)
 NOK=as.integer(NOK)
 knownGind=as.integer(knownGind)
 relGind=as.integer(relGind)
 ibdLong=as.numeric(ibdLong) 

 #t(matrix(outG1contr[1:(nAlleles[1]*nGenos[1])],nrow=nAlleles[1]))
 #datList[[1]]$outG1contr
 alleleNames = lapply(datList,function(x) x$alleles) #obtain allele names for each marker (only observed)
 names(alleleNames) = locs #insert locus names
 
 retlist = list( NOC=NOC, NOK=NOK, NOU=NOU,knownGind=knownGind, relGind=relGind,ibdLong=ibdLong, #hypotheseis
              nLocs=nLocs,nReps=length(repNames), nRepMarkers=nReps,nAlleles=nAlleles,nAlleles2=nAlleles2, #dimensions
              nGenos=nGenos, nJointCombs=nJointCombs,  freq=freq,peakHeights=peakHeights, basepairs=basepairs,#data
              repNames=repNames, repID=repID, startIndMarker_nRepMarkers=startIndMarker_nRepMarkers,  #rep info
              AT = AT,fst=fst, lambda=lambda, dropinProb=dropinProb, dropinWeight=dropinWeight, nTyped=nTyped,maTyped=maTyped, #settings data
              nStutters=nStutters,stuttFromInd=stuttFromInd,stuttToInd=stuttToInd,stuttParamInd=stuttParamInd, #stutter info
              startIndMarker_nAlleles=startIndMarker_nAlleles,startIndMarker_nAllelesReps=startIndMarker_nAllelesReps,
              startIndMarker_nStutters=startIndMarker_nStutters,startIndMarker_nGenos=startIndMarker_nGenos, startIndMarker_nJointCombs=startIndMarker_nJointCombs,#cumulative
              outG1allele=outG1allele, outG1contr=outG1contr,startIndMarker_outG1allele=startIndMarker_outG1allele, startIndMarker_outG1contr=startIndMarker_outG1contr,  #contributor matrix
              locNames=locs,Gset=Gset,alleleNames=alleleNames)  #meta info
 
 return(retlist)
} #end function

