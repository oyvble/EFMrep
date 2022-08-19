#' @title contLikMLE2
#' @author Oyvind Bleka
#' @description contLikMLE optimizes the likelihood function of the DNA mixture model 
#'
#' @param nC Number of contributors in model
#' @param samples A List with samples which for each samples has locus-list elements with list elements adata and hdata. 'adata' is a qualitative (allele) data vector and 'hdata' is a quantitative (peak heights) data vector.
#' @param popFreq A list of allele frequencies for a given population.
#' @param refData Reference objects with (2-size) allele-vector given in list element [[i]][[s]].
#' @param condOrder Specify conditioning references from refData (must be consistent order). For instance condOrder=(0,2,1,0) means that we restrict the model such that Ref2 and Ref3 are respectively conditioned as 2. contributor and 1. contributor in the model. 
#' @param knownRef Specify known non-contributing references from refData (index). For instance knownRef=(1,2) means that reference 1 and 2 is known non-contributor in the hypothesis. This affectes coancestry correction.
#' @param kit A vector of kitname for each sample name, used to model degradation. Must be one of the shortnames of kit: check getKit()
#' @param AT Analytical threshold. A vector per sample, or a list with marker elements per sample
#' @param pC Dropin prob parameter. A vector per sample, or a list with marker elements per sample
#' @param lambda DropinPH parameter. A vector per sample, or a list with marker elements per sample
#' @param fst The co-ancestry coefficient. Default is 0. A vector per marker.
#' @param mixProp Mixture proportion params: Indicate which shared parameter group each samples belongs to
#' @param PHexp PeakHeight expectation params: Indicate which shared parameter group each samples belongs to
#' @param PHvar PeakHeight variability params: Indicate which shared parameter group each samples belongs to
#' @param DEG Degradation slope params: Indicate which shared parameter group each samples belongs to
#' @param stuttBW Backward stutter params: Indicate which shared parameter group each samples belongs to
#' @param stuttFW Forward stutter params: Indicate which shared parameter group each samples belongs to
#' @param minF The freq value included for new alleles (new alleles as potential stutters will have 0). Default NULL is using min.observed in popFreq.
#' @param normalize A boolean of whether normalization should be applied or not. Default is FALSE.
#' @param maxIter Maximum number of iterations for the optimization. Default is 100.
#' @param steptol Argument used in the nlm function for faster return from the optimization (tradeoff is lower accuracy).
#' @param nDone Number of optimizations required providing equivalent results (same logLik value obtained)
#' @param maxThreads Maximum number of threads to be executed by the parallelization
#' @param delta Scaling of variation of normal distribution when drawing random startpoints. Default is 1.
#' @param seed The user can set seed if wanted
#' @param verbose Boolean whether printing optimization progress. Default is TRUE.
#' @param knownRel gives the index of the reference which the 1st unknown is related to.
#' @param ibd the identical by decent coefficients of the relationship (specifies the type of relationship)
#' @return ret A list from model fit 
#' @export

#Function fit MLE for replicates
#maxIter=1000;steptol=1e-6; nDone=3; maxThreads=128; delta=1; verbose=TRUE; seed=1;knownRel=NULL;ibd=NULL;#minF=NULL;normalize=FALSE
contLikMLE2 = function(nC,samples,popFreq,refData=NULL, condOrder=NULL,knownRef = NULL, kit=NULL, AT=50,pC=0.05,lambda=0.01,fst=0, 
                       mixProp=1:length(samples), PHexp=1:length(samples), PHvar=1:length(samples), DEG=rep(0,length(samples)), 
                       stuttBW=rep(0,length(samples)), stuttFW=rep(0,length(samples)), minF=NULL, normalize=FALSE,
                       maxIter=100, steptol=1e-6, nDone=3, maxThreads=0, delta=1, seed=NULL,verbose=TRUE,
                       knownRel=NULL, ibd=NULL) {

  if(!is.null(seed)) set.seed(seed)
  if( any(DEG>0) && length(DEG)!=length(kit)) stop("Length of kit argument must be same as length of DEG")
  if(nC < sum(condOrder>0)) stop("Number of contributors can't be specified less than number of conditional references.")
  
  #Make sure that the order is correct for the model configurations:
  fixOrder = function(x,canContainZero=FALSE) {
    if(!canContainZero && 0%in%x) stop( paste0("A model config argument contained zero!"))
    
    notZero = x>0
    if(!any(notZero)) return(x)
    grps = unique(x[notZero]) #get groups
    for(g in seq_along(grps)) {
      x[x==grps[g]] = g
    }
    return(x)
  }
  
  mixProp = fixOrder(mixProp)
  PHexp = fixOrder(PHexp)
  PHvar = fixOrder(PHvar)
  DEG = fixOrder(DEG,TRUE)
  stuttBW = fixOrder(stuttBW,TRUE)
  stuttFW = fixOrder(stuttFW,TRUE)
  
  #Prepare data:
  dat = prepareData2(samples,refData,popFreq,minF = minF, normalize = normalize, AT=AT) #obtain data to use for analysis
  
  repNames=names(samples) #obtain replicate names (this is the index order for model params)
  nReps = length(repNames)
  #Prepare data for C-object:
# incBS=any(stuttBW>0);incFS= any(stuttFW>0)
  c <- prepareC2(dat, repNames, nC, condOrder,knownRef,kit,AT,pC,lambda,fst,incBS=any(stuttBW>0), incFS= any(stuttFW>0),knownRel=knownRel, ibd=ibd)
  c$freq
#c$locNames
#c$stuttParamInd
#c$startIndMarker_nAlleles
#c$startIndMarker_nAllelesReps
  #Prefitting data based on the model for sum of the peak heights  to find proper startvalues for MLE 
  nLocs = c$nLocs #number of markers
  sumY <- meanbp <- matrix(NA,nrow=nReps,ncol=nLocs) #sum PH, bp per marker per replicates
  startindMarker = 0 #marker start index (over all reps)
  for(marker in 1:nLocs) { #for each marker
    for(rep in 1:nReps) { #for each replicate
      indR = (0:(c$nAlleles[marker]-1))*c$nRepMarkers[marker] + rep + startindMarker #obtain indices for replicate r
      heights = c$peakHeights[indR] #obtain PHs
      sumY[rep,marker]  <- sum(c$peakHeights[indR]) #take sum of the peak heights
      meanbp[rep,marker] <- mean(c$basepairs[indR]) 
    }
    startindMarker = startindMarker + c$nRepMarkers[marker]*c$nAlleles[marker] #skip to next marker
  }
  
  theta_Reps = matrix(ncol=3,nrow=nReps) #obtain start values for each replicates (PHexp,PHvar,Deg)
  for(rep in 1:nReps) { #for each replicate
    usedeg=DEG[rep]>0 #whether degradation are modelled for start model
    
    meanbp1 = meanbp[rep,]
    if(!usedeg) meanbp1 = NULL
    th0 <- fitgammamodel2(y=sumY[rep,],x=meanbp1,delta=delta,offset=0,scale=1,restrictDeg=FALSE)
    if(is.null(th0)) stop("Failed to find suitable start values in fitgammamodel!")
    theta_Reps[rep,1:length(th0)] = th0 #insert pre-fitted value
  }   
  
  #PREPEARING THE LIKELIHOOD OPTIMZATION (what parameters are provided?):
  
  #INNER HELPFUNCTIONS:
  logit = function(x) log(x/(1-x))
  drawMixprop = function(NOC) {
    #CONSIDER Mixture porpotions    
    ncond = sum(condOrder>0)
    nU = NOC - sum(condOrder>0)
    mxrnd = rgamma(NOC,1) #Draw simplex (flat)
    mxrnd = mxrnd/sum(mxrnd)
    if(nU>1) { #sort if more than 1 unknown
      ind = (ncond+1):nC #sort Mix-prop for the unknowns 
      mxrnd[ind] = sort(mxrnd[ind],decreasing = TRUE) #sort Mx in decreasing order
    }
    
    #convert Mx values to real domain (nu:
    nurnd = numeric()
    if(NOC>1) {
      cs = 0 #c( 0,cumsum(mxrnd)) #cumulative sum of mixture proportins
      for(cc in 1:(NOC-1)) { #traverse contributors (Restricted)
        nurnd = c(nurnd, logit( mxrnd[cc]/(1-cs))) 
        cs = cs + mxrnd[cc] #update sum
      }
    }
    return(nurnd)
  }
  
  drawPHvars = function(th0,delta=1) {   #CONSIDER PH prop variables
    th0 = na.omit(th0) #remove possible NAs
    sdPH = delta*0.15*th0 #obtain considered SD of PH props 
    return( log( abs( rnorm(length(th0),th0,sd=sdPH)) ) )  #Obtain random start for mu/sigma/tau, Note using the delta here (should be small)
  }
  
  drawStutterProp = function(exp=0.05,sd=0.5) {
    return( rnorm(1,logit(exp),sd=sd) )
  }
  
  convParamBack = function(phi) {
    #convert param vector back to individual param used in lik-function
    par_PHexp <- par_PHvar <- par_DEG <- rep(1,nReps) #same length
    par_mixProp <- rep(1,nReps*nC) #span out all mix-prop (easier in C-code)
    par_stutt <- rep(0,nReps*2) #we have two types of stutters (per replicate). ALWAYS
    cc = 1 #index traverser of phi vector

    if(nC>1) { #OBTAIN mixture proportions if more than 1 contr:
      for(id in 1:max(mixProp)) {
        tmpvec = rep(NA, nC-1) #values
        mxprop = phi[cc + 1:(nC-1) - 1] #temp
        cs = 0 
        for(k in 1:(nC-1)) { #for each mix props
          tmpvec[k] = (1-cs)/(1+exp(-mxprop[k])) #convert back (ilogit)
          cs = cs + tmpvec[k] #add mixture proportion to cumulative sum 
        }
        tmpvec = c(tmpvec,1-sum(tmpvec)) #inlcude last
        cc = cc + (nC - 1) #update counter when completed
        
        #INSERT Mx-VALUES
        indIDs = which(mixProp==id) #obtain replicate IDs
        for(indID in indIDs) par_mixProp[ nC*(indID-1) + 1:nC ] = tmpvec
      }
    }
    
    #OBTAIN PH params (PHexp,PHvar,DEG) after eachother:
    for(id in 1:max(PHexp)) {
      par_PHexp[PHexp==id] = exp(phi[cc])
      cc = cc + 1 #update to next
    }
    for(id in 1:max(PHvar)) {
      par_PHvar[PHvar==id] = exp(phi[cc])
      cc = cc + 1 #update to next
    }
    if(any(DEG>0)) { #draw backward stutters
      for(id in 1:max(DEG)) {
        par_DEG[DEG==id] = exp(phi[cc])
        cc = cc + 1 #update to next
      } 
    }
    
    #DRAW STUTTERS
    if(any(stuttBW>0)) { #draw backward stutters
      insIND = 2*(1:nReps)-1 #indices to insert for BW stutter (odd numbers)
      for(id in 1:max(stuttBW)) {
        par_stutt[ insIND[stuttBW==id] ] =  1/(1+exp( -phi[cc] ) )
        cc = cc + 1 #update to next
      } 
    }
    if(any(stuttFW>0)) { #draw forward stutters
      insIND = 2*(1:nReps) #indices to insert for FW stutter (even numbers)
      for(id in 1:max(stuttFW)) {
        par_stutt[ insIND[stuttFW==id] ] =  1/(1+exp( -phi[cc] ) )
        cc = cc + 1 #update to next
      } 
    } 
  
    return( list(mixProp=par_mixProp, PHexp=par_PHexp, PHvar=par_PHvar, DEG=par_DEG, stutt=par_stutt) )
  }
  
  
  
  #function for calling on C-function: Must convert "real domain" values back to model params 
  negloglik_phi <- function(phi) { #assumed order: mixprop(1:C-1),mu,sigma,beta,xi
    
    par = convParamBack(phi) #obtain list with parameters (on original scale)
    
    #Hypothesis:
    #c$nJointCombs #join genotype combinations per marker
    #c$NOC #number of contirbutors (constant)
    #c$NOK #Number of known contributors (per marker)
    
    #Other params: 
    #as.numeric(c$AT) #analytical/detection threshold (per locus)
    #as.numeric(c$fst) #theta-correction (per locus)
    #as.numeric(c$nNoiseParam) #number of noise (per locus)
    #as.numeric(c$noiseSizeWeight) #noise size weighting (per allele)
    
    #names(c)
    #data structure and observations
    #as.integer(c$nLocs) #number of loci
    #as.integer(c$nSamples) #number of samples per locus
    #as.integer(c$nAlleles) #number of observed alleles (inlcuding Q)
    #as.integer(c$nAlleles2) #number of observed alleles (inlcuding Q) + number of potential stutters
    #as.integer(c$nGenos) #number of genotype combination (1 contributor)
    #as.numeric(c$peakHeights) #coverage vector (including Q)
    #as.numeric(c$freq) #frequency vector
    #as.numeric(c$nTyped) #Number of previously typed alleles. Used for theta-correction
    #as.numeric(c$maTyped) #Number of previously typed alleles (per allele). Used for theta-correction
    #as.numeric(c$basepairs) #adjusted fragmentlength in basepairs (used for degradation)
    #as.numeric(c$repID) #Indicate which replicate a specific replicate in each marker corresponds to (IMPORTANT)
    
    #stutter part:  
    #as.integer(c$nStutters) #Number of stutters
    #as.integer(c$stuttFromInd) #index stutters from (starts from 0)
    #as.integer(c$stuttToInd) #index stutters to  (starts from 0)
    #as.integer(c$stuttParamInd) #index of stutter parameter to use  (starts from 0)
    
    #Allele outcome/contribution  
    #as.integer(c$outG1allele) #Allele index for genotype outcome
    #as.integer(c$outG1contr) #Allele contriubtion for genotype outcome
    
    #Cumulative outcome   (start index for markers)
    #as.integer(c$startIndMarker_nAlleles)  #Startindex for number of alleles
    #as.integer(c$startIndMarker_nAllelesReps)  #Startindex for number of alleles (including replicates)
    #as.integer(c$startIndMarker_nRepMarkers)  #Startindex for number of replicates per marker 
    #as.integer(c$startIndMarker_nStutters) #Startindex for number of stutters
    #as.integer(c$startIndMarker_nGenos) #Startindex for number of genotypes
    #as.integer(c$startIndMarker_outG1allele) #Startindex for number of alleles
    #as.integer(c$startIndMarker_outG1contr) #Startindex for number of alleles
    #as.integer(c$startIndMarker_nJointCombs) #Startindex for number of joint genotypes
    
    #Relatedness
    #as.integer(c$relGind) #genotype of related reference
    #as.integer(c$ibdLong) #relatedness definition of each contributors (also known)
    
    calc = .C("loglikGamma_allcomb1",as.numeric(0), c$nJointCombs, c$NOC, c$NOK,
              as.numeric(par$mixProp),  as.numeric(par$PHexp), as.numeric(par$PHvar), as.numeric(par$DEG), as.numeric(par$stutt), 
              as.numeric(c$AT),as.numeric(c$fst),as.numeric(c$dropinProb),as.numeric(c$dropinWeight),
              c$nReps, c$nLocs, c$nRepMarkers, c$nAlleles, c$nAlleles2, c$startIndMarker_nAlleles, c$startIndMarker_nAllelesReps,
              c$peakHeights, c$freq, c$nTyped, c$maTyped, c$basepairs, c$repID, c$startIndMarker_nRepMarkers, 
              c$nGenos, c$outG1allele, c$outG1contr, c$startIndMarker_outG1allele, c$startIndMarker_outG1contr, c$startIndMarker_nJointCombs,
              c$nStutters, c$stuttFromInd, c$stuttToInd, c$stuttParamInd , c$startIndMarker_nStutters,
              c$knownGind, as.integer(maxThreads), c$relGind, c$ibdLong  ) 
    
    loglik = calc[[1]] #obtain log-likelihood value
    #if(is.null(xi))  loglik <- loglik + log(pXi(xiB)) #weight with prior of xi
    #if(is.null(xiFW))  loglik <- loglik + log(pXiFW(xiF)) #weight with prior of xiFW
    return(-loglik) #weight with prior of stutter.
  }
  
  #PERFORM OPTIMIZATION:
  nparam = max(mixProp)*(nC-1) + max(PHexp) + max(PHvar) + max(DEG) + max(stuttBW) + max(stuttFW) #Obtain number of parameters (over all replicates)
  nOK <- 0 #number of times for reaching largest previously seen optimum
  maxL <- -Inf #value of maximum obtained loglik
  maxITERS <- maxIter #100 #number of possible times to be INF or not valid optimum before any acceptance
  nITER <- 0 #number of times beeing INF (invalid value)
  
  logLik_tolerance = 0.01 #tolerance of accepting similar liklihood optimization
  if(verbose) print(paste0("Number of parameters to optimize: ",nparam))
  suppressWarnings({
  while(nOK<nDone) {
    
    #FIRST: generate start values (real domain): 
    phi0 = NULL
    
    if(nC>1) { #DRAW mixture proportions if more than 1 contr:
      for(id in 1:max(mixProp)) phi0 = append(phi0, drawMixprop(nC))
    }
    #DRAW PH params (PHexp,PHvar,DEG) after eachother:
    for(id in 1:max(PHexp)) {
      indid = which(PHexp==id)
      phi0 = append(phi0, drawPHvars( mean(theta_Reps[indid,1]), delta=delta ))  #draw mixture proportions:
    }
    for(id in 1:max(PHvar)) {
      indid = which(PHvar==id)
      phi0 = append(phi0, drawPHvars( mean(theta_Reps[indid,2]), delta=delta ))  #draw mixture proportions:
    }
    #DRAW DEGRADATION:
    if(any(DEG>0)) { #draw backward stutters
      for(id in 1:max(DEG)) {
        indid = which(DEG==id)
        phi0 = append(phi0, drawPHvars( mean(theta_Reps[indid,3]), delta=delta ))  #draw mixture proportions:
      }
    }
    #DRAW STUTTERS
    if(any(stuttBW>0)) { #draw backward stutters
      for(id in 1:max(stuttBW)) phi0 = append(phi0,  drawStutterProp(0.05) ) 
    }
    if(any(stuttFW>0)) { #draw forward stutters
      for(id in 1:max(stuttFW)) phi0 = append(phi0,  drawStutterProp(0.01) ) 
    } 
# phi=phi0
    #CHECK
    if(length(phi0)!=nparam) stop("The number of proposed params was not correct!") #check only
    
    # phi=phi0
    timeOneCall = system.time({ #estimate the time for calling the likelihood function one time 
      likval <- -negloglik_phi(phi=phi0)   #check if start value was accepted
    })[3] #obtain time in seconds
    
    if( is.infinite(likval) ) { #if it was infinite (invalid)
      nITER = nITER + 1	 
      
    } else { #PERFORM OPTIMIZATION
      
      tryCatch( {
        foo <- nlm(f=negloglik_phi, p=phi0, iterlim=1000,steptol=steptol, hessian=TRUE)#,print.level=2)
        Sigma <- solve(foo$hessian)
        
        if(all(diag(Sigma)>=0) && foo$iterations>2) { #} && foo$code%in%c(1,2)) { #REQUIREMENT FOR BEING ACCEPTED
          nITER <- 0 #reset INF if accepted
          likval <- -foo$min #obtain local maximum
          
          #was the maximum (approx) equal the prev: Using decimal numbers as difference 
          isEqual = !is.infinite(maxL) && abs(likval-maxL) < logLik_tolerance #all.equal(likval,maxL, tolerance = 1e-2) # # #was the maximum (approx) equal the prev?
          
          if(isEqual) { 
            if(verbose)  print(paste0("Equal maximum found: loglik=",likval))
            nOK = nOK + 1 #add counter by 1
          } else {  #if values were different
            if(likval>maxL) { #if new value is better
              nOK = 1 #first accepted optimization found
              maxL <- likval #maximized likelihood
              maxPhi <- foo$est #set as topfoo     
              maxSigma <- Sigma 
              if(verbose) print(paste0("New maximum at loglik=",likval))
            } else {
              if(verbose)  print(paste0("Local (non-global) maximum found at logLik=",likval))
            }
          } 
          if(verbose) print(paste0(" (",nOK,"/",nDone,") optimizations done"))
          #flush.console()
        } else { #NOT ACCEPTED
          nITER <- nITER + 1 
        }
      },error=function(e) e,finally = {nITER <- nITER + 1} ) #end trycatch (update counter)
    } #end if else 
    
    if(nOK==0 && nITER>maxITERS) {
      nOK <- nDone #finish loop
      maxL <- -Inf #maximized likelihood
      maxPhi <- rep(NA,nparam) #Set as NA
      maxSigma <- matrix(NA,nparam,nparam)#Set as NA
      break #stop loop if too many iterations  
    }
  } #end while loop
  })
  
  #Convert param back:
  par = convParamBack(maxPhi) #obtain list with parameters (on original scale)
  
  #Naming parameters (from replicates etc):
  names(par$PHexp) <- names(par$PHvar) <- names(par$DEG) <- repNames #one param per replicate
  contrNames = paste0("C",1:nC)
  
  #Name stutter prop names:
  outcomb1 = expand.grid(contrNames,repNames) #expand combinations
  outcomb2 = expand.grid(c("BW","FW"),repNames) #expand combinations
  sep = "-"  #separator in name
  names(par$mixProp) = paste0(outcomb1[,2],sep,outcomb1[,1])
  names(par$stutt) = paste0(outcomb2[,2],sep,outcomb2[,1])

  #STORE PARAM ESTIMATES IN A MORE PRETTY FORMAT
  mat1 = matrix(par$mixProp,nrow=nReps,byrow = TRUE, dimnames = list(repNames, contrNames)) #create table
  mat2 = cbind(PHexp=par$PHexp,PHvar=par$PHvar,DEG=par$DEG) #also add DEG param
  mat3 = matrix(par$stutt,nrow=nReps,byrow = TRUE, dimnames = list(repNames, c("BWS","FWS"))) #create table  
  par2 = cbind(mat1,mat2,mat3)

  ret <- mget(names(formals()),sys.frame(sys.nframe())) #return all values in argument call
  ret$fit=list(par=par,par2=par2,phihat=maxPhi,phiSigma=maxSigma,loglik=maxL)
  ret$samples <- ret$popFreq <- ret$refData <- NULL #remove these
  ret$data = dat #include the processed data instead
  ret$prepareC = c
  
  #Post-calculations:
  ret$logLiki = logLiki2(mlefit=ret) #calculate marker-specific results 
  
  return( ret )
  
} #end function

