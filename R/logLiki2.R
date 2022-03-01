#' @title logLiki2
#' @author Oyvind Bleka
#' @description logLiki returns the likelihood for each markers (based on the MLE results)
#'
#' @param mlefit Fitted object using contLikMLE
#' @return ret A vector with log-likelihood-values for each locus for given model
#' @export

#Get value of likelihood for each marker given model fit.
#mlefit = mle;maxThreads=128
logLiki2 <- function(mlefit){
  par = mlefit$fit$par
  c <- mlefit$prepareC #returned from prepareC
  NOC = c$NOC #number of contrs
  nReps = c$nReps #number of repeats
  nM = c$nLocs #number of markers
  
  logLikv = rep(NA,nM) #obtain loglik per marker
  names(logLikv) = c$locNames #name marker names
  
  if(any(is.na(unlist(par)))) return(logLikv) #return if params not found
  
  #Step 1) Calculate L(E|thetahat) for each marker
  for(m in 1:nM) { #extract info in c relevant for each markers:
#    print(m)
    NOCrange = NOC*(m-1) + 1:NOC #obtain NOC indices: Is range 1:NOC per marker
    logLikv[m] = .C("loglikGamma_allcomb1",as.numeric(0), c$nJointCombs[m], c$NOC, c$NOK[m],
              as.numeric(par$mixProp),  as.numeric(par$PHexp), as.numeric(par$PHvar), as.numeric(par$DEG), as.numeric(par$stutt), 
              as.numeric(c$AT),as.numeric(c$fst[m]),as.numeric(c$dropinProb),as.numeric(c$dropinWeight),
              c$nReps, as.integer(1), c$nRepMarkers[m], c$nAlleles[m], c$nAlleles2[m], c$startIndMarker_nAlleles[m], c$startIndMarker_nAllelesReps[m],
              c$peakHeights, c$freq, c$nTyped[m], c$maTyped, c$basepairs, c$repID, c$startIndMarker_nRepMarkers[m], 
              c$nGenos[m], c$outG1allele, c$outG1contr, c$startIndMarker_outG1allele[m], c$startIndMarker_outG1contr[m], c$startIndMarker_nJointCombs[m],
              c$nStutters[m], c$stuttFromInd, c$stuttToInd, c$stuttParamInd , c$startIndMarker_nStutters[m],
              c$knownGind[NOCrange], as.integer(mlefit$maxThreads), c$relGind[NOCrange], c$ibdLong )[[1]]
  }    
        
 return( logLikv ) #simply returns  
}



