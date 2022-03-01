#' @title validMLEmodel2
#' @author Oyvind Bleka
#' @description validMLEmodel makes model validation whether the observed peak heights fits the maximum likelihood fitted gamma distribution.
#' @details The cumulative probability of the observed allele peaks are calculated and compared with a uniform distribution.
#' Function calls a 'cumulative procedure' in C++ which uses the package Boost and performs paralellisation (OpenMP).
#'
#' @param mlefit Fitted object using contLikMLE
#' @param plottitle Maintitle text used in the PP-plot
#' @param alpha The significance level used for the envelope test. Default is 0.01
#' @param createplot Boolean of whether plot should be created
#' @param verbose Boolean of whether printing out information about significant alleles
#' @return retinfo A dataframe with information from the calculated validation model
#' @export

#mlefit=mleHp; plottitle="PP-plot";alpha=0.01;createplot=TRUE;verbose=TRUE
validMLEmodel2 <- function(mlefit,plottitle="PP-plot",alpha=0.01,createplot=TRUE,verbose=TRUE) {
  cols = c("black","#df462a","#466791","#60bf37","#953ada","#5a51dc","#d49f36","#507f2d","#db37aa","#5b83db","#c76c2d","#552095","#82702d","#55baad","#62aad3","#8c3025","#417d61")
  xlab="Expected probabilities"#: Unif(0,1)"
  ylab="Observed probabilities"#: (Pr(Yj<=yj|Y_{-j}<=y_{-j},Yj>=thresh,model))"
  sz <- 1.5
  
  pchmax = 26 #number of possible point types
  c <- mlefit$prepareC #Object already stored in mlefit. returned from prepareC function
  par = mlefit$fit$par #estimated params

  if(any(is.na(unlist(par)))) return(NULL) #return if params not found
  
  #Obtaining number of alleles to obtain an upper peak height boundary to integrate up to (maxY)
  nAtotObs = sum(c$peakHeights>0)
  alphaQ <- 0.001 #ensure very far out in quantile (used for estimating probs in gamma-distribution).
  alpha2 <- alphaQ/nAtotObs #"bonferroni outlier"
  suppressWarnings({ #don't show missing allele warning
    maxYobs <- c$peakHeights #max observation
    maxYexp <- qgamma(1-alpha2,2/par$PHvar^2,scale=par$PHexp*par$PHvar^2) #max observation in theory
    maxY <- ceiling( max( na.omit(c(maxYobs,maxYexp))) ) #get max observed (scaled up)
  }) 
  
  #PArt 1/2: Consider evaluation at PH-levels  
  dropinWeight1 <- dropinWeight2 <- c$dropinWeight #modify dropin weights to account for cdf instead of pdf
  for(locind in 1:c$nLocs) { #traverse each marker
    for(rind in 1:c$nRepMarkers[locind]) { 
      lambda0 = c$lambda[ c$nRepMarkers[locind] + rind]
      AT0 =  c$AT[ c$nRepMarkers[locind] + rind]
      for(aind in 1:c$nAlleles[locind]) {
        cind = c$startIndMarker_nAllelesReps[locind] + (aind-1)*c$nRepMarkers[locind] + rind
        peak = c$peakHeights[cind]
        if(peak < AT0) next #skip if less than zero
        if(peak==AT0) peak = peak + 1 #avoid Inf in dropin model (equivalent with 'AT-1' threshold)
        freq = c$freq[ c$startIndMarker_nAlleles[locind] + aind] #obtain frequency
        
        dropinWeight1[cind] = log(freq) + pexp(peak-AT0,lambda0,log.p=TRUE) #= log( 1 - exp(- (lambda0)*(peak-AT0)))
        dropinWeight2[cind] = log(freq) + pexp(maxY-AT0,lambda0,log.p=TRUE) #= log(freq) + 1
	   }
    }
  }

  #obtain cumulative vals int_0^y 
  pvalVEC = c$peakHeights #must have a copy
  UaPH = .C("loglikGamma_cumprob",as.numeric(pvalVEC), as.numeric(0), c$nJointCombs, c$NOC, c$NOK,
                 as.numeric(par$mixProp),  as.numeric(par$PHexp), as.numeric(par$PHvar), as.numeric(par$DEG), as.numeric(par$stutt), 
                 as.numeric(c$AT),as.numeric(c$fst),as.numeric(c$dropinProb),as.numeric(c$dropinWeight), as.numeric(dropinWeight1),
                 c$nReps, c$nLocs, c$nRepMarkers, c$nAlleles, c$nAlleles2, c$startIndMarker_nAlleles, c$startIndMarker_nAllelesReps,
                 c$peakHeights, c$freq, c$nTyped, c$maTyped, c$basepairs, c$repID, c$startIndMarker_nRepMarkers, 
                 c$nGenos, c$outG1allele, c$outG1contr, c$startIndMarker_outG1allele, c$startIndMarker_outG1contr, c$startIndMarker_nJointCombs,
                 c$nStutters, c$stuttFromInd, c$stuttToInd, c$stuttParamInd , c$startIndMarker_nStutters,
                 c$knownGind, as.integer(mlefit$maxThreads), c$relGind, c$ibdLong )[[1]] #obtain 
  #print(UaPH)
  
  #PArt 2/2: Consider evaluation at PH-levels  
  #obtain cumulative vals int_0^inf 
  pvalVEC = c$peakHeights #must have a copy
  UaMAX = .C("loglikGamma_cumprob",as.numeric(pvalVEC), as.numeric(maxY), c$nJointCombs, c$NOC, c$NOK,
                as.numeric(par$mixProp),  as.numeric(par$PHexp), as.numeric(par$PHvar), as.numeric(par$DEG), as.numeric(par$stutt), 
                as.numeric(c$AT),as.numeric(c$fst),as.numeric(c$dropinProb),as.numeric(c$dropinWeight), as.numeric(dropinWeight2),
                c$nReps, c$nLocs, c$nRepMarkers, c$nAlleles, c$nAlleles2, c$startIndMarker_nAlleles, c$startIndMarker_nAllelesReps,
                c$peakHeights, c$freq, c$nTyped, c$maTyped, c$basepairs, c$repID, c$startIndMarker_nRepMarkers, 
                c$nGenos, c$outG1allele, c$outG1contr, c$startIndMarker_outG1allele, c$startIndMarker_outG1contr, c$startIndMarker_nJointCombs,
                c$nStutters, c$stuttFromInd, c$stuttToInd, c$stuttParamInd , c$startIndMarker_nStutters,
                c$knownGind, as.integer(mlefit$maxThreads), c$relGind, c$ibdLong )[[1]] #obtain 
  #print(UaMAX)
  #Combing the product from the two parts:  
  cumprobi = UaPH/UaMAX #obtaining cumulative probs
  #print(cumprobi)
  #hist(cumprobi)
  #Obtain names:
  locNames = c$locNames #obtain locus names
  alleleNames = c$alleleNames
  repNames = c$repNames
  
  tab = NULL
  for(locind in 1:c$nLocs) { #traverse each marker
    loc0 = locNames[locind] #obtain locus name
    for(rind in 1:c$nRepMarkers[locind]) { 
      rep0 =  repNames[ c$repID[ c$startIndMarker_nRepMarkers[locind] + rind]+1] #obtain replicate name
      for(aind in 1:c$nAlleles[locind]) {
        cind = c$startIndMarker_nAllelesReps[locind] + (aind-1)*c$nRepMarkers[locind] + rind
        peak = c$peakHeights[cind]
        if(peak==0) next 
        allele = c$alleleNames[[loc0]][aind] #obtain allele name
        new = c( loc0, rep0, allele, peak, cumprobi[cind])
        tab = rbind(tab, new)
      }
    }
  }
  colnames(tab) = c("Marker","Replicate","Allele","Height","cumProb")
  #nrow(tab)==sum(!is.nan(cumprobi))
  tab = as.data.frame(tab)
  tab$cumProb = as.numeric(tab$cumProb)
  
  cumprobi = as.numeric(tab[,ncol(tab)]) #obtain values again
  N <- nrow(tab) #length(cumprobi) #number of peak heights
  alpha2 <- alpha/N #0.05 significanse level with bonferroni correction
  cumunif <- ((1:N)-0.5)/N #=punif((1:N)-0.5,0,N)
  
  #hist(cumprobi)
  #sum(cumprobi<=0.5)/N
  ord <- order(cumprobi,decreasing = FALSE)
  ord2 <- match(1:length(ord),ord) #get reverse index
  pval <- 1-pbeta(cumprobi[ord],N*cumunif,N-N*cumunif+1) #one sided p-value
  
  #Must indicate the values below the line (Symmetry)
  ind <- cumprobi[ord]<cumunif #those below the line (two-sided p-value)
  pval[ind] <- pbeta(cumprobi[ord],N*cumunif,N-N*cumunif+1)[ind] #calculate pval in opposite direction
  #cumprobi<qbeta(alpha/2,N*cumunif,N-N*cumunif+1) | cumprobi<qbeta(1-alpha/2,N*cumunif,N-N*cumunif+1)
  outside <- pval[ord2] < (alpha2/2) #criterion outside region (divide by 2 to get two-sided)
  tab = cbind(tab,pvalue=pval[ord2],Significant=outside) #update table
  
  if(verbose) { #print info:
    #print points outside the bonferroni-adjusted envolopment
    print(paste0("Total number of peak height observation=",N))
    print(paste0("Significance level=",signif(alpha*100,digits=2),"%"))
    print(paste0("Bonferroni-adjusted significance level=",signif(alpha2*100,digits=2),"%"))
    print("List of observations outside the envelope (with the Bonferroni-level):") #two sided check
    print(tab[tab$Significant,],drop=FALSE)  #print list of outliers (outside envelop)
  }
  
  #plot   
  if(createplot) {
    zones <- matrix(c(1,1,1, 2,5,4, 0,3,0), ncol = 3, byrow = TRUE)
    layout(zones, widths=c(0.3,7,1), heights = c(1,7,.75))
    par(xaxt="n", yaxt="n",bty="n",  mar = c(0,0,0,0))   # for all three titles: 
    plot(x=1,y=1,type="n",ylim=c(-1,1), xlim=c(-1,1))   # fig 1 from the layout
    text(0,0,paste(plottitle), cex=2)  
    plot(x=1,y=1,type="n",ylim=c(-1,1), xlim=c(-1,1))  # fig 2
    text(0,0,paste(ylab), cex=2, srt=90)   
    plot(x=1,y=1,type="n",ylim=c(-1,1), xlim=c(-1,1)) # fig 3
    text(0,0,paste(xlab), cex=2)  
    
    # fig 4: Rigth margin shows distr per replicate
    par(mar = c(1,0,0,0))
    
    plot(0,0,xlim=0:1,ylim=0:1,ty="n")
    for(rep in repNames) {
#rep=repNames[1]
      rind = which(repNames==rep)
      ind = tab$Replicate==rep #obtain index of certain replicate
      subDat = tab[ind,]
      locind = match(subDat$Marker,locNames) #obtain loci to plot
      
      points( rep(rind/(length(repNames)+1),sum(ind)), subDat$cumProb ,pch=(locind-1)%%pchmax,col=cols[rind],cex=sz)
    }
    rect(0,0,1,1)
    
    # fig 5, finally, the scatterplot-- needs regular axes, different margins
    par(mar = c(1,2,0,.5), xaxt="s", yaxt="s", bty="n")
    
    plot(0,0,ty="n",xlim=0:1,ylim=c(0,1),cex.axis=sz,asp=1)
    segments(x0=0,y0=0,x1=1,y1=1,lwd=1.2)
    
    #Goodness of fit test  #REMOVED: pval <- ks.test(cumprobi, "punif")$p.value
    #Updated block in 0.6.2: DRAW ENVELOPE LINES FOR ORDER STATISTICS UNDER RANDOMNESS:
    xsq <- seq(0,1,l=1000) #Draw envolope lines
    ysq <- c(alpha,alpha2) #quantiles to consider
    for(qq in ysq) { #for each quanitle
#      qq=ysq[1]
      lines(xsq,qbeta(qq/2,N*xsq,N-N*xsq+1),col=which(ysq==qq),lty=2)
      lines(xsq,qbeta(1-qq/2,N*xsq,N-N*xsq+1),col=which(ysq==qq),lty=2)
    }
    legend("bottomright",legend=paste0("1-Envelope-coverage=",c("",paste0(alpha,"/",N,"=")),signif(ysq,2)),col=1:length(ysq),lty=2,cex=1.3)
    
    #  abline(0,1)
    locind = match(tab$Marker,locNames) #obtain loci to plot
    repind = match(tab$Replicate,repNames) #get replicate indices
    points(cumunif,tab$cumProb[ord],pch=(locind[ord]-1)%%pchmax,cex=sz,col=cols[repind[ord]])
    
    legend("topleft",legend=repNames,pch=19,cex=sz,col=cols[1:length(repNames)])
    dev.new()
    op <- par(no.readonly = TRUE)
    dev.off()
    par(op)
  } #end createplot
  
  #return information
  return(tab)
  
}

