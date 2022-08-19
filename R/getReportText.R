#' @title getReportText
#' @author Oyvind Bleka
#' @description Helpfunction to get report text in efm2
#' @param set fitted MLE object (includes hd and possibly also hp)
#' @param envir environment object 
#' @param sig significant numbers present for numeric
#' @param colps column separator to use
#' @return txt return text of report
#' @export 

#library(EFMreps)
#load("C:\\Users\\oyvbl\\Dropbox\\Forensic\\MixtureProj\\myDev\\EFMrepsDEV\\example3\\projCalc.Rdata")
#envir=mmTK;type="EVID";sig=2;colps="\t"
#txt=getReportText(type="EVID",envir=mmTK)
getReportText = function(set, envir,sig=2,colps="\t") {
  
  #HELPFUNCTIONS:
  getws = function(len) paste0( rep(" ",len),collapse="") #max whitespace
  
  #Helpfunction to print a table
  printTable = function(tab) {
    tab2 = signif(tab,sig)
    rowNames = rownames(tab) #obtain rownames
    colNames = colnames(tab) #obtain rownames
    
    txt = paste0(  getws(max(nchar(rowNames))+1) , paste0(colNames, collapse=getws(sig)))
    for(j in 1:nrow(tab2)) txt = paste0(txt,"\n",rowNames[j]," ", paste0(tab2[j,],collapse=" ") )
    return(txt)
  }
  
  printMLE <- function(mlefit,hyp) {
    mle <- mlefit$fit
    txt <- paste0("\n\n-------Parameter estimates under ",hyp,"---------\n")
    txt <- paste0(txt,printTable(mle$par2)) #obtain param table
    txt <- paste0(txt,"\nlogLik=",round(mle$loglik,sig)) #store loglik val
    return(txt)
  }
  
  printSET <- function(mlefit) { #print settings used
    txt <- paste0("\n#Markers=",paste0( length(mlefit$prepareC$locNames) ,collapse="/"))
    txt <- paste0(txt,"\nReplicate(s)=",paste0(mlefit$prepareC$repNames,collapse="/"))
    txt <- paste0(txt,"\nKit(s)=",paste0(mlefit$kit,collapse="/"))
    txt <- paste0(txt,"\nAnalytical threshold=",paste0(mlefit$AT,collapse="/"))
    txt <- paste0(txt,"\nProbability of drop-in=",paste0(mlefit$pC,collapse="/"))
    txt <- paste0(txt,"\nHyperparam lambda=",paste0(mlefit$lambda,collapse="/"))
    txt <- paste0(txt,"\nFst-correction (global)=",paste0(mlefit$fst,collapse="/"))
    return(txt)
  }
  
  printMOD <- function(mlefit,hyp) { #print refs
    txt <- paste0("\n\n-------Hypothesis ",hyp,"---------")
    txt <- paste0(txt,"\nNumber of contributors: ",mlefit$nC) #Number of contributors
    refNames = names(mlefit$data[[1]]$refs) #unique( unlist( lapply( mlefit$data, function(x) names(x$refs) ) ) ) #obtain reference names (use first locus only?)
    txt <- paste0(txt,"\nKnown contributors: ",paste0(refNames[which(mlefit$condOrder>0)],collapse="/")) #conditional references
    if(length(mlefit$knownRef)) txt <- paste0(txt,"\nKnown non-contributors: ",paste0(refNames[mlefit$knownRef],collapse="/")) #conditional references
    #if(length(mlefit$knownRel)) txt <- paste0(txt,"\nAssumed relationship: 1st Unknown is a ",names(mlefit$ibd)[1]," to reference ",names(mlefit$knownRel)) #Relationship
    return(txt)
  }
  
  printFREQ <- function(dat) { #print freqs
    locs = names(dat)
    txt <- paste0("\n\n-------Frequency data---------")
    for(loc in locs) {
      freq = dat[[loc]]$freq
      tmp <- paste0(names(freq),"=",freq) #get allele names with freqs
      txt <- paste0(txt,"\n",loc,": ",paste0(tmp,collapse="/"))
    }
    return(txt)
  }
  
  txt <- paste0("EFMrep version ",packageVersion("EFMrep"))
  txt <- paste0(txt,"\nR-version: ",R.version.string) #Add R-version used
  txt <- paste0(txt,"\nCreated: ",Sys.time())
  txt <- paste0(txt,"\nUser: ",Sys.getenv("USERNAME"),"\n")
  
  txt <- paste0(txt,"\n-------Settings-------")
  txt <- paste0(txt,"\nPopulation frequency file: ", basename(get("Freqfile",envir=envir)) ) #store population freq file used (only base name)
  
  txt <-  paste0(txt,printSET(set$mlefit_hd)) #Print Data and model options under Hd
  
  #Hypotheses:  
  if(!is.null(set$mlefit_hp))  txt <- paste0(txt,printMOD(mlefit=set$mlefit_hp,hyp="Hp")) #Print hypothesis Hp:
  if(!is.null(set$mlefit_hd))  txt <- paste0(txt,printMOD(mlefit=set$mlefit_hd,hyp="Hd")) #Print hypothesis Hd:
  
  
  #Obtain and show LR results
  res = set$resLR  
  if(!is.null(res) ) {
    txt0 <- paste0("LR (MLE)=",signif(res$LRmle),sig)
    txt1 <- paste0("log10LR (MLE)=",signif(log10(res$LRmle)),sig)
    txt3 <- paste0("log10LR (Upper boundary)=",signif(log10(res$LRupper),sig)) 
    
    txt5 <- paste0(paste0(names(res$LRi),colps,signif(res$LRi,sig)),collapse="\n")
    txt <- paste0(txt,"\n\n-------LR (all markers)------\n",txt0,"\n",txt1,"\n",txt3,"\n")
    txt <- paste0(txt,"\n-------LR (per marker)------\n",txt5,"\n")
    
    #Obtain and insert NUMBER OF FAILED PP-plot points outside envelope: 
    alpha <- 0.01 #as.numeric(getValueUser("Set significance level \nin model validation:",0.01))
    #checkPositive(alpha,"The significance level",strict=TRUE)
#    mlefit=set$mlefit_hp
    validHp = validMLEmodel2(set$mlefit_hp,alpha=alpha,createplot = FALSE,verbose=FALSE)
    validHd = validMLEmodel2(set$mlefit_hd,alpha=alpha,createplot = FALSE,verbose=FALSE)
    nFailedHp = sum(validHp$Significant)
    nFailedHd = sum(validHd$Significant)
    txtValid = paste0("Under H",c("p","d"),": ",c(nFailedHp,nFailedHd))
    
    txt <- paste0(txt,"\n-------Model validation------\nNumber of fails (signif level=",alpha,"):\n",txtValid[1],"\n",txtValid[2])
  }
  #store consLR - estimate
 
  #store parameter estimates
  if(!is.null(set$mlefit_hp)) txt <- paste0(txt, printMLE(set$mlefit_hp,"Hp"))
  if(!is.null(set$mlefit_hd)) txt <- paste0(txt, printMLE(set$mlefit_hd,"Hd"))
  
  txt <- paste0(txt,"\n\n-------Optimalisation setting-------")
  txt <- paste0(txt,"\nRequired number of (identical) optimizations: ",set$mlefit_hd$nDone) 
  txt <- paste0(txt,"\nAccuracy of optimisations (steptol): ",set$mlefit_hd$steptol) 
  txt <- paste0(txt,"\nSeed for optimisations: ", ifelse(is.null(set$mlefit_hd$seed),"NONE",set$mlefit_hd$seed)) 
  
  
  #Print allele freqs last: ADDED in v2.0.1 (notice that only Hd is printed)
  txt <-  paste0(txt,printFREQ(dat=set$mlefit_hd$data)) 
  optFreq = get("optFreq",envir=envir)
  minFreq = optFreq$minF
  if(optFreq$freqsize>0) {
    txt <-  paste0(txt,"\nSize of frequency database (N): ",optFreq$freqsize)  
  } else {
    txt <-  paste0(txt,"\nRare allele frequency (minFreq): ",optFreq$minF)  
  }
  txt <-  paste0(txt,"\nNormalized after impute: ", ifelse(optFreq$normalize==1,"Yes","No") ) 
  txt <- as.matrix(txt)
  
  colnames(txt) <- paste0("This is a generated report from")
  return(txt)
} #end savetableALL
  
