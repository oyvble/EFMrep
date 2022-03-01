
#HELPFUNCTION TO TEST FUNCTIONS (ONLY CHANGED HYPOTHESIS)
#knownRef=NULL;knownRel=NULL;ibd=NULL
testFuns = function(condOrder, knownRef=NULL,knownRel=NULL,ibd=NULL,compareEFM=FALSE,steptol0=1e-3,nDone0=1,seed0=1 ) {

  expect_approx = function(x,y,tol=1e-8) { #helpfunction with specified tolerance
    expect_equal(as.numeric(x),as.numeric(y),tolerance = tol)
  }
  
  mle = contLikMLE2(NOC,samples,popFreq,refData, condOrder, knownRef, kit, ATv,pCv,lambdav,fstv, mixProp, PHexp, PHvar, DEG, stuttBW, stuttFW, knownRel=knownRel,ibd=ibd,steptol=steptol0, seed=seed0,nDone=nDone0,verbose=FALSE)
#mle$fit$loglik
#  mle$fit$par = list(mixProp=thEFM[1:NOC], PHexp=thEFM[NOC+1], PHvar=thEFM[NOC+2], DEG=thEFM[NOC+3], stutt=thEFM[NOC+4:5])
  logLikv = logLiki2(mle) #obtain per marker resutls
  expect_approx(sum(logLikv),mle$fit$loglik)
  
  #MANUAL DERIVATION:
  par = mle$fit$par#obtain parameters
  likValList = getLogLiki(samples,refData,popFreq,condOrder,knownRef, par ,NOC,ATv,pCv,lambdav,fstv, knownRel=knownRel,ibd=ibd)
  logLikv2 = log(sapply( likValList,sum)) 
  expect_approx(logLikv,logLikv2 ) 
  
#  abs(logLikv-logLikv2)>1e-3
  #Calculate deconvolution (DC) and check normalized vector for each marker:
  DC = deconvolve2(mle)
  for(loc in names(likValList)) {
    #  loc=names(likValList)[1]
    tab = DC$rankGi[[loc]]
    vals1=as.numeric(tab[,ncol(tab)])
    vals1=vals1[vals1>0]    
    
    vals2 = likValList[[loc]]
    vals2 = vals2[vals2>0] #keep only non-zero values
    vals2 = sort(vals2/sum(vals2),decreasing = TRUE)

    expect_approx( log(vals1) , log(vals2))
  }
  
  #CHECK validation data (manual derived must be close)
  #validMLEmodel2(mle)
  #valid = validMLEmodel2(mle,createplot=FALSE,alpha=0.01,verbose = FALSE)
  #expect_approx(valid$cumProb,checkVals$cumProb)
  
  #Calc WIth EFM
  if(compareEFM) {
    locs = intersect(names(popFreq),names(samples[[1]]))
    dat = prepareData(samples,refData,popFreq[locs],threshT = AT,normalize = FALSE,minF=NULL)
    #plotEPG2(dat$samples,dat$refData ,kit = kit,AT = AT)
  
    #samples=dat$samples;popFreq=dat$popFreq;refData=dat$refData;knownRef=NULL;xi=NULL;prC=pC;threshT=AT;pXi=function(x)1;delta=1;kit=kit;verbose=TRUE;maxIter=100;knownRel=NULL;ibd=c(1,0,0);xiFW=NULL;pXiFW=function(x)1;seed=NULL;maxThreads=32;steptol=steptol0
    mleEFM = contLikMLE(NOC,dat$samples,dat$popFreq, dat$refData , condOrder,knownRef=NULL,xi = NULL,pC,nDone=nDone0, AT,  fst, lambda, xiFW = NULL, kit=kit,steptol = steptol0, seed=seed0,verbose=FALSE)
  #  thEFM = mleEFM$fit$thetahat2
 #   mleEFM$fit$thetahat = unlist(par)[-NOC]
    #mle$fit$par2
    logLikvEFM = logLiki(mleEFM)
    #sum(logLikvEFM)
    expect_approx( mleEFM$fit$loglik,mle$fit$loglik, 1e-6 )
    expect_approx( logLikvEFM, logLikv, 1e-5 )
  }
}
