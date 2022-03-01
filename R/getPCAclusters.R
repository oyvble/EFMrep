#' @title getPCAclusters
#' @author Oyvind Bleka
#' @description Perform PCA on samples
#' @param X summed PH per marker per sample
#' @param grps Indicate which group each sample belongs to 
#' @param txt Title name
#' @param getClusters Whether clusters should be estimated
#' @param alpha Coverage grade for ellipsoid of each cluster 
#' @return ret Cluster classification of each sample
#' @export 

getPCAclusters = function(X,grps=NULL, txt="",getClusters=TRUE, alpha=0.1) {
  if(nrow(X)<2) stop("Less than 2 samples observed. No PCA carried out!")
  cols = c("black","#df462a","#466791","#60bf37","#953ada","#5a51dc","#d49f36","#507f2d","#db37aa","#5b83db","#c76c2d","#552095","#82702d","#55baad","#62aad3","#8c3025","#417d61")
  
  #Organize groups
  if(is.null(grps)) {
    grps = rep(1,nrow(pred0)) 
  } else {
    grps[is.na(grps)] = 0 #set this group as zero
    if(any(grps==0)) grps <- grps + 1 #increment group to have index at least 1
  }
  unGrps = unique(grps)     
  
  #Performing PCA to discriminate samples
  pca0 <- prcomp(X, center = TRUE, scale. = TRUE) #features as columns (don't scale)
  #pca0 <- prcomp(X, center = FALSE, scale. = FALSE) #features as columns (don't scale)
  # plot(pca0,ty="l")
  propV <- pca0$sdev^2/sum(pca0$sdev^2)
  pred0 <- predict(pca0)[,1:2]
  
  rngX = range(pred0[,1])
  rngY = range(pred0[,2])
  deltaX = diff(rngX)*0.1
  deltaY = diff(rngY)*0.1
  limX = rngX + c(-1,1)*deltaX
  limY = rngY + c(-1,1)*deltaY
  
  giveAX <- function(a) paste0("PC ",a,"(",round(propV[a]*100,2),"%)")
  plot(pred0[,1],pred0[,2],main=txt,xlab=giveAX(1),ylab=giveAX(2),xlim=limX,ylim=limY)
  for(g in 1:length(unGrps)) { #for each group
    ind = grps==unGrps[g]
    points(pred0[ind,1],pred0[ind,2],col=cols[g]) 
    text(pred0[ind,1],pred0[ind,2],labels=rownames(pred0)[ind],col=cols[g],cex=0.8,pos = 1)#,pch=1)   
  }
  if(length(unGrps)>1) legend("topright",legend=unGrps,col=cols[1:length(unGrps)],pch=1,cex=0.8)
  
  #obtain clustered groups:
  if(!getClusters) return() #no clustering provided
  
  if(!require(mclust)) stop("Cluster method not installed.")
  
  #BIC <- mclustBIC(X)

  clust = mclust::Mclust(pred0,verbose = FALSE)
  #clust$classification
  #clust = mclust::Mclust(X)
  nClusters = clust$G #estimated clusters
  mtext(paste("Estimated number of clusters:",nClusters))
  par = clust$parameters #obtain estimated params
  
  #plot(clust, what = "classification",addEllipses = TRUE)
  
  classes = clust$classification #get predicted classes
  names(classes) = rownames(X)
  for(g in 1:nClusters) {
    reps = paste0( names(classes)[classes==g], collapse="/")
    print(paste0("Cluster ",g,": ", reps ))
    
    #Adding uncertainty elipse to plot
    MU0 = par$mean[,g] #obtain mean
    SIGMA0 = par$variance$sigma[,,g] #obtain within-group covariance
    if(require(mixtools)) {
      mixtools::ellipse(MU0, SIGMA0,alpha = alpha,col="gray")
    }
  }
  return(clust$classification) #return cluster classifications
} #nend function