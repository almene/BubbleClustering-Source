library(Rcpp)
library(RcppProgress)
library(mclust)

#' Perform hierarchical clustering using bubble clustering methods.  Produces an hclust object
#'
#' @param data A numerical matrix
#' @param level Level of desired bubble coverage generates n*10^level bubbles are generated
#' @param Brad maximum bubble radius in the range of 0-1 of maximum distance between points
#' @param verbose Boolean indicator of progress bar presence or absence
#' @return hclust object of hierarcical tree of \code{data}
#' @examples
#' bubbleClustering(myData, 4, .75)
#' bubbleClustering(myData)


bubbleClustering <-
  function(data, level=3, Brad=1, verbose=FALSE){
    ##if (RewardType=="membership") RT=0
    ##else if (RewardType=="radius") RT=1
   ## else (return(NULL))

    if (is.matrix(data)==FALSE){
      message("The data needs to be a numeric matrix")
      stop()
      }
    if (is.numeric(data)==FALSE){
      message("The data needs to be a numeric matrix")
      stop()
    }

    samples<-c(dim(data)[1])
    features<-c(dim(data)[2])
    output<-.Call("bubble", data, samples, features, level, Brad, verbose, PACKAGE = "BubbleClustering")
    test<-.Call("ATree", output, PACKAGE = "BubbleClustering")

    result<-structure(list(merge = test[[1]], height = test[[2]], order = test[[3]],
                           labels = rownames(data), method = "Bubble",
                           dist.method = "euclidean"),
                      class = "hclust")
    return(result)
  }

  #' Perform hierarchical clustering using bubble clustering methods with sequential development of the associator matrix without the use of random values.  Produces an hclust object
#'
#' @param data A numerical matrix
#' @param level Level of desired bubble coverage generates n*10^level bubbles are generated
#' @param Brad maximum bubble radius in the range of 0-1 of maximum distance between points
#' @param verbose Boolean indicator of progress bar presence or absence
#' @return hclust object of hierarcical tree of \code{data}
#' @examples
#' bubbleOrder(myData, 4, .75)
#' bubbleOrder(myData)

bubbleOrder <-
  function(data, level=3, verbose=FALSE){
    ##if (RewardType=="membership") RT=0
    ##else if (RewardType=="radius") RT=1
    ## else (return(NULL))

    if (is.matrix(data)==FALSE){
      message("The data needs to be a numeric matrix")
      stop()
    }
    if (is.numeric(data)==FALSE){
      message("The data needs to be a numeric matrix")
      stop()
    }

    samples<-c(dim(data)[1])
    features<-c(dim(data)[2])
    output<-.Call("bubbleOrder", data, samples, features, level, verbose, PACKAGE = "BubbleClustering")
    test<-.Call("ATree", output, PACKAGE = "BubbleClustering")

    result<-structure(list(merge = test[[1]], height = test[[2]], order = test[[3]],
                           labels = rownames(data), method = "Bubble",
                           dist.method = "euclidean"),
                      class = "hclust")
    return(result)
  }


#' Calculate the average snip-rebuild distance using bubble clustering
#'
#' @param data A numerical matrix
#' @param level Level of desired bubble coverage generates n*10^level bubbles are generated
#' @param Brad maximum bubble radius in the range of 0-1 of maximum distance between points
#' @return Returns the average snip-rebuild distance for each sample in \code{data} using bubble clustering
#' @examples
#' BubStab(myData, 4, .75)
#' BubStab(myData)
#'
BubStab<- function(data,level=4, Brad=1 ,verbose=FALSE, Type="original"){ ## data is the data and tree is the Bubble output
  if (is.matrix(data)==FALSE){
    message("The data needs to be a numeric matrix")
    stop()
  }
  if (is.numeric(data)==FALSE){
    message("The data needs to be a numeric matrix")
    stop()
  }
if (Type=="original"){
  tree<-bubbleClustering(data,level, Brad, verbose)
  merge<-tree$merge
  totalStab=0
  LD<-length(data[1,])
  LP<-length(data[,1])
  for (i in 1:LP){
    if (i==1) {rbldDat<-data[c(2:LP),]}
    else if(i==LP) {rbldDat<-data[c(1:(LP-1)), ]}
    else {rbldDat<-data[c(1:(i-1), (i+1):LP),]} ## make the appropriate merge tables with the target data point removed
    RBub<-bubbleClustering(rbldDat,level)

    Sdata<-.Call("snip", merge, i)
    Rdata<-.Call("rebuild", RBub$merge)

    runStab<-sum(sqrt((Sdata- Rdata)^2))
    totalStab<-totalStab+runStab
  }
  totalStab<-totalStab/LP
  return(totalStab)}

  if (Type=="ordered"){
    tree<-bubbleOrder(data,level, verbose)
    merge<-tree$merge
    totalStab=0
    LD<-length(data[1,])
    LP<-length(data[,1])
    for (i in 1:LP){
      if (i==1) {rbldDat<-data[c(2:LP),]}
      else if(i==LP) {rbldDat<-data[c(1:(LP-1)), ]}
      else {rbldDat<-data[c(1:(i-1), (i+1):LP),]} ## make the appropriate merge tables with the target data point removed
      RBub<-bubbleOrder(rbldDat,level, verbose)

      Sdata<-.Call("snip", merge, i)
      Rdata<-.Call("rebuild", RBub$merge)

      runStab<-sum(sqrt((Sdata- Rdata)^2))
      totalStab<-totalStab+runStab
    }
    totalStab<-totalStab/LP
    return(totalStab)
  }
  message("The type need to be either ordered or original")
  stop()
}


#' Calculate the snip-rebuild distance using bubble clustering for each taxa and returns the values as a vector
#'
#' @param data A numerical matrix
#' @param level Level of desired bubble coverage generates n*10^level bubbles are generated
#' @param Brad maximum bubble radius in the range of 0-1 of maximum distance between points
#' @return Returns the snip-rebuild distance for each sample in \code{data} using bubble clustering as a vector
#' @examples
#' BubStab(myData, 4, .75)
#' BubStab(myData)
#'
BubStabHist<- function(data,level=4, Brad=1 ,verbose=FALSE, Type="original"){ ## data is the data and tree is the Bubble output
  if (is.matrix(data)==FALSE){
    message("The data needs to be a numeric matrix")
    stop()
  }
  if (is.numeric(data)==FALSE){
    message("The data needs to be a numeric matrix")
    stop()
  }
  if (Type=="original"){
    tree<-bubbleClustering(data,level, Brad, verbose)
    merge<-tree$merge
    totalStab=0
    LD<-length(data[1,])
    LP<-length(data[,1])
    for (i in 1:LP){
      if (i==1) {rbldDat<-data[c(2:LP),]}
      else if(i==LP) {rbldDat<-data[c(1:(LP-1)), ]}
      else {rbldDat<-data[c(1:(i-1), (i+1):LP),]} ## make the appropriate merge tables with the target data point removed
      RBub<-bubbleClustering(rbldDat,level)

      Sdata<-.Call("snip", merge, i)
      Rdata<-.Call("rebuild", RBub$merge)

      runStab<-sum(sqrt((Sdata- Rdata)^2))
      totalStab<-c(c(totalStab),runStab)
    }
    totalStab<-totalStab[2:(LP+1)]
    return(totalStab)}

  if (Type=="ordered"){
    tree<-bubbleOrder(data,level, verbose)
    merge<-tree$merge
    totalStab=0
    LD<-length(data[1,])
    LP<-length(data[,1])
    for (i in 1:LP){
      if (i==1) {rbldDat<-data[c(2:LP),]}
      else if(i==LP) {rbldDat<-data[c(1:(LP-1)), ]}
      else {rbldDat<-data[c(1:(i-1), (i+1):LP),]} ## make the appropriate merge tables with the target data point removed
      RBub<-bubbleOrder(rbldDat,level, verbose)

      Sdata<-.Call("snip", merge, i)
      Rdata<-.Call("rebuild", RBub$merge)

      runStab<-sum(sqrt((Sdata- Rdata)^2))
      totalStab<-c(c(totalStab),runStab)
    }
    totalStab<-totalStab[2:(LP+1)]
    return(totalStab)
  }
  message("The type need to be either ordered or original")
  stop()

}


#' Calculate the average snip-rebuild distance using hclust()
#'
#' @param data A numerical matrix
#' @param method Desired methodology for joining clusters in hclust() and as such should be  one of "ward.D", "ward.D2", "single", "complete", "average" (= UPGMA), "mcquitty" (= WPGMA), "median" (= WPGMC) or "centroid" (= UPGMC).
#' @return Returns the average snip-rebuild distance for each sample in \code{data} using hclust()
#' @examples
#' BubStab(myData, "complete")
#' BubStab(myData, method= "ward.D")
#'
HClustStab<- function(data, method="complete"){ ## data is the data and tree is the Bubble output
  Distdata<-dist(data)
  tree<-hclust(Distdata, method)
  merge<-tree$merge
  totalStab=0
  LD<-length(data[1,])
  LP<-length(data[,1])
  for (i in 1:LP){
    if (i==1) {rbldDat<-data[c(2:LP), ]}
    else if(i==LP) {rbldDat<-data[c(1:(LP-1)), ]}
    else {rbldDat<-data[c(1:(i-1), (i+1):LP),]} ## make the appropriate merge tables with the target data point removed
    RDist<-dist(rbldDat)
    RClust<-hclust(RDist, method)

    Sdata<-.Call("snip", merge, i)
    Rdata<-.Call("rebuild", RClust$merge)

    runStab<-sum(sqrt((Sdata- Rdata)^2))
    totalStab<-totalStab+runStab
  }
  totalStab<-totalStab/LP
  return(totalStab)
}

#' Calculate the snip-rebuild distance using hclust() for each taxa and returns the value as a vector
#'
#' @param data A numerical matrix
#' @param method Desired methodology for joining clusters in hclust() and as such should be  one of "ward.D", "ward.D2", "single", "complete", "average" (= UPGMA), "mcquitty" (= WPGMA), "median" (= WPGMC) or "centroid" (= UPGMC).
#' @return Returns the snip-rebuild distance for each sample in \code{data} using hclust() as a vector
#' @examples
#' BubStab(myData, "complete")
#' BubStab(myData, method= "ward.D")
#'
HClustStabHist<- function(data, method="complete"){ ## data is the data and tree is the Bubble output
  Distdata<-dist(data)
  tree<-hclust(Distdata, method)
  merge<-tree$merge
  totalStab=0
  LD<-length(data[1,])
  LP<-length(data[,1])
  for (i in 1:LP){
    if (i==1) {rbldDat<-data[c(2:LP), ]}
    else if(i==LP) {rbldDat<-data[c(1:(LP-1)), ]}
    else {rbldDat<-data[c(1:(i-1), (i+1):LP),]} ## make the appropriate merge tables with the target data point removed
    RDist<-dist(rbldDat)
    RClust<-hclust(RDist, method)

    Sdata<-.Call("snip", merge, i)
    Rdata<-.Call("rebuild", RClust$merge)

    runStab<-sum(sqrt((Sdata- Rdata)^2))
    totalStab<-c(c(totalStab),runStab)
  }
  totalStab<-totalStab[2:(LP+1)]
  return(totalStab)
}


treeDist<-function(tree1, tree2){
  test<-dim(tree1$merge)-dim(tree2$merge)
  test<-abs(sum(test))
  if (test>0){
    message("The trees need to be the same size")
    stop()
  }

  MCC1<-.Call("rebuild", tree1$merge)
  MCC2<-.Call("rebuild", tree2$merge)

  distT<-sum(sqrt((MCC1- MCC2)^2))

  return(distT)
}

treeDistAdd<-function(Ltree, Stree, addP){
  test<-dim(Ltree$merge)-dim(Stree$merge)
  test<-abs(sum(test)-1)
  if (test>0){
    message("The trees need to in the correct order and only differ by one leaf")
    stop()
  }
#note make check for leaf being present in sample
  MCC1<-.Call("snip", Ltree$merge, addP)
  MCC2<-.Call("rebuild", Stree$merge)

  distT<-sum(sqrt((MCC1- MCC2)^2))

  return(distT)
}

noIOhc<-
  function(data){

    if (is.matrix(data)==FALSE){
      message("The data needs to be a numeric matrix")
      stop()
    }
    if (is.numeric(data)==FALSE){
      message("The data needs to be a numeric matrix")
      stop()
    }
    distMat<-as.matrix(dist(data))
    distMat<-1/distMat*max(distMat)
    samples<-dim(distMat)[1]
    for (i in 1:samples){distMat[i,i]=0.000001}


    test<-.Call("ATree", distMat, PACKAGE = "BubbleClustering")

    result<-structure(list(merge = test[[1]], height = test[[2]], order = test[[3]],
                           labels = rownames(data), method = "noIOhc",
                           dist.method = "euclidean"),
                      class = "hclust")
    return(result)
  }

noIOhcStab<-
  function(data){
      tree<-noIOhc(data)
      merge<-tree$merge
      totalStab=0
      LD<-length(data[1,])
      LP<-length(data[,1])
      for (i in 1:LP){
        if (i==1) {rbldDat<-data[c(2:LP), ]}
        else if(i==LP) {rbldDat<-data[c(1:(LP-1)), ]}
        else {rbldDat<-data[c(1:(i-1), (i+1):LP),]} ## make the appropriate merge tables with the target data point removed
        RClust<-noIOhc(rbldDat)

        Sdata<-.Call("snip", merge, i)
        Rdata<-.Call("rebuild", RClust$merge)

        runStab<-sum(sqrt((Sdata- Rdata)^2))
        totalStab<-totalStab+runStab
      }
      totalStab<-totalStab/LP
      return(totalStab)
    }

