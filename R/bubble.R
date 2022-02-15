library(Rcpp)
library(RcppProgress)
library(mclust)
set.seed(2020)

#' Perform hierarchical clustering using bubble clustering methods.  Produces an hclust object
#'
#' @param data A numerical matrix
#' @param level Level of desired bubble coverage generates n*10^level bubbles are generated
#' @param Brad % maximum bubble radius in the range of 0-1 of maximum distance between points
#' @param RadiiDist Selection of radii distribution, 1=random uniform 2= quadratically favour smaller, 3=log favour smaller
#' @return hclust object of hierarcical tree of \code{data}
#' @examples
#' bubbleClustering(data, 4, .75)
#' bubbleClustering(data)


bubbleClustering <-
  function(data, level=3, Brad=1, RadiiDist=1, rewardType=1, isDist=FALSE, verbose=TRUE){
    ##if (RewardType=="membership") RT=0
    ##else if (RewardType=="radius") RT=1
    ## else (return(NULL))
    if (isDist==FALSE){
      isDistance=0
    }
    if (isDist==T){
      isDistance=1
      data<-data*data
    }
    if (verbose==FALSE){
      isVerbose=0
    }
    if (verbose==TRUE){
      isVerbose=1
    }
    if (RadiiDist<1){
      message("That is not an appropriate distribution selection")
      stop()
    }
    if (RadiiDist>5){
      message("That is not an appropriate distribution selection")
      stop()
    }
    if (rewardType<1){
      message("That is not an appropriate distribution selection")
      stop()
    }
    if (rewardType>4){
      message("That is not an appropriate distribution selection")
      stop()
    }


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
    output<-.Call("bubble", data, samples, features, level, Brad, RadiiDist, rewardType, c(isDistance, isVerbose), PACKAGE = "BubbleClusteringT")
    #message("The association between samples 1 and 2 is: ", output[1,2])
    test<-.Call("ATree", output, PACKAGE = "BubbleClusteringT")
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
#' BubStab(data, 4, .75)
#' BubStab(data)
#'
BubStab<- function(data,level=3, Brad=1 ,RadDist=1, rewardT=1, is_distance=F, verbose=TRUE){ ## data is the data and tree is the Bubble output
  if (is.matrix(data)==FALSE){
    message("The data needs to be a numeric matrix")
    stop()
  }
  if (is.numeric(data)==FALSE){
    message("The data needs to be a numeric matrix")
    stop()
  }

  tree<-bubbleClustering(data,level, Brad, RadDist,rewardT, is_distance, verbose)
  merge<-tree$merge
  totalStab=0
  LD<-length(data[1,])
  LP<-length(data[,1])
  if (is_distance==F){
    for (i in 1:LP){
      if (i==1) {rbldDat<-data[c(2:LP),]}
      else if(i==LP) {rbldDat<-data[c(1:(LP-1)), ]}
      else {rbldDat<-data[c(1:(i-1), (i+1):LP),]} ## make the appropriate merge tables with the target data point removed
      RBub<-bubbleClustering(rbldDat,level, Brad, RadDist,rewardT, is_distance, verbose)
      Sdata<-.Call("snip", merge, i)
      Rdata<-.Call("rebuild", RBub$merge)
      runStab<-sqrt(sum((Sdata- Rdata)^2))
      totalStab<-totalStab+runStab
    }
  }
  if (is_distance==T){
    for (i in 1:LP){
      if (i==1) {rbldDat<-data[c(2:LP),c(2:LP)]}
      else if(i==LP) {rbldDat<-data[c(1:(LP-1)),c(1:(LP-1)) ]}
      else {rbldDat<-data[c(1:(i-1), (i+1):LP),c(1:(i-1), (i+1):LP)]} ## make the appropriate merge tables with the target data point removed
      RBub<-bubbleClustering(rbldDat,level, Brad, RadDist,rewardT, is_distance, verbose)
      Sdata<-.Call("snip", merge, i)
      Rdata<-.Call("rebuild", RBub$merge)
      runStab<-sqrt(sum((Sdata- Rdata)^2))
      totalStab<-totalStab+runStab
    }
  }
      totalStab<-totalStab/LP
      return(totalStab)
}



    #' Calculate the snip-rebuild distance using bubble clustering for each taxa and returns the values as a vector
    #'
    #' @param data A numerical matrix
    #' @param level Level of desired bubble coverage generates n*10^level bubbles are generated
    #' @param Brad maximum bubble radius in the range of 0-1 of maximum distance between points only relevant for uniform radii
    #' @return Returns the snip-rebuild distance for each sample in \code{data} using bubble clustering as a vector
    #' @examples
    #' BubStabHist(data, 4, .75)
    #' BubStabHist(data)
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

          runStab<-sqrt(sum((Sdata- Rdata)^2))
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

          runStab<-sqrt(sum((Sdata- Rdata)^2))
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
    #' HClustStab(data, "complete")
    #' HClustStab(data, method= "ward.D")
    #'
    HClustStab<- function(data, method="complete", is_dist=F){ ## data is the data and tree is the Bubble output
      if (is_dist==F){
        Distdata<-dist(data)
      }
      else{
        Distdata<-as.dist(data)
      }

      tree<-hclust(Distdata, method)
      merge<-tree$merge
      totalStab=0
      LD<-length(data[1,])
      LP<-length(data[,1])
      if (is_dist==F){
        for (i in 1:LP){
          if (i==1) {rbldDat<-data[c(2:LP), ]}
          else if(i==LP) {rbldDat<-data[c(1:(LP-1)), ]}
          else {rbldDat<-data[c(1:(i-1), (i+1):LP),]} ## make the appropriate merge tables with the target data point removed
          RDist<-dist(rbldDat)
          RClust<-hclust(RDist, method)

          Sdata<-.Call("snip", merge, i)
          Rdata<-.Call("rebuild", RClust$merge)

          runStab<-sqrt(sum((Sdata- Rdata)^2))
          totalStab<-totalStab+runStab
        }
      }
      if (is_dist==T){
        for (i in 1:LP){
          if (i==1) {rbldDat<-data[c(2:LP),c(2:LP) ]}
          else if(i==LP) {rbldDat<-data[c(1:(LP-1)),c(1:(LP-1)) ]}
          else {rbldDat<-data[c(1:(i-1), (i+1):LP),c(1:(i-1), (i+1):LP)]} ## make the appropriate merge tables with the target data point removed
          RDist<-dist(rbldDat)
          RClust<-hclust(RDist, method)

          Sdata<-.Call("snip", merge, i)
          Rdata<-.Call("rebuild", RClust$merge)

          runStab<-sqrt(sum((Sdata- Rdata)^2))
          totalStab<-totalStab+runStab
        }
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
    #' HClustStabHist(myData, "complete")
    #' HClustStabHist(myData, method= "ward.D")
    #'
    HClustStabHist<- function(data, method="complete", is_dist=F){ ## data is the data and tree is the Bubble output
      if (is_dist==F)
      { Distdata<-dist(data)}
      else{Distdata<-as.dist(data)}
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

        runStab<-sqrt(sum((Sdata- Rdata)^2))
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

      distT<-sqrt(sum((MCC1- MCC2)^2))

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

      distT<-sqrt(sum((MCC1- MCC2)^2))

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


        test<-.Call("ATree", distMat, PACKAGE = "BubbleClusteringT")

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

          runStab<-sqrt(sum((Sdata- Rdata)^2))
          totalStab<-totalStab+runStab
        }
        totalStab<-totalStab/LP
        return(totalStab)
      }

    #' Generates an associator matrix based hierarcical tree by comparing the distances between triplets of points
    #'
    #' @param data A numerical matrix
    #' @param verbose Boolean indicator of progress bar presence or absence
    #' @return hclust object of hierarcical tree of \code{data}
    #' @examples
    #'

    tripletTree<-
      function(data, verbose=FALSE){
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
        output<-.Call("tripleComp", data, samples, features, verbose, PACKAGE = "BubbleClusteringT")
        test<-.Call("ATree", output, PACKAGE = "BubbleClusteringT")

        result<-structure(list(merge = test[[1]], height = test[[2]], order = test[[3]],
                               labels = rownames(data), method = "Bubble",
                               dist.method = "euclidean"),
                          class = "hclust")
        return(result)
      }

    TripleStab<-
      function(data){
        tree<-tripletTree(data)
        merge<-tree$merge
        totalStab=0
        LD<-length(data[1,])
        LP<-length(data[,1])
        for (i in 1:LP){
          if (i==1) {rbldDat<-data[c(2:LP), ]}
          else if(i==LP) {rbldDat<-data[c(1:(LP-1)), ]}
          else {rbldDat<-data[c(1:(i-1), (i+1):LP),]} ## make the appropriate merge tables with the target data point removed
          RClust<-tripletTree(rbldDat)

          Sdata<-.Call("snip", merge, i)
          Rdata<-.Call("rebuild", RClust$merge)

          runStab<-sqrt(sum((Sdata- Rdata)^2))
          totalStab<-totalStab+runStab
        }
        totalStab<-totalStab/LP
        return(totalStab)
      }

    ################### Building Normalization function ################

    normFactor<- matrix(c(3,1.41421,4,4,5,6.63325,6,10.3923,7,14.8324,8,20.3961,9,25.8844,10,32.4962,11,39.7492,12,47.7912,13,56.2494,14,65.8635,15,75.8288,16,86.487,17,97.8673,18,110.3,19,122.882,20,136.558,21,150.725,22,165.886,23,182.242,24,198.192,25,215.606,26,232.787,27,252.103,28,270.5,29,290.83,30,311.419,31,334.275,32,354.799,33,378.74,34,401.475,35,426.818,36,449.709,37,475.615,38,501.512,39,528.927,40,558.84,41,586.58,42,614.853,43,643.654,44,674.571,45,706.504,46,738.585,47,771.793,48,804.025,49,839.382,50,874.039,51,909.521,52,947.737,53,981.967,54,1023.47,55,1061.12,56,1098.93,57,1137.33,58,1180.86,59,1222.3,60,1262.54,61,1306.08,62,1348.45,63,1394.47,64,1438.7,65,1484.13,66,1528.29,67,1573.61,68,1623.13,69,1671.08,70,1721.78,71,1772.16,72,1819.77,73,1871.61,74,1922.84,75,1978.04,76,2033.98,77,2084.63,78,2140.12,79,2195.11,80,2252.01,81,2307.26,82,2362.77,83,2430.46,84,2482.64,85,2541.92,86,2601.39,87,2665.22,88,2728.65,89,2790.52,90,2851.93,91,2916.72,92,2977.44,93,3045.41,94,3111,95,3179.22,96,3246.56,97,3319.76,98,3383.82,99,3454.77,100,3520.94,101,3594.57,102,3665.68,103,3742.26,104,3816.33,105,3886.68,106,3959.11,107,4037.33,108,4112.2,109,4190.98,110,4268.57,111,4343.45,112,4422.96,113,4508.26,114,4584.42,115,4668.13,116,4747.13,117,4827.12,118,4910.71,119,5002.29,120,5085.54,121,5167.41,122,5249.91,123,5340.49,124,5428.25,125,5511.25,126,5599.14,127,5694.72,128,5779.81,129,5874.1,130,5966.47,131,6063.89,132,6154.36,133,6259.96,134,6347.42,135,6442.91,136,6535.67,137,6632.68,138,6725.74,139,6823.96,140,6924.39,141,7025.88,142,7125.7,143,7225.17,144,7324.7,145,7429.86,146,7529.13,147,7639.71,148,7739.17,149,7849.62,150,7958.57,151,8061.03,152,8166.57,153,8275.91,154,8392.02,155,8495.73,156,8602.31,157,8717.66,158,8820.2,159,8947.32,160,9057.5,161,9166.72,162,9274.62,163,9392.66,164,9513.71,165,9629.26,166,9744.81,167,9866.3,168,9980.82,169,10091.2,170,10223.1,171,10349.2,172,10465.7,173,10591.1,174,10707.3,175,10839.1,176,10964.2,177,11088.2,178,11212.9,179,11338.5,180,11460.8,181,11586.8,182,11711.9,183,11844.9,184,11979.1,185,12112.6,186,12248.9,187,12384.6,188,12528.3,189,12642.8,190,12776.4,191,12908.4,192,13056.4,193,13175.8,194,13319.7,195,13455.6,196,13596.5,197,13742.8,198,13874.6,199,14017.8,200,14161.2),
                        byrow = T, ncol = 2)
    BubbleNorm<-function(data, n.leaves){
      if (n.leaves>200){
        message("The desired number of leaves is beyond the scope of the experimental data")
        stop()
      }
      normalized<-data/normFactor[(n.leaves-2),2]
      return(normalized)
    }
