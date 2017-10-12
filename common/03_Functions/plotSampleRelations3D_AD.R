plotSampleRelations3D_AD<-function (x, subset = NULL, cv.Th = 0.1, standardize = TRUE, 
                                 method = c("cluster", "mds"), dimension = c(1, 2,3), color = NULL, backgr = NULL,plotchar=24,spheresize=2,
                                 main = NULL, ...) 
{
  if (is(x, "ExpressionSet")) {
    dataMatrix <- exprs(x)
  }
  else if (is.matrix(x)) {
    dataMatrix <- x
  }
  else {
    stop("The class of \"x\" should be matrix or LumiBatch!")
  }
  if (standardize) 
    dataMatrix <- scale(dataMatrix)
  if (is.null(subset)) {
    probeList <- rownames(dataMatrix)
    if (is.null(probeList)) 
      probeList <- 1:nrow(dataMatrix)
    if (cv.Th > 0) {
      cv.gene <- apply(dataMatrix, 1, function(x) sd(x)/mean(x))
      subset <- probeList[abs(cv.gene) > cv.Th]
      if (is.null(main)) 
        main <- paste("Sample relations based on", length(subset), 
                      "genes with sd/mean >", cv.Th)
    }
    else {
      subset <- probeList
      if (is.null(main)) 
        main <- paste("Sample relations based on", length(subset), 
                      "genes")
    }
  }
  else {
    if (length(subset) == 1 && is.numeric(subset)) {
      subset <- sample(1:nrow(dataMatrix), min(subset, 
                                               nrow(dataMatrix)))
    }
    if (is.null(main)) 
      main <- paste("Sample relations based on", length(subset), 
                    "selected genes")
  }
  dd <- dist(t(dataMatrix[subset, ]))
  method <- match.arg(method)
  if (method == "cluster") {
    hc = hclust(dd, "ave")
    plot(hc, xlab = "Sample", main = main, ...)
    attr(hc, "geneNum") <- length(subset)
    attr(hc, "threshold") <- cv.Th
    return(invisible(hc))
  }
  else {
    mds.result <- cmdscale(dd, k = max(dimension), eig = TRUE)
    ppoints <- mds.result$points
    eig <- mds.result$eig
    percent <- round(eig/sum(eig) * 100, 1)
    if (is.null(color)) {
      color <- 1
    }
    else {
      if (!is.numeric(color)) {
        allColor <- colors()
        if (!all(is.element(color, allColor))) {
          color <- as.numeric(factor(color, levels = unique(color)))
        }
      }
    }
    plot3d(ppoints[, dimension[1]], ppoints[, dimension[2]], ppoints[,dimension[3]],
         type = "s", size=spheresize,
           xlab = paste("Principal Component ",dimension[1], " (", percent[dimension[1]], "%)",sep = ""), 
           ylab = paste("Principal Component ",dimension[2], " (", percent[dimension[2]], "%)",sep = ""), 
           zlab = paste("Principal Component ",dimension[3], " (", percent[dimension[3]], "%)",sep = ""),
           main = main,col=color)
    #points(ppoints[, dimension[1]], ppoints[, dimension[2]], ppoints[,dimension[3]],col = color,bg=backgr, cex = 1,pch=plotchar)
  
  }
}