plotSampleRelationsAD<-function (x, subset = NULL, cv.Th = 0.1, standardize = TRUE, 
          method = c("cluster", "mds"), dimension = c(1, 2), color = NULL, backgr = NULL,plotchar=24,
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
    plot(ppoints[, dimension[1]], ppoints[, dimension[2]], 
         type = "n", xlab = paste("Principal Component ", 
                                  dimension[1], " (", percent[dimension[1]], "%)", 
                                  sep = ""), ylab = paste("Principal Component ", 
                                                          dimension[2], " (", percent[dimension[2]], "%)", 
                                                          sep = ""), main = main, ...)
    points(ppoints[, dimension[1]], ppoints[, dimension[2]], 
         col = color,bg=backgr, cex = 1,pch=plotchar)
    text(ppoints[, dimension[1]], ppoints[, dimension[2]],x$sample_name, 
           col = "lightgrey",bg=backgr, cex = 0.8,pch=plotchar)
    attr(ppoints, "geneNum") <- length(subset)
    attr(ppoints, "threshold") <- cv.Th
    return(invisible(ppoints))
  }
}