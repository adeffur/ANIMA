################################################################
## Some utilities for supervised classification using KNN

################################################################
## Compute cross-validation statistics
cv.stats <- function(cv.table) {
  cv.stats <-  c("correct"=sum(diag(cv.table)),
                 "errors"=sum(cv.table) - sum(diag(cv.table)),
                 "hit.rate"=sum(diag(cv.table))/sum(cv.table))
  return(cv.stats)
}



################################################################
## Compute the hit rate as a function of the number of top variables
knn.hit.curve.m <- function(train, ## Training set, a matrix with N subjects and M variables
                            cl, ## Cluster definition
                            k=10, ## Number of neighbours
                            plot=TRUE, ## plot the hit rate curve
                            main='Hit rate curve: impact of the number of variables on KNN') {
  top.variables <- 1:ncol(train)

  cv.stats.table <- data.frame()

  for (v in top.variables) {
    classifier <- knn.cv(train=train[,1:v], cl=cl, k=k, l = 0, prob = TRUE, use.all = TRUE)
    pred.class <- as.vector(classifier)
    cv.table <- table("true"=cl, "predicted"=pred.class)
    cv.stats.table <- rbind(cv.stats.table,
                            c(v=v,cv.stats(cv.table)))
  }
  names(cv.stats.table) <- c("v", "correct", "errors", "hit.rate")


  
  ## Plot the effect of variable numbers on the hit rate
  if (plot) {
    plot(cv.stats.table$v,
         cv.stats.table$hit.rate,
         main=main,
         xlab='Number of variables',
         ylab='Hit rate',
         ylim=c(0,1),
         type='l',
         panel.first=c(
           abline(v=seq(from=0,to=400,by=50), col='#888888'),
           abline(h=(0:10)/10, col='#888888'))
         )

  }
  return(cv.stats.table)
}
