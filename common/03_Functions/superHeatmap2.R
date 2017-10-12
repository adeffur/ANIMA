superHeatmap2<-function(x,y=selected.probes,phenomatrix=phenomatrix,scale=scale,addTit=NULL,low=low,high=high){

#cluster
genes.cor <- cor(t(exprs(x)[y,]), use="pairwise.complete.obs",method="pearson")
genes.cor.dist <- as.dist(1-genes.cor)
genes.tree <- hclust(genes.cor.dist,method='average')
samples.cor.spearman <- cor(exprs(x)[y,],use="pairwise.complete.obs",method="spearman")
samples.cor.spearman.dist <- as.dist(1-samples.cor.spearman)
samples.tree <- hclust(samples.cor.spearman.dist,method='average')
#dendrogram colors

#coloured dendrogram
smpCol <- as.dendrogram(samples.tree)
local({
  colLab <<- function(n) {
    if(is.leaf(n)) {
      a <- attributes(n)
      i <<- i+1
      #         attr(n, "nodePar") <-
      #             c(a$nodePar, list(lab.col = mycols[i], lab.font= i%%3),pch=NULL)
      attr(n, "edgePar") <-
        c(a$edgePar, list(col = mycols[i]),lwd=3)
    }
    n
  }
  mycols <- as.vector(phenomatrix[,1])[order.dendrogram(smpCol)]
  i <- 0
})
dL <- dendrapply(smpCol, colLab)

max.level <- max(abs(range(exprs(x)[y,])))
min.level <- min(abs(range(exprs(x)[y,])))
low<-1.05*min.level
high<-0.6*max.level

br <- c(seq(from=min.level,to=low,length.out=30),
        seq(from=low, to=high,length.out=105),
        seq(from=high,to=max.level,length.out=120),        
        max.level
)

#heatmap
heatmap_ad(as.matrix(exprs(x)[y,]),
           scale=scale, ## AVOID RESCALING THE VALUES OF ALL COLUMNS
           col=palette.BYR(),
           #breaks=br,
           Rowv=as.dendrogram(genes.tree),
           Colv=dL,
           main=paste(nrow(exprs(x)[y,]),"probes",addTit),
           ColSideColors=phenomatrix,
           margins=c(9,5),
           #RowSideColors=rsc,
           noan=2,labRow="",verbose=TRUE   
)
}