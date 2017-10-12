superHeatmapmatrix<-function(x,y=selected.probes,phenomatrix=phenomatrix,scale="none",addTit=NULL){

#cluster
genes.cor <- cor(t(x[y,]), use="pairwise.complete.obs",method="pearson")
genes.cor.dist <- as.dist(1-genes.cor)
genes.tree <- hclust(genes.cor.dist,method='average')
samples.cor.spearman <- cor(x[y,],use="pairwise.complete.obs",method="spearman")
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
  mycols <- as.vector(phenomatrix[,2])[order.dendrogram(smpCol)]
  i <- 0
})
dL <- dendrapply(smpCol, colLab)

max.level <- max(abs(range(x[y,])))

br <- c(seq(from=-13,to=-6,by=4/72),
            seq(from=-6, to=1.1,by=4/72),
            max.level
)
br

#heatmap
heatmap_ad(x[y,],
           scale=scale, ## AVOID RESCALING THE VALUES OF ALL COLUMNS
           col=palette.BYR(),
           breaks=br,
           Rowv=as.dendrogram(genes.tree),
           Colv=dL,
           main=paste(nrow(x[y,]),"probes",addTit),
           ColSideColors=phenomatrix,
           margins=c(9,5),
           #RowSideColors=rsc,
           noan=2,labRow="",verbose=TRUE   
)
}