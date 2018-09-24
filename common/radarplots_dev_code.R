
library(ggradar2)
library(gridExtra)
#query<-paste("MATCH (n:wgcna {square:'nonLTBI.LTBI',name:'red'})-[r]-(x) WHERE (x:cellEx) RETURN DISTINCT n.name AS module,r.qvalue as qvalue,x.name as name, labels(x) as kind",sep="")
query<-paste("MATCH (n:wgcna {square:'blood.PCF.defPC',name:'brown'})-[r]-(x) WHERE (x:cellEx) RETURN DISTINCT n.name AS module,r.qvalue as qvalue,x.name as name, labels(x) as kind",sep="")
anq<-"MATCH (c:cellEx)-[r]-(c2:CELL) RETURN DISTINCT c.name as cells,c2.name AS group"
cells<-cypher(graph,anq)
cells$SIG<-rep(0,nrow(cells))

res<-cypher(graph,query)
res$SIG<--log10(res$qvalue)
res$SIG[which(res$SIG>5)]<-5
res2<-data.frame(as.numeric(res$SIG))
rownames(res2)<-res$name
colnames(res2)<-"sig"
for (row in 1:nrow(res2)){
  cells$SIG[which(cells$cells==res$name[row])]<-res2$sig[row]
 
}
res3<-data.frame(as.numeric(cells$SIG))
rownames(res3)<-cells$cells
res4<-t(res3)
plot1<-ggradar2(res4,grid.max=5,values.radar = c("0","1","5"),plot.title = "cellEx sig",axis.label.size = 4)

###################

query<-paste("MATCH (n:wgcna {square:'blood.PCF.defPC',name:'brown',edge:5})-[r]-(x) WHERE (x:cellprop) RETURN DISTINCT n.name AS module,r.qvalue as qvalue,r.weight AS pr,r.Rsq AS Rsq,x.name as name, labels(x) as kind",sep="")
anq<-"MATCH (c:cellprop)-[r]-(c2:CELL) RETURN DISTINCT c.name as cells,c2.name AS group"
cells<-cypher(graph,anq)
cells$SIG<-rep(0,nrow(cells))

res<-cypher(graph,query)
res$SIG<--log10(res$qvalue)
res$SIG[which(res$SIG>10)]<-10
res2<-data.frame(as.numeric(res$SIG))
rownames(res2)<-res$name
colnames(res2)<-"sig"
for (row in 1:nrow(res2)){
  cells$SIG[which(cells$cells==res$name[row])]<-res2$sig[row]
  
}
res3<-data.frame(as.numeric(cells$SIG))
rownames(res3)<-cells$cells
res4<-t(res3)
plot2a<-ggradar2(res4,grid.max=10,values.radar = c("0","1","10"),plot.title = "cellprop sig",axis.label.size = 4)

res3b<-res[which(res$pr!=0),]
res3c<-res3b[order(res3b$pr),]
res4<-t(res3c$pr)
data.frame(res4)
colnames(res4)<-res3c$name
plot2<-ggradar2(res4,grid.max=1,grid.min=-1, grid.mid=0,values.radar = c("-1","0","1"),plot.title = "cellprop cor",axis.label.size = 4)


###################


query<-paste("MATCH (n:wgcna {square:'blood.PCF.probPC',name:'brown',edge:5})-[r]-(x) WHERE (x:pheno) RETURN DISTINCT n.name AS module,r.qvalue as qvalue,r.weight as pr,r.Rsq AS Rsq,x.name as name, labels(x) as kind",sep="")
anq<-"MATCH (c:pheno) RETURN DISTINCT c.name as phenos"
phenos<-cypher(graph,anq)

phenos$qvalue<-rep(0,nrow(phenos))

res<-cypher(graph,query)

res$SIG<--log10(res$qvalue)
res$SIG[which(res$SIG>10)]<-10
res2<-data.frame(as.numeric(res$SIG))#
rownames(res2)<-res$name
colnames(res2)<-"qvalue"
for (row in 1:nrow(res2)){
  phenos$qvalue[which(phenos$phenos==res$name[row])]<-res2$qvalue[row]
  
}
res3<-data.frame(as.numeric(phenos$qvalue))
res3$phenos<-phenos$phenos
rownames(res3)<-phenos$phenos
colnames(res3)<-c("qvalue","phenos")

res3b<-res3[which(res3$qvalue!=0),]
res3c<-res3b[order(res3b$qvalue),]
res4<-t(res3c$qvalue)

colnames(res4)<-rownames(res3c)
plot3a<-ggradar2(res4,grid.max=10,grid.min=0, grid.mid=5,values.radar = c("0","5","10"),plot.title = "pheno sig",axis.label.size = 4)


###


phenos$pr<-rep(0,nrow(phenos))

res<-cypher(graph,query)
res2<-data.frame(as.numeric(res$pr))#
rownames(res2)<-res$name
colnames(res2)<-"pr"
for (row in 1:nrow(res2)){
  phenos$pr[which(phenos$phenos==res$name[row])]<-res2$pr[row]
  
}
res3<-data.frame(as.numeric(phenos$pr))
res3$phenos<-phenos$phenos
rownames(res3)<-phenos$phenos
colnames(res3)<-c("pr","phenos")

res3b<-res3[which(res3$pr!=0),]
res3c<-res3b[order(res3b$pr),]
res4<-t(res3c$pr)

colnames(res4)<-rownames(res3c)
plot3<-ggradar2(res4,grid.max=1,grid.min=-1, grid.mid=0,values.radar = c("-1","0","1"),plot.title = "pheno cor",axis.label.size = 4)


###################

query<-paste("MATCH (n:wgcna {square:'blood.PCF.defPC',name:'brown',edge:5})-[r]-(x) WHERE (x:ImmunePW) OR (x:PalWangPW) OR (x:reactomePW) RETURN DISTINCT n.name AS module,r.qvalue as qvalue ,x.name as name, labels(x) as kind",sep="")
anq<-"MATCH (x) WHERE (x:ImmunePW) OR (x:PalWangPW) OR (x:reactomePW) RETURN DISTINCT x.name as pathways"
pathways<-cypher(graph,anq)
pathways$SIG<-rep(0,nrow(pathways))

res<-cypher(graph,query)
res$SIG<--log10(res$qvalue)
res$SIG[which(res$SIG>10)]<-10
res2<-data.frame(as.numeric(res$SIG))
rownames(res2)<-res$name
colnames(res2)<-"SIG"
for (row in 1:nrow(res2)){
  pathways$SIG[which(pathways$pathways==res$name[row])]<-res2$SIG[row]
  
}
res3<-data.frame(as.numeric(pathways$SIG))
res3$pathways<-pathways$pathways
rownames(res3)<-pathways$pathways
colnames(res3)<-c("SIG","pathways")

res3b<-res3[which(res3$SIG!=0),]
res3c<-res3b[order(res3b$SIG),]
res4<-t(res3c$SIG)

colnames(res4)<-rownames(res3c)
plot4<-ggradar2(res4,grid.max=10,grid.min=0, grid.mid=5,values.radar = c("0","5","10"),plot.title = "pathway sig",axis.label.size = 4)

#####
grid.arrange(plot1,plot2a,plot2,plot3a,plot3,plot4,nrow=2)
