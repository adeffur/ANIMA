#module correlation plot: correlation matrix of expression data, with cell type and pathway annotations

modCorPlot<-function(study="Berry",squareM,moduleM,edgeM){
  dataM=eval(parse(text=squareM))
  if(class(datalist[[i]])=="LumiBatch"){data.v<-lumiT(dataM,simpleOutput=FALSE)}else{data.v<-lumiT(datalist[[i]],simpleOutput=FALSE,method="log2")}
  data.q<-lumiN(data.v,method="quantile")
  pd<-pData(data.q)
  if(study=="Berry"){
    subsets<-list(pd$sex=="Female",pd$sex=="Male",pd$class=="notActiveTB",pd$class=="activeTB",pd$class=="activeTB"|pd$class=="notActiveTB")
    subset=subsets[[as.numeric(edgeM)]]
  }else if (study=="IMPI" & squareM == "square1"){
    subsets<-list(pd$HIV.Status=="negative",pd$HIV.Status=="positive",pd$class3=="notActiveTB",pd$class3=="activeTB",pd$class3=="activeTB"|pd$class3=="notActiveTB")
    subset=subsets[[as.numeric(edgeM)]]
  }else if (study=="IMPI" & squareM == "square4"){
    subsets<-list(pd$HIV.Status=="negative",pd$HIV.Status=="positive",pd$Compartment=="blood",pd$Compartment=="fluid",pd$Compartment=="blood"|pd$Compartment=="fluid")
    subset=subsets[[as.numeric(edgeM)]]
  }
  
  data.qs<-data.q[,subset]
  
  #get nuIDs for module probes
  query<-paste("MATCH (n:wgcna {name:'",moduleM,"', square:'",squareM,"', edge:",edgeM,"})-[r0]-(p:PROBE) RETURN p.name AS nuID",sep="")
  nuIDs<-cypher(graph,query)
  #subset the data and make correlation matrix
  turq.data<-exprs(data.qs)[nuIDs$nuID,]
  cormat<-cor(t(turq.data),method="pearson")
  dendro.manual<-as.dendrogram(hclust(dist(cormat)))
  idvec<-order.dendrogram(dendro.manual)
  dimnames(cormat)<-list(NULL,NULL)
  
  
  #generate annotation matrices
  #CELLS
  query<-paste("MATCH (c:cellEx)-[r3]-(n:wgcna {name:'",moduleM,"', square:'",squareM,"', edge:",edgeM,"})-[r0]-(p:PROBE)-[r1]-(s:SYMBOL)-[r2]-(c:cellEx) RETURN p.name AS probe, s.name AS gene, c.name AS cell",sep="")
  cellexlist<-cypher(graph,query)
  cells<-unique(cellexlist$cell)
  numcells<-length(unique(cellexlist$cell))
  cellmatrix<-matrix(nrow=nrow(cormat),ncol=numcells,data=0)
  if (numcells>0){
    dimnames(cellmatrix)<-list(nuIDs$nuID,cells)
    for(cell in 1:numcells){
      cellmatrix[,cell]<-nuIDs$nuID%in%cellexlist$probe[which(cellexlist$cell==cells[cell])]
    }
  }else{cellmatrix<-matrix(nrow=nrow(cormat),ncol=1,data=0)
  dimnames(cellmatrix)<-list(nuIDs$nuID,"None enriched")}
  #PATHWAYS and logfc
  query<-paste("MATCH (pw)-[r3]-(n:wgcna {name:'",moduleM,"', square:'",squareM,"', edge:",edgeM,"})-[r0]-(p:PROBE)-[r1]-(s:SYMBOL)-[r2]-(pw) WHERE (pw:ImmunePW OR pw:reactomePW) RETURN p.name AS probe, s.name AS gene, pw.name AS pathway, r3.qvalue AS qvalue",sep="")
  pwlist<-cypher(graph,query)
  try(pwlist<-pwlist[order(pwlist$qvalue,decreasing=FALSE),])
  pw<-unique(pwlist$pathway)
  numpw<-length(unique(pwlist$pathway))
  if (numpw < 3){
    query<-paste("MATCH (pw)-[r3]-(n:wgcna {name:'",moduleM,"', square:'",squareM,"', edge:",edgeM,"})-[r0]-(p:PROBE)-[r1]-(s:SYMBOL)-[r2]-(pw) WHERE (pw:ImmunePW OR pw:reactomePW OR pw:PalWangPW) RETURN p.name AS probe, s.name AS gene, pw.name AS pathway, r3.qvalue AS qvalue",sep="")
    pwlist<-cypher(graph,query)
    try(pwlist<-pwlist[order(pwlist$qvalue,decreasing=FALSE),])
    pw<-unique(pwlist$pathway)
    numpw<-length(unique(pwlist$pathway))
  }
  if(numpw>20){
    numpw<-20
    pw<-pw[1:numpw]
    pwlist<-pwlist[which(pwlist$pathway%in%pw),]
  }
  pwmatrix<-matrix(nrow=nrow(cormat),ncol=numpw,data=0)
  dimnames(pwmatrix)<-list(nuIDs$nuID,pw)
  
  for(path in 1:numpw){
    pwmatrix[,path]<-nuIDs$nuID%in%pwlist$probe[which(pwlist$pathway==pw[path])]
  }
  pwmatrix<-pwmatrix[,ncol(pwmatrix):1]
  query<-paste("MATCH (n:wgcna {name:'",moduleM,"', square:'",squareM,"', edge:",edgeM,"})-[r0]-(p:PROBE) RETURN p.name AS probe, p.logfc as logfc",sep="")
  logfclist<-cypher(graph,query)
  fc<-matrix(nrow=nrow(cormat),ncol=1,data=0)
  dimnames(fc)<-list(nuIDs$nuID,"logfc")
  fc[,1]<-logfclist$logfc
  fc2<-as.numeric(fc)
  total<-cbind(pwmatrix,fc)
  
  blues<-fc2<0
  reds<-fc2>0
  mcolvec<-character(length=length(fc2))
  mcolvec[blues]<-"blue"
  mcolvec[reds]<-"red"
  
  mcolvec2<-cbind(row.names(fc),mcolvec)
  row.names(mcolvec2)<-row.names(fc)
  mcolvec2[idvec,2]
  
  res<-annHeatmap2(cormat,dendrogram=list("status"="hidden","dendro"=dendro.manual),legend=TRUE,col=blueWhiteRed,scale="none",annotation=list(Row=list("data"=cellmatrix,"asIs"=TRUE),Col=list("data"=total,"asIs"=TRUE,control=list("pch"=16,"col.pch"=mcolvec2[idvec,2],numfac=4))))
  #plot(res,widths=c(6,1),heights=c(6,1))
  res2<-res
  res2$layout$plot<-matrix(c(0,0,0,0,4,0,0,1,2,0,3,0),ncol=4,dimnames=list(c("","image","rowAnn"),c("","leg2","image","colAnn")))
  plot(res2,widths=c(1.3,1,4,1),heights=c(1.5,4,1.5*(numpw/sqrt(numpw))))
  
}#end function