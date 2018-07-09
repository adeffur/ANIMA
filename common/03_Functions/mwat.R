###MWAT FUNCTION####
#' Title
#'
#' @param study 
#' @param squareM 
#' @param moduleM 
#' @param edgeM 
#'
#' @return
#' @export
#'
#' @examples
mwat<-function(study="Berry",squareM,moduleM,edgeM){
  library(fastcluster)
  dataM=eval(parse(text=squareM))
  if(class(dataM)=="LumiBatch"){data.v<-lumiT(dataM,simpleOutput=FALSE)}else{data.v<-lumiT(dataM,simpleOutput=FALSE,method="log2")}
  data.q<-lumiN(data.v,method="quantile")
  pd<-pData(data.q)
  
  #pseudo classvar2:control, classvar2:cases, classvar1:control, classvar1:case, classvar1:control|classvar1:case
  query<-paste("MATCH (n:wgcna {square:'",squareM,"', name:'blue'}) RETURN n.edge AS edge, n.contrastvar AS contrastvar, n.contrast AS contrast",sep="")
  res<-cypher(graph,query)
  classvar1<-res$contrastvar[1]
  classvar2<-res$contrastvar[3]
  c1con<-unlist(strsplit(res$contrast[1]," - "))[2]
  c1case<-unlist(strsplit(res$contrast[1]," - "))[1]
  c2con<-unlist(strsplit(res$contrast[3]," - "))[2]
  c2case<-unlist(strsplit(res$contrast[3]," - "))[1]
  subsets<-eval(parse(text=paste("list(pd$",classvar2,"=='",c2con,"',pd$",classvar2,"=='",c2case,"',pd$",classvar1,"=='",c1con,"',pd$",classvar1,"=='",c1case,"',pd$",classvar1,"=='",c1con,"'|","pd$",classvar1,"=='",c1case,"')",sep="")))
  subset=subsets[[as.numeric(edgeM)]]
 
  data.qs<-data.q[,subset]
  
  #get nuIDs for module probes
  query<-paste("MATCH (n:wgcna {name:'",moduleM,"', square:'",squareM,"', edge:",edgeM,"})-[r0]-(p:PROBE) RETURN p.name AS nuID",sep="")
  nuIDs<-cypher(graph,query)
  #subset the data and make correlation matrix
  turq.data<-exprs(data.qs)[nuIDs$nuID,]
  if (prod(dim(turq.data)) != 0){
    print("cormat")
    print(dim(turq.data))
    cormat<-cor(t(turq.data),method="pearson")
    print("begin clustering")
    clustres<-fastcluster::hclust(dist(cormat),method='ward.D2')
    print("begin dendro")
    dendro.manual<-as.dendrogram(clustres)
    print("begin ordering")
    idvec<-order.dendrogram(dendro.manual)
    dimnames(cormat)<-list(NULL,NULL)
    print("end clustering")
    print(dim(cormat))
    
    #generate annotation matrices
    #CELLS
    print("begin cells")
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
    print("begin pw and logfc")
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
    print("numpw")
    print(numpw)
    
    if(numpw>3){
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
    reds<-fc2>=0
    mcolvec<-character(length=length(fc2))
    mcolvec[blues]<-"blue"
    mcolvec[reds]<-"red"
    
    mcolvec2<-cbind(row.names(fc),mcolvec)
    row.names(mcolvec2)<-row.names(fc)
    mcolvec2[idvec,2]
    
    print("begin wgcna")
    #WGCNA module stats for annotation
    query<-paste("MATCH (n:wgcna {name:'",moduleM,"', square:'",squareM,"', edge:",edgeM,"}) RETURN n.name AS name, n.square AS square, n.edge AS edge, n.modAUC1 as modAUC1, n.modAUC2 AS modAUC2, n.diffME AS diffME, n.sigenrich AS sigenrich",sep="")
    wa<-cypher(graph,query)
    
    print("begin plot")
    
    res<-annHeatmap2(cormat,
                     dendrogram=list("status"="hidden","dendro"=dendro.manual),
                     annotation=list(
                       Row=list("data"=cellmatrix,"asIs"=TRUE,control=list(boxh=2)),
                       Col=list("data"=total,"asIs"=TRUE,control=list("pch"=16,"col.pch"=mcolvec2[idvec,2],numfac=4))
                     ),
                     legend=TRUE,
                     col=blueWhiteRed,
                     scale="none"
                     #labels=list(Row=list(cex=cexval),Col=list(cex=cexval)),
                     
    ) 
    #res<-annHeatmap2(cormat,dendrogram=list("status"="hidden","dendro"=dendro.manual),legend=TRUE,col=blueWhiteRed,scale="none",annotation=list(Row=list("data"=cellmatrix,"asIs"=TRUE,Col=list("data"=total,"asIs"=TRUE,control=list("pch"=16,"col.pch"=mcolvec2[idvec,2],numfac=4)))))
    
    #plot(res,widths=c(6,1),heights=c(6,1))
    #labels=list(Row=list(cex=input$labelFontSize, side=4),Col=list(cex=1)),

    res2<-res
    print("layout")
    res2$layout$plot<-matrix(c(0,5,0,0,4,0,0,1,2,0,3,0),ncol=4,dimnames=list(c("","image","rowAnn"),c("","leg2","image","colAnn")))
    #res2$layout$plot<-matrix(c(0,5,0,0,4,0,0,1,2,0,3,0),ncol=4,dimnames=list(c("","image","rowAnn"),c("","leg2","image","colAnn")))
    #res2$layout$plot<-matrix(c(0,4,0,0,1,2,0,3,0),ncol=3,dimnames=list(c("","image","rowAnn"),c("leg2","image","colAnn")))
    print("actual plot")
    print(res2)
    #plot(res2,widths=c(1.8,1,4,1)/2,heights=c(1.5,4,1.5*(numpw/sqrt(numpw)))/2)
    #plot(res2,widths=c(1.5,.5,3,1)/2,heights=c(1,3,1.0*(numpw/sqrt(numpw)))/2)
    plot(res2,widths=c(1,0.5,3,.8)/2,heights=c(1,4.5,1.0*(numpw/sqrt(numpw)))/3)
    #text(0,0,"TEST",add=TRUE)
    #plot(10,5,type="n",axes=FALSE,ann=FALSE,xlim=c(0, 10),ylim = c(0,10))
    mtext(paste(squareM,"edge:",edgeM,"\n","module:",moduleM,"\n","modAUC1:",sprintf("%.3f",wa$modAUC1),"\n","modAUC2:",sprintf("%.3f",wa$modAUC2),"\n","diffME:",sprintf("%.3f",wa$diffME),"\n","sigenrich:",sprintf("%.3f",wa$sigenrich)),cex=1.0,line=0)
    #text(-5,10,paste(squareM,"edge:",edgeM),cex=2.0,pos=4,xpd=NA)
    #text(-5,9,paste("module:",moduleM),cex=2.0,pos=4,xpd=NA)
    #text(-5,8,paste("modAUC1:",sprintf("%.3f",wa$modAUC1)),cex=2,pos=4,xpd=NA)
    #text(-5,7,paste("modAUC2:",sprintf("%.3f",wa$modAUC2)),cex=2,pos=4,xpd=NA)
    #text(-5,6,paste("diffME:",sprintf("%.3f",wa$diffME)),cex=2,pos=4,xpd=NA)
    #text(-5,5,paste("sigenrich:",sprintf("%.3f",wa$sigenrich)),cex=2,pos=4,xpd=NA)
    }#end final numpw check
    }#end turq.data dimension check
}#end function
