###MWAT FUNCTION####
mwat<-function(study="Berry",squareM,moduleM,edgeM){
  library(fastcluster)
  dataM=eval(parse(text=squareM))
  if(class(dataM)=="LumiBatch"){data.v<-lumiT(dataM,simpleOutput=FALSE)}else{data.v<-lumiT(dataM,simpleOutput=FALSE,method="log2")}
  data.q<-lumiN(data.v,method="quantile")
  pd<-pData(data.q)
  # if(study=="Berry"){
  #   subsets<-list(pd$sex=="Female",pd$sex=="Male",pd$class=="notActiveTB",pd$class=="activeTB",pd$class=="activeTB"|pd$class=="notActiveTB")
  #   subset=subsets[[as.numeric(edgeM)]]
  # }else if (study=="IMPI" & squareM == "square1"){
  #   subsets<-list(pd$HIV.Status=="negative",pd$HIV.Status=="positive",pd$class3=="notActiveTB",pd$class3=="activeTB",pd$class3=="activeTB"|pd$class3=="notActiveTB")
  #   subset=subsets[[as.numeric(edgeM)]]
  # }else if (study=="IMPI" & squareM == "square2"){
  #   subsets<-list(pd$HIV.Status=="negative",pd$HIV.Status=="positive",pd$class=="CON_healthy",pd$class=="CON_LTBI",pd$class=="CON_LTBI"|pd$class=="CON_healthy")
  #   subset=subsets[[as.numeric(edgeM)]]
  # }else if (study=="IMPI" & squareM == "square3"){
  #   subsets<-list(pd$HIV.Status=="negative",pd$HIV.Status=="positive",pd$class=="CON_PTB",pd$class=="TBPC",pd$class=="TBPC"|pd$class=="CON_PTB")
  #   subset=subsets[[as.numeric(edgeM)]]
  # }else if (study=="IMPI" & squareM == "square4"){
  #   subsets<-list(pd$HIV.Status=="negative",pd$HIV.Status=="positive",pd$Compartment=="blood",pd$Compartment=="fluid",pd$Compartment=="blood"|pd$Compartment=="fluid")
  #   subset=subsets[[as.numeric(edgeM)]]
  # }else if (study=="IMPI" & squareM == "square5"){
  #   subsets<-list(pd$ECP=="Eff",pd$ECP=="EC",pd$Compartment=="blood",pd$Compartment=="fluid",pd$Compartment=="blood"|pd$Compartment=="fluid")
  #   subset=subsets[[as.numeric(edgeM)]]
  # }else if (study=="IMPI" & squareM == "square6"){
  #   subsets<-list(pd$HIV.Status=="negative",pd$HIV.Status=="positive",pd$class=="CON_LTBI",pd$class=="CON_PTB",pd$class=="CON_PTB"|pd$class=="CON_LTBI")
  #   subset=subsets[[as.numeric(edgeM)]]
  # }else if (study=="IMPI" & squareM == "rec1"){
  #   subsets<-list(pd$HIV.Status=="negative",pd$HIV.Status=="positive",pd$class=="CON_healthy",pd$class=="CON_PTB",pd$class=="CON_PTB"|pd$class=="CON_healthy")
  #   subset=subsets[[as.numeric(edgeM)]]
  # }else if (study=="IMPI" & squareM == "rec2"){
  #   subsets<-list(pd$HIV.Status=="negative",pd$HIV.Status=="positive",pd$class=="CON_healthy",pd$class=="TBPC",pd$class=="TBPC"|pd$class=="CON_healthy")
  #   subset=subsets[[as.numeric(edgeM)]]
  # }else if (study=="Chong" & squareM == "SLE5"){
  #   subsets<-list(pd$set=="train",pd$set=="test",pd$class=="control",pd$class=="DLE",pd$class=="control"|pd$class=="DLE")
  #   subset=subsets[[as.numeric(edgeM)]]
  # }else if (study=="Chiche" & squareM == "SLE1"){
  #   subsets<-list(pd$sex=="female",pd$sex=="male",pd$stageName=="control",pd$stageName=="inactiveSLE",pd$stageName=="control"|pd$stageName=="inactiveSLE")
  #   subset=subsets[[as.numeric(edgeM)]]
  # }else if (study=="Chiche" & squareM == "SLE6"){
  #   subsets<-list(pd$sex=="female",pd$sex=="male",pd$stageName=="control",pd$stageName=="activeSLE",pd$stageName=="control"|pd$stageName=="activeSLE")
  #   subset=subsets[[as.numeric(edgeM)]]
  # }else if (study=="Chiche" & squareM == "SLE2"){
  #   subsets<-list(pd$sex=="female",pd$sex=="male",pd$stageName=="inactiveSLE",pd$stageName=="activeSLE",pd$stageName=="control"|pd$stageName=="activeSLE")
  #   subset=subsets[[as.numeric(edgeM)]]
  # }else if (study=="Chiche" & squareM == "SLE7"){
  #   subsets<-list(pd$sex=="female",pd$sex=="male",pd$group=="control",pd$group=="SLE",pd$group=="control"|pd$group=="SLE")
  #   subset=subsets[[as.numeric(edgeM)]]
  #   
  # }else if (study=="HE" & squareM == "HErec1"){
  #   subsets<-list(pd$Sex=="F",pd$Sex=="M",pd$CLASS2=="NIL",pd$CLASS2=="ACT",pd$CLASS2=="NIL"|pd$CLASS2=="ACT")
  #   subset=subsets[[as.numeric(edgeM)]]
  #   
  # }else if (study=="EU" & squareM == "EUrec3"){
  #   subsets<-list(pd$HIV.Status=="negative",pd$HIV.Status=="positive",pd$class=="LTBI",pd$class=="OD",pd$class=="LTBI"|pd$class=="OD")
  #   subset=subsets[[as.numeric(edgeM)]]
  # }else if (study=="EU" & squareM == "EUrec4"){
  #   subsets<-list(pd$HIV.Status=="negative",pd$HIV.Status=="positive",pd$class=="OD",pd$class=="PTB",pd$class=="OD"|pd$class=="PTB")
  #   subset=subsets[[as.numeric(edgeM)]]
  # }else if (study=="EU" & squareM == "EUsquare6"){
  #   subsets<-list(pd$HIV.Status=="negative",pd$HIV.Status=="positive",pd$class=="LTBI",pd$class=="PTB",pd$class=="LTBI"|pd$class=="PTB")
  #   subset=subsets[[as.numeric(edgeM)]]
  # }else if (study=="EU" & squareM == "EUsquare6MF_neg"){
  #   subsets<-list(pd$gend=="F",pd$gend=="M",pd$class=="LTBI",pd$class=="PTB",pd$class=="LTBI"|pd$class=="PTB")
  #   subset=subsets[[as.numeric(edgeM)]]
  # }else if (study=="EU" & squareM == "EUsquare6MF_pos"){
  #   subsets<-list(pd$gend=="F",pd$gend=="M",pd$class=="LTBI",pd$class=="PTB",pd$class=="LTBI"|pd$class=="PTB")
  #   subset=subsets[[as.numeric(edgeM)]]
  # 
  # }else if (study=="HIV"){
  #   subsets<-list(pd$sex=="female",pd$sex=="male",pd$class=="Healthy",pd$class=="HIV",pd$class=="Healthy"|pd$class=="HIV")
  #   subset=subsets[[as.numeric(edgeM)]]
  # }else if (study=="malaria"){
  #   subsets<-list(pd$sex=="F",pd$sex=="M",pd$malaria.class=="asymptomatic",pd$malaria.class=="symptomatic",pd$malaria.class=="asymptomatic"|pd$malaria.class=="symptomatic")
  #   subset=subsets[[as.numeric(edgeM)]]
  # }else if (study=="respVir"){
  #   subsets<-list(pd$sex=="female",pd$sex=="male",pd$time=="Baseline",pd$time=="Day0",pd$time=="Baseline"|pd$time=="Day0")
  #   subset=subsets[[as.numeric(edgeM)]]
  # }
  
  ##New subsets
  
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
  ##end new subsets
  
  data.qs<-data.q[,subset]
  
  #get nuIDs for module probes
  query<-paste("MATCH (n:wgcna {name:'",moduleM,"', square:'",squareM,"', edge:",edgeM,"})-[r0]-(p:PROBE) RETURN p.name AS nuID",sep="")
  nuIDs<-cypher(graph,query)
  #subset the data and make correlation matrix
  turq.data<-exprs(data.qs)[nuIDs$nuID,]
  if (prod(dim(turq.data)) != 0){
    print("cormat")
    cormat<-cor(t(turq.data),method="pearson")
    print("begin clustering")
    clustres<-fastcluster::hclust(dist(cormat))
    print("begin dendro")
    dendro.manual<-as.dendrogram(clustres)
    print("begin ordering")
    idvec<-order.dendrogram(dendro.manual)
    dimnames(cormat)<-list(NULL,NULL)
    print("end clustering")
    
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
                       Row=list("data"=cellmatrix,"asIs"=TRUE,control=list(cex=2,boxh=2)),
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
    print("actual plot")
    #plot(res2,widths=c(1.8,1,4,1)/2,heights=c(1.5,4,1.5*(numpw/sqrt(numpw)))/2)
    plot(res2,widths=c(1.5,.5,4,1)/2,heights=c(2,4,1.5*(numpw/sqrt(numpw)))/2)
    plot(10,5,type="n",axes=FALSE,ann=FALSE,xlim=c(0, 10),ylim = c(0,10))
    text(-9,10,paste(squareM,"edge:",edgeM),cex=2.0,pos=4,xpd=NA)
    text(-9,9,paste("module:",moduleM),cex=2.0,pos=4,xpd=NA)
    text(-9,8,paste("modAUC1:",sprintf("%.3f",wa$modAUC1)),cex=2,pos=4,xpd=NA)
    text(-9,7,paste("modAUC2:",sprintf("%.3f",wa$modAUC2)),cex=2,pos=4,xpd=NA)
    text(-9,6,paste("diffME:",sprintf("%.3f",wa$diffME)),cex=2,pos=4,xpd=NA)
    text(-9,5,paste("sigenrich:",sprintf("%.3f",wa$sigenrich)),cex=2,pos=4,xpd=NA)
  }#end turq.data dimension check
}#end function
