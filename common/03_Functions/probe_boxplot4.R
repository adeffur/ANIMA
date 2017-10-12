#DEDICATED boxplot 4 function

probe_boxplot4<-function(squareC,generegex,probetype=FALSE,miny=6,maxy=14,orderP=c(1,2,3,4)){
  dataC=eval(parse(text=squareC))
  if(class(dataC)=="LumiBatch"){
    data.v<-lumiT(dataC,simpleOutput=FALSE,verbose=TRUE,ifPlot=T)
  }else{
      data.v<-lumiT(dataC,simpleOutput=FALSE,method="log2")
      }
  print("quantile normalisation")
  data.q<-lumiN(data.v,method="quantile")
  pd<-pData(data.q)
  print("orderP")
  orderP<-orderP
  #pseudo classvar2:control, classvar2:cases, classvar1:control, classvar1:case, classvar1:control|classvar1:case
  query<-paste("MATCH (n:wgcna {square:'",squareC,"', name:'blue'}) RETURN n.edge AS edge, n.contrastvar AS contrastvar, n.contrast AS contrast",sep="")
  res<-cypher(graph,query)
  classvar1<-res$contrastvar[1]
  classvar2<-res$contrastvar[3]
  c1con<-unlist(strsplit(res$contrast[1]," - "))[2]
  c1case<-unlist(strsplit(res$contrast[1]," - "))[1]
  c2con<-unlist(strsplit(res$contrast[3]," - "))[2]
  c2case<-unlist(strsplit(res$contrast[3]," - "))[1]
  subsets<-eval(parse(text=paste("list(pd$",classvar2,"=='",c2con,"',pd$",classvar2,"=='",c2case,"',pd$",classvar1,"=='",c1con,"',pd$",classvar1,"=='",c1case,"',pd$",classvar1,"=='",c1con,"'|","pd$",classvar1,"=='",c1case,"')",sep="")))
  
  ##end new subsets
  
  #loop over edges to make boxplot data structure
  boxplotdatalist<-list()
  caseslist<-list()
  controlslist<-list()
  groupnamelist<-list()
  
  for (edgeC in 1:4){ 
  subset=subsets[[as.numeric(edgeC)]]
  data.qs<-data.q[,subset]
  pd<-pData(data.q)[subset,]
  #get probes for subset (re-use subsetting function from cellcor) and boxplot annotation information
  query0<-paste("MATCH (s:SYMBOL)-[r]-(p:PROBE) WHERE s.name =~ '",generegex,"' AND p.square= '",squareC,"' AND p.edge= ",edgeC," RETURN s.name AS Symbol, p.edge as Edge, p.name AS probe, p.logfc AS logfc, p.adjPVAL AS pval, p.contrastvar AS contrastvar, p.contrast AS contrast",sep="")
  res0<-cypher(graph,query0)
  
  if(probetype==TRUE){
  query<-paste("MATCH (p0:PROBETYPE)-[r0]-(s:SYMBOL) WHERE s.name =~ '",generegex,"' RETURN p0.name as probeID, s.name AS Symbol",sep="")
  res<-cypher(graph,query)
  res<-res[order(res$probeID),]
  }else{query<-paste("MATCH (p0:PROBE)-[r0]-(s:SYMBOL) WHERE s.name =~ '",generegex,"' RETURN p0.name as probeID, s.name AS Symbol",sep="")
  res<-cypher(graph,query)
  res<-res[order(res$probeID),]}
  res<-res[which(res$probeID%in%res0$probe),]
  res<-unique(res[order(res$Symbol),])
  if(is.null(res)){print(paste("no result for",squareC))}
  if(!is.null(res)){
  classifier<-res0$contrastvar[1]
  case<-unlist(strsplit(res0$contrast,"-"))[1]
  control<-unlist(strsplit(res0$contrast[1],"-"))[2]
  boxplotdatalist[[edgeC]]<-exprs(data.qs)[unique(res$probeID),]
  caseslist[[edgeC]]<-which(eval(parse(text=paste("pd$",classifier,sep="")))==case)
  controlslist[[edgeC]]<-which(eval(parse(text=paste("pd$",classifier,sep="")))==control)
  groupnamelist[[edgeC]]<-list(control,case)
  }#endif
  }#end edgewise
  
  if(!is.null(res)&nrow(res)>1){
   #par(mfrow=c(max(3,ceiling(nrow(res)/3)),3),mar=c(12,10,6,3))
    #par(mfrow=c(10,3),mar=c(10,9,4,2))
    #par(mfrow=c(min(2,ceiling(nrow(res)/3)),3),mar=c(10,8,4,2),mgp=c(3,1,0))
    #par(mfrow=c(nrow(res),1),mar=c(4,3,2,1))
    par(mfrow=c(ceiling(nrow(res)/3),3),mar=c(10,8,4,2),mgp=c(3,1,0))
  for (row in 1:nrow(boxplotdatalist[[1]])){
    #group1data
    #[order(rownames(boxplotdatalist[[1]])),] [order(rownames(boxplotdatalist[[4]])),] [order(rownames(boxplotdatalist[[3]])),] [order(rownames(boxplotdatalist[[2]])),]
    g1d<-boxplotdatalist[[1]][row,controlslist[[1]]]
    g1dname<-paste(groupnamelist[[1]][[1]],"\n",groupnamelist[[3]][[1]])
    g2d<-boxplotdatalist[[4]][row,controlslist[[4]]]
    g2dname<-paste(groupnamelist[[1]][[2]],"\n",groupnamelist[[4]][[1]])
    g3d<-boxplotdatalist[[3]][row,caseslist[[3]]]
    g3dname<-paste(groupnamelist[[2]][[1]],"\n",groupnamelist[[3]][[2]])
    g4d<-boxplotdatalist[[2]][row,caseslist[[2]]]
    g4dname<-paste(groupnamelist[[2]][[2]],"\n",groupnamelist[[4]][[2]])
    
    gdlist<-list(g1d,g2d,g3d,g4d)
    stats<-kruskal.test(gdlist)
 
    boxplot(gdlist[orderP],ylim=c(miny,maxy),col=c("white","white","white","white"),ylab="expression:\nlog2 intensity\n",names=c(g1dname,g2dname,g3dname,g4dname)[orderP],las=2,main=paste(res$Symbol[row],"(",res$probeID[row],")\n",stats$method,"\n P value =",sprintf("%.2f",stats$p.value)),outpch=NA,boxwex=.5,border="gray50",cex.axis=1.5,cex.lab=1.5)
    stripchart(gdlist[orderP],vertical=TRUE,add=TRUE,pch=21,bg=c("red"),col=c("black"),method = "jitter",cex=1.3)
    
    #boxplot(g1d,g2d,g3d,g4d,ylim=c(miny,maxy),col=c("white","white","white","white"),ylab="expression\n",names=c(g1dname,g2dname,g3dname,g4dname),las=2,main=paste(res$Symbol[row],"(",res$probeID[row],")\n",stats$method,"\n P value =",sprintf("%.2f",stats$p.value)),outpch=NA,boxwex=.5,border="gray50",cex.axis=1.5,cex.lab=1.5)
    #stripchart(list(g1d,g2d,g3d,g4d),vertical=TRUE,add=TRUE,pch=21,bg=c("red"),col=c("black"),method = "jitter",cex=1.3)
  }#real
 
  par(mfrow=c(1,1))
  }else if (nrow(res)==1){
    par(mfrow=c(1,1))
    #group1data
    g1d<-boxplotdatalist[[1]][controlslist[[1]]]
    g1dname<-paste(groupnamelist[[1]][[1]],"\n",groupnamelist[[3]][[1]])
    g2d<-boxplotdatalist[[4]][controlslist[[4]]]
    g2dname<-paste(groupnamelist[[1]][[2]],"\n",groupnamelist[[4]][[1]])
    g3d<-boxplotdatalist[[3]][caseslist[[3]]]
    g3dname<-paste(groupnamelist[[2]][[1]],"\n",groupnamelist[[3]][[2]])
    g4d<-boxplotdatalist[[2]][caseslist[[2]]]
    g4dname<-paste(groupnamelist[[2]][[2]],"\n",groupnamelist[[4]][[2]])
    
    gdlist<-list(g1d,g2d,g3d,g4d)
    stats<-kruskal.test(gdlist)
    
    boxplot(g1d,g2d,g3d,g4d,ylim=c(miny,maxy),col=c("white","white","white","white"),ylab="expression",names=c(g1dname,g2dname,g3dname,g4dname),las=2,main=paste(res$Symbol[row],"(",res$probe[order(res$probe)][row],")\n",stats$method,"\n P value =",sprintf("%.2f",stats$p.value)),outpch=NA,boxwex=.5,border="gray50")
    
  }#end else for special case of nrow(res==1)
  #end if 2
  return(nrow(res))
}#end function
