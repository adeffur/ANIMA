#DEDICATED boxplot 4 function

probe_boxplot5<-function(squareC,generegex,probetype=FALSE,miny=6,maxy=14,orderP=c(1,2,3,4)){
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
  print("nrow(res)")
  print(nrow(res))
  print("res")
  print(res)
  res<-res[order(res$edge),]#key!
  print(res)
  classvar1<-res$contrastvar[1]
  classvar2<-res$contrastvar[3]
  c1con<-unlist(strsplit(res$contrast[1]," - "))[2]
  c1case<-unlist(strsplit(res$contrast[1]," - "))[1]
  c2con<-unlist(strsplit(res$contrast[3]," - "))[2]
  c2case<-unlist(strsplit(res$contrast[3]," - "))[1]
  
  subsets<-eval(parse(text=paste("list(pd$",classvar2,"=='",c2con,"',pd$",classvar2,"=='",c2case,"',pd$",classvar1,"=='",c1con,"',pd$",classvar1,"=='",c1case,"',pd$",classvar1,"=='",c1con,"'|","pd$",classvar1,"=='",c1case,"')",sep="")))
  #print("subsets")
  #print(subsets)
  ##end new subsets
  
  #loop over edges to make boxplot data structure
  boxplotdatalist<-list()
  caseslist<-list()
  controlslist<-list()
  groupnamelist<-list()
  
  for (edgeC in 1:4){
    
  print(paste("edge:",edgeC))
  subset=subsets[[as.numeric(edgeC)]]
  #print(subset)
  data.qs<-data.q[,subset]
  pd<-pData(data.q)[subset,]
  #get probes for subset (re-use subsetting function from cellcor) and boxplot annotation information
  query0<-paste("MATCH (s:SYMBOL)-[r]-(p:PROBE) WHERE s.name =~ '",generegex,"' AND p.square= '",squareC,"' AND p.edge= '",edgeC,"' RETURN s.name AS Symbol, p.edge as Edge, p.name AS probe, p.logfc AS logfc, p.adjPVAL AS pval, p.contrastvar AS contrastvar, p.contrast AS contrast",sep="")
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
  #print("res")
  #print(res)
  
  if(is.null(res)){print(paste("no result for",squareC))}
  
  if(!is.null(res)){
  classifier<-res0$contrastvar[1]
  print("classifier")
  print(classifier)
  case<-unlist(strsplit(res0$contrast,"-"))[1]
  print("case")
  print(case)
  control<-unlist(strsplit(res0$contrast[1],"-"))[2]
  print("control")
  print(control)
  boxplotdatalist[[edgeC]]<-exprs(data.qs)[unique(res$probeID),]
  
  caseslist[[edgeC]]<-which(eval(parse(text=paste("pd$",classifier,sep="")))==case)
  print("caseslist[[edgeC]]")
  print(caseslist[[edgeC]])
  
  controlslist[[edgeC]]<-which(eval(parse(text=paste("pd$",classifier,sep="")))==control)
  print("controlslist[[edgeC]]")
  print(controlslist[[edgeC]])
  
  groupnamelist[[edgeC]]<-list(control,case)
  
  
  }#endif
  }#end edgewise
  
  print("intermediate debug")
  
  print("caseslist")
  print(caseslist)
  print("controlslist")
  print(controlslist)
  print("groupnamelist")
  print(groupnamelist)

  if(!is.null(res)&nrow(res)>1){
   #par(mfrow=c(max(3,ceiling(nrow(res)/3)),3),mar=c(12,10,6,3))
    #par(mfrow=c(10,3),mar=c(10,9,4,2))
    #par(mfrow=c(min(2,ceiling(nrow(res)/3)),3),mar=c(10,8,4,2),mgp=c(3,1,0))
    #par(mfrow=c(nrow(res),1),mar=c(4,3,2,1))
    #par(mfrow=c(ceiling(nrow(res)/3),3),mar=c(10,8,4,2),mgp=c(3,1,0))
    
    allgd<-data.frame(probe = character(0), symbol = character(0), df = character(0), y = numeric(0) , x = character(0))
    
    
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
    
    gdlist<-list(
      "group1"=data.frame(y=g1d,x=g1dname),
      "group2"=data.frame(y=g2d,x=g2dname),
      "group3"=data.frame(y=g3d,x=g3dname),
      "group4"=data.frame(y=g4d,x=g4dname))
    
    
    #print(gdlist)
    #stats<-kruskal.test(gdlist)
 
    #boxplot(gdlist[orderP],ylim=c(miny,maxy),col=c("white","white","white","white"),ylab="expression:\nlog2 intensity\n",names=c(g1dname,g2dname,g3dname,g4dname)[orderP],las=2,main=paste(res$Symbol[row],"(",res$probeID[row],")\n",stats$method,"\n P value =",sprintf("%.2f",stats$p.value)),outpch=NA,boxwex=.5,border="gray50",cex.axis=1.5,cex.lab=1.5)
    #stripchart(gdlist[orderP],vertical=TRUE,add=TRUE,pch=21,bg=c("red"),col=c("black"),method = "jitter",cex=1.3)
    print("class(gdlist[[1]])")
    print(class(gdlist[[1]]))
    
    newgd<-bind_rows(gdlist[orderP],.id="df")
    newgd2<-cbind(rep(res$Symbol[row],nrow(newgd)),newgd)
    newgd3<-cbind(rep(res$probe[row],nrow(newgd2)),newgd2)
    print(head(newgd3))
    
    allgd<-rbind(allgd,newgd3)
    
    #boxplot(g1d,g2d,g3d,g4d,ylim=c(miny,maxy),col=c("white","white","white","white"),ylab="expression\n",names=c(g1dname,g2dname,g3dname,g4dname),las=2,main=paste(res$Symbol[row],"(",res$probeID[row],")\n",stats$method,"\n P value =",sprintf("%.2f",stats$p.value)),outpch=NA,boxwex=.5,border="gray50",cex.axis=1.5,cex.lab=1.5)
    #stripchart(list(g1d,g2d,g3d,g4d),vertical=TRUE,add=TRUE,pch=21,bg=c("red"),col=c("black"),method = "jitter",cex=1.3)
     }#real row-wise
    
    print("final debug")
    colnames(allgd)<-c("probe","symbol","df","y","Category")
    print(head(allgd))
    print("str(allgd)")
    print("class(allgd$Category)")
    print(str(allgd))
    print(class(allgd$Category))
    levelvec<-c(g1dname,g2dname,g3dname,g4dname)
    levelvec<-levelvec[orderP]
    allgd$Category<-factor(allgd$Category,levels=levelvec)
    #levels(allgd$Category)<-levelvec
    print("levelvec")
    print(levelvec)
    
    #boxplot comes here
    nplot<-ggplot(allgd, aes(Category,y)) +
      geom_boxplot(fill="lightyellow") +
      geom_jitter(colour="blue",alpha=.4,size=2) +
      ylim(miny,maxy) +
      labs(y="expression:\nlog2intensity") +
      stat_compare_means() +
      facet_wrap(c("symbol","probe"),ncol=4) +
      theme_bw() +
      theme(strip.text = element_text(size = rel(1.5))) +
      theme(axis.text = element_text(size = rel(1.5))) +
      theme(axis.title.y = element_text(size = rel(2.5), angle = 90, color = "slategrey",margin = margin(t = 0, r = 50, b = 50, l = 50))) +
      theme(axis.title.x = element_text(size = rel(2.5), color = "slategrey",margin = margin(t = 50, r = 50, b = 0, l = 50))) + 
      theme(legend.position="none")
      # theme(legend.text = element_text(size = rel(1.2))) +
      # theme(legend.title = element_text(size = rel(1.2))) +
      # theme(legend.text=element_text(size=rel(1.4))) 
      
    
    
    print(nplot)
    
   
 
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
