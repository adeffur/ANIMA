#DEDICATED boxplot function

probe_boxplot<-function(study,squareC,edgeC,generegex){
  dataC=eval(parse(text=squareC))
  data.v<-lumiT(dataC,simpleOutput=FALSE)
  data.q<-lumiN(data.v,method="quantile")
  pd<-pData(data.q)
  # #subsetting####
  # if(study=="Berry"){
  #   subsets<-list(pd$sex=="Female",pd$sex=="Male",pd$class=="notActiveTB",pd$class=="activeTB",pd$class=="activeTB"|pd$class=="notActiveTB")
  #   subset=subsets[[as.numeric(edgeC)]]
  # }else if (study=="IMPI" & squareC == "square1"){
  #   subsets<-list(pd$HIV.Status=="negative",pd$HIV.Status=="positive",pd$class3=="notActiveTB",pd$class3=="activeTB",pd$class3=="activeTB"|pd$class3=="notActiveTB")
  #   subset=subsets[[as.numeric(edgeC)]]
  # }else if (study=="IMPI" & squareC == "square2"){
  #   subsets<-list(pd$HIV.Status=="negative",pd$HIV.Status=="positive",pd$class=="CON_healthy",pd$class=="CON_LTBI",pd$class=="CON_healthy"|pd$class=="CON_LTBI")
  #   subset=subsets[[as.numeric(edgeC)]]
  # }else if (study=="IMPI" & squareC == "square3"){
  #   subsets<-list(pd$HIV.Status=="negative",pd$HIV.Status=="positive",pd$class=="CON_PTB",pd$class=="TBPC",pd$class=="CON_PTB"|pd$class=="TBPC")
  #   subset=subsets[[as.numeric(edgeC)]]
  # }else if (study=="IMPI" & squareC == "square4"){
  #   subsets<-list(pd$HIV.Status=="negative",pd$HIV.Status=="positive",pd$Compartment=="blood",pd$Compartment=="fluid",pd$Compartment=="blood"|pd$Compartment=="fluid")
  #   subset=subsets[[as.numeric(edgeC)]]
  # }else if (study=="IMPI" & squareC == "square5"){
  #   subsets<-list(pd$EC=="Eff",pd$EC=="EC",pd$Compartment=="blood",pd$Compartment=="fluid",pd$Compartment=="blood"|pd$Compartment=="fluid")
  #   subset=subsets[[as.numeric(edgeC)]]
  # }else if (study=="IMPI" & squareC == "square6"){
  #   subsets<-list(pd$HIV.Status=="negative",pd$HIV.Status=="positive",pd$class=="CON_LTBI",pd$class=="CON_PTB",pd$class=="CON_LTBI"|pd$class=="CON_PTB")
  #   subset=subsets[[as.numeric(edgeC)]]
  # }else if (study=="IMPI" & squareC == "rec1"){
  #   subsets<-list(pd$HIV.Status=="negative",pd$HIV.Status=="positive",pd$class=="CON_healthy",pd$class=="CON_PTB",pd$class=="CON_healthy"|pd$class=="CON_PTB")
  #   subset=subsets[[as.numeric(edgeC)]]
  # }else if (study=="IMPI" & squareC == "rec2"){
  #   subsets<-list(pd$HIV.Status=="negative",pd$HIV.Status=="positive",pd$class=="CON_healthy",pd$class=="TBPC",pd$class=="CON_healthy"|pd$class=="TBPC")
  #   subset=subsets[[as.numeric(edgeC)]]
  # }
  
  ##New subsets
  
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
  subset=subsets[[as.numeric(edgeC)]]
  ##end new subsets
  
  data.qs<-data.q[,subset]
  pd<-pData(data.q)[subset,]
  #get probes for subset (re-use subsetting function from cellcor) and boxplot annotation information
  query<-paste("MATCH (s:SYMBOL)-[r]-(p:PROBE) WHERE s.name =~ '",generegex,"' AND p.square= '",squareC,"' AND p.edge= ",edgeC," RETURN s.name AS Symbol, p.name AS probe, p.logfc AS logfc, p.adjPVAL AS pval, p.contrastvar AS contrastvar, p.contrast AS contrast",sep="")
  res<-cypher(graph,query)
  res<-res[order(res$probe),]
  classifier<-res$contrastvar[1]
  case<-unlist(strsplit(res$contrast," - "))[1]
  control<-unlist(strsplit(res$contrast," - "))[2]
  
  boxplotdata<-exprs(data.qs[res$probe,])
  cases<-which(eval(parse(text=paste("pd$",classifier,sep="")))==case)
  controls<-which(eval(parse(text=paste("pd$",classifier,sep="")))==control)
  
  par(mfrow=c(min(3,ceiling(nrow(boxplotdata)/3)),4))
  for (row in 1:nrow(boxplotdata)){
    casedata<-boxplotdata[row,cases]
    controldata<-boxplotdata[row,controls]
    boxplot(controldata,casedata,col=c("white","white"),ylab="expression",xlab=paste(res$Symbol[row]," (",res$probe[order(res$probe)][row],")",sep=""),names=c(control,case),main=paste(squareC,"edge",edgeC,"\n",res$contrast[1],"\n","(Log2FC:",sprintf("%.2f",res$logfc[row]),"q-value:",sprintf("%.3f",res$pval[row]),")"),outpch=NA,boxwex=.4,border="gray50")
    stripchart(list(controldata,casedata),vertical=TRUE,add=TRUE,pch=21,bg=c("black"),col=c("black"),method = "jitter",cex=1.4)
  }
  par(mfrow=c(1,1))
  
}#end function
