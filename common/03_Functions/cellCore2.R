######
cellCore2<-function(study="Berry",squareC,edgeC,cellGroup,pwm,PalWang,plotCell=TRUE){
  library(fastcluster)
  dataC=eval(parse(text=squareC))
  if(class(dataC)=="LumiBatch"){data.v<-lumiT(dataC,simpleOutput=FALSE)}else{data.v<-lumiT(dataC,simpleOutput=FALSE,method="log2")}
  data.q<-lumiN(data.v,method="quantile")
  pd<-pData(data.q)
  #subsetting####
 
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
  
  #get the individual celltypes####
  print("query1: get cells for group")
  query<-paste("MATCH (c:CELL {name:'",cellGroup,"'})-[r0]-(c1:cellEx) RETURN DISTINCT c1.name AS cellgrouplist",sep="")
  cellgrouplist<-cypher(graph,query)
  cellgrouplist<-cellgrouplist$cellgrouplist
  cellscorelist<-list()
  for (cellC in cellgrouplist){
    print(cellC)
    query<-paste("MATCH q1=(c:cellEx {name:'",cellC,"'})-[r0]-(s:SYMBOL)-[r1]-(p:PROBE {square:'",squareC,"',edge:",edgeC,"})-[r2]-(n:wgcna {square:'",squareC,"',edge:",edgeC,"}) RETURN DISTINCT p.name as cellprobename",sep="")
    cellprobes<-cypher(graph,query)
    cellprobes<-unique(cellprobes$cellprobename)
    
    #get nuIDs for cell probes and tightly correlated ones, regardless which module####
    query<-paste("MATCH q1=(c:cellEx {name:'",cellC,"'})-[r0]-(s:SYMBOL)-[r1]-(p:PROBE {square:'",squareC,"',edge:",edgeC,"})-[r2]-(n:wgcna {square:'",squareC,"',edge:",edgeC,"}) WITH p OPTIONAL MATCH (p:PROBE {square:'",squareC,"',edge:",edgeC,"})-[r2]-(p2:PROBE {square:'",squareC,"',edge:",edgeC,"}) WHERE r2.TOMweight > .1 WITH collect(p) + collect(p2) AS x  UNWIND x AS y WITH y MATCH (y:PROBE  {square:'",squareC,"',edge:",edgeC,"})-[rx]-(n:wgcna) RETURN DISTINCT y.name AS probename",sep="")
    probes<-cypher(graph,query)
    nuID<-unique(probes$probename)
    nuID<-nuID[which(!is.na(nuID))]
    nuID<-unique(c(nuID,cellprobes))
    print(paste("query2: nuID",length(nuID)))
    
    if(length(nuID)>5){
      #subset the data and make correlation matrix####
      turq.data<-exprs(data.qs)[nuID,]
      cormat<-cor(t(turq.data),method="pearson")
      dendro.manual<-as.dendrogram(fastcluster::hclust(dist(cormat)))
      idvec<-order.dendrogram(dendro.manual)
      dimnames(cormat)<-list(NULL,NULL)
      if(nrow(cormat)>5){
        #wgcna cols####
        print("start query 3")
        
        query<-paste("MATCH q1=(c:cellEx {name:'",cellC,"'})-[r0]-(s:SYMBOL)-[r1]-(p:PROBE {square:'",squareC,"',edge:",edgeC,"})-[r2]-(n:wgcna {square:'",squareC,"',edge:",edgeC,"}) WITH p OPTIONAL MATCH (p:PROBE {square:'",squareC,"',edge:",edgeC,"})-[r2]-(p2:PROBE {square:'",squareC,"',edge:",edgeC,"}) WHERE r2.TOMweight > .1 WITH collect(p) + collect(p2) AS x  UNWIND x AS y WITH y MATCH (y:PROBE {square:'",squareC,"',edge:",edgeC,"})-[r3]-(n:wgcna {square:'",squareC,"',edge:",edgeC,"}) RETURN DISTINCT y.name AS probe, n.name AS wgcnaname",sep="")
        wgcnalist0<-unique(cypher(graph,query))
        addnames<-cbind(nuID[!nuID%in%wgcnalist0$probe],rep("NONE",length(nuID[!nuID%in%wgcnalist0$probe])))
        colnames(addnames)<-c("probe","wgcnaname")
        wgcnalist<-rbind(wgcnalist0,addnames)
        print(paste("wgcnalist dim:",dim(wgcnalist)))
        
        wgcnacolmap<-as.matrix(wgcnalist)
        #wgcnacolmap<-wgcnacolmap[idvec,]
        module.names<-unique(wgcnacolmap[,2])
        nmod<-length(module.names)
        wcmat<-matrix(nrow=nrow(cormat),ncol=nmod,data=0)
        dimnames(wcmat)<-list(wgcnacolmap[,1],module.names)
        for (mod in 1:nmod){
          tv<-which(wgcnacolmap[,2]==module.names[mod])
          wcmat[tv,mod]<-1
        }
        print(paste("wcmat dim:",dim(wcmat)))
        
        #PATHWAYS####
        print(paste("start query 4;pwm=",pwm,";PalWang=",PalWang,sep=""))
        if (pwm=="wgcna"){
          query<-paste("MATCH (c:CELL)-[r0]-(c1:cellEx {name:'",cellC,"'})-[r1]-(s:SYMBOL)-[r2]-(p:PROBE {square:'",squareC,"',edge:",edgeC,"})-[r3]-(n:wgcna {square:'",squareC,"',edge:",edgeC,"}) WITH n MATCH (n:wgcna {square:'",squareC,"',edge:",edgeC,"})-[r4]-(pw) WHERE (pw:reactomePW OR pw:ImmunePW) RETURN DISTINCT n.name as module, pw.name as pathway,r4.qvalue as qvalue",sep="")
          
          pwlist<-cypher(graph,query)
          try(pwlist<-pwlist[order(pwlist$qvalue,decreasing=FALSE),])
          module.names<-unique(pwlist[,1])
          nmod<-length(module.names)
          pwx<-data.frame()
          for (mod in 1:nmod){
            pw0<-pwlist[which(pwlist$module==module.names[mod]),][1:2,]
            pwx<-rbind(pwx,pw0)
          }
          pw<-unique(pwx$pathway)
          numpw<-length(pw)
          
          pwchar<-paste("'",pw,"'",collapse=",",sep="")
          query<-paste("MATCH (pw)-[r1]-(s:SYMBOL)-[r2]-(p:PROBE {square:'",squareC,"',edge:",edgeC,"}) WHERE pw.name IN [",pwchar,"] RETURN DISTINCT p.name AS probe, pw.name as pathway",sep="")
          probelist<-cypher(graph,query)
          
          pwmatrix<-matrix(nrow=nrow(cormat),ncol=numpw,data=0)
          dimnames(pwmatrix)<-list(nuID,pw)
          
          for(path in 1:numpw){
            pwmatrix[,path]<-nuID%in%probelist$probe[which(probelist$pathway==pw[path])]
          }
          pwmatrix<-pwmatrix[,ncol(pwmatrix):1]
          
          #end wgcna method
          
          
        }else{
          if(PalWang==TRUE){
            query<-paste("MATCH q1=(c:cellEx {name:'",cellC,"'})-[r0]-(s:SYMBOL)-[r1]-(p:PROBE {square:'",squareC,"',edge:",edgeC,"})-[r2a]-(n:wgcna {square:'",squareC,"',edge:",edgeC,"}) WITH p OPTIONAL MATCH (p:PROBE)-[r2]-(p2:PROBE) WHERE r2.TOMweight > .1 WITH collect(p) + collect(p2) AS x  UNWIND x AS y WITH y MATCH (y:PROBE {square:'",squareC,"',edge:",edgeC,"})-[r5]-(s:SYMBOL)-[r6]-(pw) WHERE (pw:ImmunePW OR pw:reactomePW OR pw:PalWangPW) RETURN y.name as probe, pw.name as pw",sep="")
          }else{
            query<-paste("MATCH q1=(c:cellEx {name:'",cellC,"'})-[r0]-(s:SYMBOL)-[r1]-(p:PROBE {square:'",squareC,"',edge:",edgeC,"})-[r2a]-(n:wgcna {square:'",squareC,"',edge:",edgeC,"}) WITH p OPTIONAL MATCH (p:PROBE)-[r2]-(p2:PROBE) WHERE r2.TOMweight > .1 WITH collect(p) + collect(p2) AS x  UNWIND x AS y WITH y MATCH (y:PROBE {square:'",squareC,"',edge:",edgeC,"})-[r5]-(s:SYMBOL)-[r6]-(pw) WHERE (pw:ImmunePW OR pw:reactomePW) RETURN y.name as probe, pw.name as pw",sep="")
          }
          pwlist<-cypher(graph,query)
          pwlist<-unique(pwlist)
          print(paste("dim pwlist",dim(pwlist)))
          
          if(is.null(pwlist)){
            query<-paste("MATCH q1=(c:cellEx {name:'",cellC,"'})-[r0]-(s:SYMBOL)-[r1]-(p:PROBE {square:'",squareC,"',edge:",edgeC,"})-[r2a]-(n:wgcna {square:'",squareC,"',edge:",edgeC,"}) WITH p OPTIONAL MATCH (p:PROBE)-[r2]-(p2:PROBE) WHERE r2.TOMweight > .1 WITH collect(p) + collect(p2) AS x  UNWIND x AS y WITH y MATCH (y:PROBE {square:'",squareC,"',edge:",edgeC,"})-[r5]-(s:SYMBOL)-[r6]-(pw) WHERE (pw:ImmunePW OR pw:reactomePW OR pw:PalWangPW) RETURN y.name as probe, pw.name as pw",sep="")
            pwlist<-cypher(graph,query)
            pwlist<-unique(pwlist)
          }#add PalWang if
          #now we should have a valid pwlist (after 2 tries)
          
          
          if(!is.null(pwlist)&dim(pwlist)[1]>2){
            
            pwnames<-unique(pwlist$pw)
            print("dim(pwlist")
            print(dim(pwlist))
            
            print("debug pwnames")
            print(pwnames[1:5])
            npw<-length(pwnames)
            pwdf<-data.frame()
            for (path in 1:npw){
              pwdf<-rbind(pwdf,c(pwnames[path],nrow(pwlist[which(pwlist$pw==pwnames[path]),])))
            }#short for loop
            colnames(pwdf)<-c("pathway","nprobes")
            pwdf$nprobes<-as.numeric(pwdf$nprobes)
            rowchoose<-min(20,nrow(pwdf))
            pwdf<-pwdf[order(pwdf[,2],decreasing = TRUE),][1:rowchoose,]
            pw<-unique(pwdf$pathway)
            numpw<-length(pw)
            
            pwchar<-paste("'",pw,"'",collapse=",",sep="")
            #query<-paste("MATCH (pw)-[r1]-(s:SYMBOL)-[r2]-(p:PROBE {square:'",squareC,"',edge:",edgeC,"}) WHERE pw.name = 'Pyruvate metabolism and Citric Acid *TCA* cycle' RETURN DISTINCT p.name AS probe, pw.name as pathway",sep="")
            
            #query<-paste("MATCH (pw)-[r1]-(s:SYMBOL)-[r2]-(p:PROBE {square:'",squareC,"',edge:",edgeC,"}) WHERE pw.name IN [",pwchar,"] RETURN DISTINCT p.name AS probe, pw.name as pathway",sep="")
            query<-paste("MATCH (pw)-[r1]-(s:SYMBOL)-[r2]-(p:PROBE {square:'",squareC,"',edge:",edgeC,"})-[r2a]-(n:wgcna {square:'",squareC,"',edge:",edgeC,"}) WHERE pw.name =~ ",gsub("\\*","[\\*]",gsub("','","|",pwchar))," RETURN DISTINCT p.name AS probe, pw.name as pathway",sep="")
            
            probelist<-cypher(graph,query)
            probelist<-unique(probelist)
            
            pwmatrix<-matrix(nrow=nrow(cormat),ncol=numpw,data=0)
            dimnames(pwmatrix)<-list(nuID,pw)
            
            for(path in 1:numpw){
              unique.probes<-unique(probelist$probe)
              pwmatrix[,path]<-nuID%in%probelist$probe[which(probelist$pathway==pw[path])]
            }
            pwmatrix<-pwmatrix[,ncol(pwmatrix):1]
          }else{
            numpw=1
            pwmatrix<-matrix(nrow=nrow(cormat),ncol=1,data=0)
            dimnames(pwmatrix)<-list(nuID,"NONE")
          }#end second pwlist check after adding PalWang
          
        }#end non-wgna else
        
        #logfc####
        query<-paste("MATCH q1=(c:cellEx {name:'",cellC,"'})-[r0]-(s:SYMBOL)-[r1]-(p:PROBE {square:'",squareC,"',edge:",edgeC,"})-[r2a]-(n:wgcna {square:'",squareC,"',edge:",edgeC,"}) WITH p OPTIONAL MATCH (p:PROBE {square:'",squareC,"',edge:",edgeC,"})-[r2]-(p2:PROBE {square:'",squareC,"',edge:",edgeC,"}) WHERE r2.TOMweight > .1 WITH collect(p) + collect(p2) AS x  UNWIND x AS y WITH y MATCH (y:PROBE  {square:'",squareC,"',edge:",edgeC,"}) RETURN DISTINCT y.name as probe, y.logfc AS logfc ",sep="")
        
        logfclist0<-cypher(graph,query)
        addlogfc<-paste("'",paste(cellprobes[!cellprobes%in%logfclist0$probe],collapse="','"),"'",sep="")
        query<-paste("MATCH (p:PROBE {square:'",squareC,"',edge:",edgeC,"}) WHERE p.name IN [",addlogfc,"] RETURN DISTINCT p.name as probe, p.logfc AS logfc ",sep="")
        addlogfclist<-cypher(graph,query)
        logfclist<-rbind(logfclist0,addlogfclist)
        
        fc<-matrix(nrow=nrow(cormat),ncol=1,data=0)
        dimnames(fc)<-list(logfclist$probe,"logfc")
        fc[,1]<-logfclist$logfc
        fc2<-as.numeric(fc)
        total<-cbind(pwmatrix,fc)
        
        #the original cell-specific probes:
        
        cellprobematrix<-matrix(nrow=nrow(cormat),ncol=1,data=0)
        dimnames(cellprobematrix)<-list(nuID,"cellprobes")
        tv<-which(nuID%in%cellprobes)
        cellprobematrix[tv,]<-1
        total<-cbind(total,cellprobematrix)
        
        blues<-fc2<0
        reds<-fc2>=0
        mcolvec<-character(length=length(fc2))
        mcolvec[blues]<-"blue"
        mcolvec[reds]<-"red"
        
        mcolvec2<-cbind(row.names(fc),mcolvec)
        row.names(mcolvec2)<-row.names(fc)
        #mcolvec2[idvec,2]
        
        dimnames(cormat)<-list(NULL,NULL)
        #plot(annHeatmap2(cormat,dendrogram=list("status"="hidden","dendro"=dendro.manual),legend=TRUE,col=blueWhiteRed,scale="none",annotation=list(Row=list("data"=wcmat,"asIs"=TRUE))))
        if(plotCell==TRUE){
        res<-annHeatmap2(cormat,dendrogram=list("status"="hidden","dendro"=dendro.manual),legend=TRUE,col=blueWhiteRed,scale="none",annotation=list(Row=list("data"=wcmat,"asIs"=TRUE),Col=list("data"=total,"asIs"=TRUE,control=list("pch"=16,"col.pch"=mcolvec2[idvec,2],numfac=4))))
        #plot(res,widths=c(6,1),heights=c(6,1))
        #res2<-res
        res$layout$plot<-matrix(c(0,5,0,0,4,0,0,1,2,0,3,0),ncol=4,dimnames=list(c("","image","rowAnn"),c("","leg2","image","colAnn")))
        plot(res,widths=c(1.8,1,4,1.3*(ncol(wcmat)/sqrt(ncol(wcmat)))),heights=c(1.5,4,1.5*(numpw/sqrt(numpw))))
        plot(10,5,type="n",axes=FALSE,ann=FALSE,xlim=c(0, 10),ylim = c(0,10))
        text(-9,10,paste(squareC,"edge:",edgeC),cex=2.0,pos=4,xpd=NA)
        text(-9,9,paste("cell:",cellC),cex=2.0,pos=4,xpd=NA)
        }
        #end####    
        print("debug1: dim pwmatrix")
        
        scorepw<-colnames(pwmatrix)
        print(scorepw)
        cellscore<-vector()
        for (spw in scorepw){
          tv<-as.logical(total[,spw])
          pwlfc<-total[,"logfc"][tv]
          pwscore<-mean(pwlfc)
          print("debug2")
          print(pwscore)
          cellscore[spw]<-pwscore
        }
        if(!is.null(cellscore)&!is.na(cellscore)&!is.nan(cellscore)){
          cellscorelist[[cellC]]<-cellscore
        }#check NA, NaN, NULL
      }#end if length nuID check
    }#end nuid length check
    #Export cell column
    
  }#end cellgroup loop
  return(cellscorelist)
}#end function
