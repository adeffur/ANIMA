#############
#META FUNCTION####
#########

moduleMeta<-function(type="sets",extent=111,setlist=c("berry.train","berry.test","berry.val"),edgelist=c(1,2,5,3,4),plotsimple=TRUE){
  n.edges<-length(edgelist)
  n.squares<-length(setlist)
  
  if(type=="sets"){
    #loop over datasets
    par(mfrow=c(1,1))
    hmm.mat<-matrix(nrow=n.edges,ncol=258,data=NA)
    hmm.mat.list<-list()
    hmm.mat.col<-matrix(nrow=n.edges,ncol=258,data=NA)
    hmm.mat.col.list<-list()
    for (set in setlist) {
      if(plotsimple==TRUE){
      #plot scaffold
      par(bg="white",mar=c(4,6,2,2))
      finallab<-character()
      baselabel<-paste("M",9:1,sep=".")
      for(lab in baselabel){finallab<-c(finallab,paste(lab,edgelist[n.edges:1],sep="_edge_"))}
      plot(c(0, extent), c(1, 9*n.edges), type= "n", xlab = "", ylab = "",frame=FALSE,axes=FALSE,main=set,xpd=NA)
      axis(1, at=1:extent, labels=1:extent,las=2,cex.axis=.6)
      axis(2, at=1:(9*n.edges), labels=finallab,cex.axis=.6,las=2)
      
      abline(h=seq(0.5,9*n.edges+0.5,n.edges),col="gray",lwd=1)
      abline(v=seq(0.5,extent+0.5,10),col="gray",lwd=1)
      abline(v=seq(5.5,extent+0.5,10),col="gray",lwd=1,lty=2)
      #query.base<-query
      #edgewise loop, with same M adjacent!
      }
      
      for (edge in 1:n.edges){
        query.base<-paste("MATCH (b:baylor) WHERE b.square = '",set,"' AND b.edge =",edgelist[edge]," OPTIONAL MATCH (b:baylor)-[r]-(n:wgcna) WHERE b.square = '",set,"' AND b.edge =",edgelist[edge]," RETURN b.name AS module, b.diffEX as diffEX, n.name AS wgcnacol, r.qvalue as qval",sep="")
        #print(query.base)
        res<-cypher(graph,query.base)
        res<-res[order(res$qval),]
        res<-res[!duplicated(res$module),]
      
        if(nrow(res)>258){res<-res[which(res$module!='M8.90'),]}#hardcoded fix for M8.90 bug
        resA<-gsub("M", "", res$module)
        resAA<-unlist(strsplit(resA,"[.]"))[seq(2,length(resA)*2,2)]
        resAAA<-sprintf("%03s",resAA)
        resAAAA<-unlist(strsplit(resA,"[.]"))[seq(1,length(resA)*2,2)]
        resAAAAA<-paste("M",resAAAA,".",resAAA,sep="")
        res<-res[order(resAAAAA),]
        res$wgcnacol[which(is.na(res$wgcnacol))]<-"gray"
        
        hmdata.col<-res$wgcnacol
        hmdata<-res$diffEX
        hmm.mat[edge,]<-hmdata
        hmm.mat.col[edge,]<-hmdata.col
        dimensions<-strsplit(res$module,".",fixed=T)
        coords.rows=as.numeric(unlist(lapply(dimensions,function(x){strsplit(x[1],"M",fixed=T)[[1]][2]})))*n.edges-(2+n.edges-edge)
        #
        coords.cols=unlist(lapply(dimensions,function(x){x[2]}))
        coords.all=cbind(as.numeric(coords.rows),as.numeric(coords.cols))
        coords.all[,1]<-(coords.all[,1]-9*n.edges+1)*-1
        for (cgc in 1:nrow(coords.all)){
          cord<-c(coords.all[cgc,2],coords.all[cgc,1])
          #print(cord)
          modc<-numbers2colors(hmdata,centered=TRUE)[cgc]
          #print(modc)
          if(plotsimple==TRUE){rect(cord[1]-0.3,cord[2]-0.3,cord[1]+0.3,cord[2]+0.3,col=modc,border=NA)}
        }
        
      }#end edgewise
      #   labfin<-character()
      #   for (ju in 1:9){
      #     int2<-grep(paste("M",ju,".",sep=""),res$module)
      #     int3<-1:length(int2)
      #     int4<-paste("M",ju,".",int3,sep="")
      #     labfin<-c(labfin,int4)
      #   }
      #   
      #   dimnames(hmm.mat)<-list(paste("edge",edgelist,sep=""),labfin)
      dimnames(hmm.mat)<-list(unlist(edgelist),res$module)
      dimnames(hmm.mat.col)<-list(unlist(edgelist),res$module)
      hmm.mat.list[[set]]<-hmm.mat
      hmm.mat.col.list[[set]]<-hmm.mat.col
    }#end sets
    
  }else if(type=="edges"){
    par(mfrow=c(1,1))
    hmm.mat<-matrix(nrow=n.squares,ncol=258,data=NA)
    hmm.mat.list<-list()
    for (edge in 1:n.edges){
      #plot scaffold
      par(bg="white",mar=c(4,8,2,2))
      finallab<-character()
      baselabel<-paste("M",9:1,sep=".")
      for(lab in baselabel){finallab<-c(finallab,paste(lab,setlist[n.squares:1],sep="_square_"))}
      plot(c(0, extent), c(1, 9*n.squares), type= "n", xlab = "", ylab = "",frame=FALSE,axes=FALSE,main=edgelist[edge],xpd=NA)
      axis(1, at=1:extent, labels=1:extent,las=2,cex.axis=.6)
      axis(2, at=1:(9*n.squares), labels=finallab,cex.axis=.6,las=2)
      
      abline(h=seq(0.5,9*n.squares+0.5,n.squares),col="gray",lwd=1)
      abline(v=seq(0.5,extent+0.5,10),col="gray",lwd=1)
      abline(v=seq(5.5,extent+0.5,10),col="gray",lwd=1,lty=2)
      
      for (set in 1:n.squares){
        query.base<-paste("MATCH (b:baylor) WHERE b.square = '",setlist[[set]],"' AND b.edge =",edgelist[edge]," RETURN b.name AS module, b.diffEX as diffEX",sep="")
        
        res<-cypher(graph,query.base)
        resA<-gsub("M", "", res$module)
        resAA<-unlist(strsplit(resA,"[.]"))[seq(2,length(resA)*2,2)]
        resAAA<-sprintf("%03s",resAA)
        resAAAA<-unlist(strsplit(resA,"[.]"))[seq(1,length(resA)*2,2)]
        resAAAAA<-paste("M",resAAAA,".",resAAA,sep="")
        res<-res[order(resAAAAA),]
        
        hmdata<-rep(NA,258)
        hmdata<-res$diffEX
        hmm.mat[set,]<-hmdata
        dimensions<-strsplit(res$module,".",fixed=T)
        coords.rows=as.numeric(unlist(lapply(dimensions,function(x){strsplit(x[1],"M",fixed=T)[[1]][2]})))*n.squares-(2+n.squares-set)
        #
        coords.cols=unlist(lapply(dimensions,function(x){x[2]}))
        coords.all=cbind(as.numeric(coords.rows),as.numeric(coords.cols))
        coords.all[,1]<-(coords.all[,1]-9*n.squares+1)*-1
        for (cgc in 1:nrow(coords.all)){
          cord<-c(coords.all[cgc,2],coords.all[cgc,1])
          #print(cord)
          modc<-numbers2colors(hmdata,centered=TRUE)[cgc]
          #print(modc)
          rect(cord[1]-0.3,cord[2]-0.3,cord[1]+0.3,cord[2]+0.3,col=modc,border=NA)
        }
        
      }#end sets
      #     labfin<-character()
      #     for (ju in 1:9){
      #       int2<-grep(paste("M",ju,".",sep=""),res$module)
      #       int3<-1:length(int2)
      #       int4<-paste("M",ju,".",int3,sep="")
      #       labfin<-c(labfin,int4)
      #     }
      #     dimnames(hmm.mat)<-list(unlist(setlist),labfin)
      dimnames(hmm.mat)<-list(unlist(setlist),res$module)
      hmm.mat.list[[edge]]<-hmm.mat
    }#end edgewise
  }
  moduleMetaList<-list(hmm.mat.list,hmm.mat.col.list)
  return(moduleMetaList)
}

