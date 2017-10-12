tbhivplot<-function(square){
  hivquery=paste("MATCH (p:PROBE {square:'",square,"',edge:3})-[r]-(s:SYMBOL) WHERE s.name =~ '(?i)GBP.*' RETURN s.name, p.name AS e3probe, p.logfc AS e3logfc, p.adjPVAL AS e3PVAL",sep="")
  hiv<-cypher(graph,hivquery)
  tbquery=paste("MATCH (p:PROBE {square:'",square,"',edge:1})-[r]-(s:SYMBOL) WHERE s.name =~ '(?i)GBP.*' RETURN s.name, p.name AS e1probe, p.logfc AS e1logfc, p.adjPVAL AS e1PVAL",sep="")
  tb<-cypher(graph,tbquery)
  tbhivquery=paste("MATCH (p:PROBE {square:'",square,"',edge:2})-[r]-(s:SYMBOL) WHERE s.name =~ '(?i)GBP.*' RETURN s.name, p.name AS e2probe, p.logfc AS e2logfc, p.adjPVAL AS e2PVAL",sep="")
  tbhiv<-cypher(graph,tbhivquery)
  
  merged<-merge(merge(tb,hiv,by.x="e1probe",by.y="e3probe",all.x=TRUE,all.y=TRUE),tbhiv,by.x="e1probe",by.y="e2probe",all.x=TRUE,all.y=TRUE)[,c(1,2,3,6,9)]
  
  total<-cbind(merged,merged$e3logfc+merged$e2logfc)
  colnames(total)<-c("Probe","Symbol","TB logfc","HIV logfc","TB logfc (HIVpos)", "TB-HIV logfc")
  
  par(mfrow=c(3,1))
  barplot(merged$e1logfc,names.arg=merged$s.name.x,col=c("red"),legend.text=c("TB only"),ylim=c(floor(min(tb$e1logfc)),ceiling(max(tb$e1logfc))),ylab="log2FC")
  barplot(merged$e3logfc,names.arg=merged$s.name.x,col=c("blue"),legend.text=c("HIV only"),ylim=c(floor(min(hiv$e3logfc)),ceiling(max(hiv$e3logfc))),ylab="log2FC")
  barplot(t(matrix(c(merged$e3logfc,merged$e2logfc,merged$e3logfc+merged$e2logfc),ncol=3)),beside=T,names.arg=merged$s.name.x,col=c("blue","red","yellow"),legend.text=c("HIV","TB","TB-HIV"),ylim=c(floor(min(cbind(merged$e3logfc,merged$e2logfc))),ceiling(max(merged$e3logfc+merged$e2logfc))),ylab="log2FC")
  
  #return total
  par(mfrow=c(1,1),mar=c(4,5,5,12),xpd=T)
  matplot(t(total[,3:6]),type=c("b"),pch=1:13,col=rainbow(7),lty=1,ylab="Log2 fold change",xaxt="n",main=paste(square,": TB-HIV GBP",sep=""))
  axis(1,at=1:4,labels=c("TB (E1)","HIV (E3)","TB in HIV (E2)","TB-HIV (additive effect)"))
  legend(4.2,max(total[,3:6]),total$Symbol,pch=1:13,col=rainbow(7),cex=1)
  
  
}#end function