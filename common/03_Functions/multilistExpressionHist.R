multilistExpressionHist<-function(cellC,squareC,edgeC){
  
  
    query1<-paste("MATCH q1=(c1:CELL {name:'",cellC,"'})-[r1]-(c2:cellEx)-[r2]-(s:SYMBOL)-[r3]-(p:PROBE {square:'",squareC,"',edge:",edgeC,"}) RETURN DISTINCT s.name AS SYMBOL, p.name AS probe, p.logfc AS logfc, p.adjPVAL AS qvalue",sep="")
    res1<-cypher(graph,query1)
    
    query2<-paste("MATCH q1=(c1:CELL {name:'",cellC,"'})-[r1]-(c2:cellEx)-[r2]-(s:SYMBOL)-[r3]-(p:PROBE {square:'",squareC,"',edge:",edgeC,"}) WITH p OPTIONAL MATCH (p:PROBE {square:'",squareC,"',edge:",edgeC,"})-[r4]-(p2:PROBE {square:'",squareC,"',edge:",edgeC,"})-[r5]-(s2:SYMBOL) RETURN DISTINCT s2.name AS SYMBOL, p2.name AS probe, p2.logfc AS logfc, p2.adjPVAL AS qvalue",sep="")
    res2<-cypher(graph,query2)
    notInOne<-(!res2$probe%in%res1$probe)
    df1<-data.frame(res1$logfc)
    df2<-data.frame(res2$logfc[notInOne])
    df1$type <- 'cellist'
    df2$type <- 'probecor'
    colnames(df1)<-c("logfc","listType")
    colnames(df2)<-c("logfc","listType")
    #and combine into your new data frame vegLengths
    combo <- rbind(df1,df2)
    #now make your lovely plot
    gpr<-ggplot(combo, aes(logfc, fill = listType)) + geom_density(alpha = 0.2) + ggtitle(paste(cellC,squareC,"edge:",edgeC))
    print(gpr)
    #ggplot(combo, aes(logfc, fill = listType)) + geom_histogram(alpha = 0.5, aes(y = ..density..), position = 'identity')
    
  return(list(res1,res2[notInOne,]))
  
}