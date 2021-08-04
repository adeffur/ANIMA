
install.packages("treemap")
library(treemap)

#wgcna
query<-"MATCH (n:wgcna {square:'blood.PCF.defPC',edge:'5'})-[r]-(p:PROBE) RETURN n.name AS colour, COUNT(r) as size"
query<-"MATCH (n:wgcna {square:'blood.PCF.probPC',edge:'5'})-[r]-(p:PROBE) RETURN n.name AS colour, COUNT(r) as size"
query<-"MATCH (n:wgcna {square:'prob.def.blood',edge:'5'})-[r]-(p:PROBE) RETURN n.name AS colour, COUNT(r) as size"
query<-"MATCH (n:wgcna {square:'prob.def.fluid',edge:'5'})-[r]-(p:PROBE) RETURN n.name AS colour, COUNT(r) as size"


res<-cypher(graph,query)

treemap(res,index="colour",vColor="colour",vSize="size",type="color")



#baylormod
query<-"MATCH (b:baylor {square:'prob.def.fluid',edge:'5'})-[r]-(p:PROBE) RETURN b.name AS name, COUNT(r) as size, toFloat(b.diffEX) as value"
res<-cypher(graph,query)
treemap(res,index="name",vSize="size",vColor="value",type="value",palette="-RdYlBu")
head(res)


#Sankey
#query builder
square<-"blood.PCF.defPC"
query1=paste("MATCH (n:wgcna {square:'",square,"',edge:'5'})-[r0]-(p:PROBE)-[r1]-(b:baylor {square:'",square,"',edge:'5'}) ",sep="")

query2<-" RETURN n.name AS node1, COUNT(r1) AS ngenes, b.name AS node2"

query<-paste(query1,query2,sep="")
print(query)
res<-cypher(graph,query)
print(head(res))

links<-res
colnames(links)<-c("source","value","target")
nodes<-unique(c(links$source,links$target))
nodes<-data.frame(nodes)
colnames(nodes)<-"name"
nodes<-cbind(seq(0,nrow(nodes)-1),nodes)
colnames(nodes)<-c("index","name")

for(index in 1:nrow(nodes)){
  links[links==nodes[index,2]]<-nodes[index,1]
}
links<-data.frame(apply(links,2,as.numeric))
class(links)
print("start plot")
print("print(head(links))")
print(head(links))
links2<-cbind(links,res)
print(head(links2))

#links2<-links2[links2$value>10,]
#nodes<-nodes[which(nodes$name%in%links2$node1|nodes$name%in%links2$node2),]

#Scolors <- paste(links2$node1, collapse = '", "')
#colorJS <- paste('d3.scaleOrdinal(["', Scolors, '"])')
#colourScale = colorJS, 

sankeyNetwork(Links = links2, Nodes = nodes, Source = 'source',Value='value',Target = 'target',NodeID = 'name',NodeGroup="name",LinkGroup="node1", iterations = 128)
