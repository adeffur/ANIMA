#While the igraph package is very flexible, we require some additional functionality; specifically the ability encode information in the border colour of network nodes. As the default border width is fixed, we here define a function that allows us to flexibly specify the width and colour of the node borders.
#Code attribution: Gabor Csardi at \url{http://lists.gnu.org/archive/html/igraph-help/2013-03/msg00030.html}

mycircle <- function(coords, v=NULL, params) {
  vertex.color <- params("vertex", "color")
  if (length(vertex.color) != 1 && !is.null(v)) {
    vertex.color <- vertex.color[v]
  }
  vertex.size  <- 1/200 * params("vertex", "size")
  if (length(vertex.size) != 1 && !is.null(v)) {
    vertex.size <- vertex.size[v]
  }
  vertex.frame.color <- params("vertex", "frame.color")
  if (length(vertex.frame.color) != 1 && !is.null(v)) {
    vertex.frame.color <- vertex.frame.color[v]
  }
  vertex.frame.width <- params("vertex", "frame.width")
  if (length(vertex.frame.width) != 1 && !is.null(v)) {
    vertex.frame.width <- vertex.frame.width[v]
  }
  
  mapply(coords[,1], coords[,2], vertex.color, vertex.frame.color,
         vertex.size, vertex.frame.width,
         FUN=function(x, y, bg, fg, size, lwd) {
           symbols(x=x, y=y, bg=bg, fg=fg, lwd=lwd,
                   circles=size, add=TRUE, inches=FALSE)
         })
}

add.vertex.shape("fcircle", clip=igraph.shape.noclip,plot=mycircle, parameters=list(vertex.frame.color=1,vertex.frame.width=1))

###################################
#igraph_plotter function definition
###################################

igraph_plotter<-function(query.base,nodelist,edgetrips,rimpar,plot=TRUE,csv=FALSE,prefix=cytodir,filename="",vertexsize=15,lay_out=layout_with_dh(ig,maxiter=200,cool.fact=0.95,weight.node.dist=.65),optlabel="",optvalue="",optchar="",make.vertex.label=TRUE,return_graph=FALSE,plot_bipartite=FALSE,vertex.label.cex=.55,legendcex=1,plotd3=FALSE,plotwhich=c(1,2,3),main=NULL){
 
  #Step 1: return the nodes
  allnodes<-data.frame(node_id=character(),nodename=character(),kind=character(),square=character(),edge=numeric(),logFC=numeric(),diffME=numeric())
  for (nodetype in nodelist){
    query.nodes<-paste("RETURN DISTINCT ID(",nodetype,") AS node_id, ",nodetype,".name AS nodename, labels(",nodetype,") AS kind, ",nodetype,".square AS square, ",nodetype,".edge AS edge, ",nodetype,".ratio AS ratio, ",nodetype,".diffP AS diffP, ",nodetype,".diffQ AS diffQ, ",nodetype,".logfc AS logfc, ",nodetype,".adjPVAL AS adjPVAL, ",nodetype,".sigenrich AS sigenrich, ",nodetype,".modAUC1 AS modAUC1, ",nodetype,".modAUC2 AS modAUC2, ",nodetype,".diffME AS diffME, ",nodetype,".diffEX AS diffEX",sep="")
    query<-paste(query.base,query.nodes,sep="")
    #print(query)
    r<-cypher(graph,query)
    allnodes<-rbind(allnodes,r)
  }
  allnodes$node_id2<-paste("id",allnodes$node_id,sep="")
  allnodes[is.na(allnodes)]<-""
  allnodes<-unique(allnodes)
  allnodes<-allnodes[which(allnodes$node_id!=""),]
  allnodes2<-allnodes
  if ("PROBE" %in% allnodes2$kind) {
    allnodes2[which(allnodes2$kind=="PROBE"),]$nodename<-""
  }
  
  #Step 2: return the edges
  alledges<-data.frame(source=character(),target=character(),qvalue=numeric())
  for (edgetrip in edgetrips){
    query.edges<-paste("RETURN DISTINCT ID(",edgetrip[1],") AS source, ID(",edgetrip[3],") AS target, ",edgetrip[2],".qvalue AS qvalue, ","TYPE(",edgetrip[2],") AS reltype, ",edgetrip[2],".Rsq AS RSQ, ",edgetrip[2],".weight AS PearsonR, ",edgetrip[2],".TOMweight AS TOMweight",sep="")
    query<-paste(query.base,query.edges,sep="")
    r<-cypher(graph,query)
    alledges<-rbind(alledges,r)
  }
  alledges$source2<-paste("id",alledges$source,sep="")
  alledges$target2<-paste("id",alledges$target,sep="")
  alledges[is.na(alledges)]<-""
  alledges<-alledges[which(alledges$source2!="idNA"|alledges$target2!="idNA"),]
  alledges<-unique(alledges)
  
  #make sure all nodes in edgelist are also in nodelist
  print(nrow(alledges))
  alledges<-alledges[which(alledges$source%in%allnodes2$node_id&alledges$target%in%allnodes2$node_id),]
  print(nrow(alledges))
  
  #Step 3: optionally export the nodes and edges in a format that can easily be imported into Cytoscape
  if (csv==TRUE){
    write.csv(allnodes,file.path(prefix,paste(filename,"_nodes.csv",sep="")),quote=F)
    write.csv(alledges,file.path(prefix,paste(filename,"_edges.csv",sep="")),quote=F)
  }
  
  #Step 3.5 Check if code needs to continue
  #if(plot==FALSE){stop("No plot generated! Look for csv files!")}
  
  #Step 4. Create the igraph object
  ig<-graph.data.frame(alledges,directed = FALSE,vertices=allnodes2)
  kindcol<-factor(V(ig)$kind)
  
  #Step 5. Assign standard colours to the various node types (colours match the figures in the main manuscript)
  ty<-c("baylor","cellEx","cellprop","ImmunePW","PalWangPW","pheno","PROBE","PROBETYPE","reactomePW","SYMBOL","wgcna","CELL","flowcellprop","cellpropFprop","bigc")
  tycol<-c("darkgoldenrod1","chartreuse1","chartreuse3","cadetblue1","cadetblue2","plum1","azure","gray88","cadetblue2","yellow","cornsilk2","seagreen2","springgreen2","springgreen3","pink")
  dicti<-cbind(ty,tycol)
  
  start<-V(ig)$kind
  legstart<-levels(as.factor(V(ig)$kind))
  for (i in 1:nrow(dicti)){
    start[which(start==as.character(dicti[i,1]))]<-dicti[i,2]
    legstart[which(legstart==as.character(dicti[i,1]))]<-dicti[i,2]
  }
  kindcolvector<-start
  
  #Step 6. Calculate rimcolours based on a specific node property. Only nodes with that property will have their rims coloured this way.
  rimvals<-eval(parse(text=paste("as.numeric(V(ig)$",rimpar,")",sep="")))
  rimvals[is.na(rimvals)]<-0
  extreme<-max(abs(rimvals))
  rimcol<-as.character(numbers2colors(rimvals,signed=TRUE,lim=c(-extreme,extreme),centered=T))
  ig<-graph.data.frame(alledges,vertices=allnodes2)
  framewidth<-rep(1,length(rimvals))
  framewidth[which(rimvals!=0)]<-4
  
  #Step 7. Calculate the node positions according to a standard graph layout algorithm. The default can be overriden if lay_out is passed as an argument to the function
  ig.lay<-lay_out
  
  #Step 8. Generate node (vertex) labels. The default uses the node name, but additional information can be added
  if(make.vertex.label==TRUE){
    vertexlabel<-V(ig)$nodename
    addvalue<-sprintf("%.4s",eval(parse(text=optvalue))) #this hast to be based on V(ig)$something
    addvalue2<-eval(parse(text=optchar))
    addlabel<-paste(optlabel,addvalue2[which(addvalue!="")],addvalue[which(addvalue!="")])
    vertexlabel[which(addvalue!="")]<-paste(vertexlabel[which(addvalue!="")],addlabel)
  }else{vertexlabel=""}
  #Step 9. Generate types for bipartite graphs
  
  if (length(unique(V(ig)$kind))==2){
    ig.b<-bipartite.mapping(ig)
    V(ig)$type<-ig.b$type
    #     nodetypes<-as.factor(V(ig)$kind)
    #     levels(nodetypes)<-c(0,1)
    #     V(ig)$types<-as.numeric(nodetypes)
  }
  
  #Step 10. Generate the plot and add two legends (one identifying the node types, the other showing a scale for the rim colour)
  if (plotd3 != TRUE&1%in%plotwhich){
  #plot(ig,vertex.size=vertexsize,vertex.shape="fcircle",vertex.label=vertexlabel,vertex.color=as.character(kindcolvector),layout=ig.lay,vertex.label.cex=vertex.label.cex,vertex.label.color="dimgray",vertex.frame.color=rimcol,vertex.frame.width=framewidth,edge.label=sprintf("%.5s",E(ig)$PearsonR),edge.label.cex=.7,xlim = c(-1.0, 1.0),ylim = c(-1.7, 1.0),margin=c(0,0,0,0),edge.arrow.mode=0,main=main)
  plot(ig,vertex.size=vertexsize,vertex.shape="fcircle",vertex.label=vertexlabel,vertex.color=as.character(kindcolvector),layout=ig.lay,vertex.label.cex=vertex.label.cex,vertex.label.color="dimgray",vertex.frame.color=rimcol,vertex.frame.width=framewidth,edge.label=sprintf("%.5s",E(ig)$PearsonR),edge.label.cex=.7,edge.arrow.mode=0,main=main)#,xlim = c(-1.0, 1.0),ylim = c(-1.7, 1.0),margin=c(0,0,0,0)
    
  legend("bottomleft",levels(factor(V(ig)$kind)),fill=legstart,cex=legendcex,title="Vertex Type")
  legend("bottomright",sprintf("%.3f",seq(-extreme,extreme,length.out=11)),fill=as.character(numbers2colors(seq(-extreme,extreme,length.out=11),signed=TRUE,lim=c(-extreme,extreme),centered=T)),cex=legendcex,title=rimpar,ncol=2)
  }
  #Step 11. Optionally also plot two bipartite projections
  #make projections
  
  if(plot_bipartite==TRUE){
    bp<-bipartite.projection(ig)
    b1<-bp[[1]]
    b2<-bp[[2]]
    
    if (csv==TRUE){
      
      bp1e<-as_data_frame(bp[[1]],"edges")
      bp2e<-as_data_frame(bp[[2]],"edges")
      bp1v<-as_data_frame(bp[[1]],"vertices")
      bp2v<-as_data_frame(bp[[2]],"vertices")
      
      write.csv(bp1e,file.path(prefix,paste(filename,"_bp1_edges.csv",sep="")),quote=F)
      write.csv(bp2e,file.path(prefix,paste(filename,"_bp2_edges.csv",sep="")),quote=F)
      write.csv(bp1v,file.path(prefix,paste(filename,"_bp1_vertices.csv",sep="")),quote=F)
      write.csv(bp2v,file.path(prefix,paste(filename,"_bp2_vertices.csv",sep="")),quote=F)
      
       }
    
    #image 2: projection1
    b1.lay<-layout.kamada.kawai(b1,dim=2)
    start<-V(b1)$kind
    legstart<-levels(as.factor(V(b1)$kind))
    for (i in 1:nrow(dicti)){
      start[which(start==as.character(dicti[i,1]))]<-dicti[i,2]
      legstart[which(legstart==as.character(dicti[i,1]))]<-dicti[i,2]
    }
    kindcolvector<-start
    if(2%in%plotwhich){
    plot(b1,vertex.label=V(b1)$nodename,vertex.color=as.character(kindcolvector),layout=b1.lay,vertex.label.cex=vertex.label.cex)
    }
    
    #image 3: projection2
    #b2.lay<-layout.fruchterman.reingold(b2,dim=2)
    b2.lay<-layout.kamada.kawai(b2,dim=2)
    start<-V(b2)$kind
    legstart<-levels(as.factor(V(b2)$kind))
    for (i in 1:nrow(dicti)){
      start[which(start==as.character(dicti[i,1]))]<-dicti[i,2]
      legstart[which(legstart==as.character(dicti[i,1]))]<-dicti[i,2]
    }
    kindcolvector<-start
    if(3%in%plotwhich){
    plot(b2,vertex.label=V(b2)$nodename,vertex.color=as.character(kindcolvector),layout=b2.lay,vertex.label.cex=vertex.label.cex)
    }
  }
  #Step 12. Return a graph object
  
  if(return_graph==TRUE){
    return(ig)
  }
  
  #Step 13. D3 output
  if (plotd3 == TRUE){
    # Convert to object suitable for networkD3
    ig_d3 <- igraph_to_networkD3(ig, group = V(ig)$kind)
    
    # Create force directed network plot
    forceNetwork(Links = ig_d3$links, Nodes = ig_d3$nodes, 
                 Source = 'source', Target = 'target', 
                 NodeID = 'name', Group = 'group')
  }
  
  }