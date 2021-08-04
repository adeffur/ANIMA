# 0. Initialize
library(dplyr)
source("~/scripts/ANIMA_setup.R")
source(file.path(dir.data_root,"questions.R")) #this is specific to the project; each project will get its own set of questions tied to specific data.
source("~/source_data/setlist.R")
graph<-startGraph(graphstring)#using docker for Mac!


######
#query<-"MATCH (p:personCell {square:'blood.PCF.defPC'}) WHERE p.personName =~ 'PC124.*' RETURN p.personName as person,p.name AS cell, p.value AS value"
query<-"MATCH (p:personCell {square:'blood.PCF.defPC'}) RETURN p.personName as person,p.name AS cell, p.value AS value"


res<-cypher(graph,query)

persons<-unique(res$person)
cells<-unique(res$cell)
df<-data.frame()

for (person in persons) {
  res2<-res[which(res$person==person),]
  for (cell in cells){
   df[person,cell]<-res2[which(res2$cell==cell),"value"]
  }
  
}

##########

setwd("~/output/build/")
folders<-list.files()
folders<-folders[grep("Q.*",folders)]
for (folder in folders){
  setwd(paste("~/output/build/",folder,sep=""))
  setwd("5_WGCNA_4k_tv/results")
  data<-read.csv("ModuleEigengenes_all_samples.csv")
  square<-strsplit(folder,"_")[[1]][3]
  print(square)
  print(nrow(data))
  
  graphdata<-data
  rownames(graphdata)<-data[,1]
  graphdata<-graphdata[,2:ncol(graphdata)]
  MEcol<-unlist(strsplit(colnames(graphdata),"ME"))[seq(2,ncol(graphdata)*2,2)]
  colnames(graphdata)<-MEcol
  MElist<-MEcol
  t = suppressMessages(newTransaction(graph))
  
  for (therow in 1:nrow(graphdata)) {
    personName = as.character(rownames(graphdata)[therow])
    
    for (theME in MElist){
      persME = graphdata[therow,theME]
      if(is.na(persME)){persME <- "not done"}
      #class1g = graphdata[therow,contrast.variable1]
      #class2g = graphdata[therow,contrast.variable2]
      #query = "MERGE (person:person {name:{personName},square:{squareg},class1:{class1g},class2:{class2g}}) MERGE (personCell:personCell {name:{cellname},personName:{personName},value:{persCell},square:{squareg},class1:{class1g},class2:{class2g}}) MERGE (person)-[:has]->(personCell)"
      query = "MERGE (person:person {name:{personName},square:{squareg}}) MERGE (personME:personME {name:{MEname},personName:{personName},value:{persME},square:{squareg}}) MERGE (person)-[:has]->(personME)"
      
      #if(!is.na(persPheno)){
      suppressMessages(appendCypher(t, 
                                    query,
                                    personName = personName,
                                    persME = persME,
                                    MEname = theME,
                                    squareg = square
                                    
                                    
      ))
      #}#end if
      
    }
    
    
    
  }#end row
  
  suppressMessages(commit(t))
  
  
  
}
