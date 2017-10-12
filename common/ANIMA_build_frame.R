#initial text#### 
  #Run checklist
    #fresh session
    #previous neo4j database backed up
    #no debug statements 
    #neoinsert is on

#fixed: ratio order in flat cellprop ratios
#fixed: probetype mapping for question 1 and 4
#fixed: double mapsTo relationships for (probe)-[:mapsTo]-(gene): CREATE replaced with MERGE

##TO DO####
  #URGENT:  
  #fix levs to always have cases first, then controls, in BOTH dimensions (TB and HIV/sex), otherwise cell ratios are weird (factor levels are usually alphabetic, maybe simply use "pericardial fluid" and "whole blood"; this will result in better labels too)
  #database index to speed up writing. Use pathway1 and pathway2 variables!
  #make PBMC subset and make barplot comparisons by PBMC subset with stats (also send to neo4j as PBMC_ratio)
  #implement better ha2 cell-type enrichment analysis (based on nuID rather than gene symbols; code in rcypher)
  #fix cypher query for interprobe netowrks; way too expensive; also fix index if needed
  #put into WGCNA all cell-specific probes regardless of fc and variability (i.e. make a list of probes that always have to be there; cellprobes)
  #!!!!SOME PROBES DO NOT MAP TO WGCNA MODULES (it is because of tv at the wgcna step)

  #sort phenodata csv file correctly so that wgcna plots are more meaningful - DONE
  #fix platform when mapping PROBETYPE-SYMBOL relationships - DONE
  #new phenodata to be used with new version of impi data, where ldh is now in a single column - DONE
  #use revised matrixPD
    #use same phenodata for square1 and square4 - Applied
    #add extra datastructure for immunology data, and do separate correlation analyis (this will add another node type and network)
  #neo4j code: check where to use match and create vs merge, may speed things up - Applied
  #need to allocate more heap space to Cytoscape java rpc, and overall to VM - Applied
  
  #implement fix for cell to cellex mapping as current one is incomplete (check that new csv file works) -- seems to...
  #link pathway genes (starting with reactome) using graphite
    
  #WGCNA do a matplot and ME superimposed - DONE
  #CEGS
  #Map CEGS probes to CEGS
  #CEGS correlation
  #CEGS ME
  #CEGS ME correlation
  #gr
  #Add flow data to test set...this could be 2 networks: flow-cellprop and wgcna-flow - DONE

##TO DO MAJOR##
  #Main loop: add flattened squares as extra edges - DONE
  #Merge in IMPI-MA data - DONE
  #Harmonise phenodata handling between datasets
  #add extra pheno data for certain subsets (e.g. SET1)
  #Implement non-redunant orthogonal annotaion space strategy (e.g. GO, KEGG, InterPRO)
    #Annotation pipelines in R
      #Singular
        #DAVID:DAVIDQuery
      #Gene-set enrichment
        #GSEABase
        #gage
      #Modular
  #sigenrich with phyper

#Setup####

  #Current date string used for versioning output
  source("scripts/ANIMA_setup.R")
  par("xpd"=FALSE)
  #WGCNA multithreading (optional)
  enableWGCNAThreads()
  #Time Zone
  Sys.setenv(TZ="Africa/Johannesburg")
  #Output folders
  ds<-gsub(":","_",gsub(" ","_",Sys.time(),fixed=TRUE),fixed=TRUE)
  #versioned output
  dir.output.version<-file.path(dir.main,paste("build/version_",ds,sep=""))
  rsc="none"
  for (dir in c(dir.main,dir.output.version)) {
    if (!file.exists(dir)) {
      dir.create(dir, recursive=T,showWarnings=T)
    }
  }
  sink(type="output",file=file.path(dir.output.version,"session_output.txt"))
  sessionInfo()
  print(paste("Begin time:",ds))

#switch NEO4J and Cytoscape on or off independently ####
  neoinsert<-"off" #debug
  if (neoinsert=="on"){
  #initialise neo4j db
    graph<-startGraph("http://192.168.65.2:8474/db/data/")#using docker for Mac!
    graph$version
    print("clearing graph")
    clear(graph,input=FALSE)
    print("graph clear, making indexes")
    addIndex(graph,"PROBE","name")
    addIndex(graph,"PROBE","square")
    addIndex(graph,"PROBE","edge")
    addIndex(graph,"SYMBOL","name")
    addIndex(graph,"PROBETYPE","name")
    addIndex(graph,"baylor","name")
    addIndex(graph,"PalWangPW","name")
    addIndex(graph,"ImmunePW","name")
    addIndex(graph,"reactomePW","name")
    print("These are the indexes:")
    getIndex(graph)
  }

#Switch Cytoscape on/off
  cytoscapelink<-"off"
#Reduce pheno
  reducePheno<-FALSE

#Default wgcna_override
  #wgcna_override<-FALSE
      
#import Baylor module data for use later####
  modgenes<-read.csv(file.path(dir.data_root,"modgenes.csv"),header=T)
  mod.names<-as.character(unique(modgenes$moduleID))
  modlist<-list()
  #Gene names for each module
  for (mod in mod.names){
    data<-modgenes[modgenes[,1]==mod,4]
    modlist[[mod]]<-data
  }
  modlist2<-lapply(modlist,function(x){as.character(x)})
  modlist3<-lapply(modlist2,function(x){unique(unlist(strsplit(x,"///")))})
  modlist4<-lapply(modlist3,function (x) gsub("^\\s+|\\s+$", "", x))
  
 #new Baylor modules!
  modgenes2<-read.csv(file.path(dir.data_root,"modgenes2.csv"),header=T)
  mod.names2<-as.character(unique(modgenes2$Module))
  modlist2<-list()
  #Gene names for each module
  for (mod2 in mod.names2){
    data2<-modgenes2[modgenes2[,1]==mod2,8]
    modlist2[[mod2]]<-data2
  }
  modlist2.2<-lapply(modlist2,function(x){as.character(x)})
  modlist3.2<-lapply(modlist2.2,function(x){unique(unlist(strsplit(x,"///")))})
  modlist4.2<-lapply(modlist3.2,function (x) gsub("^\\s+|\\s+$", "", x))

#write modlist4.2 to neo4j for mappings to baylor modules here...
  
  modentrez<-list()
  #Gene names for each module
  for (mod in mod.names){
    data<-modgenes[modgenes[,1]==mod,3]
    modentrez[[mod]]<-data
  }
  modentrez2<-lapply(modentrez,function(x){as.character(x)})
  modentrez3<-lapply(modentrez2,function(x){unique(unlist(strsplit(x,"///")))})
  modentrez4<-lapply(modentrez3,function (x) gsub("^\\s+|\\s+$", "", x))

  modNU<-list()
  #Gene names for each module
  for (mod in mod.names2){
    #data<-modgenes2[modgenes2[,1]==mod,7]
    dataBase<-illuminaHumanv2NUID[mappedkeys(illuminaHumanv2NUID)]
    data<-as.character(unlist(mget(as.character(modgenes2[modgenes2[,1]==mod,7]),dataBase)))
    modNU[[mod]]<-data
  }
  
  modEN<-lapply(lapply(lapply(modNU,nuID2EntrezID,lib.mapping="lumiHumanIDMapping"),as.character),function(x){x[which(x!="")]})

#WGCNA Module preservation container####
  #make empty list and put in RData file for use later
  allmodules<-list()
  save(allmodules,file=file.path(dir.output.version,"allmodules.RData"))
  rm(allmodules)

#Make questions####
source(file.path(dir.data_root,"questions.R")) #this is specific to the project; each project will get its own set of questions tied to specific data.
  
mastervar<-c(ls())
#################################################################################################  
#open outer for loop####DEBUG Off####
#debug
#wgcna_override<-FALSE
#real outer loop 
for (q.i in 1:length(questions)){
  print("before_source_mem_used")
  print(mem_used())
  print(system("free"))
source("~/common/bodybuild.R",local=TRUE)
  print("after_source_mem_used")
  print(mem_used())
  print(system("free"))
}#exit outer loop

#################################################################################################  
  
#11_ADDITIONAL MAPPINGS######

if(q.i==length(questions)&edge==5){
    
print("now for the final mappings")

if (neoinsert=="on"){
  nodes<-RNeo4j::nodes

#NETWORK 16: Map SYMBOL to PalWang####
data(PWLists)

#Make compatible names
newPWNames<-sapply(as.character(PWLists[,2]),function(x){as.character(strsplit(x,"__")[[1]][1])})

names(newPWNames)<-NULL
newPWNames2<-sapply(newPWNames,function(x){
  splits<-unlist(strsplit(x,"_",fixed = TRUE))
  newst<-paste(splits[3:length(splits)],sep=" ",collapse=" ")
  return(newst)
})
names(newPWNames2)<-NULL

newPWNames3<-sapply(newPWNames,function(x){
  splits<-unlist(strsplit(x,"_",fixed = TRUE))
  newst<-paste(splits[2:length(splits)],sep=" ",collapse=" ")
  return(newst)
})
names(newPWNames3)<-NULL

newPWNames4<-gsub("'","prime",gsub(",",";;",gsub(")","*",gsub("(","*",newPWNames3,fixed=TRUE),fixed=TRUE),fixed=TRUE),fixed=TRUE)

#first get only the PalWang pathways in current DB
cPW<-suppressMessages(cypher(graph,"MATCH (n:PalWangPW) RETURN DISTINCT n.name"))
cPW$n.name<-gsub("'","prime",cPW$n.name)

#get the dataframe to use for relationships and nodes
graphdata.raw<-data.frame(cbind(PWLists[,1],newPWNames4))
graphdata.sub<-graphdata.raw[graphdata.raw$newPWNames4%in%cPW$n.name,]
graphdata<-graphdata.sub
colnames(graphdata)<-c("SYMBOL","PalWangPW")

##new
write.csv(graphdata,file.path(dir.figures,"tmp.csv"))

query = paste("LOAD CSV WITH HEADERS FROM 'file://",file.path(dir.figures,"tmp.csv"),"' AS csvLine 
              MATCH (gene:SYMBOL {name:csvLine.SYMBOL}) MATCH (palwang:PalWangPW {name:csvLine.PalWangPW}) MERGE (gene)-[:mapsTo]->(palwang)",sep="")
cypher(graph,query)

##end new


#initial transaction
# tx = suppressMessages(newTransaction(graph))
# 
# for (therow in 1:nrow(graphdata)) {
#   genename = as.character(graphdata[therow, ]$SYMBOL)
#   #print(genename)
#   palwangname = as.character(graphdata[therow, ]$PalWangPW)
#   #print(palwangname)
#   query = paste("MATCH (gene:SYMBOL {name:'",genename,"'}) MATCH (palwang:PalWangPW {name:'",palwangname,"'}) MERGE (gene)-[:mapsTo]->(palwang)",sep="")
#   #print(query)
#   # Upload in blocks of 1000.
#   if(therow %% 1000 == 0) {
#     # Commit current transaction.
#     suppressMessages(commit(tx))
#     print(paste("Batch", therow / 1000, "committed."))
#     # Open new transaction.
#     tx = suppressMessages(newTransaction(graph))
#   }
#   
#   suppressMessages(appendCypher(tx,query)) 
#   
# }
# 
# suppressMessages(commit(tx))

#NETWORK 17: Map SYMBOL to reactomePW####
#reactome
library(ReactomePA)
library(lumi)
library(reactome.db)
library(lumiHumanIDMapping)

nuIDs <- featureNames(datalist[[1]])
mappingInfo <- nuID2EntrezID(nuIDs, lib.mapping='lumiHumanIDMapping')
head(mappingInfo)
map1<-data.frame(cbind(names(mappingInfo),as.character(mappingInfo)))
colnames(map1)<-c("nuID","ENTREZID")

map2ids<-keys(reactome.db,keytype="PATHNAME")
map2<-select(reactome.db,keys=map2ids,columns=c("ENTREZID"),keytype="PATHNAME")

mappingInfo2 <- nuID2targetID(nuIDs, lib.mapping='lumiHumanIDMapping')
head(mappingInfo2)
map3<-data.frame(cbind(names(mappingInfo2),as.character(mappingInfo2)))
colnames(map3)<-c("nuID","SYMBOL")

#merge 3 maps
merged<-merge(merge(map1,map2,by.x="ENTREZID",by.y="ENTREZID",all.x=T,all.y=T),map3,by.x="nuID",by.y="nuID",all.x=T,all.y=T)
sum(!is.na(merged$PATHNAME))
merged2<-merged[!is.na(merged$PATHNAME),]#this has incompatible characters!

theOld<-merged2$PATHNAME
theNew<-iconv(theOld,from="LATIN-9",to="UTF-8")

merged2$PATHNAME2<-gsub("Homo sapiens: ","",theNew,fixed=TRUE)
merged2$PATHNAME3<-gsub("'","",gsub(",",";;",gsub("~","approx.",gsub(")","*",gsub("(","*",merged2$PATHNAME2,fixed=TRUE),fixed=TRUE),fixed=TRUE),fixed=TRUE),fixed=TRUE)

graphdata.raw<-merged2[,c("SYMBOL","PATHNAME3")]

indb<-suppressMessages(cypher(graph,"MATCH (p:reactomePW) RETURN DISTINCT p.name"))

graphdata<-graphdata.raw[graphdata.raw$PATHNAME3%in%indb$p.name,]
colnames(graphdata)<-c("SYMBOL","reactomePW")

##new
write.csv(graphdata,file.path(dir.figures,"tmp.csv"))

query = paste("LOAD CSV WITH HEADERS FROM 'file://",file.path(dir.figures,"tmp.csv"),"' AS csvLine 
              MATCH (gene:SYMBOL {name:csvLine.SYMBOL}) MATCH (reactome:reactomePW {name:csvLine.reactomePW}) MERGE (gene)-[:mapsTo]->(reactome)",sep="")
cypher(graph,query)

##end new



#initial transaction
# tx = suppressMessages(newTransaction(graph))
# 
# for (therow in 1:nrow(graphdata)) {
#   #print(i)
#   genename = as.character(graphdata[therow, ]$SYMBOL)
#   #print(genename)
#   reactomename = as.character(graphdata[therow, ]$reactomePW)
#   #print(palwangname)
#   query = paste("MATCH (gene:SYMBOL {name:'",genename,"'}) MATCH (reactome:reactomePW {name:'",reactomename,"'}) MERGE (gene)-[:mapsTo]->(reactome)",sep="")
#   #print(query)
#   # Upload in blocks of 1000.
#   if(therow %% 1000 == 0) {
#     # Commit current transaction.
#     suppressMessages(commit(tx))
#     print(paste("Batch", therow / 1000,"of",nrow(graphdata), "committed."))
#     # Open new transaction.
#     tx = suppressMessages(newTransaction(graph))
#   }
#   
#   suppressMessages(appendCypher(tx,query)) 
#   
# }
# 
# suppressMessages(commit(tx))

#NETWORK 18: Map SYMBOL to ImmunePW####

data(ImmunePathwayLists)
head(ImmunePathwayLists)

#Make compatible names
newIpwNames<-gsub(",",";;",gsub(")","*",gsub("(","*",ImmunePathwayLists[,2],fixed=TRUE),fixed=TRUE),fixed=TRUE)

#first get only the PalWang pathways in current DB
cIPW<-suppressMessages(cypher(graph,"MATCH (n:ImmunePW) RETURN DISTINCT n.name"))
cIPW$n.name<-gsub("'","prime",cIPW$n.name)

#get the dataframe to use for relationships and nodes
graphdata.raw<-data.frame(cbind(ImmunePathwayLists[,1],newIpwNames))
graphdata.sub<-graphdata.raw[graphdata.raw$newIpwNames%in%cIPW$n.name,]
graphdata<-graphdata.sub
colnames(graphdata)<-c("SYMBOL","ImmunePW")

##new
write.csv(graphdata,file.path(dir.figures,"tmp.csv"))

query = paste("LOAD CSV WITH HEADERS FROM 'file://",file.path(dir.figures,"tmp.csv"),"' AS csvLine 
              MATCH (gene:SYMBOL {name:csvLine.SYMBOL}) MATCH (immune:ImmunePW {name:csvLine.ImmunePW}) MERGE (gene)-[:mapsTo]->(immune)",sep="")
cypher(graph,query)

##end new


#initial transaction
# tx = suppressMessages(newTransaction(graph))
# 
# for (therow in 1:nrow(graphdata)) {
#   #print(i)
#   genename = as.character(graphdata[therow, ]$SYMBOL)
#   #print(genename)
#   immunePWname = as.character(graphdata[therow, ]$ImmunePW)
#   #print(palwangname)
#   query = paste("MATCH (gene:SYMBOL {name:'",genename,"'}) MATCH (immune:ImmunePW {name:'",immunePWname,"'}) MERGE (gene)-[:mapsTo]->(immune)",sep="")
#   #print(query)
#   # Upload in blocks of 1000.
#   if(therow %% 1000 == 0) {
#     # Commit current transaction.
#     suppressMessages(commit(tx))
#     print(paste("Batch", therow / 1000, "committed."))
#     # Open new transaction.
#     tx = suppressMessages(newTransaction(graph))
#   }
#   
#   suppressMessages(appendCypher(tx,query)) 
#   
# }
# 
# suppressMessages(commit(tx))

#END NEO4J


#NETWORK 20: Map SYMBOL to cellEx####

data(BloodLists)
head(BloodLists)

#Make compatible datasets for all three methods (BloodLists, HaemAtlas and Abbas)

#fromBloodLists
newCellNames1<-paste(unlist(lapply(strsplit(unlist(lapply(strsplit(BloodLists[,2],"__"),function(x){x[[1]]})),"_"),function(x){x[[1]]})),"_xp_1",sep="")
cell.sub1<-cbind(BloodLists[,1],newCellNames1)
#From HaemAtlas
newCellNames2<-paste(ha.2[,2],"_xp_2",sep="")
cell.sub2<-cbind(ha.2[,1],newCellNames2)
#From Abbas
newCellNames3<-paste(abbasvals[,2],"_xp_3",sep="")
cell.sub3<-cbind(abbasvals[,1],newCellNames3)

threecell<-rbind(cell.sub1,cell.sub2,cell.sub3)

#first get only the celltypes in current DB
cCell<-suppressMessages(cypher(graph,"MATCH (n:cellEx) RETURN DISTINCT n.name"))

#get the dataframe to use for relationships and nodes
graphdata.raw<-data.frame(threecell)
graphdata.sub<-graphdata.raw[graphdata.raw$newCellNames%in%cCell$n.name,]
graphdata<-graphdata.sub
colnames(graphdata)<-c("SYMBOL","cellEx")

##new
write.csv(graphdata,file.path(dir.figures,"tmp.csv"))

query = paste("LOAD CSV WITH HEADERS FROM 'file://",file.path(dir.figures,"tmp.csv"),"' AS csvLine 
              MATCH (gene:SYMBOL {name:csvLine.SYMBOL}) MATCH (cell:cellEx {name:csvLine.cellEx}) MERGE (gene)-[:mapsTo]->(cell)",sep="")
cypher(graph,query)

##end new


#initial transaction
# tx = suppressMessages(newTransaction(graph))
# 
# for (therow in 1:nrow(graphdata)) {
#   #print(i)
#   genename = as.character(graphdata[therow, ]$SYMBOL)
#   #print(genename)
#   cellname = as.character(graphdata[therow, ]$cellEx)
#   #print(palwangname)
#   ##proposed fix: replace match with merge if i want to include cellex that is nowhere significantly enriched in any wgcna module
#   query = paste("MERGE (gene:SYMBOL {name:'",genename,"'}) MERGE (cell:cellEx {name:'",cellname,"'}) MERGE (gene)-[:mapsTo]->(cell)",sep="")
#   #print(query)
#   # Upload in blocks of 1000.
#   if(therow %% 1000 == 0) {
#     # Commit current transaction.
#     suppressMessages(commit(tx))
#     print(paste("Batch", therow / 1000, "committed."))
#     # Open new transaction.
#     tx = suppressMessages(newTransaction(graph))
#   }
#   
#   suppressMessages(appendCypher(tx,query)) 
#   
# }
# 
# suppressMessages(commit(tx))

#END NEO4J

#NETWORK 21 and NETWORK 22: Map CELL to cellEx and cellprop #NEED TO COMPLETE THIS####

#get the dataframe to use for relationships and nodes
graphdata<-read.csv(file.path(dir.annot,"cellmapEX.csv"),header=F)
colnames(graphdata)<-c("cellEx","cell")

##new
write.csv(graphdata,file.path(dir.figures,"tmp.csv"))

query = paste("LOAD CSV WITH HEADERS FROM 'file://",file.path(dir.figures,"tmp.csv"),"' AS csvLine 
              MATCH (c:cellEx {name:csvLine.cellEx}) MERGE (cell:CELL {name:csvLine.cell}) CREATE (c)-[:mapsTo]->(cell)",sep="")
cypher(graph,query)

##end new
#get the dataframe to use for relationships and nodes
graphdata<-read.csv(file.path(dir.annot,"cellmapPROP.csv"),header=F)
colnames(graphdata)<-c("cellprop","cell")

##new
write.csv(graphdata,file.path(dir.figures,"tmp.csv"))

query = paste("LOAD CSV WITH HEADERS FROM 'file://",file.path(dir.figures,"tmp.csv"),"' AS csvLine 
              MATCH (c2:cellprop {name:csvLine.cellprop}) MERGE (cell:CELL {name:csvLine.cell}) CREATE (c2)-[:mapsTo]->(cell)",sep="")
cypher(graph,query)

##end new


#initial transaction
# tx = suppressMessages(newTransaction(graph))
# 
# for (therow in 1:nrow(graphdata)) {
#   #print(i)
#   name1 = as.character(graphdata[therow, ]$exprop)
#   name2 = as.character(graphdata[therow, ]$cell)
#   query = paste("MERGE (exprop {name:'",name1,"'}) MERGE (cell:CELL {name:'",name2,"'}) MERGE (exprop)-[:mapsTo]->(cell)",sep="")
#   #print(query)
#   # Upload in blocks of 1000.
#   if(therow %% 1000 == 0) {
#     # Commit current transaction.
#     suppressMessages(commit(tx))
#     print(paste("Batch", therow / 1000, "committed."))
#     # Open new transaction.
#     tx = suppressMessages(newTransaction(graph))
#   }
#   
#   suppressMessages(appendCypher(tx,query)) 
#   
# }
# 
# suppressMessages(commit(tx))

#END NEO4J

source("~/scripts/ANIMA_MC.R")

#end final mappings####
try(sendmail("armin.deffur@icloud.com", subject="Final mappings",message="Final mappings for reactome, PalWang, ImmunePW and cellEx are done", password="rmail"))

}#end neoinsert if
  }else{try(sendmail("armin.deffur@icloud.com", subject="Oops!",message=print(paste("I think something went wrong at",Sys.time(),"at question",q.i,"and analysis",analysis,"\n","error:\n",geterrmessage())), password="rmail"))}

sink()