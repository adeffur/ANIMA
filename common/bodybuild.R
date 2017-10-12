#protected variables
mastervar2<-c(mastervar,"q.i","mastervar")

print("start_loop_mem_used")
print(mem_used())
print(system("free"))
#Definitions for the current question####
#Question (c.q: "current question")
c.q=questions[[q.i]]
datalist.filenames=c.q$datalist_filenames
dataset.variables=c.q$dataset_variables
dataset.names=c.q$dataset_names
contrast.variable=c.q$contrast_variable
contrast.variable1=c.q$contrast_variable[[1]]
contrast.variable2=c.q$contrast_variable[[2]]
colour.variable1=c.q$colour_variable[[1]]
colour.variable2=c.q$colour_variable[[2]]
colourmap1=c.q$colourmap1
colourmap2=c.q$colourmap2
pqlist=c.q$pqlist
validation=c.q$validation
pathway1<-c.q$pathway1
pathway2<-c.q$pathway2
pathway1COL<-c.q$pathway1COL
pathway2COL<-c.q$pathway2COL
pqlist<-c.q$pqlist
matrixPDname<-c.q$matrixPDname
filtervar=c.q$filter
#New additions
#wgcna_override=c.q$wgcna_override
wgcnaSTP=c.q$wgcnaSTP
wgcnaDS=c.q$wgcnaDS

print(paste("Results for Question ",q.i,": ",names(questions)[q.i],sep=""))

#Define folders for storing output####
# Define folder for saving results
dir.results <- file.path(dir.output.version,paste("Q_",q.i,"_",dataset.variables,sep=""), "results")
# Define folder for saving figures
dir.figures <- file.path(dir.output.version,paste("Q_",q.i,"_",dataset.variables,sep=""), "figures")
# create folders
for (dir in c(dir.main, dir.figures, dir.results)) {
  if (!file.exists(dir)) {
    dir.create(dir, recursive=T,showWarnings=T)
  }
}

#load data and make datalist####

for (filename in datalist.filenames){
  load(file.path(dir.rdata,filename))
}

datalist<-list()
for (dataindex in 1:length(datalist.filenames)){
  datalist[[dataindex]]<-eval(parse(text=dataset.variables[[dataindex]]))
  names(datalist)[[dataindex]]<-dataset.names[[dataindex]]
}

#properly sort samples by classifiers
classif.1<-eval(parse(text=paste("datalist[[1]]$",contrast.variable1,sep="")))
classif.2<-eval(parse(text=paste("datalist[[1]]$",contrast.variable2,sep="")))

datalist[[1]]<-datalist[[1]][,order(classif.2,classif.1)]

print(paste(c("Data to be analysed:",names(datalist))))

setwd(dir.figures)

#open inner for loop####DBUG ON####
i<-1#DEBUG

#for (i in 1:1){

#real
#for (i in 1:length(datalist)){

print(paste("current dataset: ",names(datalist)[[i]]))
figure<-1

#QC#####
setwd(dir.figures)
#Outlier plot
pdf(paste('fig',q.i,'.',figure,'_',dataset.variables[[i]],'_outlierPlot.pdf',sep=""))
plot(datalist[[i]],what="outlier",main=names(datalist)[i])
dev.off()
figure<-figure+1

#Definitions for current dataset####
classifier1<-eval(parse(text=paste("datalist[[i]]$",contrast.variable1,sep="")))
classifier2<-eval(parse(text=paste("datalist[[i]]$",contrast.variable2,sep="")))
#Make class levels and colours (used to be called tbstat)
classfactor1=list()
classfactor2=list()
classfactor1[[i]]<-as.factor(eval(parse(text=paste("datalist[[i]]$",colour.variable1,sep=""))))
classfactor2[[i]]<-as.factor(eval(parse(text=paste("datalist[[i]]$",colour.variable2,sep=""))))

classcolours1=list()
classcolours2=list()
classcolours1[[i]]<-classfactor1[[i]]
classcolours2[[i]]<-classfactor2[[i]]
levels(classcolours1[[i]])<-colourmap1
levels(classcolours2[[i]])<-colourmap2

pathways=list(
  "edge1"=c.q$pathway1,
  "edge2"=c.q$pathway1,
  "edge3"=c.q$pathway2,
  "edge4"=c.q$pathway2,
  "edge5"=c.q$pathway1
)

pathwaysCOL=list(
  "edge1"=c.q$pathway1COL,
  "edge2"=c.q$pathway1COL,
  "edge3"=c.q$pathway2COL,
  "edge4"=c.q$pathway2COL,
  "edge5"=c.q$pathway1COL
)

edges=list(
  "edge1"=which(classifier2==c.q$pathway2[[1]]),#women
  "edge2"=which(classifier2==c.q$pathway2[[2]]),#men
  "edge3"=which(classifier1==c.q$pathway1[[1]]),#control
  "edge4"=which(classifier1==c.q$pathway1[[2]]),#case
  "edge5"=which(classifier1==c.q$pathway1[[1]]|classifier1==c.q$pathway1[[2]])#all 
)

classifiers=list(
  "edge1"=classifier1[edges[[1]]],
  "edge2"=classifier1[edges[[2]]],
  "edge3"=classifier2[edges[[3]]],
  "edge4"=classifier2[edges[[4]]],
  "edge5"=classifier1[edges[[5]]]
)

classes=list(
  "edge1"=classfactor1[[i]][edges[[1]]],
  "edge2"=classfactor1[[i]][edges[[2]]],
  "edge3"=classfactor2[[i]][edges[[3]]],
  "edge4"=classfactor2[[i]][edges[[4]]],
  "edge5"=classfactor1[[i]][edges[[5]]]
)

colourfactors=list(
  "edge1"=classcolours1[[i]][edges[[1]]],
  "edge2"=classcolours1[[i]][edges[[2]]],
  "edge3"=classcolours2[[i]][edges[[3]]],
  "edge4"=classcolours2[[i]][edges[[4]]],
  "edge5"=classcolours1[[i]][edges[[5]]]
)

colourvars=list(
  "edge1"=colour.variable1,
  "edge2"=colour.variable1,
  "edge3"=colour.variable2,
  "edge4"=colour.variable2,
  "edge5"=colour.variable1
)

contrastvars=list(
  "edge1"=contrast.variable1,
  "edge2"=contrast.variable1,
  "edge3"=contrast.variable2,
  "edge4"=contrast.variable2,
  "edge5"=contrast.variable1
)

#NETWORK 0A: Send probe-gene mappings to Neo4j once per platform####
allfeatures<-rownames(fData(datalist[[i]]))
annotSYM_all<-as.character(getSYMBOL(allfeatures,"lumiHumanAll.db"))#need new annotSYM that matches targetID annotation. targetID is full of aliases!!!
#first name: annotSYM_all (these map to HGNC names)
#second choice:if annotSYM_all entry is NA, then use "targetID" from LumiBatch object (better than nothing)
if(ncol(fData(datalist[[i]]))>1){
  probemap<-data.frame(cbind(allfeatures,fData(datalist[[i]])[,2],annotSYM_all),stringsAsFactors=F)
  replace<-as.character(probemap[which(is.na(probemap[,3])),2])
  faulty<-as.numeric(which(is.na(as.character(probemap[,3]))))
  probemap$annotSYM_all[faulty]<-replace
  colnames(probemap)<-c("nuID","targetID","SYMBOL")
}else{
  probemap<-data.frame(cbind(allfeatures,fData(datalist[[i]])[,1],annotSYM_all),stringsAsFactors=F)
  replace<-as.character(probemap[which(is.na(probemap[,3])),2])
  faulty<-as.numeric(which(is.na(as.character(probemap[,3]))))
  probemap$annotSYM_all[faulty]<-replace
  colnames(probemap)<-c("nuID","targetID","SYMBOL")#this is a hack to accommodate incomplete data
}

#new probemap
GEO16<-read.table(file.path(dir.data_root,"Illumina_HT-12v4_CDF_noHead.txt"),header=TRUE,sep="\t",fill=TRUE,quote = "")
newnuID<-IlluminaID2nuID(GEO16$Probes,species="Human")
probemap2<-merge(probemap,newnuID,by.x="nuID",by.y="nuID",all.x=TRUE,all.y=TRUE)
probemap3<-merge(probemap2,GEO16,by.x="Probe_Id",by.y="Probes",all.x=TRUE,all.y=TRUE)

probemap4<-probemap3[which(probemap3$SYMBOL!=probemap3$Symbol&probemap3$Symbol!=""&probemap3$Symbol!=probemap3$targetID),]
probemap5<-probemap4[which(probemap4$SYMBOL!=probemap4$Gene_symbol),]
probemap6<-probemap4[which(probemap4$SYMBOL!=probemap4$Gene_symbol&probemap4$Gene_symbol!=""),]

raData<-read.table(file.path(dir.data_root,"humanHt12v4.txt"),header=TRUE,sep="\t")
probemap7<-merge(probemap3,raData,by.x="Probe_Id",by.y="X.PROBE_ID",all.x=TRUE,all.y=TRUE)
probemap8<-probemap7[which(probemap7$SYMBOL!=probemap7$Gene_symbol.x&probemap7$Gene_symbol.x!=""&probemap7$Gene_symbol.y!=""),]

#DEBUG: commented out
if (neoinsert=="on"){
  # if ((q.i==1 & i==1)|(q.i==4 & i==1)) {#debug; hardcoded
  # 
  #   nodes<-RNeo4j::nodes
  #   #network-specific cypher query
  #   query = "
  #   CREATE (probe:PROBETYPE {name:{nuID},platform:{plat}})
  #   MERGE (gene:SYMBOL {name:{targetID}})
  #   CREATE (probe)-[:mapsTo]->(gene)
  #   "
  #   blocks<-split(1:nrow(probemap),ceiling(seq_along(1:nrow(probemap))/1000))
  #   platform = getChipInfo(datalist[[1]])$chipVersion[1]
  #   
  #   for (block in blocks){
  #   t = suppressMessages(newTransaction(graph))
  #   graphdata=probemap[block,]
  #   for (rowd in 1:nrow(graphdata)) {
  #     probe = as.character(graphdata[rowd, ]$nuID)
  #     gene = as.character(graphdata[rowd, ]$SYMBOL)
  #     
  #       
  #     suppressMessages(appendCypher(t, 
  #                  query, 
  #                  nuID = probe, 
  #                  targetID = gene,
  #                  plat = platform
  #     ))
  #   }
  #   
  #   print("committing :)")
  #   suppressMessages(commit(t))
  #   }#end blocks
  #   #END NEO4J
  # }#end if
  
  ##new
  #if ((q.i==1 & i==1)|(q.i==4 & i==1)) {#debug; hardcoded
  if ((q.i==1 & i==1)|(q.i==2 & i==1)) {#debug; hardcoded
    
    nodes<-RNeo4j::nodes
    
    platform = getChipInfo(datalist[[i]])$chipVersion[1]
    write.csv(probemap,"./tmp.csv")
    
    if (class(datalist[[i]])=="LumiBatch"){
      query = paste("LOAD CSV WITH HEADERS FROM 'file://",file.path(dir.figures,"tmp.csv"),"' AS csvLine CREATE (probe:PROBETYPE {name:csvLine.nuID,platform:'",platform,"'}) MERGE (gene:SYMBOL {name:csvLine.targetID}) CREATE (probe)-[:mapsTo]->(gene)",sep="")
      cypher(graph,query)
    }else{
      query = paste("LOAD CSV WITH HEADERS FROM 'file://",file.path(dir.figures,"tmp.csv"),"' AS csvLine CREATE (probe:PROBETYPE {name:csvLine.nuID,platform:'",platform,"'}) MERGE (gene:SYMBOL {name:csvLine.SYMBOL}) CREATE (probe)-[:mapsTo]->(gene)",sep="")
      cypher(graph,query)
    }
    
    
    #END NEO4J
  }#end if
  #endnew
  
  
}#end neoinsert if

#1_Phenodata analysis####
analysis="1_PhenoData"
an.count<-1
print("Analysis 0: Sample phenotype summaries")

#make relevant output folders
## Define folder for storing results NB do this for each analysis
dir.results <- file.path(dir.output.version,paste("Q_",q.i,"_",dataset.variables,sep=""),analysis, "results")
## Define folder for saving figures
dir.figures <- file.path(dir.output.version,paste("Q_",q.i,"_",dataset.variables,sep=""),analysis, "figures")

for (dir in c(dir.figures, dir.results)) {
  if (!file.exists(dir)) {
    dir.create(dir, recursive=T,showWarnings=T)
  }
}



#Get phenodata
pd2<-pData(datalist[[i]])
#Save phenodata for the given dataset
write.csv(pd2,file=file.path(dir.results,paste("PhenoData",dataset.names[[i]],".csv",sep="")))

#TO DO: project-specific varclass encoded in questions, or separate file (matrixPD)
#TO DO: project-specific matrixPD.csv encodes which variables are to be used in which question. There could be a default.

pdm2<-read.csv(file.path(dir.data_root,matrixPDname),row.names=1,stringsAsFactors = FALSE)

pd3<-pd2[,pdm2[dataset.variables[[i]],]==1]
write.csv(pd3,file=file.path(dir.results,paste("PhenoData_short",dataset.names[[i]],".csv",sep="")))
#need to make a vector of variable type
names(pd2)
#this depends on the structure of the csv file
varclass<-as.character(pdm2["varclass",])
#analyse phenodata

#categorical
pd3cat<-pd2[,pdm2[dataset.variables[[i]],]==1&varclass=="c"]

if(!is.null(pd2$sample_name)){
  rownames(pd3cat)<-pd2$sample_name
}else{
  rownames(pd3cat)<-pd2$sampleID 
}
#numerical
pd3num<-pd2[,pdm2[dataset.variables[[i]],]==1&varclass=="n"]
pd3num<-data.frame(lapply(pd3num,as.numeric))
if(!is.null(pd2$sample_name)){
  rownames(pd3num)<-pd2$sample_name
}else{
  rownames(pd3num)<-pd2$sampleID 
}
pd3num[,contrast.variable[[1]]]<-as.factor(classifier1)
pd3num[,contrast.variable[[2]]]<-as.factor(classifier2)

plistcat<-list()
for (colu in 1:ncol(pd3cat)){
  
  print(paste("The results for:",names(pd3cat)[colu],"   ***********************************"))
  
  theTable_c1<-ftable(xtabs(
    eval(parse(text=paste("~",contrast.variable[[1]],"+",names(pd3cat)[colu],sep="")))
    ,data=pd3cat[complete.cases(pd3cat),]))
  print(theTable_c1)
  
  
  theTable_c12<-as.data.frame(theTable_c1)
  write.csv(theTable_c12,file=file.path(dir.results,paste("counts_",contrast.variable1,names(pd3cat)[colu],".csv",sep="")))
  if (nrow(theTable_c12)>1){
    testres_c1<-chisq.test(theTable_c1)
    print(testres_c1)
    testres_c1_df<-data.frame(cbind(names(testres_c1),as.character(unlist(lapply(testres_c1, paste, collapse=" ")))))
    write.csv(testres_c1_df,file=file.path(dir.results,paste("testres_c1_colu",colu,contrast.variable1,names(pd3cat)[colu],".csv",sep="")))
  }
  
  theTable_c2<-ftable(xtabs(
    eval(parse(text=paste("~",contrast.variable[[2]],"+",names(pd3cat)[colu],sep="")))
    ,data=pd3cat[complete.cases(pd3cat),]))
  print(theTable_c2)
  
  
  theTable_c22<-as.data.frame(theTable_c2)
  write.csv(theTable_c22,file=file.path(dir.results,paste("counts_",contrast.variable2,names(pd3cat)[colu],".csv",sep="")))
  if (nrow(theTable_c22)>1){
    testres_c2<-chisq.test(theTable_c2)
    print(testres_c2)
    testres_c2_df<-data.frame(cbind(names(testres_c2),as.character(unlist(lapply(testres_c2, paste, collapse=" ")))))
    write.csv(testres_c2_df,file=file.path(dir.results,paste("testres_c2_coul_",colu,contrast.variable1,names(pd3cat)[colu],".csv",sep="")))
  }
  plistcat[[colu]]<-NULL
  if (nrow(theTable_c12)>1&nrow(theTable_c22)>1){
    plistcat[[colu]]<-ggplot(data=pd3cat,aes_string(x=names(pd3cat)[colu],fill=names(pd3cat)[colu]))+geom_bar(aes(y = (..count..)/sum(..count..)))+facet_grid(eval(parse(text=paste(contrast.variable[[1]],"~",contrast.variable[[2]],sep=""))))+theme_bw()+theme(axis.text.x = element_text(angle = 45, hjust = 1))+scale_x_discrete(name=paste("Proportion by",contrast.variable2,"\nP=",format(testres_c2$p.value,digits=3)))+scale_y_continuous(name=paste("Proportion by",contrast.variable1,"\nP=",format(testres_c1$p.value,digits=3)))+ggtitle(names(pd3cat)[[colu]])
  }#endif
}

pdf(file.path(dir.figures,"cat.pdf"),height=7*ceiling(ncol(pd3cat)/4),width=1.7*7,onefile=TRUE)
multiplot(plotlist=plistcat,cols=2)
dev.off()

#numeric
#across category subsets
c1cats<-levels(factor(classifier1))
c2cats<-levels(factor(classifier2))

#Describe data sets and subsets
byc1<-describeBy(pd3num[,1:(ncol(pd3num)-2)],classifier1,mat=T)
byc2<-describeBy(pd3num[,1:(ncol(pd3num)-2)],classifier2,mat=T)

#all active TB
data.xs.contrast.variable1.2.all<-pd3num[eval(parse(text=paste("pd3num$",contrast.variable1,sep="")))==c1cats[[1]],1:(ncol(pd3num)-2)]
byc1.2<-describeBy(data.xs.contrast.variable1.2.all,classifier2[eval(parse(text=paste("pd3num$",contrast.variable1,sep="")))==c1cats[[1]]],mat=T)

#all not active TB
data.xs.contrast.variable1.1.all<-pd3num[eval(parse(text=paste("pd3num$",contrast.variable1,sep="")))==c1cats[[2]],1:(ncol(pd3num)-2)]
byc1.1<-describeBy(data.xs.contrast.variable1.1.all,classifier2[eval(parse(text=paste("pd3num$",contrast.variable1,sep="")))==c1cats[[2]]],mat=T)

#all HIV negative
data.xs.contrast.variable2.1.all<-pd3num[eval(parse(text=paste("pd3num$",contrast.variable2,sep="")))==c2cats[[1]],1:(ncol(pd3num)-2)]
byc2.1<-describeBy(data.xs.contrast.variable2.1.all,classifier1[eval(parse(text=paste("pd3num$",contrast.variable2,sep="")))==c2cats[[1]]],mat=T)

#all HIV positive
data.xs.contrast.variable2.2.all<-pd3num[eval(parse(text=paste("pd3num$",contrast.variable2,sep="")))==c2cats[[2]],1:(ncol(pd3num)-2)]
byc2.2<-describeBy(data.xs.contrast.variable2.2.all,classifier1[eval(parse(text=paste("pd3num$",contrast.variable2,sep="")))==c2cats[[2]]],mat=T)

write.csv(byc1,file=file.path(dir.results,"byc1.csv"))
write.csv(byc2,file=file.path(dir.results,"byc2.csv"))
write.csv(byc1.2,file=file.path(dir.results,"byc1.2.csv"))
write.csv(byc1.1,file=file.path(dir.results,"byc1.1.csv"))
write.csv(byc2.1,file=file.path(dir.results,"byc2.1.csv"))
write.csv(byc2.2,file=file.path(dir.results,"byc2.2.csv"))

plistnum<-list()
kwt1=list()
kwt2=list()
kwt.contrast.variable1.1=list()
kwt.contrast.variable1.2=list()
kwt.contrast.variable2.1=list()
kwt.contrast.variable2.2=list()

for (colu in 1:(ncol(pd3num)-2)){
  #stats
  kwt1[[colu]]<-"NA"
  kwt.contrast.variable1.1[[colu]]<-"NA"
  kwt.contrast.variable1.2[[colu]]<-"NA"
  kwt.contrast.variable2.1[[colu]]<-"NA"
  kwt.contrast.variable2.2[[colu]]<-"NA"
  name<-names(pd3num)[colu]
  #across categories
  data.xs<-as.numeric(pd3num[,colu])
  xs1<-try(summary(data.xs~as.factor(classifier1),data=pd3num))
  xs2<-try(summary(data.xs~as.factor(classifier2),data=pd3num))
  kwt1[[colu]]<-NA
  kwt2[[colu]]<-NA
  try(kwt1[[colu]]<-kruskal.test(data.xs~as.factor(classifier1),data=pd3num)$p.value)
  try(kwt2[[colu]]<-kruskal.test(data.xs~as.factor(classifier2),data=pd3num)$p.value)
  
  #all HIV negative
  data.xs.contrast.variable2.1<-as.numeric(pd3num[eval(parse(text=paste("pd3num$",contrast.variable2,sep="")))==c2cats[[1]],colu])
  try(kwt.contrast.variable2.1[[colu]]<-kruskal.test(data.xs.contrast.variable2.1~as.factor(classifier1[eval(parse(text=paste("pd3num$",contrast.variable2,sep="")))==c2cats[[1]]]),data=pd3num)$p.value)
  #all HIV positive
  data.xs.contrast.variable2.2<-as.numeric(pd3num[eval(parse(text=paste("pd3num$",contrast.variable2,sep="")))==c2cats[[2]],colu])
  try(kwt.contrast.variable2.2[[colu]]<-kruskal.test(data.xs.contrast.variable2.2~as.factor(classifier1[eval(parse(text=paste("pd3num$",contrast.variable2,sep="")))==c2cats[[2]]]),data=pd3num)$p.value)
  #all active TB
  data.xs.contrast.variable1.2<-as.numeric(pd3num[eval(parse(text=paste("pd3num$",contrast.variable1,sep="")))==c1cats[[1]],colu])
  try(kwt.contrast.variable1.2[[colu]]<-kruskal.test(data.xs.contrast.variable1.2~as.factor(classifier2[eval(parse(text=paste("pd3num$",contrast.variable1,sep="")))==c1cats[[1]]]),data=pd3num)$p.value)
  #all not active TB
  data.xs.contrast.variable1.1<-as.numeric(pd3num[eval(parse(text=paste("pd3num$",contrast.variable1,sep="")))==c1cats[[2]],colu])
  try(kwt.contrast.variable1.1[[colu]]<-kruskal.test(data.xs.contrast.variable1.1~as.factor(classifier2[eval(parse(text=paste("pd3num$",contrast.variable1,sep="")))==c1cats[[2]]]),data=pd3num)$p.value)
  
  
  #plots  
  plistnum[[colu]]<-ggplot(data=pd3num, aes_string(x=contrast.variable2,fill=contrast.variable2, y=names(pd3num[colu])))+theme_bw()+theme(legend.position="none") + geom_boxplot(outlier.colour = NULL,colour="#333333")+ geom_jitter(position=position_jitter(w=0.1,h=0.1),size=2,colour="#000000")+facet_grid(paste(contrast.variable1,"~.",sep=""))+scale_y_continuous(name=paste("Median difference by",contrast.variable1,"\nP=",format(kwt1[[colu]],digits=3)))+scale_x_discrete(name=paste("Median difference by",contrast.variable2,"\nP=",format(kwt2[[colu]],digits=3)))+ggtitle(names(pd3num[colu]))
  
}

pdf(file.path(dir.figures,"num.pdf"),height=7*ceiling(ncol(pd3num)/4),width=1.7*7,onefile=TRUE)
multiplot(plotlist=plistnum,cols=2)
dev.off()


# table of inner p-vals
resultmat<-matrix(as.numeric(cbind(unlist(kwt.contrast.variable2.1),unlist(kwt.contrast.variable2.2),unlist(kwt.contrast.variable1.2),unlist(kwt.contrast.variable1.1))),ncol=4)
row.names(resultmat)<-names(pd3num[1:(ncol(pd3num)-2)])


colnamecounter<-0
namevector<-vector()
for (c.1 in 1:2){
  for (c.2 in 1:2){
    colnamecounter<-colnamecounter+1
    namevector[colnamecounter]<-paste("by",contrast.variable[[c.1]],"in",contrast.variable[[which(c(1,2)!=c.1)]],"is",list(c2cats,c1cats)[[c.1]][c.2])
  }
}
colnames(resultmat)<-namevector

write.csv(resultmat,file=file.path(dir.results,"pd3num_inner_pvals.csv"))


##2_Expression data: preprocessing ************************####

#TO DO: make platform-specifc preprocessing
#IlluminaHT12v3
#IlluminaHT12v4
#Affy
#RT-PCR

#in loop, source this pre-processing pipeline (specify in questions.R which preprocessing and DE pipelines to use; the output is ranked gene lists written to Neo4j and used to make inputs for WGCNA)
#raw - transform - normalize - filter

analysis="2_Preprocessing"
an.count<-2
figure<-1
print("Analysis 2: Preprocessing")

#make relevant output folders
## Define folder for storing results NB do this for each analysis
dir.results <- file.path(dir.output.version,paste("Q_",q.i,"_",dataset.variables,sep=""),analysis, "results")
## Define folder for saving figures
dir.figures <- file.path(dir.output.version,paste("Q_",q.i,"_",dataset.variables,sep=""),analysis, "figures")

for (dir in c(dir.figures, dir.results)) {
  if (!file.exists(dir)) {
    dir.create(dir, recursive=T,showWarnings=T)
  }
}
setwd(dir.figures)
#Show raw data####

#Outlier plot by edge after removing outliers

pdf(paste('fig_',q.i,'.',an.count,'.',figure,'__outlierPlot.pdf',sep=""))
par(mfrow=c(2,2))
for (edge in 1:length(edges)){
  edgeset<-edges[[edge]]
  plot(datalist[[i]][,edgeset],what="outlier",main=names(edges)[edge])
}
par(mfrow=c(1,1))
dev.off()
figure<-figure+1

legendcolours1<-unique(cbind(as.character(classfactor1[[1]]),as.character(classcolours1[[1]])))
legendcolours2<-unique(cbind(as.character(classfactor2[[1]]),as.character(classcolours2[[1]])))
legendcolourlist<-list(
  unique(cbind(as.character(classfactor1[[1]][edges[[1]]]),as.character(classcolours1[[1]][edges[[1]]]))),
  unique(cbind(as.character(classfactor1[[1]][edges[[2]]]),as.character(classcolours1[[1]][edges[[2]]]))),
  unique(cbind(as.character(classfactor2[[1]][edges[[3]]]),as.character(classcolours2[[1]][edges[[3]]]))),
  unique(cbind(as.character(classfactor2[[1]][edges[[4]]]),as.character(classcolours2[[1]][edges[[4]]]))),
  unique(cbind(as.character(classfactor1[[1]][edges[[5]]]),as.character(classcolours1[[1]][edges[[5]]])))
)

#boxplots and stripcharts

#change identifiers to "prettier" sample name if available in phenodata; useful for plots, etc. Requires sample_name column
identifier<-eval(parse(text=paste("datalist[[i]]$","sample_name",sep="")))
if(!is.null(identifier)){sampleNames(datalist[[i]])<-identifier}#new in R3.3.3, used to be colnames(exprs(X))


# the order of this depends completely on the order in the csv file - it should not impact on the actual analysis. Just need to make sure that the class assignments for DE are correct

pdf(paste('fig_',q.i,'.',an.count,'.',figure,'_boxplot_intensities.pdf',sep=""),height=pdf.options()$height*2)
par(mfrow=c(2,1))
boxplot(datalist[[i]],main=paste(dataset.variables[[i]],contrast.variable1,sep=":"),col=as.character(classcolours1[[i]]))
boxplot(datalist[[i]],main=paste(dataset.variables[[i]],contrast.variable2,sep=":"),col=as.character(classcolours2[[i]]))
par(mfrow=c(1,1))
dev.off()
figure<-figure+1

pdf(paste('fig_',q.i,'.',an.count,'.',figure,'_boxplot_intensities_by_edge.pdf',sep=""),height=pdf.options()$height*1.5,width=pdf.options()$width*1.5)
par(mfrow=c(2,4))
for (edge in 1:length(edges)){
  edgeset<-edges[[edge]]
  boxplot(datalist[[i]][,edgeset],main=paste(names(edges)[edge],contrast.variable1,sep=":"),col=as.character(classcolours1[[i]])[edgeset])
  boxplot(datalist[[i]][,edgeset],main=paste(names(edges)[edge],contrast.variable2,sep=":"),col=as.character(classcolours2[[i]])[edgeset])
}
par(mfrow=c(1,1))
par(parbackup)
dev.off()
figure<-figure+1

#MDS all probes

#FIXREQ
pdf(paste('fig_',q.i,'.',an.count,'.',figure,'_MDS_all_samples_by_',contrast.variable1,'.pdf',sep=""),height=pdf.options()$height*1.1,width=pdf.options()$width*1.1)
plotSampleRelationsAD(datalist[[i]],method="mds",color=as.character(classcolours1[[i]]),plotchar=16)
# legend("bottomleft",
#        levels(classfactor1[[i]]),
#        pch=16,
#        col=levels(classcolours1[[i]]),
#        cex=.9
# )
legend("bottomleft",
       legendcolours1[,1],
       pch=16,
       col=legendcolours1[,2],
       cex=.9
)
dev.off()
figure<-figure+1

# ##removeBatchEffect
# batch=factor(pData(datalist[[i]])$study)
# pdf(paste('fig_',q.i,'.',an.count,'.',figure,'_MDS_all_samples_by_',contrast.variable1,'.pdf',sep=""),height=pdf.options()$height*1.1,width=pdf.options()$width*1.1)
# plotSampleRelationsAD(removeBatchEffect(datalist[[i]],batch),method="mds",color=as.character(classcolours1[[i]]),plotchar=16)
# # legend("bottomleft",
# #        levels(classfactor1[[i]]),
# #        pch=16,
# #        col=levels(classcolours1[[i]]),
# #        cex=.9
# # )
# legend("bottomleft",
#        legendcolours1[,1],
#        pch=16,
#        col=legendcolours1[,2],
#        cex=.9
# )
# dev.off()
# figure<-figure+1
# 
# 
# ##end

pdf(paste('fig_',q.i,'.',an.count,'.',figure,'_MDS_all_samples_by_',contrast.variable2,'.pdf',sep=""),height=pdf.options()$height*1.1,width=pdf.options()$width*1.1)
plotSampleRelationsAD(datalist[[i]],method="mds",color=as.character(classcolours2[[i]]),plotchar=16)
# legend("bottomleft",
#        levels(classfactor2[[i]]),
#        pch=16,
#        col=levels(classcolours2[[i]]),
#        cex=.8
# )
legend("bottomleft",
       legendcolours2[,1],
       pch=16,
       col=legendcolours2[,2],
       cex=.9
)
dev.off()
figure<-figure+1

pdf(paste('fig_',q.i,'.',an.count,'.',figure,'_MDS_all_samples_by_edge.pdf',sep=""),height=pdf.options()$height*1.1,width=pdf.options()$width*1.1)
par(mfrow=c(2,2))
for (edge in 1:length(edges)){
  edgeset<-edges[[edge]]
  plotSampleRelationsAD(datalist[[i]][,edgeset],method="mds",color=as.character(colourfactors[[edge]]),plotchar=16)
  legend("bottomleft",
         legendcolourlist[[edge]][,1],
         pch=16,
         col=legendcolourlist[[edge]][,2],
         cex=.8
  )
}
par(mfrow=c(1,1))
dev.off()
figure<-figure+1


#PCA all probes

#FIXREQ
pdf(paste('fig_',q.i,'.',an.count,'.',figure,'_PCA_all_samples_by_',contrast.variable1,'.pdf',sep=""),height=pdf.options()$height*1.1,width=pdf.options()$width*1.1)
spca <- SamplePCA(exprs(datalist[[i]]), classfactor1[[i]])
plot(spca, col=legendcolours1[,2],main=paste(nrow(exprs(datalist[[i]]))," probes"))#,cex=2,cex.axis=2,cex.main=3,cex.lab=2.5,pch=16)
mtext(sprintf("%.2f",spca@variances[1]/sum(spca@variances)),side=1,line=3,adj=1,cex=1.5)
mtext(sprintf("%.2f",spca@variances[2]/sum(spca@variances)),side=2,line=3,adj=1,cex=1.5)
text(spca@scores[,1],spca@scores[,2],labels=row.names(spca@scores),col="black",cex=0.8)
legend("bottomleft",
       legendcolours1[,1],
       pch=16,
       col=legendcolours1[,2],
       cex=.9
)
# mark the group centers
for (type in 1:length(levels(classfactor1[[i]]))){
  x1 <- predict(spca, matrix(apply(t(exprs(datalist[[i]])[,grep(levels(classfactor1[[i]])[[type]],eval(parse(text=paste("datalist[[i]]$",colour.variable1,sep=""))))]), 2, mean), ncol=1))
  points(x1[1], x1[2], col=legendcolours1[,2][[type]], cex=6,pch=9)}
dev.off()
figure<-figure+1 

pdf(paste('fig_',q.i,'.',an.count,'.',figure,'_PCA_all_samples_by_',contrast.variable2,'.pdf',sep=""),height=pdf.options()$height*1.1,width=pdf.options()$width*1.1)
spca <- SamplePCA(exprs(datalist[[i]]), classfactor2[[i]])
plot(spca, col=legendcolours2[,2],main=paste(nrow(exprs(datalist[[i]]))," probes"))#,cex=2,cex.axis=2,cex.main=3,cex.lab=2.5,pch=16)
mtext(sprintf("%.2f",spca@variances[1]/sum(spca@variances)),side=1,line=3,adj=1,cex=1.5)
mtext(sprintf("%.2f",spca@variances[2]/sum(spca@variances)),side=2,line=3,adj=1,cex=1.5)
text(spca@scores[,1],spca@scores[,2],labels=row.names(spca@scores),col="grey",cex=0.8)
legend("bottomleft",
       legendcolours2[,1],
       pch=16,
       col=legendcolours2[,2],
       cex=.9
)
# mark the group centers
for (type in 1:length(levels(classfactor2[[i]]))){
  x1 <- predict(spca, matrix(apply(t(exprs(datalist[[i]])[,grep(levels(classfactor2[[i]])[[type]],eval(parse(text=paste("datalist[[i]]$",colour.variable2,sep=""))))]), 2, mean), ncol=1))
  points(x1[1], x1[2], col=legendcolours2[,2][[type]], cex=6,pch=9)}
dev.off()
figure<-figure+1 

pdf(paste('fig_',q.i,'.',an.count,'.',figure,'_PCA_all_samples_by_edge.pdf',sep=""),height=pdf.options()$height*1.1,width=pdf.options()$width*1.1)
par(mfrow=c(2,2))
for (edge in 1:length(edges)){
  edgeset<-edges[[edge]]
  spca <- SamplePCA(exprs(datalist[[i]][,edgeset]), classes[[edge]])
  plot(spca, col=legendcolourlist[[edge]][,2],main=paste(nrow(exprs(datalist[[i]][,edgeset]))," probes"))#,cex=2,cex.axis=2,cex.main=3,cex.lab=2.5,pch=16)
  mtext(sprintf("%.2f",spca@variances[1]/sum(spca@variances)),side=1,line=3,adj=1,cex=1.5)
  mtext(sprintf("%.2f",spca@variances[2]/sum(spca@variances)),side=2,line=3,adj=1,cex=1.5)
  text(spca@scores[,1],spca@scores[,2],labels=row.names(spca@scores),col="grey",cex=0.8)
  legend("bottomleft",
         legendcolourlist[[edge]][,1],
         pch=16,
         col=legendcolourlist[[edge]][,2],
         cex=.8
  )
  # mark the group centers
  for (type in 1:length(levels(classes[[edge]]))){
    x1 <- predict(spca, matrix(apply(t(exprs(datalist[[i]][,edgeset])[,grep(levels(classes[[edge]])[[type]],eval(parse(text=paste("datalist[[i]][,edgeset]$",colourvars[[edge]],sep=""))))]), 2, mean), ncol=1))
    points(x1[1], x1[2], col=legendcolourlist[[edge]][,2][[type]], cex=6,pch=9)}
}
par(mfrow=c(1,1))
dev.off()
figure<-figure+1 
#preprocess####

#Variance stabilising transformation or log2 transformation
if(class(datalist[[i]])=="LumiBatch"){
  data.v<-lumiT(datalist[[i]],simpleOutput=FALSE)
  data.v.edges<-list()
  for (edge in 1:length(edges)){
    edgeset=edges[edge]
    data.v.edges[[edge]]<-lumiT(datalist[[i]][,edges[[edge]]],simpleOutput=FALSE)
  }
}else{
  data.v<-lumiT(datalist[[i]],simpleOutput=FALSE,method="log2")
  data.v.edges<-list()
  for (edge in 1:length(edges)){
    edgeset=edges[edge]
    data.v.edges[[edge]]<-lumiT(datalist[[i]][,edges[[edge]]],simpleOutput=FALSE,method="log2")
  }
}
#Quantile normalise
data.q<-lumiN(data.v,method="quantile")
data.q.edges<-list()
#by edge; this is not used anywhere
for (edge in 1:length(edges)){
  data.q.edges[[edge]]<-lumiN(data.v.edges[[edge]],method="quantile")
}
#show preprocessed data####
#boxplots of quantile normalised data

pdf(paste('fig_',q.i,'.',an.count,'.',figure,'_boxplot_intensities_normalised.pdf',sep=""),height=pdf.options()$height*2)
par(mfrow=c(2,1))
boxplot(data.q,main=paste(dataset.variables[[i]],contrast.variable1,sep=":"),col=as.character(classcolours1[[i]]))
boxplot(data.q,main=paste(dataset.variables[[i]],contrast.variable2,sep=":"),col=as.character(classcolours2[[i]]))
par(mfrow=c(1,1))
dev.off()
figure<-figure+1

pdf(paste('fig_',q.i,'.',an.count,'.',figure,'_boxplot_intensities_normalised_by_edge.pdf',sep=""),height=pdf.options()$height*1.3,width=pdf.options()$width*1.3)
par(mfrow=c(2,4))
for (edge in 1:length(edges)){
  edgeset<-edges[[edge]]
  boxplot(data.q.edges[[edge]],main=paste(names(edges)[edge],contrast.variable1,sep=":"),col=as.character(classcolours1[[i]])[edgeset])
  boxplot(data.q.edges[[edge]],main=paste(names(edges)[edge],contrast.variable2,sep=":"),col=as.character(classcolours2[[i]])[edgeset])
}
par(mfrow=c(1,1))
par(parbackup)
dev.off()
figure<-figure+1

#density plots of quantile normalised data

pdf(paste('fig_',q.i,'.',an.count,'.',figure,'_denisty_all_arrays_normalised.pdf',sep=""))
density(data.q,main=dataset.variables[[i]],col=as.character(classcolours1[[i]]),addLegend=FALSE)
abline(h=0,col="red",lwd=2)
dev.off()
figure<-figure+1

pdf(paste('fig_',q.i,'.',an.count,'.',figure,'_denisty_all_arrays_by_edge_normalised.pdf',sep=""))
par(mfrow=c(2,3))
for (edge in 1:length(edges)){
  edgeset<-edges[[edge]]
  density(data.q.edges[[edge]],main=names(edges)[edge],col=list("deeppink","dodgerblue3","forestgreen","goldenrod2","red")[[edge]],addLegend=FALSE)
}
par(mfrow=c(1,1))
par(parbackup)
dev.off()
figure<-figure+1

#MDS all probes (normalised)
#FIXREQ
pdf(paste('fig_',q.i,'.',an.count,'.',figure,'_MDS_all_samples_by_',contrast.variable1,'_normalised.pdf',sep=""),height=pdf.options()$height*1.1,width=pdf.options()$width*1.1)
plotSampleRelationsAD(data.q,method="mds",color=as.character(classcolours1[[i]]),plotchar=16)
legend("bottomleft",
       legendcolours1[,1],
       pch=16,
       col=legendcolours1[,2],
       cex=.9
)
dev.off()  
figure<-figure+1

pdf(paste('fig_',q.i,'.',an.count,'.',figure,'_MDS_all_samples_by_',contrast.variable2,'_normalised.pdf',sep=""),height=pdf.options()$height*1.1,width=pdf.options()$width*1.1)
plotSampleRelationsAD(data.q,method="mds",color=as.character(classcolours2[[i]]),plotchar=16)
legend("bottomleft",
       legendcolours2[,1],
       pch=16,
       col=legendcolours2[,2],
       cex=.9
)
dev.off()  
figure<-figure+1

pdf(paste('fig_',q.i,'.',an.count,'.',figure,'_MDS_all_samples_by_edge_normalised.pdf',sep=""),height=pdf.options()$height*1.1,width=pdf.options()$width*1.1)
par(mfrow=c(2,2))
for (edge in 1:length(edges)){
  edgeset<-edges[[edge]]
  plotSampleRelationsAD(data.q.edges[[edge]],method="mds",color=as.character(colourfactors[[edge]]),plotchar=16)
  legend("bottomleft",
         legendcolourlist[[edge]][,1],
         pch=16,
         col=legendcolourlist[[edge]][,2],
         cex=.9
  )
}
par(mfrow=c(1,1))
dev.off()

#PCA all probes (normalised)

#FIXREQ
pdf(paste('fig_',q.i,'.',an.count,'.',figure,'_PCA_all_samples_by_',contrast.variable1,'_normalised.pdf',sep=""),height=pdf.options()$height*1.1,width=pdf.options()$width*1.1)
spca <- SamplePCA(exprs(data.q), classfactor1[[i]])
plot(spca, col=legendcolours1[,2],main=paste(nrow(exprs(data.q))," probes (normalised)"))#,cex=2,cex.axis=2,cex.main=3,cex.lab=2.5,pch=16)
mtext(sprintf("%.2f",spca@variances[1]/sum(spca@variances)),side=1,line=3,adj=1,cex=1.5)
mtext(sprintf("%.2f",spca@variances[2]/sum(spca@variances)),side=2,line=3,adj=1,cex=1.5)
text(spca@scores[,1],spca@scores[,2],labels=row.names(spca@scores),col="grey",cex=0.8)
legend("bottomleft",
       legendcolours1[,1],
       pch=16,
       col=legendcolours1[,2],
       cex=.9
)
# mark the group centers
for (type in 1:length(levels(classfactor1[[i]]))){
  x1 <- predict(spca, matrix(apply(t(exprs(data.q)[,grep(levels(classfactor1[[i]])[[type]],eval(parse(text=paste("datalist[[i]]$",colour.variable1,sep=""))))]), 2, mean), ncol=1))
  points(x1[1], x1[2], col=legendcolours1[,2][[type]], cex=6,pch=9)}
dev.off()
figure<-figure+1

pdf(paste('fig_',q.i,'.',an.count,'.',figure,'_PCA_all_samples_by_',contrast.variable2,'_normalised.pdf',sep=""),height=pdf.options()$height*1.1,width=pdf.options()$width*1.1)
spca <- SamplePCA(exprs(data.q), classfactor2[[i]])
plot(spca, col=legendcolours2[,2],main=paste(nrow(exprs(data.q))," probes (normalised)"))#,cex=2,cex.axis=2,cex.main=3,cex.lab=2.5,pch=16)
mtext(sprintf("%.2f",spca@variances[1]/sum(spca@variances)),side=1,line=3,adj=1,cex=1.5)
mtext(sprintf("%.2f",spca@variances[2]/sum(spca@variances)),side=2,line=3,adj=1,cex=1.5)
text(spca@scores[,1],spca@scores[,2],labels=row.names(spca@scores),col="grey",cex=0.8)
legend("bottomleft",
       legendcolours2[,1],
       pch=16,
       col=legendcolours2[,2],
       cex=.9
)
# mark the group centers
for (type in 1:length(levels(classfactor2[[i]]))){
  x1 <- predict(spca, matrix(apply(t(exprs(data.q)[,grep(levels(classfactor2[[i]])[[type]],eval(parse(text=paste("datalist[[i]]$",colour.variable2,sep=""))))]), 2, mean), ncol=1))
  points(x1[1], x1[2], col=legendcolours2[,2][[type]], cex=6,pch=9)}
dev.off()
figure<-figure+1


pdf(paste('fig_',q.i,'.',an.count,'.',figure,'_PCA_all_samples_by_edge_normalised.pdf',sep=""),height=pdf.options()$height*1.1,width=pdf.options()$width*1.1)
par(mfrow=c(2,2))
for (edge in 1:length(edges)){
  edgeset<-edges[[edge]]
  spca <- SamplePCA(exprs(data.q.edges[[edge]]),classes[[edge]])
  plot(spca, col=legendcolourlist[[edge]][,2],main=paste(nrow(exprs(data.q.edges[[edge]]))," probes"))#,cex=2,cex.axis=2,cex.main=3,cex.lab=2.5,pch=16)
  mtext(sprintf("%.2f",spca@variances[1]/sum(spca@variances)),side=1,line=3,adj=1,cex=1.5)
  mtext(sprintf("%.2f",spca@variances[2]/sum(spca@variances)),side=2,line=3,adj=1,cex=1.5)
  text(spca@scores[,1],spca@scores[,2],labels=row.names(spca@scores),col="grey",cex=0.8)
  legend("bottomleft",
         legendcolourlist[[edge]][,1],
         pch=16,
         col=legendcolourlist[[edge]][,2],
         cex=.9
  )
  # mark the group centers
  for (type in 1:length(levels(classes[[edge]]))){
    x1 <- predict(spca, matrix(apply(t(exprs(data.q.edges[[edge]])[,grep(levels(classes[[edge]])[[type]],eval(parse(text=paste("data.q.edges[[edge]]$",colourvars[[edge]],sep=""))))]), 2, mean), ncol=1))
    points(x1[1], x1[2], col=legendcolourlist[[edge]][,2][[type]], cex=6,pch=9)}
}
par(mfrow=c(1,1))
dev.off()
figure<-figure+1
#filter####

#Filter 1 for HT12v3
#generate filter based on probe quality
x <- pqlist[[i]]
mapped_probes <- mappedkeys(x)
# Convert to a list
xx <- as.list(x[mapped_probes])
ids <- as.character(nuID2IlluminaID(as.character(featureNames(data.q)),chipVersion=getChipInfo(data.q)$chipVersion[[1]]))
is<-intersect(names(unlist(xx)),ids)#remove probes not in database, then the mget call will not fail
qual <- unlist(mget(is, pqlist[[i]]))
table(qual)

pdf(paste('fig_',q.i,'.',an.count,'.',figure,'_BarplotProbequality.pdf',sep=""),height=pdf.options()$height*1.0,width=pdf.options()$width*1.0)
barplot(table(qual)[c(5,1,2,3,4,6,7,8)],col=rainbow(8))
dev.off()
figure<-figure+1

pdf(paste('fig_',q.i,'.',an.count,'.',figure,'_PIE_Probequality.pdf',sep=""),height=pdf.options()$height*1.0,width=pdf.options()$width*1.0)
pie(table(qual)[c(5,1,2,3,4,6,7,8)],col=rainbow(8),border=FALSE)
dev.off()
figure<-figure+1

rem <- qual == "No match" | qual == "Bad" | is.na(qual)
#filter raw, vst and qnorm data

data.raw.fil<- datalist[[i]][!rem, ]  
data.q.fil<- data.q[!rem, ]
data.v.fil<- data.v[!rem, ]

#########NEW FILTER for HT12v4
if(filtervar=="new"){
  #other pipeline
  #NB! this file is only 44K rows, not 47K!!! (control probes missing?)
  rad<-read.table(file.path(dir.data_root,"humanHt12v4.txt"),header=TRUE,sep="\t")
  #rad[grep("^GBP.*",rad$Gene_symbol),]
  #hist(rad$BWA_NUM_MISMATCHES)
  bad<-as.character(rad[which(rad$uniq!=1),]$X.PROBE_ID)
  bad.nuID<-IlluminaID2nuID(bad)[,"nuID"]
  good<-as.character(rad[which(rad$uniq==1),]$X.PROBE_ID)
  good.nuID<-IlluminaID2nuID(good)[,"nuID"]
  rem2<-rad$uniq!=1
  names(rem2)<-rad$X.PROBE_ID[which(rad$uniq!=1)]
  rad2<-rad[which(rad$uniq!=1),]
  
  # data.raw.fil<- datalist[[i]][!rem2, ]  
  # data.q.fil<- data.q[!rem2, ]
  # data.v.fil<- data.v[!rem2, ]
  
  data.raw.fil<- datalist[[i]][which(!rownames(datalist[[i]])%in%bad.nuID),]  
  data.q.fil<- data.q[which(!rownames(data.q)%in%bad.nuID),]
  data.v.fil<- data.v[which(!rownames(data.v)%in%bad.nuID),]
}  
#########END NEW FILTER

data.raw.fil.edges<-list()
for (edge in 1:length(edges)){
  edgeset=edges[[edge]]
  data.raw.fil.edges[[edge]]<-data.raw.fil[,edgeset]
}

data.q.fil.edges<-list()
for (edge in 1:length(edges)){
  edgeset=edges[[edge]]
  data.q.fil.edges[[edge]]<-data.q.fil[,edgeset]
}

data.v.fil.edges<-list()
for (edge in 1:length(edges)){
  edgeset=edges[[edge]]
  data.v.fil.edges[[edge]]<-data.v.fil[,edgeset]
}
#show filtered probes ####

#filtered probes (PROBEQUALITY)

pdf(paste('fig_',q.i,'.',an.count,'.',figure,'_boxplot_intensities_filtered.pdf',sep=""),height=pdf.options()$height*2,width=pdf.options()$width*1.5)
par(mfrow=c(2,1),mar=c(7,5,1,1))
boxplot(data.q.fil,las=2,col=as.character(classcolours1[[i]]),main=paste(dataset.variables[[i]],contrast.variable1,"filtered",sep=":"))
boxplot(data.q.fil,las=2,col=as.character(classcolours2[[i]]),main=paste(dataset.variables[[i]],contrast.variable2,"filtered",sep=":"))
par(parbackup)
dev.off()
figure<-figure+1

pdf(paste('fig_',q.i,'.',an.count,'.',figure,'_boxplot_intensities_filtered_by_edge.pdf',sep=""),height=pdf.options()$height*1.2,width=pdf.options()$width*1.5)
par(mfrow=c(2,4))
for (edge in 1:length(edges)){
  edgeset<-edges[[edge]]
  boxplot(data.q.fil.edges[[edge]],main=paste(names(edges)[edge],contrast.variable1,sep=":"),col=as.character(classcolours1[[i]])[edgeset])
  boxplot(data.q.fil.edges[[edge]],main=paste(names(edges)[edge],contrast.variable2,sep=":"),col=as.character(classcolours2[[i]])[edgeset])
}
par(mfrow=c(1,1))
par(parbackup)
dev.off()
figure<-figure+1


#MDS filtered probes
#FIXREQ
pdf(paste('fig_',q.i,'.',an.count,'.',figure,'_MDS_all_samples_by_',contrast.variable1,'_filtered.pdf',sep=""),height=pdf.options()$height*1.0,width=pdf.options()$width*1.0)
plotSampleRelationsAD(data.q.fil,method="mds",color=as.character(classcolours1[[i]]),plotchar=16)
legend("bottomleft",
       legendcolours1[,1],
       pch=16,
       col=legendcolours1[,2],
       cex=.9
)
dev.off()
figure<-figure+1

pdf(paste('fig_',q.i,'.',an.count,'.',figure,'_MDS_all_samples_by_',contrast.variable2,'_filtered.pdf',sep=""),height=pdf.options()$height*1.0,width=pdf.options()$width*1.0)
plotSampleRelationsAD(data.q.fil,method="mds",color=as.character(classcolours2[[i]]),plotchar=16)
legend("bottomleft",
       legendcolours2[,1],
       pch=16,
       col=legendcolours2[,2],
       cex=.9
)
dev.off()
figure<-figure+1

pdf(paste('fig_',q.i,'.',an.count,'.',figure,'_MDS_all_samples_by_edge_filtered.pdf',sep=""),height=pdf.options()$height*1.2,width=pdf.options()$width*1.2)
par(mfrow=c(2,2))
for (edge in 1:length(edges)){
  edgeset<-edges[[edge]]
  plotSampleRelationsAD(data.q.fil.edges[[edge]],method="mds",color=as.character(colourfactors[[edge]]),plotchar=16)
  legend("bottomleft",
         legendcolourlist[[edge]][,1],
         pch=16,
         col=legendcolourlist[[edge]][,2],
         cex=.9
  )
}
par(mfrow=c(1,1))
dev.off()
figure<-figure+1


#PCA filtered probes

#FIXREQ
pdf(paste('fig_',q.i,'.',an.count,'.',figure,'_PCA_all_samples_by_',contrast.variable1,'_filtered.pdf',sep=""),height=pdf.options()$height*1.0,width=pdf.options()$width*1.0)
spca <- SamplePCA(exprs(data.q.fil), classfactor1[[i]])
plot(spca, col=legendcolours1[,2],main=paste(nrow(exprs(data.q.fil))," probes (filtered)"))#,cex=2,cex.axis=2,cex.main=3,cex.lab=2.5,pch=16)
mtext(sprintf("%.2f",spca@variances[1]/sum(spca@variances)),side=1,line=3,adj=1,cex=1.5)
mtext(sprintf("%.2f",spca@variances[2]/sum(spca@variances)),side=2,line=3,adj=1,cex=1.5)
text(spca@scores[,1],spca@scores[,2],labels=row.names(spca@scores),col="grey",cex=0.8)
legend("bottomleft",
       legendcolours1[,1],
       pch=16,
       col=legendcolours1[,2],
       cex=.9
)
# mark the group centers
for (type in 1:length(levels(classfactor1[[i]]))){
  x1 <- predict(spca, matrix(apply(t(exprs(data.q.fil)[,grep(levels(classfactor1[[i]])[[type]],eval(parse(text=paste("datalist[[i]]$",colour.variable1,sep=""))))]), 2, mean), ncol=1))
  points(x1[1], x1[2], col=legendcolours1[,2][[type]], cex=6,pch=9)}
dev.off()
figure<-figure+1

pdf(paste('fig_',q.i,'.',an.count,'.',figure,'_PCA_all_samples_by_',contrast.variable2,'_filtered.pdf',sep=""),height=pdf.options()$height*1.0,width=pdf.options()$width*1.0)
spca <- SamplePCA(exprs(data.q.fil), classfactor2[[i]])
plot(spca, col=legendcolours2[,2],main=paste(nrow(exprs(data.q.fil))," probes (filtered)"))#,cex=2,cex.axis=2,cex.main=3,cex.lab=2.5,pch=16)
mtext(sprintf("%.2f",spca@variances[1]/sum(spca@variances)),side=1,line=3,adj=1,cex=1.5)
mtext(sprintf("%.2f",spca@variances[2]/sum(spca@variances)),side=2,line=3,adj=1,cex=1.5)
text(spca@scores[,1],spca@scores[,2],labels=row.names(spca@scores),col="grey",cex=0.8)
legend("bottomleft",
       legendcolours2[,1],
       pch=16,
       col=legendcolours2[,2],
       cex=.9
)
# mark the group centers
for (type in 1:length(levels(classfactor2[[i]]))){
  x1 <- predict(spca, matrix(apply(t(exprs(data.q.fil)[,grep(levels(classfactor2[[i]])[[type]],eval(parse(text=paste("datalist[[i]]$",colour.variable2,sep=""))))]), 2, mean), ncol=1))
  points(x1[1], x1[2], col=legendcolours2[,2][[type]], cex=6,pch=9)}
dev.off()
figure<-figure+1

pdf(paste('fig_',q.i,'.',an.count,'.',figure,'_PCA_all_samples_by_edge_filtered.pdf',sep=""),height=pdf.options()$height*1.2,width=pdf.options()$width*1.2)
spca <- SamplePCA(exprs(data.q.fil), classfactor2[[i]])
par(mfrow=c(2,2))
for (edge in 1:length(edges)){
  edgeset<-edges[[edge]]
  spca <- SamplePCA(exprs(data.q.fil.edges[[edge]]),classes[[edge]])
  plot(spca, col=legendcolourlist[[edge]][,2],main=paste(nrow(exprs(data.q.fil.edges[[edge]]))," probes"))#,cex=2,cex.axis=2,cex.main=3,cex.lab=2.5,pch=16)
  mtext(sprintf("%.2f",spca@variances[1]/sum(spca@variances)),side=1,line=3,adj=1,cex=1.5)
  mtext(sprintf("%.2f",spca@variances[2]/sum(spca@variances)),side=2,line=3,adj=1,cex=1.5)
  text(spca@scores[,1],spca@scores[,2],labels=row.names(spca@scores),col="lightgrey",cex=0.6)
  legend("bottomleft",
         legendcolourlist[[edge]][,1],
         pch=16,
         col=legendcolourlist[[edge]][,2],
         cex=.9
  )
  # mark the group centers
  for (type in 1:length(levels(classes[[edge]]))){
    x1 <- predict(spca, matrix(apply(t(exprs(data.q.fil.edges[[edge]])[,grep(levels(classes[[edge]])[[type]],eval(parse(text=paste("data.q.fil.edges[[edge]]$",colourvars[[edge]],sep=""))))]), 2, mean), ncol=1))
    points(x1[1], x1[2], col=legendcolourlist[[edge]][,2][[type]], cex=6,pch=9)}
}
par(mfrow=c(1,1))
dev.off()
figure<-figure+1

#PCA of all square 1 in 4 categories

a<-levels(as.factor(paste(classifier1,classifier2)))
b<-paste(classifier1,classifier2)

#FIXREQ
pdf(paste('fig_',q.i,'.',an.count,'.',figure,'_PCA_all_samples_by_category_filtered.pdf',sep=""),height=pdf.options()$height*1.0,width=pdf.options()$width*1.0)
spca <- SamplePCA(exprs(data.q.fil), as.factor(paste(classifier1,classifier2)))
plot(spca, col=brewer.pal(4,"Set2"),main=paste(nrow(exprs(data.q.fil))," probes (filtered)"))#,cex=2,cex.axis=2,cex.main=3,cex.lab=2.5,pch=16)
mtext(sprintf("%.2f",spca@variances[1]/sum(spca@variances)),side=1,line=3,adj=1,cex=1.5)
mtext(sprintf("%.2f",spca@variances[2]/sum(spca@variances)),side=2,line=3,adj=1,cex=1.5)
text(spca@scores[,1],spca@scores[,2],labels=row.names(spca@scores),col="lightgrey",cex=0.6)
legend("bottomleft",
       levels(as.factor(paste(classifier1,classifier2))),
       pch=c(9,10,12,13),
       col=brewer.pal(4,"Set2"),
       cex=.8)
# mark the group centers
x1 <- predict(spca, matrix(apply(t(exprs(data.q.fil)[,which(b==a[1])]), 2, mean), ncol=1))
points(x1[1], x1[2], col=brewer.pal(4,"Set2")[1], cex=6,pch=9)
x2 <- predict(spca, matrix(apply(t(exprs(data.q.fil)[,which(b==a[2])]), 2, mean), ncol=1))
points(x2[1], x2[2], col=brewer.pal(4,"Set2")[2], cex=6,pch=10)
x3 <- predict(spca, matrix(apply(t(exprs(data.q.fil)[,which(b==a[3])]), 2, mean), ncol=1))
points(x3[1], x3[2], col=brewer.pal(4,"Set2")[3], cex=6,pch=12)
x4 <- predict(spca, matrix(apply(t(exprs(data.q.fil)[,which(b==a[4])]), 2, mean), ncol=1))
points(x4[1], x4[2], col=brewer.pal(4,"Set2")[4], cex=6,pch=13)
dev.off()
figure<-figure+1

#PCA of all square 1 in 8 categories

#FIXREQ
pdf(paste('fig_',q.i,'.',an.count,'.',figure,'_PCA_all_samples_by_6_categories_filtered.pdf',sep=""),height=pdf.options()$height*1.0,width=pdf.options()$width*1.0)
spca <- SamplePCA(exprs(data.q.fil), as.factor(paste(classfactor1[[i]],classfactor2[[i]])
))
plot(spca, col=brewer.pal(8,"Paired"),main=paste(nrow(exprs(data.q.fil))," probes (filtered)"))#,cex=2,cex.axis=2,cex.main=3,cex.lab=2.5,pch=16)
mtext(sprintf("%.2f",spca@variances[1]/sum(spca@variances)),side=1,line=3,adj=1,cex=1.5)
mtext(sprintf("%.2f",spca@variances[2]/sum(spca@variances)),side=2,line=3,adj=1,cex=1.5)
text(spca@scores[,1],spca@scores[,2],labels=row.names(spca@scores),col="grey",cex=0.5)
legend("bottomleft",
       levels(as.factor(paste(classfactor1[[i]],classfactor2[[i]])
       )),
       pch=c(9,9,10,10,12,12,13,13),
       col=brewer.pal(8,"Paired"),
       cex=.65)
# mark the group centers
for (cat1 in 1:(length(levels(as.factor(paste(classfactor1[[i]])))))){
  for (cat2 in 1:(length(levels(as.factor(paste(classfactor2[[i]])))))){
    selector=which(eval(parse(text=paste("data.q.fil$",colour.variable1,sep="")))==levels(classfactor1[[i]])[cat1]&
                     eval(parse(text=paste("data.q.fil$",colour.variable2,sep="")))==levels(classfactor2[[i]])[cat2])
    x1 <- predict(spca, matrix(apply(t(exprs(data.q.fil)[,selector]), 2, mean), ncol=1))
    points(x1[1], x1[2], col=matrix(brewer.pal(8,"Paired"),nrow=2)[cat2,cat1], cex=6,pch=c(9,10,12,13)[cat1])
  }
}
dev.off()
figure<-figure+1


##3_Class Description: DE analysis************************####
analysis="3_DE"
an.count<-3
figure<-1
print("Analysis 3: Differential expression")


## Define folder for storing results
dir.results <- file.path(dir.output.version,paste("Q_",q.i,"_",dataset.variables,sep=""),analysis, "results")
## Define folder for saving figures
dir.figures <- file.path(dir.output.version,paste("Q_",q.i,"_",dataset.variables,sep=""),analysis, "figures")

for (dir in c(dir.figures, dir.results)) {
  if (!file.exists(dir)) {
    dir.create(dir, recursive=T,showWarnings=T)
  }
}
setwd(dir.figures)
#Flat limma####
flat<-classifier1
flat<-factor(flat)
design<-model.matrix(~0+flat)
colnames(design)<-levels(flat)
# batch=NULL
# if("study"%in%colnames(pData(datalist[[i]]))){
#   batch=factor(pData(datalist[[i]])$study)
#   if(length(levels(batch))>1){
#     design<-model.matrix(~0+flat+batch)
#     colnames(design)<-c(levels(flat),paste("batch",1:(ncol(design)-length(levels(flat))),sep=""))
#     
#   }
#   }



flatfit<-lmFit(data.q.fil,design) 

v0a<-pathway1[[2]]
v0b<-pathway1[[1]]

cont.matrix.flat<-makeContrasts(
  edge5=eval(parse(text=v0a))-eval(parse(text=v0b)),
  levels=design
)
flatfit2<-contrasts.fit(flatfit,cont.matrix.flat)
flatfit2<-eBayes(flatfit2)
#flatfit2$genes$ID<-featureData(data.q.fil)$TargetID
if(filtervar=="new"){flatfit2$genes$ID<-probemap$SYMBOL[which(!rownames(data.v)%in%bad.nuID)]
}else{flatfit2$genes$ID<-probemap$SYMBOL[!rem]}



#three different methods of multiple testing adjustment are illustrated here. See limma userguide for details
flat.results.hierarchical.1<-decideTests(flatfit2,method="hierarchical")
flat.results.global.1<-decideTests(flatfit2,method="global")
flat.results.separate.1<-decideTests(flatfit2,method="separate")

pdf(paste('fig_',q.i,'.',an.count,'.',figure,'_Venn_Diagram_3BH_methods_flat.pdf',sep=""),height=pdf.options()$height*0.7,width=pdf.options()$width*1.2)
par(mfrow=c(1,3))
vennDiagram(flat.results.hierarchical.1,main=paste(contrast.variable1,"(hierarchical P)"),include=c("up","down"),counts.col=c("red","darkgreen"))
vennDiagram(flat.results.global.1,main=paste(contrast.variable1,"(global P)"),include=c("up","down"),counts.col=c("red","darkgreen"))
vennDiagram(flat.results.separate.1,main=paste(contrast.variable1,"(separate P)"),include=c("up","down"),counts.col=c("red","darkgreen"))
par(mfrow=c(1,1))
dev.off()
figure<-figure+1

#Top tables for various numbers of top probes for edge1, edge2 and the interaction i.e. all significant (P<0.05, separate P), top 500 and top 12000
edge5.flat.TT.sig=topTable(flatfit2,coef=1,sort.by="B",resort.by="logFC",number=Inf,adjust.method="BH",p.value=0.05)#TT uses separate BH adjustment
edge5.flat.TT.500=topTable(flatfit2,coef=1,sort.by="B",resort.by="logFC",number=500)
edge5.flat.TT.12K=topTable(flatfit2,coef=1,sort.by="B",resort.by="logFC",number=12000)
edge5.flat.TT.all=topTable(flatfit2,coef=1,sort.by="B",resort.by="logFC",number=Inf)

#volcanoplot
pdf(paste('fig_',q.i,'.',an.count,'.',figure,'_Volcanoplot_3BH_methods_flat.pdf',sep=""),height=pdf.options()$height*1,width=pdf.options()$width*1)
par(mfrow=c(1,1))
volcanoplot(flatfit2,coef=1,highlight=30,col="grey",names=flatfit2$genes$ID,main="Edge 5")

par(mfrow=c(1,1))
dev.off()
figure<-figure+1
#Factorial_limma and all results exported (edges1-5)####

#DE statistics:Factorial design: normalised and filtered data stratified by TB class and sex.
TS<-paste(classifier1,classifier2,sep=".")
TS<-factor(TS)
design<-model.matrix(~0+TS)
colnames(design)<-levels(TS)
# batch=NULL
# if("study"%in%colnames(pData(datalist[[i]]))){
#   batch=factor(pData(datalist[[i]])$study)
#   if(length(levels(batch))>1){
#     design<-model.matrix(~0+TS+batch)
#     colnames(design)<-c(levels(TS),paste("batch",1:(ncol(design)-length(levels(TS))),sep=""))
#     
#   }
# }



fit<-lmFit(data.q.fil,design) 

#Defining individual edges for use in the contrast matrix
#the pathway variable is used as it always differentiates cases from controls, and is not affected by lexicographical order. v1-4 represent the 4 edges of the square
v1<-paste(pathway1[[2]],pathway2[[1]],sep=".")
v2<-paste(pathway1[[2]],pathway2[[2]],sep=".")
v3<-paste(pathway1[[1]],pathway2[[1]],sep=".")
v4<-paste(pathway1[[1]],pathway2[[2]],sep=".")

#Contrast matrix 1: to be used to extract results for edge1, edge2 and the interaction
#The contrasts and edges are synonymous
cont.matrix1<-makeContrasts(
  edge1=eval(parse(text=v1))-eval(parse(text=v3)),
  edge2=eval(parse(text=v2))-eval(parse(text=v4)),
  Diff1=(eval(parse(text=v1))-eval(parse(text=v3)))-(eval(parse(text=v2))-eval(parse(text=v4))),
  levels=design
)

#extract the results for the 4 edges/ contrasts
fit2<-contrasts.fit(fit,cont.matrix1)
fit2<-eBayes(fit2)
if(filtervar=="new"){fit2$genes$ID<-probemap$SYMBOL[which(!rownames(data.v)%in%bad.nuID)]
}else{fit2$genes$ID<-probemap$SYMBOL[!rem]}


#three different methods of multiple testing adjustment are illustrated here. See limma userguide for details
results.hierarchical.1<-decideTests(fit2,method="hierarchical")
results.global.1<-decideTests(fit2,method="global")
results.separate.1<-decideTests(fit2,method="separate")

pdf(paste('fig_',q.i,'.',an.count,'.',figure,'_Venn_Diagram_3BH_methods.pdf',sep=""),height=pdf.options()$height*0.7,width=pdf.options()$width*1.2)
par(mfrow=c(1,3))
vennDiagram(results.hierarchical.1,main=paste(contrast.variable1,"by",contrast.variable2,"(hierarchical P)"),include=c("up","down"),counts.col=c("red","darkgreen"))
vennDiagram(results.global.1,main=paste(contrast.variable1,"by",contrast.variable2,"(global P)"),include=c("up","down"),counts.col=c("red","darkgreen"))
vennDiagram(results.separate.1,main=paste(contrast.variable1,"by",contrast.variable2,"(separate P)"),include=c("up","down"),counts.col=c("red","darkgreen"))
par(mfrow=c(1,1))
dev.off()
figure<-figure+1

#Top tables for various numbers of top probes for edge1, edge2 and the interaction i.e. all significant (P<0.05, separate P), top 500 and top 12000
edge1.fac.TT.sig=topTable(fit2,coef=1,sort.by="B",resort.by="logFC",number=Inf,adjust.method="BH",p.value=0.05)#TT uses separate BH adjustment
edge1.fac.TT.500=topTable(fit2,coef=1,sort.by="B",resort.by="logFC",number=500)
edge1.fac.TT.12K=topTable(fit2,coef=1,sort.by="B",resort.by="logFC",number=12000)

edge2.fac.TT.sig=topTable(fit2,coef=2,sort.by="B",resort.by="logFC",number=Inf,adjust.method="BH",p.value=0.05)#TT uses separate BH adjustment
edge2.fac.TT.500=topTable(fit2,coef=2,sort.by="B",resort.by="logFC",number=500)
edge2.fac.TT.12K=topTable(fit2,coef=2,sort.by="B",resort.by="logFC",number=12000)

diff1.fac.TT.sig=topTable(fit2,coef=3,sort.by="B",resort.by="logFC",number=Inf,adjust.method="BH",p.value=0.05)#TT uses separate BH adjustment
diff1.fac.TT.500=topTable(fit2,coef=3,sort.by="B",resort.by="logFC",number=500)
diff1.fac.TT.12K=topTable(fit2,coef=3,sort.by="B",resort.by="logFC",number=12000)

#volcanoplot
pdf(paste('fig_',q.i,'.',an.count,'.',figure,'_Volcanoplot_3BH_methods.pdf',sep=""),height=pdf.options()$height*0.7,width=pdf.options()$width*1.2)
par(mfrow=c(1,3))
volcanoplot(fit2,coef=1,highlight=30,col="grey",names=fit2$genes$ID,main="Edge 1")
volcanoplot(fit2,coef=2,highlight=30,col="grey",names=fit2$genes$ID,main="Edge 2")
volcanoplot(fit2,coef=3,highlight=30,col="grey",names=fit2$genes$ID,main="Interaction")
par(mfrow=c(1,1))
dev.off()
figure<-figure+1

#Contrast matrix 2: to be used to extract results for edge3, edge4 and the interaction
#The contrasts and edges are synonymous
cont.matrix2<-makeContrasts(
  edge3=eval(parse(text=v4))-eval(parse(text=v3)),
  edge4=eval(parse(text=v2))-eval(parse(text=v1)),
  Diff2=(eval(parse(text=v4))-eval(parse(text=v3)))-(eval(parse(text=v2))-eval(parse(text=v1))),
  levels=design
)

fit3<-contrasts.fit(fit,cont.matrix2)
fit3<-eBayes(fit3)
if(filtervar=="new"){fit3$genes$ID<-probemap$SYMBOL[which(!rownames(data.v)%in%bad.nuID)]
}else{fit3$genes$ID<-probemap$SYMBOL[!rem]}


results.hierarchical.2<-decideTests(fit3,method="hierarchical",p.value=0.05)
results.global.2<-decideTests(fit3,method="global",p.value=0.05)
results.separate.2<-decideTests(fit3,method="separate",p.value=0.05)

pdf(paste('fig_',q.i,'.',an.count,'.',figure,'_Venn_Diagram_3BH_methods_contrastMatrix2.pdf',sep=""),height=pdf.options()$height*0.7,width=pdf.options()$width*1.2)
par(mfrow=c(1,3))
vennDiagram(results.hierarchical.2,main=paste(contrast.variable1,"by",contrast.variable2,"(hierarchical P)"),include=c("up","down"),counts.col=c("red","darkgreen"))
vennDiagram(results.global.2,main=paste(contrast.variable1,"by",contrast.variable2,"(global P)"),include=c("up","down"),counts.col=c("red","darkgreen"))
vennDiagram(results.separate.2,main=paste(contrast.variable1,"by",contrast.variable2,"(separate P)"),include=c("up","down"),counts.col=c("red","darkgreen"))
par(mfrow=c(1,1))
dev.off()
figure<-figure+1

edge3.fac.TT.sig=topTable(fit3,coef=1,sort.by="B",resort.by="logFC",number=Inf,adjust.method="BH",p.value=0.05)#TT uses separate BH adjustment
edge3.fac.TT.500=topTable(fit3,coef=1,sort.by="B",resort.by="logFC",number=500)
edge3.fac.TT.12K=topTable(fit3,coef=1,sort.by="B",resort.by="logFC",number=12000)

edge4.fac.TT.sig=topTable(fit3,coef=2,sort.by="B",resort.by="logFC",number=Inf,adjust.method="BH",p.value=0.05)#TT uses separate BH adjustment
edge4.fac.TT.500=topTable(fit3,coef=2,sort.by="B",resort.by="logFC",number=500)
edge4.fac.TT.12K=topTable(fit3,coef=2,sort.by="B",resort.by="logFC",number=12000)

diff2.fac.TT.sig=topTable(fit3,coef=3,sort.by="B",resort.by="logFC",number=Inf,adjust.method="BH",p.value=0.05)#TT uses separate BH adjustment
diff2.fac.TT.500=topTable(fit3,coef=3,sort.by="B",resort.by="logFC",number=500)
diff2.fac.TT.12K=topTable(fit3,coef=3,sort.by="B",resort.by="logFC",number=12000)

#volcanoplot
pdf(paste('fig_',q.i,'.',an.count,'.',figure,'_Volcanoplot_3BH_methods_CM2.pdf',sep=""),height=pdf.options()$height*0.7,width=pdf.options()$width*1.2)
par(mfrow=c(1,3))
volcanoplot(fit3,coef=1,highlight=30,col="grey",names=fit3$genes$ID,main="Edge 3")
volcanoplot(fit3,coef=2,highlight=30,col="grey",names=fit3$genes$ID,main="Edge 4")
volcanoplot(fit3,coef=3,highlight=30,col="grey",names=fit3$genes$ID,main="Interaction")
par(mfrow=c(1,1))
dev.off()
figure<-figure+1

#Contrast matrix 3: to be used to extract results for all 4 edges separately
#The contrasts and edges are synonymous
cont.matrix3<-makeContrasts(
  edge1=eval(parse(text=v1))-eval(parse(text=v3)),
  edge2=eval(parse(text=v2))-eval(parse(text=v4)),
  edge3=eval(parse(text=v4))-eval(parse(text=v3)),
  edge4=eval(parse(text=v2))-eval(parse(text=v1)),
  levels=design
)
fit4<-contrasts.fit(fit,cont.matrix3)
fit4<-eBayes(fit4)
if(filtervar=="new"){fit4$genes$ID<-probemap$SYMBOL[which(!rownames(data.v)%in%bad.nuID)]
}else{fit4$genes$ID<-probemap$SYMBOL[!rem]}

##this needs a complete lumiHumanID mapping to get at HGNC symbols

results.hierarchical.3<-decideTests(fit4,method="hierarchical")
results.global.3<-decideTests(fit4,method="global")
results.separate.3<-decideTests(fit4,method="separate")

pdf(paste('fig_',q.i,'.',an.count,'.',figure,'_Venn_Diagram_3BH_methods_contrastMatrix3.pdf',sep=""),height=pdf.options()$height*0.7,width=pdf.options()$width*1.2)
par(mfrow=c(1,3))
vennDiagram(results.hierarchical.3,main=paste("4 edges","(hierarchical P)"),include=c("up","down"),counts.col=c("red","darkgreen"))
vennDiagram(results.global.3,main=paste("4 edges","(global P)"),include=c("up","down"),counts.col=c("red","darkgreen"))
vennDiagram(results.separate.3,main=paste("4 edges","(separate P)"),include=c("up","down"),counts.col=c("red","darkgreen"))
par(mfrow=c(1,1))
dev.off()
figure<-figure+1

#Extract various numbers of probes sorted by log-likelyhood of differential abundance
edge1.fac.sep.TT.sig=topTable(fit4,coef=1,sort.by="B",resort.by="logFC",number=Inf,adjust.method="BH",p.value=0.05)#TT uses separate BH adjustment
edge1.fac.sep.TT.500=topTable(fit4,coef=1,sort.by="B",resort.by="logFC",number=500)
edge1.fac.sep.TT.4K=topTable(fit4,coef=1,sort.by="B",resort.by="logFC",number=4000)
edge1.fac.sep.TT.6K=topTable(fit4,coef=1,sort.by="B",resort.by="logFC",number=6000)
edge1.fac.sep.TT.8K=topTable(fit4,coef=1,sort.by="B",resort.by="logFC",number=8000)
edge1.fac.sep.TT.12K=topTable(fit4,coef=1,sort.by="B",resort.by="logFC",number=12000)
edge1.fac.sep.TT.all=topTable(fit4,coef=1,sort.by="B",resort.by="logFC",number=Inf)#all probes

edge2.fac.sep.TT.sig=topTable(fit4,coef=2,sort.by="B",resort.by="logFC",number=Inf,adjust.method="BH",p.value=0.05)#TT uses separate BH adjustment
edge2.fac.sep.TT.500=topTable(fit4,coef=2,sort.by="B",resort.by="logFC",number=500)
edge2.fac.sep.TT.4K=topTable(fit4,coef=2,sort.by="B",resort.by="logFC",number=4000)
edge2.fac.sep.TT.6K=topTable(fit4,coef=2,sort.by="B",resort.by="logFC",number=6000)
edge2.fac.sep.TT.8K=topTable(fit4,coef=2,sort.by="B",resort.by="logFC",number=8000)
edge2.fac.sep.TT.12K=topTable(fit4,coef=2,sort.by="B",resort.by="logFC",number=12000)
edge2.fac.sep.TT.all=topTable(fit4,coef=2,sort.by="B",resort.by="logFC",number=Inf)

edge3.fac.sep.TT.sig=topTable(fit4,coef=3,sort.by="B",resort.by="logFC",number=Inf,adjust.method="BH",p.value=0.05)#TT uses separate BH adjustment
edge3.fac.sep.TT.500=topTable(fit4,coef=3,sort.by="B",resort.by="logFC",number=500)
edge3.fac.sep.TT.4K=topTable(fit4,coef=3,sort.by="B",resort.by="logFC",number=4000)
edge3.fac.sep.TT.6K=topTable(fit4,coef=3,sort.by="B",resort.by="logFC",number=6000)
edge3.fac.sep.TT.8K=topTable(fit4,coef=3,sort.by="B",resort.by="logFC",number=8000)
edge3.fac.sep.TT.12K=topTable(fit4,coef=3,sort.by="B",resort.by="logFC",number=12000)
edge3.fac.sep.TT.all=topTable(fit4,coef=3,sort.by="B",resort.by="logFC",number=Inf)

edge4.fac.sep.TT.sig=topTable(fit4,coef=4,sort.by="B",resort.by="logFC",number=Inf,adjust.method="BH",p.value=0.05)#TT uses separate BH adjustment
edge4.fac.sep.TT.500=topTable(fit4,coef=4,sort.by="B",resort.by="logFC",number=500)
edge4.fac.sep.TT.4K=topTable(fit4,coef=4,sort.by="B",resort.by="logFC",number=4000)
edge4.fac.sep.TT.6K=topTable(fit4,coef=4,sort.by="B",resort.by="logFC",number=6000)
edge4.fac.sep.TT.8K=topTable(fit4,coef=4,sort.by="B",resort.by="logFC",number=8000)
edge4.fac.sep.TT.12K=topTable(fit4,coef=4,sort.by="B",resort.by="logFC",number=12000)
edge4.fac.sep.TT.all=topTable(fit4,coef=4,sort.by="B",resort.by="logFC",number=Inf)

#volcanoplot
pdf(paste('fig_',q.i,'.',an.count,'.',figure,'_Volcanoplot_3BH_methods_CM3.pdf',sep=""),height=pdf.options()$height*0.7,width=pdf.options()$width*1.5)
par(mfrow=c(2,2))
volcanoplot(fit4,coef=1,highlight=30,col="grey",names=fit4$genes$ID)
volcanoplot(fit4,coef=2,highlight=30,col="grey",names=fit4$genes$ID)
volcanoplot(fit4,coef=3,highlight=30,col="grey",names=fit4$genes$ID)
volcanoplot(fit4,coef=4,highlight=30,col="grey",names=fit4$genes$ID)
par(mfrow=c(1,1))
dev.off()
figure<-figure+1

#HEATMAPS

#heatmap by edge
phenomatrix<-matrix(cbind(as.character(classcolours1[[i]]),as.character(classcolours2[[i]])),ncol=2)

pdf(paste('fig_',q.i,'.',an.count,'.',figure,'_heatmap_top500_edge1.pdf',sep=""),height=pdf.options()$height*1,width=pdf.options()$width*1)
superHeatmap2(x=data.q.fil[,edges[[1]]],y=rownames(edge1.fac.sep.TT.500),phenomatrix=phenomatrix[edges[[1]],],scale="row")
dev.off()
figure<-figure+1

pdf(paste('fig_',q.i,'.',an.count,'.',figure,'_heatmap_top500_edge2.pdf',sep=""),height=pdf.options()$height*1,width=pdf.options()$width*1)
superHeatmap2(x=data.q.fil[,edges[[2]]],y=rownames(edge2.fac.sep.TT.500),phenomatrix=phenomatrix[edges[[2]],],scale="row")
dev.off()
figure<-figure+1

pdf(paste('fig_',q.i,'.',an.count,'.',figure,'_heatmap_top500_edge3.pdf',sep=""),height=pdf.options()$height*1,width=pdf.options()$width*1)
superHeatmap2(x=data.q.fil[,edges[[3]]],y=rownames(edge3.fac.sep.TT.500),phenomatrix=phenomatrix[edges[[3]],],scale="row")
dev.off()
figure<-figure+1

pdf(paste('fig_',q.i,'.',an.count,'.',figure,'_heatmap_top500_edge4.pdf',sep=""),height=pdf.options()$height*1,width=pdf.options()$width*1)
superHeatmap2(x=data.q.fil[,edges[[4]]],y=rownames(edge4.fac.sep.TT.500),phenomatrix=phenomatrix[edges[[4]],],scale="row")
dev.off()
figure<-figure+1

#heatmaps for 4 signatures on all data
pdf(paste('fig_',q.i,'.',an.count,'.',figure,'_heatmap_top500_edge1_allData.pdf',sep=""),height=pdf.options()$height*1,width=pdf.options()$width*1)
phenomatrix<-matrix(cbind(as.character(classcolours1[[i]]),as.character(classcolours2[[i]])),ncol=2)
superHeatmap2(x=data.q.fil,y=rownames(edge1.fac.TT.500),phenomatrix=phenomatrix,scale="row")
dev.off()
figure<-figure+1

pdf(paste('fig_',q.i,'.',an.count,'.',figure,'_heatmap_top500_edge2_allData.pdf',sep=""),height=pdf.options()$height*1,width=pdf.options()$width*1)
superHeatmap2(x=data.q.fil,y=rownames(edge2.fac.TT.500),phenomatrix=phenomatrix,scale="row")
dev.off()
figure<-figure+1

pdf(paste('fig_',q.i,'.',an.count,'.',figure,'_heatmap_top500_edge3_allData.pdf',sep=""),height=pdf.options()$height*1,width=pdf.options()$width*1)
superHeatmap2(x=data.q.fil,y=rownames(edge3.fac.TT.500),phenomatrix=phenomatrix,scale="row")
dev.off()
figure<-figure+1

pdf(paste('fig_',q.i,'.',an.count,'.',figure,'_heatmap_top500_edge4_allData.pdf',sep=""),height=pdf.options()$height*1,width=pdf.options()$width*1)
superHeatmap2(x=data.q.fil,y=rownames(edge4.fac.TT.500),phenomatrix=phenomatrix,scale="row")
dev.off()
figure<-figure+1

pdf(paste('fig_',q.i,'.',an.count,'.',figure,'_heatmap_top500_edge5_allData.pdf',sep=""),height=pdf.options()$height*1.5,width=pdf.options()$width*1)
superHeatmap2(x=data.q.fil,y=rownames(edge5.flat.TT.500),phenomatrix=phenomatrix,scale="row")
dev.off()
figure<-figure+1

#PCA plots for top 500 probes for each of the 4 edges

#FIXREQ
pdf(paste('fig_',q.i,'.',an.count,'.',figure,'_PCA_top500_allEdges.pdf',sep=""),height=pdf.options()$height1.8,width=pdf.options()$width*1)
par(mfrow=c(3,2))
r500list<-lapply(list(edge1.fac.TT.500,edge2.fac.TT.500,edge3.fac.TT.500,edge4.fac.TT.500,edge5.flat.TT.500),rownames)
for (edge in 1:length(edges)){
  edgeset<-edges[[edge]]
  spca <- SamplePCA(exprs(data.q.fil.edges[[edge]][r500list[[edge]],]),classes[[edge]])
  plot(spca, col=legendcolourlist[[edge]][,2],main=paste("top 500"," probes for edge ",edge))#,cex=2,cex.axis=2,cex.main=3,cex.lab=2.5,pch=16)
  mtext(sprintf("%.2f",spca@variances[1]/sum(spca@variances)),side=1,line=3,adj=1,cex=1.5)
  mtext(sprintf("%.2f",spca@variances[2]/sum(spca@variances)),side=2,line=3,adj=1,cex=1.5)
  text(spca@scores[,1],spca@scores[,2],labels=row.names(spca@scores),col="lightgrey",cex=0.6)
  legend("bottomleft",
         legendcolourlist[[edge]][,1],
         pch=16,
         col=legendcolourlist[[edge]][,2],
         cex=.9
  )
  # mark the group centers
  for (type in 1:length(levels(classes[[edge]]))){
    x1 <- predict(spca, matrix(apply(t(exprs(data.q.fil.edges[[edge]][r500list[[edge]],])[,grep(levels(classes[[edge]])[[type]],eval(parse(text=paste("data.q.fil.edges[[edge]]$",colourvars[[edge]],sep=""))))]), 2, mean), ncol=1))
    points(x1[1], x1[2], col=legendcolourlist[[edge]][,2][[type]], cex=6,pch=9)}
}
dev.off()
figure<-figure+1



#Selection of probes to take forward to WGCNA. Several options are provided.
interestingProbeID_try<-union(
  rownames(edge1.fac.sep.TT.6K),
  rownames(edge3.fac.sep.TT.6K))

interestingProbeID_4K<-union(
  rownames(edge1.fac.sep.TT.4K),union(
    rownames(edge2.fac.sep.TT.4K),union(
      rownames(edge3.fac.sep.TT.4K),rownames(edge4.fac.sep.TT.4K))))

interestingProbeID_6K<-union(
  rownames(edge1.fac.sep.TT.6K),union(
    rownames(edge2.fac.sep.TT.6K),union(
      rownames(edge3.fac.sep.TT.6K),rownames(edge4.fac.sep.TT.6K))))

interestingProbeID_8K<-union(
  rownames(edge1.fac.sep.TT.8K),union(
    rownames(edge2.fac.sep.TT.8K),union(
      rownames(edge3.fac.sep.TT.8K),rownames(edge4.fac.sep.TT.8K))))

interestingProbeID_12K<-union(
  rownames(edge1.fac.sep.TT.12K),union(
    rownames(edge2.fac.sep.TT.12K),union(
      rownames(edge3.fac.sep.TT.12K),rownames(edge4.fac.sep.TT.12K))))


#variability filter: take only the probes where sd/mean > 0.02
p12k.1<-exprs(data.q.fil[,edges[[1]]])[rownames(edge1.fac.sep.TT.12K),]
tv12.1<-apply(p12k.1,1,function(xx){sd(xx)/mean(xx)}>0.02)
sum(tv12.1)

p12k.2<-exprs(data.q.fil[,edges[[2]]])[rownames(edge2.fac.sep.TT.12K),]
tv12.2<-apply(p12k.2,1,function(xx){sd(xx)/mean(xx)}>0.02)
sum(tv12.2)

p12k.3<-exprs(data.q.fil[,edges[[3]]])[rownames(edge3.fac.sep.TT.12K),]
tv12.3<-apply(p12k.3,1,function(xx){sd(xx)/mean(xx)}>0.02)
sum(tv12.3)

p12k.4<-exprs(data.q.fil[,edges[[4]]])[rownames(edge4.fac.sep.TT.12K),]
tv12.4<-apply(p12k.4,1,function(xx){sd(xx)/mean(xx)}>0.02)
sum(tv12.4)

interestingProbeID_12K_tv<-union(
  names(which(tv12.1==TRUE)),union(
    names(which(tv12.2==TRUE)),union(
      names(which(tv12.3==TRUE)),names(which(tv12.4==TRUE)))))

#export data
resultlist<-list(
  "edge1.sig"=edge1.fac.TT.sig,
  "edge2.sig"=edge2.fac.TT.sig,
  "edge3.sig"=edge3.fac.TT.sig,
  "edge4.sig"=edge4.fac.TT.sig,
  "edge5.sig"=edge5.flat.TT.sig,
  "edge1.500"=edge1.fac.TT.500,
  "edge2.500"=edge2.fac.TT.500,
  "edge3.500"=edge3.fac.TT.500,
  "edge4.500"=edge4.fac.TT.500,
  "edge5.500"=edge5.flat.TT.500,
  "edge1.12K"=edge1.fac.TT.12K,
  "edge2.12K"=edge2.fac.TT.12K,
  "edge3.12K"=edge3.fac.TT.12K,
  "edge4.12K"=edge4.fac.TT.12K,
  "edge5.12K"=edge5.flat.TT.12K
)

entrez.resultlist<-list()
for (result in 1:length(resultlist)){
  dict.entrez<-nuID2EntrezID(as.character(rownames(resultlist[[result]])),filterTh = NULL,lib.mapping='lumiHumanIDMapping', returnAllInfo = TRUE)
  newresult<-cbind(dict.entrez,resultlist[[result]])
  entrez.resultlist[[result]]<-newresult
  write.csv(newresult,file=file.path(dir.results,paste(analysis,dataset.names[[i]],names(resultlist)[result],".csv",sep="")))
}  
#NETWORK 0B: Move union of top4K (factorial, separate) DE genes for 4 edges to Neo4j, per edge####

##data for insertion (duplicate probemap code!!!!!)
# allfeatures<-rownames(fData(datalist[[1]]))
# annotSYM_all<-as.character(getSYMBOL(allfeatures,"lumiHumanAll.db"))#need new annotSYM that matches targetID annotation. targetID is full of aliases!!!
# #first name: annotSYM_all (these map to HGNC names)
# #second choice:if annotSYM_all entry is NA, then use "targetID" from LumiBatch object (better than nothing)
# probemap<-data.frame(cbind(allfeatures,fData(datalist[[1]])[,2],annotSYM_all),stringsAsFactors=F)
# replace<-as.character(probemap[which(is.na(probemap[,3])),2])
# faulty<-as.numeric(which(is.na(as.character(probemap[,3]))))
# probemap$annotSYM_all[faulty]<-replace
# colnames(probemap)<-c("nuID","targetID","SYMBOL") # here, the (HGNC) SYMBOL column has been 'fixed' by inserting illumina annotation information where SYMBOL is missing

if (neoinsert=="on"){
  #toneolist<-list(edge1.fac.sep.TT.all[interestingProbeID_4K,],edge2.fac.sep.TT.all[interestingProbeID_4K,],edge3.fac.sep.TT.all[interestingProbeID_4K,],edge4.fac.sep.TT.all[interestingProbeID_4K,],edge5.flat.TT.all[interestingProbeID_4K,])
  toneolist<-list(edge1.fac.sep.TT.all,edge2.fac.sep.TT.all,edge3.fac.sep.TT.all,edge4.fac.sep.TT.all,edge5.flat.TT.all)
  
  for (edge in 1:5){
    top4Kdata<-toneolist[[edge]]
    top4Kdata$nuID<-rownames(top4Kdata)
    top4Kdata.merge<-merge(top4Kdata,probemap,by.x=0,by.y="nuID",all.x=T,all.y=F)
    colnames(top4Kdata.merge)[length(colnames(top4Kdata.merge))]<-"SYMBOL2"
    
    nodes<-RNeo4j::nodes
    #network-specific cypher query
    write.csv(top4Kdata.merge,file.path(dir.figures,"tmp.csv"))
    
    query = paste("USING PERIODIC COMMIT 10000 LOAD CSV WITH HEADERS FROM 'file://",file.path(dir.figures,"tmp.csv"),"' AS csvLine 
                  CREATE (probe:PROBE {name:csvLine.nuID,square:'",dataset.variables[[1]],"',edge:toInt(",edge,"),contrastvar:'",contrastvars[[edge]],"',contrast:'",pathways[[edge]][[2]],"-",pathways[[edge]][[1]],"',logfc:toFloat(csvLine.logFC),aveEXPR:toFloat(csvLine.AveExpr),adjPVAL:toFloat(csvLine['adj.P.Val'])})
                  MERGE (gene:SYMBOL {name:csvLine.SYMBOL2})
                  CREATE (probe)-[:mapsTo]->(gene)",sep="")
    cypher(graph,query)
    queryclean<-"MATCH (s:SYMBOL {name:'NA'}) OPTIONAL MATCH (s)-[r]-() DELETE s, r"
    cypher(graph,queryclean)
    
    # blocks<-split(1:nrow(top4Kdata.merge),ceiling(seq_along(1:nrow(top4Kdata.merge))/1000))
    # for (block in blocks){
    # t = suppressMessages(newTransaction(graph))
    # graphdata=top4Kdata.merge[block,]
    # for (therow in 1:nrow(graphdata)) { #this can get extremely slow
    #   probe = as.character(graphdata[therow, ]$nuID)
    #   if (!is.null(graphdata[therow, ]$SYMBOL)){
    #   gene = as.character(graphdata[therow, ]$SYMBOL)
    #  } else if (is.na(graphdata[therow, ]$SYMBOL.y)){#this should not be necessary as the replacements were made in probemap
    #     gene = as.character(graphdata[therow, ]$targetID)#use secondary name derived from chip annotation rather than illumina.db lookup. this will be non HGNC in many cases
    #   } else {gene = as.character(graphdata[therow, ]$SYMBOL.y)}
    #   #endif      
    #   suppressMessages(appendCypher(t, 
    #                query, 
    #                nuID = probe, 
    #                targetID = gene,
    #                #invariant node properties
    #                squareg = dataset.variables[[1]],
    #                edgeg = edge,
    #                contrastvarsg = contrastvars[[edge]],
    #                contrastg = paste(pathways[[edge]][[2]],"-",pathways[[edge]][[1]]),
    #                logfcg = graphdata[therow, ]$logFC,
    #                aveEXPRg = graphdata[therow, ]$AveExpr,
    #                adjPVALg = graphdata[therow, ]$adj.P.Val
    #                
    #   ))
    # }
    # 
    # print("committing :)")
    # suppressMessages(commit(t))
    # }#end blocks
    #END NEO4J
    
  }#end edgewise neo4j loop
}#end neoinsert if

###
#Paired analysis of blood-fluid: output additional results for compartment####
#only apples to edges 1 and 2
# if (q.i==5){#compartment conditional####
#   
#   for (edge in c(1,2)){#start edgewise
#     edgeset<-edges[[edge]]
#     #subset data
#     edgedata<-data.q.fil[,edgeset]
#     #design matrix
#     #create person id from sampleIDs
#     persons<-factor(pData(datalist[[1]][,edgeset])$ID_KERRYN)
#     samples<-factor(pData(datalist[[1]][,edgeset])$sample_name)
#     compartment<-factor(pData(datalist[[1]][,edgeset])$Compartment)
#     #design.edge<-model.matrix(~persons+compartment)
#     design.edge<-model.matrix(~0+compartment)
#     colnames(design.edge)<-levels(compartment)
#     
#     #   colnames(design.edge)<-levels(factor(classifiers[[edge]]))
#     design.edge
#     
#     corfit <- duplicateCorrelation(edgedata,design.edge,block=persons)
#     corfit$consensus
#     
#     #lmFit
#     fit.edge<-lmFit(edgedata,design.edge,block=persons,correlation=corfit$consensus)
#     #contrast matrix
#     contrast.edge<-makeContrasts(
#       paste(pathways[[edge]][2],"-",pathways[[edge]][1]),
#       levels=design.edge
#     )
#     contrast.edge
#     #contrasts.fit
#     fit2.edge<-contrasts.fit(fit.edge,contrast.edge)
#     #eBayes
#     fit2.edge<-eBayes(fit2.edge)
#     fit2.edge$genes$ID<-featureData(edgedata)$TargetID
#     #decideTests
#     results.separate.edge<-decideTests(fit2.edge,method="separate")
#     #venn diagrams
#     vennDiagram(results.separate.edge,main=paste(names(edges)[edge],"(separate P)"),include=c("up","down"),counts.col=c("red","darkgreen"))
#     #topTable
#     edge.TT.sig.2=topTable(fit2.edge,coef=1,sort.by="B",resort.by="logFC",number=Inf,adjust.method="BH",p.value=0.05)#TT uses separate BH adjustment
#     dict.entrez<-nuID2EntrezID(as.character(rownames(edge.TT.sig.2)),filterTh = NULL,lib.mapping='lumiHumanIDMapping', returnAllInfo = TRUE)
#     newresult<-cbind(dict.entrez,edge.TT.sig.2)
#     write.csv(newresult,file=file.path(dir.results,paste(analysis,dataset.names[[i]],"_edge_",edge,"_paired_SIG.csv",sep="")))
#     
#     edge.TT.500.2=topTable(fit2.edge,coef=1,sort.by="B",resort.by="logFC",number=500)
#     dict.entrez<-nuID2EntrezID(as.character(rownames(edge.TT.500.2)),filterTh = NULL,lib.mapping='lumiHumanIDMapping', returnAllInfo = TRUE)
#     newresult<-cbind(dict.entrez,edge.TT.500.2)
#     write.csv(newresult,file=file.path(dir.results,paste(analysis,dataset.names[[i]],"_edge_",edge,"_paired_500.csv",sep="")))
#     
#     edge.TT.12K.2=topTable(fit2.edge,coef=1,sort.by="B",resort.by="logFC",number=12000)
#     dict.entrez<-nuID2EntrezID(as.character(rownames(edge.TT.12K.2)),filterTh = NULL,lib.mapping='lumiHumanIDMapping', returnAllInfo = TRUE)
#     newresult<-cbind(dict.entrez,edge.TT.12K.2)
#     write.csv(newresult,file=file.path(dir.results,paste(analysis,dataset.names[[i]],"_edge_",edge,"_paired_12K.csv",sep="")))
#     
#     #added this to compare to csde lists made later
#     edge.TT.1K.2=topTable(fit2.edge,coef=1,sort.by="B",resort.by="logFC",number=1000)
#     dict.entrez<-nuID2EntrezID(as.character(rownames(edge.TT.1K.2)),filterTh = NULL,lib.mapping='lumiHumanIDMapping', returnAllInfo = TRUE)
#     newresult<-cbind(dict.entrez,edge.TT.1K.2)
#     write.csv(newresult,file=file.path(dir.results,paste(analysis,dataset.names[[i]],"_edge_",edge,"_paired_1K.csv",sep="")))
#     
#     
#     
#     #volcanoplot
#     volcanoplot(fit2.edge,coef=1,highlight=30,col="grey",names=fit2.edge$genes$ID)
#     export.plot(file.prefix=paste('fig',q.i,'.',an.count,'.',i,'.',figure,names(questions)[[q.i]],'.',analysis,names(datalist)[i],"volc_paired","_edge",edge,sep=""),export.formats=export.formats.plots,height=1*height,width=(1 + sqrt(5))/1*height)
#     figure<-figure+1
#     
#     #heatmap
#     phenomatrix_c<-matrix(cbind(as.character(colourfactors[[edge]]),as.character(colourfactors[[edge]])),ncol=2)
#     superHeatmap2(x=edgedata,y=rownames(edge.TT.500.2),phenomatrix=phenomatrix_c,scale="row")
#     
#     export.plot(file.prefix=paste('fig',q.i,'.',an.count,'.',i,'.',figure,names(questions)[[q.i]],'.',analysis,names(datalist)[i],"heat_paired","_edge",edge,sep=""),export.formats=export.formats.plots,height=1*height,width=(1 + sqrt(5))/1*height)
#     figure<-figure+1
#     
#   }#end edgewise
#   
#   
# }#end compartment conditional

###



#4_Deconvolution****************************************####
analysis="4_Deconvolution"
an.count<-4
print("Analysis 2: Deconvolution and csDE")
#make relevant output folders
## Define folder for storing results NB do this for each analysis
dir.results <- file.path(dir.output.version,paste("Q_",q.i,"_",dataset.variables,sep=""),analysis, "results")
## Define folder for saving figures
dir.figures <- file.path(dir.output.version,paste("Q_",q.i,"_",dataset.variables,sep=""),analysis, "figures")
for (dir in c(dir.figures, dir.results)) {
  if (!file.exists(dir)) {
    dir.create(dir, recursive=T,showWarnings=T)
  }
}
setwd(dir.figures)
#setup

figure<-1
#VST already done not quantile normalised data
dtd<-as(data.v,"ExpressionSet")

data(Abbas)
data(HaemAtlas)

par(mar=c(6,6,6,12))
pheatmap(exprs(Abbas),scale="row",filename=paste('fig_',q.i,'.',an.count,'.',figure,'_359_Affy_probe_matrix_heatmap.pdf',sep=""),height=pdf.options()$height*1.5,width=pdf.options()$width*1.5)
par(parbackup)
figure<-figure+1

#make Marker lists
am<-MarkerList(Abbas)
am2<-MarkerList("HaemAtlas")
am2NU<-convertIDs(am2, 'lumiHumanAll.db', 'illuminaHumanv2.db')#nuID

AbbasRS<-convertIDs(am,'REFSEQ')
Abbassym<-convertIDs(am,'SYMBOL')
AbbasNU<-convertIDs(am, 'lumiHumanAll.db', 'hgu133a.db')#nuID
Abbasv3<-convertIDs(am, 'illuminaHumanv3.db', 'hgu133a.db',verbose=2,link='SYMBOL')
Abbasv4<-convertIDs(am, 'illuminaHumanv4.db', 'hgu133a.db',verbose=2)

#currently only one deconvolution method is implemented, but this can be extended
meth=list("Abbas"=gedBlood(dtd,verbose=T)#,
          #"SSKL" <- ged(dtd, am2NU, "ssFrobenius")#,#totally crashes the system
          #, log=FALSE,rng = 1234, nrun = 14 
          #, rng = 1234, nrun = 10
          #"deconf" <- ged(dtd, am2NU,"deconf")
)

#if more methods are used, this for loop needs to be used, otherwise h=1
#  for (h in 1:1){
h=1  
#1. deconvolution results####
decdat1<-meth[[h]]
colname<-contrast.variable
#colnames(coef(decdat1))<-as.character(TS)
decdat2<-as.matrix(coef(asCBC(decdat1,drop=FALSE)))
decmat<-as.matrix(coef(decdat1))
colnames(decmat)<-as.character(TS)
colnames(decdat2)<-as.character(TS)
decmat.red<-decmat[apply(decmat,1,sum)>0,]
write.csv(decmat.red,file.path(dir.results,"cellprop_table.csv"))
numbers<-as.vector(table(colnames(decmat.red)))
names<-names(table(colnames(decmat.red)))
nstring=paste("\nN=(",names[1]," ",numbers[1],",",names[2]," ",numbers[2],",",names[3]," ",numbers[3],",",names[4]," ",numbers[4],")",sep="")  

#stacked bar plot all types
pdf(paste('fig_',q.i,'.',an.count,'.',figure,'_Prop_matrix_all_stacked_barplot.pdf',sep=""),height=pdf.options()$height*1,width=pdf.options()$width*1)
par(parbackup)
par(mar=c(8,4.5,4.5,8))
barplot(decmat[,order(colnames(decmat))],las=2,col=brewer.pal(12,"Set3"),legend.text=rownames(coef(decdat1)),main=paste("Proportions of" ,nrow(coef(decdat1)),"cell types",names(meth)[[h]],names(datalist)[[i]],nstring),beside=F,ylab="Proportion",
        args.legend=c(x=ncol(decmat.red)*1.35,y=1,cex=.7),
        cex.names=.8)
par(parbackup)
dev.off()
figure<-figure+1

#stacked barplot detected types
pdf(paste('fig_',q.i,'.',an.count,'.',figure,'_Prop_matrix_detected_stacked_barplot.pdf',sep=""),height=pdf.options()$height*1,width=pdf.options()$width*1)
par(mar=c(8,4.5,4.5,8))
barplot(decmat.red[,order(colnames(decmat.red))],las=2,col=brewer.pal(12,"Set3"),legend.text=rownames(decmat.red),main=paste("Proportions of" ,nrow(decmat.red),"cell types",names(meth)[[h]],names(datalist)[[i]],nstring),beside=F,ylab="Proportion",args.legend=c(x=ncol(decmat.red)*1.35,y=1,cex=.7),cex.names=.8)
par(parbackup)
dev.off()
figure<-figure+1

#stacked barplot CBC
pdf(paste('fig_',q.i,'.',an.count,'.',figure,'_Prop_matrix_CBC_stacked_barplot.pdf',sep=""),height=pdf.options()$height*1,width=pdf.options()$width*1)
par(mar=c(8,4.5,4.5,8))
barplot(decdat2[,order(colnames(decdat2))],las=2,col=brewer.pal(12,"Set3"),legend.text=rownames(decdat2),main=paste("Proportions of CBC cell types",names(meth)[[h]],names(datalist)[[i]],nstring),beside=F,ylab="Proportion",args.legend=c(x=ncol(decdat2)*1.35,y=1,cex=.7),cex.names=.8)
par(parbackup)
dev.off()
figure<-figure+1

#stacked barplot detected types (PBMC only) 
decmat.red.PBMC<-matrix(data=NA,nrow=nrow(decmat.red)-1,ncol=ncol(decmat.red))
rownames(decmat.red.PBMC)<-rownames(decmat.red[-nrow(decmat.red),])
colnames(decmat.red.PBMC)<-colnames(decmat.red)
for (l in 1:ncol(decmat.red)){
  decmat.red.PBMC[,l]<-decmat.red[-nrow(decmat.red),l]/sum(decmat.red[-nrow(decmat.red),l])
}

pdf(paste('fig_',q.i,'.',an.count,'.',figure,'_Prop_matrix_PBMC_stacked_barplot.pdf',sep=""),height=pdf.options()$height*1,width=pdf.options()$width*1)
par(mar=c(8,4.5,4.5,8))
barplot(decmat.red.PBMC,las=2,col=brewer.pal(12,"Set3"),legend.text=rownames(decmat.red.PBMC),main=paste("Proportions of" ,nrow(decmat.red.PBMC),"PBMC cell types",names(meth)[[h]],names(datalist)[[i]],nstring),beside=F,ylab="Proportion",args.legend=c(x=ncol(decmat.red.PBMC)*1.35,y=1,cex=.7),cex.names=.8)
par(parbackup)
dev.off()
figure<-figure+1


##
#eigencell (embryonic code; the idea is generate the 1st PC of the cell-prop matrix)
cellmat1<-decmat.red
eigencell<-WGCNA::moduleEigengenes(t(decmat.red),rep("cellME",nrow(decmat.red)))

pdf(paste('fig_',q.i,'.',an.count,'.',figure,'_eigencell.pdf',sep=""),height=pdf.options()$height*1,width=pdf.options()$width*1.2)
par(mfrow=c(2,1), mar=c(0.3, 5.8, 10, 2.7))
plotMat(t(scale(t(decmat.red))),nrgcols=150,rlabels=rownames(decmat.red),clabels=colnames(decmat.red),rcols=1,ccols=as.character(classcolours1[[i]]))
#par(mar=c(4, 3.8, 0, 0.7))
par(mar=c(4.0, 4.5, 0.0, 1.3))
barplot(as.numeric(eigencell$eigengenes$MEcellME), col="orange", main="", cex.main=2,ylab="eigencell expression")
par(parbackup)
dev.off()
figure<-figure+1

pdf(paste('fig_',q.i,'.',an.count,'.',figure,'_ec_celltype_correlation.pdf',sep=""),height=pdf.options()$height*1.3,width=pdf.options()$width*1.2)
corres<-list()
par(mfrow=c(4,3),mar=c(4,4,3,3))
for (row in 1:nrow(decmat.red)){
  ec<-eigencell$eigengenes$MEcellME
  cp<-as.numeric(decmat.red[row,])
  pval<-cor.test(cp,ec,method = "pearson")$p.value
  est<-cor.test(cp,ec,method = "pearson")$estimate
  plot(cp~ec,col=as.character(classcolours1[[i]]),pch=16,main=paste(rownames(decmat.red)[row],"(Pval = ",sprintf("%.2f",pval)," Pearson R = ", sprintf("%.2f",est),")"),ylab="prop",xlab="ec")
  abline(lm(cp~ec), col="red")
  corres[[rownames(decmat.red)[row]]]<-pval
}
dev.off()
figure<-figure+1
##
#boxplot by condition and stats (edge 5)####

#convert and calculate####  
#clN<-flat
clN<-classifier1
clvN<-rep(clN,nrow(decmat.red))
dN<-t(decmat.red)
vN<-as.vector(dN)
typesN<-rownames(decmat.red)
typevN<-list()
for (j in 1:nrow(decmat.red)){typevN[[j]]=rep(typesN[[j]],ncol(decmat.red))}
typev2N<-unlist(typevN)
ndN<-data.frame(vN,clvN,typev2N)
flevs<-sort(unlist(pathway1))

#plot####
pdf(paste('fig_',q.i,'.',an.count,'.',figure,'_cellprop_boxplot_by_condition.pdf',sep=""),height=pdf.options()$height*1,width=pdf.options()$width*1.2)
par(mar=c(12.5,4.5,4.5,8))
boxplot(vN~clvN*typev2N,data=ndN,las=2,col=c("red","yellow","green","blue")[1:length(flevs)],main=paste("Relative abundance of cell types in",names(datalist)[[i]],names(meth)[[h]]),varwidth=FALSE)
legend("topleft",legend=(flevs),fill=c("red","yellow","green","blue")[1:length(flevs)])



par(parbackup)
dev.off()
figure<-figure+1
#stats for comparisons (edge 5)####
tests<-data.frame()
#for (the in list(c(1,2))){
statlist<-list()
ratiolist<-list()
for (k in 1:nrow(decmat.red)) {
  
  statlist[[k]]<-
    wilcox.test(decmat.red[k,flat==pathway1[[2]]],decmat.red[k,flat==pathway1[[1]]])$p.value
  
  ratiolist[[k]]<-
    median(decmat.red[k,flat==pathway1[[2]]])/median(decmat.red[k,flat==pathway1[[1]]])
  
}
tests<-rbind(tests,data.frame("type"="flat","names"=rownames(decmat.red),"ratio"=as.vector(unlist(ratiolist)),"pval"=as.vector(unlist(statlist)),"man_B"=as.vector(unlist(statlist))*nrow(decmat.red),"bonf"=p.adjust(as.vector(unlist(statlist)),method="bonferroni"),"BH"=p.adjust(as.vector(unlist(statlist)),method="BH")))
#}

tests<-tests[with(tests,order(names)),]
tests$edge<-rep(5,nrow(tests))
cellproptests.flat<-tests
write.csv(cellproptests.flat,file=file.path(dir.results,paste("Stats_",q.i,".",an.count,".",i,"by_",contrast.variable1,"_flat_wilcox.csv",sep="")))
#boxplot by condition and stats (4 edges)####
#convert and calculate#### 
clN<-colnames(decmat.red)
dN<-t(decmat.red)
vN<-as.vector(dN)
clvN<-rep(clN,nrow(decmat.red))
typesN<-rownames(decmat.red)
typevN<-list()
for (j in 1:nrow(decmat.red)){typevN[[j]]=rep(typesN[[j]],ncol(decmat.red))}
typev2N<-unlist(typevN)
ndN<-data.frame(vN,clvN,typev2N)

levs2<-levels(as.factor(colnames(decmat.red)))

levs<-character()
for (p1 in pathway1) {
  for (p2 in pathway2) {
    levs<-c(levs,paste(p1,p2,sep="."))
  }
}

#plot####
pdf(paste('fig_',q.i,'.',an.count,'.',figure,'_cellprop_boxplot_by_condition.pdf',sep=""),height=pdf.options()$height*1,width=pdf.options()$width*1.2)
par(mar=c(12.5,4.5,4.5,8))
boxplot(vN~clvN*typev2N,data=ndN,las=2,col=c("red","yellow","green","blue")[1:length(levs)],main=paste("Relative abundance of cell types in",names(datalist)[[i]],names(meth)[[h]]),varwidth=FALSE)
legend("topleft",legend=(levs2),fill=c("red","yellow","green","blue"))
par(parbackup)
dev.off()
figure<-figure+1
#stats for comparisons####
tests<-data.frame()
for (the in list(c(3,1),c(4,2),c(2,1),c(4,3))){
  statlist<-list()
  ratiolist<-list()
  print(the)
  print(levs[the[[1]]])
  print(levs[the[[2]]])
  for (k in 1:nrow(decmat.red)) {
    
    statlist[[k]]<-
      wilcox.test(decmat.red[k,colnames(decmat.red)==levs[the[[1]]]],decmat.red[k,colnames(decmat.red)==levs[the[[2]]]])$p.value
    
    ratiolist[[k]]<-
      median(decmat.red[k,colnames(decmat.red)==levs[the[[1]]]])/median(decmat.red[k,colnames(decmat.red)==levs[the[[2]]]])
    
  }
  tests<-rbind(tests,data.frame("type"=paste(as.character(the),collapse="."),"names"=rownames(decmat.red),"ratio"=as.vector(unlist(ratiolist)),"pval"=as.vector(unlist(statlist)),"man_B"=as.vector(unlist(statlist))*nrow(decmat.red),"bonf"=p.adjust(as.vector(unlist(statlist)),method="bonferroni"),"BH"=p.adjust(as.vector(unlist(statlist)),method="BH")))
}

tests<-tests[with(tests,order(names)),]
tests$edge<-rep(1:4,nrow(tests)/4)
cellproptests.edge<-tests
write.csv(cellproptests.edge,file=file.path(dir.results,paste("Stats_",q.i,".",an.count,".",i,"by_",contrast.variable1,"_edges_wilcox.csv",sep="")))

cellproptests<-rbind(cellproptests.edge,cellproptests.flat)
write.csv(cellproptests,file=file.path(dir.results,paste("Stats_",q.i,".",an.count,".",i,"by_",contrast.variable1,"_wilcox.csv",sep="")))  





###########PBMC##############    
#boxplot by condition and stats (edge 5:PBMC)####
#convert and calculate####  
clN<-flat
clvN<-rep(clN,nrow(decmat.red.PBMC))
dN<-t(decmat.red.PBMC)
vN<-as.vector(dN)
typesN<-rownames(decmat.red.PBMC)
typevN<-list()
for (j in 1:nrow(decmat.red.PBMC)){typevN[[j]]=rep(typesN[[j]],ncol(decmat.red.PBMC))}
typev2N<-unlist(typevN)
ndN<-data.frame(vN,clvN,typev2N)
flevs<-sort(unlist(pathway1))
#plot####
pdf(paste('fig_',q.i,'.',an.count,'.',figure,'_cellprop_boxplot_by_condition_PBMC.pdf',sep=""),height=pdf.options()$height*1,width=pdf.options()$width*1.2)
par(mar=c(12.5,4.5,4.5,8))
boxplot(vN~clvN*typev2N,data=ndN,las=2,col=c("red","yellow","green","blue")[1:length(flevs)],main=paste("Relative abundance of cell types in",names(datalist)[[i]],names(meth)[[h]]),varwidth=FALSE)
legend("topleft",legend=(flevs),fill=c("red","yellow","green","blue"))



par(parbackup)
dev.off()
figure<-figure+1
#stats for comparisons (edge 5)####
tests<-data.frame()
for (the in list(c(1,2))){
  statlist<-list()
  ratiolist<-list()
  for (k in 1:nrow(decmat.red.PBMC)) {
    
    statlist[[k]]<-
      wilcox.test(decmat.red.PBMC[k,flat==flevs[the[[2]]]],decmat.red.PBMC[k,flat==flevs[the[[1]]]])$p.value
    
    ratiolist[[k]]<-
      median(decmat.red.PBMC[k,flat==flevs[the[[2]]]])/median(decmat.red.PBMC[k,flat==flevs[the[[1]]]])
    
  }
  tests<-rbind(tests,data.frame("type"="flat","names"=rownames(decmat.red.PBMC),"ratio"=as.vector(unlist(ratiolist)),"pval"=as.vector(unlist(statlist)),"man_B"=as.vector(unlist(statlist))*nrow(decmat.red.PBMC),"bonf"=p.adjust(as.vector(unlist(statlist)),method="bonferroni"),"BH"=p.adjust(as.vector(unlist(statlist)),method="BH")))
}

tests<-tests[with(tests,order(names)),]
tests$edge<-rep(5,nrow(tests))
cellproptests.flat.PBMC<-tests
write.csv(cellproptests.flat.PBMC,file=file.path(dir.results,paste("Stats_",q.i,".",an.count,".",i,"by_",contrast.variable1,"_flat_PBMC_wilcox.csv",sep="")))
#boxplot by condition and stats (4 edges)####
#convert and calculate#### 
clN<-colnames(decmat.red.PBMC)
dN<-t(decmat.red.PBMC)
vN<-as.vector(dN)
clvN<-rep(clN,nrow(decmat.red.PBMC))
typesN<-rownames(decmat.red.PBMC)
typevN<-list()
for (j in 1:nrow(decmat.red.PBMC)){typevN[[j]]=rep(typesN[[j]],ncol(decmat.red.PBMC))}
typev2N<-unlist(typevN)
ndN<-data.frame(vN,clvN,typev2N)

levs2<-levels(as.factor(colnames(decmat.red.PBMC)))

levs<-character()
for (p1 in pathway1) {
  for (p2 in pathway2) {
    levs<-c(levs,paste(p1,p2,sep="."))
  }
}
#plot####
pdf(paste('fig_',q.i,'.',an.count,'.',figure,'_cellprop_boxplot_by_condition_PBMC.pdf',sep=""),height=pdf.options()$height*1,width=pdf.options()$width*1.2)
par(mar=c(12.5,4.5,4.5,8))
boxplot(vN~clvN*typev2N,data=ndN,las=2,col=c("red","yellow","green","blue")[1:length(levs)],main=paste("Relative abundance of cell types in",names(datalist)[[i]],names(meth)[[h]]),varwidth=FALSE)
legend("topleft",legend=(levs2),fill=c("red","yellow","green","blue"))
par(parbackup)
dev.off()
figure<-figure+1
#stats for comparisons####
tests<-data.frame()
for (the in list(c(3,1),c(4,2),c(2,1),c(4,3))){
  statlist<-list()
  ratiolist<-list()
  #print(the)
  #print(levs[the[[1]]])
  #print(levs[the[[2]]])
  for (k in 1:nrow(decmat.red.PBMC)) {
    
    statlist[[k]]<-
      wilcox.test(decmat.red.PBMC[k,colnames(decmat.red.PBMC)==levs[the[[1]]]],decmat.red.PBMC[k,colnames(decmat.red.PBMC)==levs[the[[2]]]])$p.value
    
    ratiolist[[k]]<-
      median(decmat.red.PBMC[k,colnames(decmat.red.PBMC)==levs[the[[1]]]])/median(decmat.red.PBMC[k,colnames(decmat.red.PBMC)==levs[the[[2]]]])
    
  }
  tests<-rbind(tests,data.frame("type"=paste(as.character(the),collapse="."),"names"=rownames(decmat.red.PBMC),"ratio"=as.vector(unlist(ratiolist)),"pval"=as.vector(unlist(statlist)),"man_B"=as.vector(unlist(statlist))*nrow(decmat.red.PBMC),"bonf"=p.adjust(as.vector(unlist(statlist)),method="bonferroni"),"BH"=p.adjust(as.vector(unlist(statlist)),method="BH")))
}

tests<-tests[with(tests,order(names)),]
tests$edge<-rep(1:4,nrow(tests)/4)
cellproptests.edge.PBMC<-tests
write.csv(cellproptests.edge.PBMC,file=file.path(dir.results,paste("Stats_",q.i,".",an.count,".",i,"by_",contrast.variable1,"_edges_PBMC_wilcox.csv",sep="")))

cellproptests.PBMC<-rbind(cellproptests.edge.PBMC,cellproptests.flat.PBMC)
write.csv(cellproptests.PBMC,file=file.path(dir.results,paste("Stats_",q.i,".",an.count,".",i,"by_",contrast.variable1,"_PBMC_wilcox.csv",sep="")))  

############END PBMC########


#2. csde####
# 
# #PC added to B cells
# th<-apply(coef(decdat1)[1:2,],2,sum)
# tc<-apply(coef(decdat1)[3:4,],2,sum)
# b<-apply(coef(decdat1)[5:10,],2,sum)
# #pc<-coef(decdat1)[10,]
# nk<-apply(coef(decdat1)[11:12,],2,sum)
# mo<-apply(coef(decdat1)[13:14,],2,sum)
# dc<-apply(coef(decdat1)[15:16,],2,sum)
# neut<-coef(decdat1)[17,]
# 
# bigprops<-matrix(
#   rbind(th,tc,b,nk,mo,dc,neut    
#   ),nrow=7)
# rownames(bigprops)<-c("th","tc","b","nk","mo","dc","neut")
# colnames(bigprops)<-colnames(coef(decdat1))
# 
# #Smallprops
# data(Abbas)
# l<-apply(coef(decdat1)[1:10,],2,sum)
# nk<-apply(coef(decdat1)[11:12,],2,sum)
# apc<-apply(coef(decdat1)[13:16,],2,sum)
# neut<-coef(decdat1)[17,]
# 
# smallprops<-matrix(
#   rbind(l,nk,apc,neut   
#   ),nrow=4)
# rownames(smallprops)<-c("lym","nk","apc","neut")
# colnames(smallprops)<-colnames(coef(decdat1))
# 
# #edgewise csDE
# for (edge in 1:length(edges)){
# subset<-order(esApply(dtd,1,sd,na.rm=T))[1:5000]
# edgeset<-dtd[subset,edges[[edge]]]
# if(min(table(classifiers[[edge]]))<nrow(bigprops)) {
# props<-smallprops
# print("smallprops used")} else {
#   props<-bigprops
#   print("bigprops used")}
# 
# #ggplot and csplot don't play nice anymore.
# # #make a "default" plot as the csDEplot may fail
# # cs.de<-data.frame(c(1,2,3),c(1,2,3))
# # gp<-ggplot(data=cs.de,aes(x=cs.de[,1],y=cs.de[,2]))+geom_point()+ggtitle("fail")
# # 
# # #this may fail
# try(
#   cs.de<-ged(
#     edgeset,
#     props[,edges[[edge]]],
#     data=as.factor(classifiers[[edge]]),
#     verbose=TRUE,nperm=1000))
# par(parbackup);
# # try(gp<-csplot(cs.de));
# # 
# # pdf(paste('fig_',q.i,'.',an.count,'.',figure,'_edge_',edge,'_Cell-specific FDR plots.pdf',sep=""),height=pdf.options()$height*1.5,width=pdf.options()$width*1.5)
# # print(gp)
# # dev.off()
# # figure<-figure+1
# 
# tt<-"no csDE result"
# try(tt<-csTopTable(cs.de,decreasing=F,n=1000,sort.by="FDR"))
# if(tt!="no csDE result"){
# dict.th<-data.frame(cbind(names(tt$th),as.numeric(tt$th),as.character(nuID2targetID(as.character(names(tt$th)),species="Human"))));try(colnames(dict.th)<-c("nuID","FDR","SYMBOL"))
# write.csv(dict.th,file=file.path(dir.results,paste("csDE_Th_",q.i,'.',an.count,'.',i,"_",names(datalist)[[i]],names(meth)[[h]],"_","_edge_",edge,".csv",sep="")))
# dict.tc<-data.frame(cbind(names(tt$tc),as.numeric(tt$tc),as.character(nuID2targetID(as.character(names(tt$tc)),species="Human"))));try(colnames(dict.tc)<-c("nuID","FDR","SYMBOL"))
# write.csv(dict.tc,file=file.path(dir.results,paste("csDE_Tc_",q.i,'.',an.count,'.',i,"_",names(datalist)[[i]],names(meth)[[h]],"_","_edge_",edge,".csv",sep="")))
# dict.b<-data.frame(cbind(names(tt$b),as.numeric(tt$b),as.character(nuID2targetID(as.character(names(tt$b)),species="Human"))));try(colnames(dict.b)<-c("nuID","FDR","SYMBOL"))
# write.csv(dict.b,file=file.path(dir.results,paste("csDE_B_",q.i,'.',an.count,'.',i,"_",names(datalist)[[i]],names(meth)[[h]],"_","_edge_",edge,".csv",sep="")))
# dict.nk<-data.frame(cbind(names(tt$nk),as.numeric(tt$nk),as.character(nuID2targetID(as.character(names(tt$nk)),species="Human"))));try(colnames(dict.nk)<-c("nuID","FDR","SYMBOL"))
# write.csv(dict.nk,file=file.path(dir.results,paste("csDE_NK_",q.i,'.',an.count,'.',i,"_",names(datalist)[[i]],names(meth)[[h]],"_","_edge_",edge,".csv",sep="")))
# dict.mo<-data.frame(cbind(names(tt$mo),as.numeric(tt$mo),as.character(nuID2targetID(as.character(names(tt$mo)),species="Human"))));try(colnames(dict.mo)<-c("nuID","FDR","SYMBOL"))
# write.csv(dict.mo,file=file.path(dir.results,paste("csDE_MO_",q.i,'.',an.count,'.',i,"_",names(datalist)[[i]],names(meth)[[h]],"_","_edge_",edge,".csv",sep="")))
# dict.dc<-data.frame(cbind(names(tt$dc),as.numeric(tt$dc),as.character(nuID2targetID(as.character(names(tt$dc)),species="Human"))));try(colnames(dict.dc)<-c("nuID","FDR","SYMBOL"))
# write.csv(dict.dc,file=file.path(dir.results,paste("csDE_DC_",q.i,'.',an.count,'.',i,"_",names(datalist)[[i]],names(meth)[[h]],"_","_edge_",edge,".csv",sep="")))
# dict.neut<-data.frame(cbind(names(tt$neut),as.numeric(tt$neut),as.character(nuID2targetID(as.character(names(tt$neut)),species="Human"))));try(colnames(dict.neut)<-c("nuID","FDR","SYMBOL"))
# write.csv(dict.neut,file=file.path(dir.results,paste("csDE_neut_",q.i,'.',an.count,'.',i,"_",names(datalist)[[i]],names(meth)[[h]],"_","_edge_",edge,".csv",sep="")))
# dict.lym<-data.frame(cbind(names(tt$lym),as.numeric(tt$lym),as.character(nuID2targetID(as.character(names(tt$lym)),species="Human"))));try(colnames(dict.lym)<-c("nuID","FDR","SYMBOL"))
# write.csv(dict.lym,file=file.path(dir.results,paste("csDE_lym_",q.i,'.',an.count,'.',i,"_",names(datalist)[[i]],names(meth)[[h]],"_","_edge_",edge,".csv",sep="")))
# dict.apc<-data.frame(cbind(names(tt$apc),as.numeric(tt$apc),as.character(nuID2targetID(as.character(names(tt$apc)),species="Human"))));try(colnames(dict.apc)<-c("nuID","FDR","SYMBOL"))
# write.csv(dict.apc,file=file.path(dir.results,paste("csDE_apc_",q.i,'.',an.count,'.',i,"_",names(datalist)[[i]],names(meth)[[h]],"_","_edge_",edge,".csv",sep="")))} else print(tt)
# 
# if(tt!="no csDE result"){ 
# #pathway reactomePA
# for (celltype in 1:length(tt)){
# rgenesENT<-as.character(nuID2EntrezID(names(tt[[celltype]]),lib.mapping='lumiHumanIDMapping'))
# rgenesE2<-rgenesENT[which(rgenesENT!="")]
# 
# rpa<-NA
# rpa <- try(enrichPathway(unique(rgenesE2), organism="human",pvalueCutoff = 0.4,readable = T))
# #head(summary(rpa))
# summary(rpa)
# 
# #catch error conditions (NULL and NA content fpr rpa) and export result. This is subtle, and distinguishes between NULL and NA conditions
# #edit: added another error condition: rpa may succeed, but summary(rpa) is empty, I have no idea why
# 
# succ<-TRUE
# if (length(rpa)==0){succ<-FALSE;print(paste(names(tt)[celltype],"rpa is of zero length"))} else if (is.na(rpa)){succ<-FALSE;print(paste(names(tt)[celltype],"rpa is NA"))}else if (class(rpa)=="try-error"){succ<-FALSE;print(paste(names(tt)[celltype],"rpa gives a try error"))} else if (nrow(summary(rpa))==0){succ<-FALSE;print(paste(names(tt)[celltype],"summary(rpa) is of zero length"))} else if (nrow(summary(rpa))==1){succ<-TRUE;print(paste(names(tt)[celltype],"nrow summary rpa is of length 1"))}#edited last condition
# 
# if (succ) {
# write.csv(summary(rpa),file=file.path(dir.results,paste('table',q.i,'.',an.count,'.',i,"_",names(questions)[[q.i]],'_',analysis,names(datalist)[i],"Cell_",names(tt[celltype]),"_edge_",edge,"ReactomePA.csv",sep="")))
# if(nrow(summary(rpa))>1){
#   rpaplot<-barplot(rpa,showCategory=8)
#   
#   pdf(paste('fig_',q.i,'.',an.count,'.',figure,'_edge_',edge,names(tt[celltype]),"_cell_",'_rpaPlot.pdf',sep=""),height=pdf.options()$height*1,width=pdf.options()$width*1)      
#   rpaplot<-barplot(rpa,showCategory=8)
#   print(rpaplot)
#   dev.off()
#   figure<-figure+1
#   
# }
# #put iconv here
# #the code below fixes an error due to an invlid multibyte string in some data subsets. Apparently, the database used by rpa is encoded in latin-9, and cnetplot does not recognise all characters for all possible labels. Therefore, we convert the encoding for the labels to UTF-8 so that things don't break
# 
# theOld<-rpa@result$Description
# theNew<-iconv(theOld,from="LATIN-9",to="UTF-8")
# rpa2<-rpa
# rpa2@result$Description<-theNew
# 
# pdf(paste('fig_',q.i,'.',an.count,'.',figure,'_edge_',edge,names(tt[celltype]),"_cell_",'_cnetPlot.pdf',sep=""),height=pdf.options()$height*1.8,width=pdf.options()$width*1.8)    
# cnetplot(rpa2, categorySize = "pvalue")
# dev.off()
# figure<-figure+1
# 
# } else {print(paste("no result for rpa: ","Cell_",names(tt[celltype])))}
# 
# }#end for loop celltype
# 
# #get EntrezIDs on a per-cell basis
# cellList<-list()
# for (cell in 1:length(tt)){
# rgenesENT<-as.character(nuID2EntrezID(names(tt[[cell]]),lib.mapping='lumiHumanIDMapping'))
# cellList[[cell]]<-rgenesENT[which(rgenesENT!="")]
# }
# names(cellList)<-names(tt)
# res<-"fail"
# try(res <- compareCluster(cellList, fun = "enrichPathway"))
# if (class(res)!="character"){
# resplot<-plot(res)
# resplotb<-resplot+theme(axis.text.x = element_text(angle = 270, hjust = 0))
# pdf(paste('fig_',q.i,'.',an.count,'.',figure,'_edge_',edge,'_compareCluster_plot.pdf',sep=""),height=pdf.options()$height*1,width=pdf.options()$width*1)    
# print(resplotb)
# dev.off()
# figure<-figure+1}
# 
# #KEGG pathway enrichment
# res2<-"fail"
# try(res2 <- compareCluster(cellList, fun = "enrichKEGG"))
# if (class(res2)!="character"){
# if(nrow(summary(res2))>4){
# resplot2<-plot(res2)
# resplot2b<-resplot2+theme(axis.text.x = element_text(angle = 270, hjust = 0))
# 
#  pdf(paste('fig_',q.i,'.',an.count,'.',figure,'_edge_',edge,'_KEGG_enrich.pdf',sep=""),height=pdf.options()$height*1,width=pdf.options()$width*1)
#  print(resplot2b)
# dev.off()
# figure<-figure+1
# }#end inner if
# }#end if
# 
# 
# }#end conditional pathway analysis
# 
# }#end edgewise csde

#5_WGCNA********************************************####
analysis="5_WGCNA_4k_tv"
an.count<-5
print("Analysis 5: WGCNA")
#1. Setup####
#make relevant output folders
## Define folder for storing results NB do this for each analysis
dir.results <- file.path(dir.output.version,paste("Q_",q.i,"_",dataset.variables,sep=""),analysis, "results")
## Define folder for saving figures
dir.figures <- file.path(dir.output.version,paste("Q_",q.i,"_",dataset.variables,sep=""),analysis, "figures")

for (dir in c(dir.figures, dir.results)) {
  if (!file.exists(dir)) {
    dir.create(dir, recursive=T,showWarnings=T)
  }
}
setwd(dir.figures)
figure<-1
#make user list for later
data("HaemAtlas")
ha.2<-HaemAtlas[,c("SYMBOL","Type")]
colnames(ha.2)<-c("GENES","LISTS")
write.csv(ha.2,file.path(dir.results,"haemAtlas_old.csv"),row.names=F)
write.csv(ha.2,file.path(dir.results,"haemAtlas.csv"),row.names=F)#debug

if(neoinsert=="on"){
  ##new
  data(HaemAtlas)
  ha2<-HaemAtlas
  hanu<-HaemAtlas$nuID
  news<-nuID2targetID(hanu)
  ha2$news<-news
  ha2<-unique(ha2)
  
  # hanu2<-paste(hanu,collapse="','")
  # #need to add in platform specification
  # q<-paste("MATCH (p:PROBETYPE)-[r]-(s:SYMBOL) WHERE p.name IN ['",hanu2,"'] RETURN p.name AS nuID2, s.name AS SYMBOL2")
  # try(res<-cypher(graph,q))
  # ha3<-merge(ha2,res,by.x="nuID",by.y="nuID2",all.x=FALSE)
  # ha4<-ha3[,c("nuID","Type")]
  
  ha4<-ha2[,c("nuID","Type")]
  colnames(ha4)<-c("GENES","LISTS")
  ha4<-unique(ha4[order(ha4$LISTS),])
  write.csv(ha4,file.path(dir.results,"haemAtlas.csv"),row.names=F)
  ##\new
}

abbasvals<-cbind(as.character(marknames(Abbassym)),names(marknames(Abbassym)))
colnames(abbasvals)<-c("GENES","LISTS")
write.csv(abbasvals,file.path(dir.results,"abbasvals.csv"),row.names=F)

#2. Expression data prep####

# #batch correction
# ##removeBatchEffect
# batch=NULL
# if("study"%in%colnames(pData(datalist[[i]]))){
#   batch=factor(pData(datalist[[i]])$study)
#   if(length(levels(batch))>1){
#     exprs(data.q.fil)<-removeBatchEffect(data.q.fil,batch)
#     pdf(paste('fig_',q.i,'.',an.count,'.',figure,'_MDS_q_fil_BC_samples_by_',contrast.variable1,'.pdf',sep=""),height=pdf.options()$height*1.1,width=pdf.options()$width*1.1)
#     plotSampleRelationsAD(data.q.fil,method="mds",color=as.character(classcolours1[[i]]),plotchar=16)
#     legend("bottomleft",
#            legendcolours1[,1],
#            pch=16,
#            col=legendcolours1[,2],
#            cex=.9
#     )
#     dev.off()
#     figure<-figure+1
#   }#end test for length of study i.e number of batches
# }#end test for existence of study

#select data for WGCNA 1: use top x DE genes
subset<-is.element(featureNames(data.q.fil),interestingProbeID_4K)
tv0<-featureNames(data.q.fil)[subset]
subset2<-!subset
remain<-data.q.fil
datremain<-data.frame(t(exprs(remain)))

#look at variability
tvmads<-apply(datremain,2,function(xx){mad(xx)})
hist(tvmads)
datremain2<-datremain[,order(tvmads,decreasing=TRUE)]
tv2<-colnames(datremain2[1:5000])

tvx<-which(featureNames(data.q.fil)%in%unique(c(tv0,tv2)))
#usedata<-data.q.fil[subset|bloodGeneOverlap|immuneGeneOverlap,]
usedata<-data.q.fil[tvx,]

#if(wgcna_override==TRUE){usedata<-data.q.fil}

datExpr0 = data.frame(t(exprs(usedata)))#there is a weird bug here; about 20 percent of nuIDs get an X prepended, changing their name...
names(datExpr0)<-featureNames(usedata)

# Take a quick look at what is in the data set:
paste("Dimensions of data set:",nrow(datExpr0),"samples and",ncol(datExpr0),"probes")
#3. Data QC####
#check for missing values
gsg=goodSamplesGenes(datExpr0,verbose=3);
gsg$allOK
#remove offending samples
if (!gsg$allOK)
{
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes)>0)
    printFlush(paste("Removing genes:", paste(names(datExpr0)[!gsg$goodGenes], collapse = ", ")));
  if (sum(!gsg$goodSamples)>0)
    printFlush(paste("Removing samples:", paste(rownames(datExpr0)[!gsg$goodSamples], collapse = ", ")));
  # Remove the offending genes and samples from the data:
  datExpr0 = datExpr0[gsg$goodSamples, gsg$goodGenes]
}

#cluster samples to look for outliers
sampleTree = flashClust(dist(datExpr0), method = "average");

# Plot the sample tree
pdf(paste('fig_',q.i,'.',an.count,'.',figure,'_Flash_clust_samples.pdf',sep=""),height=pdf.options()$height*1,width=pdf.options()$width*1)
par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2)
# Plot a line to show the cut
abline(h = 80, col = "red")
par(parbackup)
dev.off()
figure<-figure+1

# Determine cluster under the line
clust = cutreeStatic(sampleTree, cutHeight = 80, minSize = 10)
table(clust)
# clust 1 contains the samples we want to keep.
keepSamples = (clust==1)
datExpr = datExpr0[keepSamples, ]
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)

#removing outliers breaks the code further down. For now, don't remove outliers at all
datExpr=datExpr0

#Now exclude probes that are less variable
tv<-apply(datExpr,2,function(xx){sd(xx)/mean(xx)}>0.01)
datExpr=datExpr[,tv]

#if(wgcna_override==TRUE){datExpr=datExpr0}

print(paste("EXPR data before tv:",dim(datExpr0)))
print(paste("EXPR data after tv:",dim(datExpr)))
#4. Prepare clinical trait data#####

cat<-data.frame(apply(pd3cat,2,function(x){replace(x,x=="",NA)}))
catnum<-lapply(cat,function(x){as.numeric(as.factor(x))})

num<-apply(pd3num[,1:(ncol(pd3num)-2)],2,as.numeric)

datTraits<-cbind(data.frame(catnum),num)
if(!is.null(pd2$P1.sampleID)){
  rownames(datTraits)<-pd2$P1.sampleID
}else if (!is.null(pd2$sample_name)) {
  rownames(datTraits)<-pd2$sample_name
}else if (!is.null(pd2$otherID)) {
  rownames(datTraits)<-pd2$otherID
}else{rownames(datTraits)<-pd2$sampleID} 

collectGarbage();
datTraits <- datTraits[,colSums(is.na(datTraits))<nrow(datTraits)]#this removes all columns where all values are NA
uniquelength <- apply(datTraits,2,function(x) length(unique(x[!is.na(x)])))

datTraits <- subset(datTraits, select=uniquelength>1)#this removes all columns where all values are the same after removing NAs
#following now commented out. Can't lose important information this way!
#if(reducePheno==TRUE){datTraits<-datTraits[,which(apply(datTraits,2,function(x){sum(is.na(x))/nrow(datTraits)<.6}))]}
datTraits<-data.frame(datTraits)

write.csv(datTraits,file=file.path(dir.results,"WGCNA_TRAIT_DATA.csv"))

#new: remove highly correlated traits

#debug
#datTraits<-read.csv("~/output/build/version_2016-08-10_17_15_48/Q_5_square4/5_WGCNA_4k_tv/results/WGCNA_TRAIT_DATA.csv",row.names=1)
#datTraits<-read.csv("~/output/build/version_2016-08-10_17_15_48/Q_10_rec1/5_WGCNA_4k_tv/results/WGCNA_TRAIT_DATA.csv",row.names=1)


# #just for display purposes; modify datTraits later #fails occasionally due to NAs?
# selectvector<-apply(datTraits,2,function(x){!sum(is.na(x))/length(x)>.25})
# vars2<-datTraits[,selectvector]
# pdf(paste('fig_',q.i,'.',an.count,'.',figure,'_mix_heatmap_reduced.pdf',sep=""),height=pdf.options()$height*1,width=pdf.options()$width*1)
# mix.heatmap(vars2, rowmar=10, legend.mat=F)
# dev.off()
# figure<-figure+1
# 
# pdf(paste('fig_',q.i,'.',an.count,'.',figure,'_distmap_variables.pdf',sep=""),height=pdf.options()$height*1,width=pdf.options()$width*1)
# distmap(vars2, what="variables", margins=c(6,6))
# dev.off()
# figure<-figure+1
# 
# pdf(paste('fig_',q.i,'.',an.count,'.',figure,'_distmap_subjects.pdf',sep=""),height=pdf.options()$height*1,width=pdf.options()$width*1)
# distmap(vars2, what="subjects", margins=c(6,6))
# dev.off()
# figure<-figure+1
# 
# pdf(paste('fig_',q.i,'.',an.count,'.',figure,'_confounder_plot.pdf',sep=""),height=pdf.options()$height*1,width=pdf.options()$width*1)
# #confounderPlot(vars2,x="Compartment",y="HIV.Status")
# confounderPlot(vars2,x=contrast.variable1,y=contrast.variable2)
# dev.off()
# figure<-figure+1

#correlate traits

tmp <- cor(datTraits,method="spearman",use="p")
#pheatmap(tmp)
pheatmap(tmp,filename=paste('fig_',q.i,'.',an.count,'.',figure,'_spearman_correlated_traits.pdf',sep=""))
tmp[upper.tri(tmp)] <- 0
diag(tmp) <- 0
selectvector<-!apply(tmp,2,function(x) any(abs(x) > 0.85))#debug. empirical value!
selectvector[contrast.variable1]<-TRUE
selectvector[contrast.variable2]<-TRUE

if(reducePheno){
  
  selectvector[which(is.na(selectvector))]<-TRUE
  datTraits2 <- datTraits[,selectvector]
  head(datTraits2)
  tmp <- cor(datTraits,method="spearman",use="p")
  #pheatmap(tmp)
  pheatmap(tmp,filename=paste('fig_',q.i,'.',an.count,'.',figure,'_spearman_correlated_traits_reduced.pdf',sep=""))
  
  #this needs to be fixed
  # #check with CluMix plots
  # selectvector2<-apply(datTraits,2,function(x){!sum(is.na(x))/length(x)>.25})
  # vars3<-datTraits[,selectvector2]
  # pdf(paste('fig_',q.i,'.',an.count,'.',figure,'_mix_heatmap_reduced.pdf',sep=""),height=pdf.options()$height*1,width=pdf.options()$width*1)
  # mix.heatmap(vars3, rowmar=10, legend.mat=F)
  # dev.off()
  # figure<-figure+1
  # 
  # pdf(paste('fig_',q.i,'.',an.count,'.',figure,'_distmap_variables_reduced.pdf',sep=""),height=pdf.options()$height*1,width=pdf.options()$width*1)
  # distmap(vars3, what="variables", margins=c(6,6))
  # dev.off()
  # figure<-figure+1
  # 
  # pdf(paste('fig_',q.i,'.',an.count,'.',figure,'_distmap_subjects_reduced.pdf',sep=""),height=pdf.options()$height*1,width=pdf.options()$width*1)
  # distmap(vars3, what="subjects", margins=c(6,6))
  # dev.off()
  # figure<-figure+1
  # 
  # pdf(paste('fig_',q.i,'.',an.count,'.',figure,'_confounder_plot_reduced.pdf',sep=""),height=pdf.options()$height*1,width=pdf.options()$width*1)
  # #confounderPlot(vars3,x="Compartment",y="HIV.Status")
  # confounderPlot(vars3,x=contrast.variable1,y=contrast.variable2)
  # dev.off()
  # figure<-figure+1
  
  write.csv(datTraits2,file=file.path(dir.results,"WGCNA_TRAIT_DATA_REDUCED.csv"))
  
  ##end new datTraits
}#end if

#expression data and phenotype data are now in analagous data frames

#correlation matrix of traits
pheatmap(cor(as.matrix(datTraits),use='pairwise.complete.obs'),filename=paste('fig_',q.i,'.',an.count,'.',figure,'_Correlation of phenotype variables.pdf',sep=""),height=8,width=8)
figure<-figure+1
#todo:send to neo4j after making edgewise

#correlation matrix of traits with cellprop
cellprops<-t(decmat.red)
names(cellprops)<-rownames(decmat.red)
pheatmap(cor(as.matrix(datTraits),cellprops,use='pairwise.complete.obs'),filename=paste('fig_',q.i,'.',an.count,'.',figure,'_Correlation of phenotype variables_with cellprop.pdf',sep=""),height=8,width=8)
figure<-figure+1
#todo:send to neo4j after making edgewise
#5. Visualisation of how the traits relate to the sample dendrogram####
# Re-cluster samples
sampleTree2 = flashClust(dist(datExpr), method = "average")

# Convert traits to a color representation: white means low, red means high, grey means missing entry
traitColors = numbers2colors(as.matrix(datTraits), signed = FALSE);

# Plot the sample dendrogram and the colors underneath.
pdf(paste('fig_',q.i,'.',an.count,'.',figure,'_Dendrogram_and_pheno.pdf',sep=""),height=pdf.options()$height*1,width=pdf.options()$width*1)
plotDendroAndColors(sampleTree2, traitColors,cex.dendroLabels = 0.5,cex.colorLabels=.5,cex.lab=.7,cex.axis=.7,groupLabels = colnames(datTraits),main = "Sample dendrogram and trait heatmap")
par(parbackup)
dev.off()
figure<-figure+1
#6. Choose the soft-thresholding power: analysis of network topology####
# Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to=20, by=2))
# Call the network topology analysis function
sft<-"fail"
sft <- try(pickSoftThreshold(datExpr, powerVector = powers, RsquaredCut = 0.80, verbose = 5))
sft <- try(pickSoftThreshold(datExpr, powerVector = powers, RsquaredCut = 0.80, verbose = 5))
if(class(sft)!="list"){
  disableWGCNAThreads()
  sft = try(pickSoftThreshold(datExpr, powerVector = powers, RsquaredCut = 0.80, verbose = 5))
}


print(paste("Soft-thresholding power:",sft$powerEstimate ))

# Plot the results:
pdf(paste('fig_',q.i,'.',an.count,'.',figure,'_scale_independence_and_mean_connectivity.pdf',sep=""),height=pdf.options()$height*1,width=pdf.options()$width*1)
par(mfrow = c(1,2));par(mar=c(4,4,2,2))
cex1 = 0.9;
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.80,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
par(parbackup)
dev.off()
figure<-figure+1
#7. Construct network and detect modules####
thispower = sft$powerEstimate
deepSplitParameter<-2

#for really small samples avoid overcalling modules
if(nrow(datExpr0)<20){deepSplitParameter=1}

if (is.na(thispower)){
  thispower=min(sft$fitIndices[sft$fitIndices$SFT.R.sq>0.75,]$Power)
  print("No power estimate for R squared of .80, using .75")
  print(paste("New soft-thresholding power:",thispower ))}

if(thispower<4){
  if(nSamples<20){thispower<-10}else if(nSamples<30){thispower<-9}else if(nSamples<40){thispower<-8}else if(nSamples<60){thispower<-9}else if(nSamples>=60){thispower<-6}
}

#if(wgcna_override==TRUE){
#  thispower=wgcnaSTP
#  deepSplitParameter=wgcnaDS
#}

#Module detection (blockwise algorithm, but this is irrelevant as max block size is set to really large number)
bwnet = blockwiseModules(datExpr, maxBlockSize = 44000,#edited
                         power = thispower,#networkType='signed',TOMType="signed",
                         minModuleSize = 20,#was 20
                         deepSplit=deepSplitParameter,#added, default is 2
                         reassignThreshold = 0, mergeCutHeight = 0.25,
                         numericLabels = FALSE,#was true, but I need colours, so false
                         saveTOMs = TRUE,
                         saveTOMFileBase = file.path(dir.results,"TOM_blockwise"),#change filename based on no probes
                         verbose = 3)

bwmoduleColors = bwnet$colors#this is only kept for q.i=1 and i=1, all subsequent ones are relabeled to match these
modsizes<-sort(table(bwmoduleColors)[which(names(table(bwmoduleColors))!="grey")],decreasing=TRUE)
print("modsizes")
print(modsizes)


# Plot the dendrogram and the module colors underneath
pdf(paste('fig_',q.i,'.',an.count,'.',figure,'_scale_independence_and_mean_connectivity.pdf',sep=""),height=pdf.options()$height*1,width=pdf.options()$width*1)
plotDendroAndColors(bwnet$dendrograms[[1]], bwmoduleColors[bwnet$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideAll=FALSE,guideCount=50,guideHang = 0.05,
                    main=paste(names(datalist[[1]]),nrow(datExpr),"samples",ncol(datExpr),"probes")
)
par(parbackup)
dev.off()
figure<-figure+1

if(length(bwnet$MEs)<10){
  thispower = thispower
  deepSplitParameter<-3
  
  #for really small samples avoid overcalling modules
  if(nrow(datExpr0)<20){deepSplitParameter=1}
  
  if (is.na(thispower)){
    thispower=min(sft$fitIndices[sft$fitIndices$SFT.R.sq>0.75,]$Power)
    print("No power estimate for R squared of .80, using .75")
    print(paste("New soft-thresholding power:",thispower ))}
  
  #Module detection (blockwise algorithm, but this is irrelevant as max block size is set to really large number)
  bwnet = blockwiseModules(datExpr, maxBlockSize = 44000,#edited
                           power = thispower,#networkType='signed',TOMType="signed",
                           minModuleSize = 20,#was 20
                           deepSplit=deepSplitParameter,#added, default is 2
                           reassignThreshold = 0, mergeCutHeight = 0.25,
                           numericLabels = FALSE,#was true, but I need colours, so false
                           saveTOMs = TRUE,
                           saveTOMFileBase = file.path(dir.results,"TOM_blockwise"),#change filename based on no probes
                           verbose = 3)
  
  bwmoduleColors = bwnet$colors#this is only kept for q.i=1 and i=1, all subsequent ones are relabeled to match these
  
  # Plot the dendrogram and the module colors underneath
  pdf(paste('fig_',q.i,'.',an.count,'.',figure,'_scale_independence_and_mean_connectivity.pdf',sep=""),height=pdf.options()$height*1,width=pdf.options()$width*1)
  plotDendroAndColors(bwnet$dendrograms[[1]], bwmoduleColors[bwnet$blockGenes[[1]]],
                      "Module colors",
                      dendroLabels = FALSE, hang = 0.03,
                      addGuide = TRUE, guideAll=FALSE,guideCount=50,guideHang = 0.05,
                      main=paste(names(datalist[[1]]),nrow(datExpr),"samples",ncol(datExpr),"probes")
  )
  par(parbackup)
  dev.off()
  figure<-figure+1
  
}

if(length(bwnet$MEs)<10){
  thispower = thispower
  deepSplitParameter<-4
  
  #for really small samples avoid overcalling modules
  if(nrow(datExpr0)<20){deepSplitParameter=1}
  
  if (is.na(thispower)){
    thispower=min(sft$fitIndices[sft$fitIndices$SFT.R.sq>0.75,]$Power)
    print("No power estimate for R squared of .80, using .75")
    print(paste("New soft-thresholding power:",thispower ))}
  
  #Module detection (blockwise algorithm, but this is irrelevant as max block size is set to really large number)
  bwnet = blockwiseModules(datExpr, maxBlockSize = 44000,#edited
                           power = thispower,#networkType='signed',TOMType="signed",
                           minModuleSize = 20,#was 20
                           deepSplit=deepSplitParameter,#added, default is 2
                           reassignThreshold = 0, mergeCutHeight = 0.25,
                           numericLabels = FALSE,#was true, but I need colours, so false
                           saveTOMs = TRUE,
                           saveTOMFileBase = file.path(dir.results,"TOM_blockwise"),#change filename based on no probes
                           verbose = 3)
  
  bwmoduleColors = bwnet$colors#this is only kept for q.i=1 and i=1, all subsequent ones are relabeled to match these
  
  # Plot the dendrogram and the module colors underneath
  pdf(paste('fig_',q.i,'.',an.count,'.',figure,'_scale_independence_and_mean_connectivity.pdf',sep=""),height=pdf.options()$height*1,width=pdf.options()$width*1)
  plotDendroAndColors(bwnet$dendrograms[[1]], bwmoduleColors[bwnet$blockGenes[[1]]],
                      "Module colors",
                      dendroLabels = FALSE, hang = 0.03,
                      addGuide = TRUE, guideAll=FALSE,guideCount=50,guideHang = 0.05,
                      main=paste(names(datalist[[1]]),nrow(datExpr),"samples",ncol(datExpr),"probes")
  )
  par(parbackup)
  dev.off()
  figure<-figure+1
  
}


#save the tree (module assignment and module eigengene information)

bwMEs = bwnet$MEs;
bwgeneTree = bwnet$dendrograms[[1]];

#calculate KMEs
sigKMEs<-signedKME(datExpr,bwMEs)

# ###NEW
# dim(datExpr)
# dimnames(datExpr)[[1]]
# 
# printFlush("Calculating connectivities...");
# # Human network:
# AdjMatHuman = abs(cor(datExpr ,use="p"))^thispower
# diag(AdjMatHuman)=0
# 
# ###END new

#load modules
load(file.path(dir.output.version,"allmodules.RData"))
#save new modules but keep in memory for now
allmodules[[paste("mDat",names(questions)[q.i],sep="")]]<-list("data"<-datExpr,"colours"<-bwmoduleColors)
save(allmodules,file=file.path(dir.output.version,"allmodules.RData"))

#Compare modules to reference module set (square1)
genes1<-colnames(allmodules[[1]][[1]])
mods1<-allmodules[[1]][[2]]
gm1<-as.matrix(cbind(genes1,mods1))
modNames1<-unique(mods1)
mlist1<-list()
for(name1 in modNames1){
  mlist1[[name1]]<-gm1[gm1[,2]==name1,1]
}

genes2<-colnames(allmodules[[paste("mDat",names(questions)[q.i],sep="")]][[1]])
mods2<-allmodules[[paste("mDat",names(questions)[q.i],sep="")]][[2]]
gm2<-as.matrix(cbind(genes2,mods2))
modNames2<-unique(mods2)
mlist2<-list()
for(name2 in modNames2){
  mlist2[[name2]]<-gm2[gm2[,2]==name2,1]
}

overmatrix<-matrix(data=NA,nrow=length(mlist1),ncol=length(mlist2))
dimnames(overmatrix)<-list(names(mlist1),names(mlist2))
for (m1 in 1:length(mlist1)){
  for (m2 in 1:length(mlist2)){
    #denom<-length(unique(c(mlist1[[m1]],mlist2[[m2]])))
    denom<-length(mlist1[[m1]])
    numer<-sum(mlist1[[m1]]%in%mlist2[[m2]])
    res<-numer/denom
    overmatrix[m1,m2]<-res
  }
}

pdf(paste('fig_',q.i,'.',an.count,'.',figure,'_overlap_with_reference_modules.pdf',sep=""),height=pdf.options()$height*1,width=pdf.options()$width*1)
heatmap.2(overmatrix,col=blueWhiteRed(11),scale="none",main="Overlap of modules with reference modules")
par(parbackup)
dev.off()
figure<-figure+1

rm(allmodules)
#8. Edgewise eigengenes####
# Define numbers of genes and samples
nGenes = ncol(datExpr);
nSamples = nrow(datExpr);

#define module eigengenes and order
MEs0=bwMEs
MEs = orderMEs(MEs0)

# DON'T Recalculate MEs with color labels for each edge, just subset them
ME.edges=list()
for (edge in 1:length(edges)){
  edgeset<-edges[[edge]]
  MEs0.edge = MEs[edgeset,]
  MEs.edge = orderMEs(MEs0.edge)
  ME.edges[[names(edges)[edge]]]<-MEs.edge
}
ME_export<-MEs
row.names(ME_export)<-row.names(datExpr)
ME_combo<-merge(ME_export,datTraits,by.x="row.names",by.y="row.names")
write.csv(ME_export,file=file.path(dir.results,"ModuleEigengenes_all_samples.csv"))

write.csv(ME_combo,file=file.path(dir.results,"ModuleEigengenes_all_samples_withPheno.csv"))

#8.b eigengene correlation with eigencell - which modules drive eigencell; or are particularly influenced by cell type composition?
#depends on MEs for current question; run in wgcna section. not edgewise, as we are looking for overall drivers
sampleids<-colnames(coef(decdat1))
MES2<-ME_export[order(row.names(ME_export)),]
ec2<-ec[order(sampleids)]

pdf(paste('fig_',q.i,'.',an.count,'.',figure,'_eigencell_cor_with_MEs.pdf',sep=""),height=pdf.options()$height*1,width=pdf.options()$width*1)
corME<-list()
par(mfrow=c(5,3),mar=c(4,4,3,3))
for (col in 1:ncol(MES2)){
  me<-as.numeric(MES2[,col])
  pval<-cor.test(me,ec2,method = "pearson")$p.value
  est<-cor.test(me,ec2,method = "pearson")$estimate
  plot(me~ec2,col=as.character(classcolours1[[i]])[order(sampleids)],pch=16,main=paste(colnames(MES2)[col],"(Pval = ",sprintf("%.2f",pval)," Pearson R = ", sprintf("%.2f",est),")"),ylab="ME",xlab="ec2")
  abline(lm(me~ec2), col="red")
  corME[[rownames(decmat.red)[row]]]<-pval
}
par(parbackup)
dev.off()
figure<-figure+1
#9. Heatmap of all module eigengenes for all samples####

if (ncol(MEs)>0){#this should prevent errors if no modules found
  allME<-t(as.matrix(MEs))
  colnames(allME)<-colnames(exprs(datalist[[i]]))
  MEcol<-unlist(strsplit(rownames(allME),"ME"))[seq(2,nrow(allME)*2,2)]
  hm<-heatmap_ad(allME,keep.dendro=T,col=blueWhiteRed(100),ColSideColors=phenomatrix,RowSideColors=matrix(c(MEcol,MEcol),ncol=2),noan=2,margins=c(6,7),main=paste("ME heatmap",names(datalist)[[i]]))
  
  pdf(paste('fig_',q.i,'.',an.count,'.',figure,'_heatmap_all_ME.pdf',sep=""),height=pdf.options()$height*1,width=pdf.options()$width*1) 
  heatmap_ad(allME,keep.dendro=T,col=blueWhiteRed(100),ColSideColors=phenomatrix,RowSideColors=matrix(c(MEcol,MEcol),ncol=2),noan=2,margins=c(6,7),main=paste("ME heatmap",names(datalist)[[i]]))
  dev.off()
  figure<-figure+1
  den0<-unlist(hm$Rowv)
  den<-den0[length(den0):1]
  par(parbackup)
}#end if

#all modules
for(edge in 1:length(edges)){
  #module subsets
  par(parbackup)
  allME.raw<-t(ME.edges[[edge]])#transpose necessary
  allME<-allME.raw[which(rownames(allME.raw)!="MEgrey"),]
  colnames(allME)<-colnames(exprs(datalist[[i]]))[edges[[edge]]]
  #boxplot(allME,las=2)
  
  MEcol<-unlist(strsplit(rownames(allME),"ME"))[seq(2,nrow(allME)*2,2)]
  
  #hm_unclust<-heatmap_ad(allME,Rowv=NA,Colv=NA,col=blueWhiteRed(100),ColSideColors=phenomatrix[edges[[edge]],],RowSideColors=matrix(c(MEcol,MEcol),ncol=2),noan=2,margins=c(6,7),main=paste("ME heatmap","edge",edge,names(datalist)[[i]]))
  
  pdf(paste('fig_',q.i,'.',an.count,'.',figure,'_heatmap_ME_edge_',edge,'.pdf',sep=""),height=pdf.options()$height*1,width=pdf.options()$width*1) 
  heatmap_ad(allME,keep.dendro=T,col=blueWhiteRed(100),ColSideColors=phenomatrix[edges[[edge]],],RowSideColors=matrix(c(MEcol,MEcol),ncol=2),noan=2,margins=c(6,7),main=paste("ME heatmap","edge",edge,names(datalist)[[i]]),scale="none")
  dev.off()
  figure<-figure+1
}
#10. NETWORK 1:Quantifying module trait associations####
#TODO: implement 2 additional methods of correlation
useMultiplePheno=FALSE
if(useMultiplePheno==TRUE){}#

MEs0=bwMEs
MEs = orderMEs(MEs0)
#NETWORK 1 raw data
#Correlate MEs and traits (all edges)
moduleTraitCor = cor(MEs, datTraits, use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);
moduleTraitQvalue = matrix(p.adjust(moduleTraitPvalue,method="BH"),nrow=nrow(moduleTraitPvalue),dimnames=list(row.names(moduleTraitPvalue),colnames(moduleTraitPvalue)))
#edges
write.csv(moduleTraitCor,file=file.path(dir.results,paste('table',q.i,'.',an.count,'.',i,'.',names(questions)[[q.i]],'.',analysis,names(datalist)[i],"ModuleTraitCor.csv",sep="")))
#edge annotation
write.csv(moduleTraitQvalue,file=file.path(dir.results,paste('table',q.i,'.',an.count,'.',i,'.',names(questions)[[q.i]],'.',analysis,names(datalist)[i],"ModuleTraitCor_qvalue.csv",sep="")))

#Display correlation matrix

# Will display correlations and their BH-adjusted p-values
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitQvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)

pdf(paste('fig_',q.i,'.',an.count,'.',figure,'_ME_pheno_correlation.pdf',sep=""),,height=nrow(moduleTraitCor)/2,width=ncol(moduleTraitCor)*0.75) 
par(mar = c(8,10,2,3));par(mfrow=c(1,1))
# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor[den,],
               xLabels = names(datTraits),
               yLabels = names(MEs)[den],
               ySymbols = names(MEs)[den],
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix[den,],
               setStdMargins = FALSE,
               cex.text = 0.9,
               zlim = c(-1,1),
               main = paste("Module-pheno correlation"))
par(parbackup)
dev.off()
figure<-figure+1

# clustered version, with NA columns removed
pdf(paste('fig_',q.i,'.',an.count,'.',figure,'_ME_pheno_correlation_clust.pdf',sep=""),height=pdf.options()$height*1,width=pdf.options()$width*1) 
par(cex.main=.7)
heatmap(moduleTraitCor[,!apply(moduleTraitCor,2,is.na)[1,]],col=blueWhiteRed(50),main="Clustered module-trait relationships",mar=c(8,1),cexRow=.8,cexCol=.8)
dev.off()
par(parbackup)
figure<-figure+1

#NETWORK 1 raw data edgewise
for (edge in 1:length(edges)){
  nSamples=length(edges[[edge]])
  #Correlate MEs and traits
  moduleTraitCor = cor(ME.edges[[edge]], datTraits[edges[[edge]],], use = "p");
  retain<-!apply(moduleTraitCor,2,is.na)[1,]
  if(sum(retain)>2){
    moduleTraitCor = moduleTraitCor[,retain] #need to now remove NA cols!
    moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);
    moduleTraitQvalue = matrix(p.adjust(moduleTraitPvalue,method="BH"),nrow=nrow(moduleTraitPvalue),dimnames=list(row.names(moduleTraitPvalue),colnames(moduleTraitPvalue)))
    #edges
    write.csv(moduleTraitCor,file=file.path(dir.results,paste('table',q.i,'.',an.count,'.',i,'.',names(edges)[edge],".",names(questions)[[q.i]],'.',analysis,names(datalist)[i],"ModuleTraitCor_edge",edge,".csv",sep="")))
    #edge annotation
    write.csv(moduleTraitQvalue,file=file.path(dir.results,paste('table',q.i,'.',an.count,'.',i,names(edges)[edge],".",'.',names(questions)[[q.i]],'.',analysis,names(datalist)[i],"ModuleTraitCor_qvalue_edge",edge,".csv",sep="")))
    
    #Display correlation matrix
    # Will display correlations and their BH-adjusted p-values
    textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                       signif(moduleTraitQvalue, 1), ")", sep = "");
    dim(textMatrix) = dim(moduleTraitCor)
    pdf(paste('fig_',q.i,'.',an.count,'.',figure,'_ME_pheno_correlation_edge_',edge,'.pdf',sep=""),height=nrow(moduleTraitCor)/2,width=ncol(moduleTraitCor)*0.75) 
    par(mar = c(8,10,2,3));par(mfrow=c(1,1))
    if(ncol(moduleTraitCor)<6){par(mar = c(4,5,1,2));par(mfrow=c(1,1))}
    
    # Display the correlation values within a heatmap plot
    labeledHeatmap(Matrix = moduleTraitCor,
                   xLabels = names(datTraits)[retain],
                   yLabels = colnames(ME.edges[[edge]]),
                   ySymbols = colnames(ME.edges[[edge]]),
                   colorLabels = FALSE,
                   colors = blueWhiteRed(50),
                   textMatrix = textMatrix,
                   setStdMargins = FALSE,
                   cex.text = 0.9,
                   zlim = c(-1,1),
                   main = paste("ME-pheno correlation",names(edges)[edge]))
    dev.off()
    par(parbackup)
    figure<-figure+1
    
    # clustered version, with NA columns removed
    pdf(paste('fig_',q.i,'.',an.count,'.',figure,'_ME_pheno_correlation_clust_',edge,'.pdf',sep=""),height=pdf.options()$height*1,width=pdf.options()$width*1)   
    par(mar = c(2,2,2,5));par(cex.main=.5);
    heatmap(moduleTraitCor[,!apply(moduleTraitCor,2,is.na)[1,]],col=blueWhiteRed(50),main=paste("Clustered module-trait relationships",names(edges)[edge]),mar=c(8,1),cexRow=.8,cexCol=.8,scale="none")
    dev.off()
    par(parbackup)
    figure<-figure+1
  }else{print("matrix is missing")}#end if
}#end edgewise NETWORK 1
#11. NETWORK 2: Quantifying module cell-proportion associations####

#NETWORK 2 raw data
# Define numbers of genes and samples
nGenes = ncol(datExpr);
nSamples = nrow(datExpr);
# Recalculate MEs with color labels
#MEs0 = moduleEigengenes(datExpr, bwmoduleColors)$eigengenes#not required

cellprops<-t(decmat.red)
names(cellprops)<-rownames(decmat.red)
MEs0=bwMEs
MEs = orderMEs(MEs0)
moduleCellCor = cor(MEs, cellprops, use = "p");
retain<-!apply(moduleCellCor,2,is.na)[1,]
moduleCellCor = moduleCellCor[,retain] #NA columns removed!

moduleCellPvalue = corPvalueStudent(moduleCellCor, nSamples);
moduleCellQvalue = matrix(p.adjust(moduleCellPvalue,method="BH"),nrow=nrow(moduleCellPvalue),dimnames=list(row.names(moduleCellPvalue),colnames(moduleCellPvalue)))
#edges
write.csv(moduleCellCor,file=file.path(dir.results,paste('table',q.i,'.',an.count,'.',i,'.',names(questions)[[q.i]],'.',analysis,names(datalist)[i],"ModuleCellCor.csv",sep="")))
#edge annotation
write.csv(moduleCellQvalue,file=file.path(dir.results,paste('table',q.i,'.',an.count,'.',i,'.',names(questions)[[q.i]],'.',analysis,names(datalist)[i],"ModuleCellCor_qvalue.csv",sep="")))
#Will display correlations and their BH adjusted q-values
textMatrix = paste(signif(moduleCellCor, 2), "\n(",
                   signif(moduleCellQvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleCellCor)

pdf(paste('fig_',q.i,'.',an.count,'.',figure,'_ME_cellprop_correlation.pdf',sep=""),height=nrow(moduleCellCor)/2,width=ncol(moduleCellCor)*0.75) 
par(mar = c(6,10,3,2));par(mfrow=c(1,1))
# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleCellCor[den,],
               xLabels = names(cellprops),
               yLabels = names(MEs)[den],
               ySymbols = names(MEs)[den],
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix[den,],
               setStdMargins = FALSE,
               cex.text = 0.9,
               zlim = c(-1,1),
               main = paste("ME-celltype correlation"))
par(parbackup)
dev.off()
figure<-figure+1

# clustered version, with NA columns removed
pdf(paste('fig_',q.i,'.',an.count,'.',figure,'_ME_cellprop_correlation_clust.pdf',sep=""),height=pdf.options()$height*1,width=pdf.options()$width*1) 
par(mar = c(2,2,2,5));par(cex.main=.5);
heatmap(moduleCellCor[,!apply(moduleCellCor,2,is.na)[1,]],col=blueWhiteRed(50),main="Clustered module-cell relationships",mar=c(5,1),cexRow=.8,cexCol=.8)
par(parbackup)
dev.off()
figure<-figure+1

#edgewise module cell associations

for (edge in 1:length(edges)){
  nSamples=length(edges[[edge]])
  #Correlate MEs and cellprops
  moduleCellCor = cor(ME.edges[[edge]], cellprops[edges[[edge]],], use = "p");
  retain<-!apply(moduleCellCor,2,is.na)[1,]
  moduleCellCor = moduleCellCor[,retain] #NA columns removed!
  moduleCellPvalue = corPvalueStudent(moduleCellCor, nSamples);
  moduleCellQvalue = matrix(p.adjust(moduleCellPvalue,method="BH"),nrow=nrow(moduleCellPvalue),dimnames=list(row.names(moduleCellPvalue),colnames(moduleCellPvalue)))
  #edges
  write.csv(moduleCellCor,file=file.path(dir.results,paste('table',q.i,'.',an.count,'.',i,'.',names(edges)[edge],".",names(questions)[[q.i]],'.',analysis,names(datalist)[i],"ModuleCellCor_edge",edge,".csv",sep="")))
  #edge annotation
  write.csv(moduleCellQvalue,file=file.path(dir.results,paste('table',q.i,'.',an.count,'.',i,names(edges)[edge],".",'.',names(questions)[[q.i]],'.',analysis,names(datalist)[i],"ModuleCellCor_qvalue_edge",edge,".csv",sep="")))
  
  #Display correlation matrix
  # Will display correlations and their BH-adjusted p-values
  textMatrix = paste(signif(moduleCellCor, 2), "\n(",
                     signif(moduleCellQvalue, 1), ")", sep = "");
  dim(textMatrix) = dim(moduleCellCor)
  
  pdf(paste('fig_',q.i,'.',an.count,'.',figure,'_ME_cellprop_correlation_edge',edge,'.pdf',sep=""),height=nrow(moduleCellCor)/2,width=ncol(moduleCellCor)*0.75) 
  par(mar = c(8,10,2,3));par(mfrow=c(1,1))
  # Display the correlation values within a heatmap plot
  labeledHeatmap(Matrix = moduleCellCor,
                 xLabels = colnames(cellprops)[retain],
                 yLabels = names(ME.edges[[edge]]),
                 ySymbols = names(ME.edges[[edge]]),
                 colorLabels = FALSE,
                 colors = blueWhiteRed(50),
                 textMatrix = textMatrix,
                 setStdMargins = FALSE,
                 cex.text = 0.9,
                 zlim = c(-1,1),
                 main = paste("ME-celltype correlation",names(edges)[edge]))
  par(parbackup)
  figure<-figure+1
  dev.off()
  
  # clustered version, with NA columns removed
  
  par(mar = c(2,2,2,5));par(cex.main=.5);
  pdf(paste('fig_',q.i,'.',an.count,'.',figure,'_ME_cellprop_correlation_edge_clust',edge,'.pdf',sep=""),height=pdf.options()$height*1,width=pdf.options()$width*1) 
  heatmap(moduleCellCor[,!apply(moduleCellCor,2,is.na)[1,]],col=blueWhiteRed(50),main=paste("Clustered module-Cell relationships",names(edges)[edge]),mar=c(5,1),cexRow=.8,cexCol=.8)
  par(parbackup)
  dev.off()
  figure<-figure+1
}
#11.b NETWORK 2b: Quantifying module flow cell-proportion associations####
if (dataset.variables[[i]]=="berry.test"){
  #NETWORK 2 raw data
  # Define numbers of genes and samples
  nGenes = ncol(datExpr);
  nSamples = nrow(datExpr);
  # Recalculate MEs with color labels
  #MEs0 = moduleEigengenes(datExpr, bwmoduleColors)$eigengenes#not required
  
  cellprops2_raw<-read.csv("/home/rstudio/source_data/Flow_props_test.csv")
  idmapper<-pd2[,c(4,2,7,12)]
  cellprops2_2<-merge(idmapper,cellprops2_raw,by.x="site_donor_id",by.y="site_donor_id",all.x=TRUE,all.y=TRUE)
  
  #remove 
  #reorder
  cellprops2_3<-cellprops2_2[which(!is.na(cellprops2_2$sample_name)),]
  cellprops2_3$comboname<-paste(cellprops2_3$class,cellprops2_3$sex,sep=".")
  cellprops2<-cellprops2_3[match(pd2$sample_name,cellprops2_3$sample_name),]
  cellprops2m<-as.matrix(cellprops2[,5:(ncol(cellprops2)-1)])
  dimnames(cellprops2m)[[1]]<-cellprops2$comboname
  
  MEs0=bwMEs
  MEs = orderMEs(MEs0)
  moduleCellCor = cor(MEs, cellprops2m, use = "p");
  retain<-!apply(moduleCellCor,2,is.na)[1,]
  moduleCellCor = moduleCellCor[,retain] #NA columns removed!
  
  moduleCellPvalue = corPvalueStudent(moduleCellCor, nSamples);
  moduleCellQvalue = matrix(p.adjust(moduleCellPvalue,method="BH"),nrow=nrow(moduleCellPvalue),dimnames=list(row.names(moduleCellPvalue),colnames(moduleCellPvalue)))
  #edges
  write.csv(moduleCellCor,file=file.path(dir.results,paste('table',q.i,'.',an.count,'.',i,'.',names(questions)[[q.i]],'.',analysis,names(datalist)[i],"ModuleFlowCellCor.csv",sep="")))
  #edge annotation
  write.csv(moduleCellQvalue,file=file.path(dir.results,paste('table',q.i,'.',an.count,'.',i,'.',names(questions)[[q.i]],'.',analysis,names(datalist)[i],"ModuleFlowCellCor_qvalue.csv",sep="")))
  #Will display correlations and their BH adjusted q-values
  textMatrix = paste(signif(moduleCellCor, 2), "\n(",
                     signif(moduleCellQvalue, 1), ")", sep = "");
  dim(textMatrix) = dim(moduleCellCor)
  
  pdf(paste('fig_',q.i,'.',an.count,'.',figure,'_ME_flowcellprop_correlation.pdf',sep=""),height=nrow(moduleCellCor)/2,width=ncol(moduleCellCor)*0.75) 
  par(mar = c(6,10,3,2));par(mfrow=c(1,1))
  # Display the correlation values within a heatmap plot
  labeledHeatmap(Matrix = moduleCellCor[den,],
                 xLabels = colnames(cellprops2m),
                 yLabels = names(MEs)[den],
                 ySymbols = names(MEs)[den],
                 colorLabels = FALSE,
                 colors = blueWhiteRed(50),
                 textMatrix = textMatrix[den,],
                 setStdMargins = FALSE,
                 cex.text = 0.9,
                 zlim = c(-1,1),
                 main = paste("ME-celltype correlation"))
  par(parbackup)
  dev.off()
  figure<-figure+1
  
  # clustered version, with NA columns removed
  pdf(paste('fig_',q.i,'.',an.count,'.',figure,'_ME_flowcellprop_correlation_clust.pdf',sep=""),height=pdf.options()$height*1,width=pdf.options()$width*1) 
  par(mar = c(2,2,2,5));par(cex.main=.5);
  heatmap(moduleCellCor[,!apply(moduleCellCor,2,is.na)[1,]],col=blueWhiteRed(50),main="Clustered module-cell relationships",mar=c(5,1),cexRow=.8,cexCol=.8)
  par(parbackup)
  dev.off()
  figure<-figure+1
  
  #edgewise module cell associations
  
  for (edge in 1:length(edges)){
    nSamples=length(edges[[edge]])
    
    #Correlate MEs and cellprops2m
    narow<-!apply(cellprops2m[edges[[edge]],],2,is.na)[,1]
    datapoints<-nrow(cellprops2m[edges[[edge]],][narow,])
    moduleCellCor = cor(ME.edges[[edge]], cellprops2m[edges[[edge]],], use = "p");
    retain<-!apply(moduleCellCor,2,is.na)[1,]
    moduleCellCor = moduleCellCor[,retain] #NA columns removed!
    moduleCellPvalue = corPvalueStudent(moduleCellCor, nSamples);
    moduleCellQvalue = matrix(p.adjust(moduleCellPvalue,method="BH"),nrow=nrow(moduleCellPvalue),dimnames=list(row.names(moduleCellPvalue),colnames(moduleCellPvalue)))
    #edges
    write.csv(moduleCellCor,file=file.path(dir.results,paste('table',q.i,'.',an.count,'.',i,'.',names(edges)[edge],".",names(questions)[[q.i]],'.',analysis,names(datalist)[i],"ModuleCellCor_edge",edge,".csv",sep="")))
    #edge annotation
    write.csv(moduleCellQvalue,file=file.path(dir.results,paste('table',q.i,'.',an.count,'.',i,names(edges)[edge],".",'.',names(questions)[[q.i]],'.',analysis,names(datalist)[i],"ModuleCellCor_qvalue_edge",edge,".csv",sep="")))
    
    #Display correlation matrix
    # Will display correlations and their BH-adjusted p-values
    textMatrix = paste(signif(moduleCellCor, 2), "\n(",
                       signif(moduleCellQvalue, 1), ")", sep = "");
    dim(textMatrix) = dim(moduleCellCor)
    
    pdf(paste('fig_',q.i,'.',an.count,'.',figure,'_ME_flowcellprop_correlation_edge',edge,'.pdf',sep=""),height=nrow(moduleCellCor)/2,width=ncol(moduleCellCor)*0.75) 
    par(mar = c(8,10,2,3));par(mfrow=c(1,1))
    # Display the correlation values within a heatmap plot
    labeledHeatmap(Matrix = moduleCellCor,
                   xLabels = colnames(cellprops2m)[retain],
                   yLabels = names(ME.edges[[edge]]),
                   ySymbols = names(ME.edges[[edge]]),
                   colorLabels = FALSE,
                   colors = blueWhiteRed(50),
                   textMatrix = textMatrix,
                   setStdMargins = FALSE,
                   cex.text = 0.9,
                   zlim = c(-1,1),
                   main = paste("ME-celltype correlation",names(edges)[edge],"with",datapoints,"datapoints"))
    par(parbackup)
    figure<-figure+1
    dev.off()
    
    # clustered version, with NA columns removed
    
    par(mar = c(2,2,2,5));par(cex.main=.5);
    pdf(paste('fig_',q.i,'.',an.count,'.',figure,'_ME_flowcellprop_correlation_edge_clust',edge,'.pdf',sep=""),height=pdf.options()$height*1,width=pdf.options()$width*1) 
    heatmap(moduleCellCor[,!apply(moduleCellCor,2,is.na)[1,]],col=blueWhiteRed(50),main=paste("Clustered module-Cell relationships",names(edges)[edge]),mar=c(5,1),cexRow=.8,cexCol=.8)
    par(parbackup)
    dev.off()
    figure<-figure+1
  }
  
  
  
  ##test ratio
  #stats for comparisons (edge 5)
  tests<-data.frame()
  testsdata<-t(cellprops2m)
  
  statlist<-list()
  ratiolist<-list()
  for (k in 1:nrow(testsdata)) {
    
    statlist[[k]]<-
      wilcox.test(testsdata[k,colnames(testsdata)==levs[3]|colnames(testsdata)==levs[4]],testsdata[k,colnames(testsdata)==levs[1]|colnames(testsdata)==levs[2]])$p.value
    
    ratiolist[[k]]<-
      median(na.omit(testsdata[k,colnames(testsdata)==levs[3]|colnames(testsdata)==levs[4]]))/median(na.omit(testsdata[k,colnames(testsdata)==levs[1]|colnames(testsdata)==levs[2]]))
    
  }
  tests<-rbind(tests,data.frame("type"="flat","names"=rownames(testsdata),"ratio"=as.vector(unlist(ratiolist)),"pval"=as.vector(unlist(statlist)),"man_B"=as.vector(unlist(statlist))*nrow(decmat.red),"bonf"=p.adjust(as.vector(unlist(statlist)),method="bonferroni"),"BH"=p.adjust(as.vector(unlist(statlist)),method="BH")))
  
  
  tests<-tests[with(tests,order(names)),]
  tests$edge<-rep(5,nrow(tests))
  flowcellproptests.flat<-tests
  write.csv(flowcellproptests.flat,file=file.path(dir.results,paste("FLOWStats_",q.i,".",an.count,".",i,"by_",contrast.variable1,"_flat_wilcox.csv",sep="")))
  
  
  
  #stats for comparisons 4 edges
  tests<-data.frame()
  
  testsdata<-t(cellprops2m)
  for (the in list(c(3,1),c(4,2),c(2,1),c(4,3))){
    statlist<-list()
    ratiolist<-list()
    for (k in 1:nrow(testsdata)) {
      
      statlist[[k]]<-
        wilcox.test(testsdata[k,colnames(testsdata)==levs[the[[1]]]],testsdata[k,colnames(testsdata)==levs[the[[2]]]])$p.value
      
      ratiolist[[k]]<-
        median(na.omit(testsdata[k,colnames(testsdata)==levs[the[[1]]]]))/median(na.omit(testsdata[k,colnames(testsdata)==levs[the[[2]]]]))
      
    }
    tests<-rbind(tests,data.frame("type"=paste(as.character(the),collapse="."),"names"=rownames(testsdata),"ratio"=as.vector(unlist(ratiolist)),"pval"=as.vector(unlist(statlist)),"man_B"=as.vector(unlist(statlist))*nrow(testsdata),"bonf"=p.adjust(as.vector(unlist(statlist)),method="bonferroni"),"BH"=p.adjust(as.vector(unlist(statlist)),method="BH")))
  }
  
  tests<-tests[with(tests,order(names)),]
  tests$edge<-rep(1:4,nrow(tests)/4)
  flowcellproptests.edge<-tests
  write.csv(flowcellproptests.edge,file=file.path(dir.results,paste("FLOW_Stats_",q.i,'.',an.count,'.',i,"by_",contrast.variable1,"_edge_wilcox.csv",sep="")))
  
  #combine stats and save
  flowcellproptests<-rbind(flowcellproptests.edge,flowcellproptests.flat)
  write.csv(flowcellproptests,file=file.path(dir.results,paste("FLOW_Stats_",q.i,'.',an.count,'.',i,"by_",contrast.variable1,"_wilcox.csv",sep="")))
  
  ##end ratio
  
  
}#end if
#11.c NETWORK 2c module cellprop (for which we have flow data)####
if (dataset.variables[[i]]=="berry.test"){
  #NETWORK 2c raw data
  # Define numbers of genes and samples
  nGenes = ncol(datExpr);
  nSamples = nrow(datExpr);
  # Recalculate edgewise MEs with color labels and select those with flow info only
  MEs0=bwMEs
  MEs = orderMEs(MEs0)
  ME.edges.flow=list()
  for (edge in 1:length(edges)){
    edgeset<-edges[[edge]]
    flowset<-which(!is.na(cellprops2$Th))
    edgeflow<-intersect(edgeset,flowset)
    MEs0.edge = MEs[edgeflow,]
    MEs.edge = orderMEs(MEs0.edge)
    ME.edges.flow[[names(edges)[edge]]]<-MEs.edge
  }
  cellprops3<-t(as.matrix(coef(decdat1)))[which(!is.na(cellprops2$Th)),]
  
  #names(cellprops)<-rownames(decmat.red)
  MEs0=bwMEs[which(!is.na(cellprops2$Th)),]#all deconvolution results where flow data exists, regardless of edge
  MEs = orderMEs(MEs0)
  moduleCellCor = cor(MEs, cellprops3, use = "p");
  retain<-!apply(moduleCellCor,2,is.na)[1,]
  moduleCellCor = moduleCellCor[,retain] #NA columns removed!
  
  moduleCellPvalue = corPvalueStudent(moduleCellCor, nSamples);
  moduleCellQvalue = matrix(p.adjust(moduleCellPvalue,method="BH"),nrow=nrow(moduleCellPvalue),dimnames=list(row.names(moduleCellPvalue),colnames(moduleCellPvalue)))
  #edges
  write.csv(moduleCellCor,file=file.path(dir.results,paste('table',q.i,'.',an.count,'.',i,'.',names(questions)[[q.i]],'.',analysis,names(datalist)[i],"ModuleCellCor_subset.csv",sep="")))
  #edge annotation
  write.csv(moduleCellQvalue,file=file.path(dir.results,paste('table',q.i,'.',an.count,'.',i,'.',names(questions)[[q.i]],'.',analysis,names(datalist)[i],"ModuleCellCor_qvalue_subset.csv",sep="")))
  #Will display correlations and their BH adjusted q-values
  textMatrix = paste(signif(moduleCellCor, 2), "\n(",
                     signif(moduleCellQvalue, 1), ")", sep = "");
  dim(textMatrix) = dim(moduleCellCor)
  
  pdf(paste('fig_',q.i,'.',an.count,'.',figure,'_ME_cellprop_subset_correlation.pdf',sep=""),height=nrow(moduleCellCor)/2,width=ncol(moduleCellCor)*0.75) 
  par(mar = c(6,10,3,2));par(mfrow=c(1,1))
  # Display the correlation values within a heatmap plot
  labeledHeatmap(Matrix = moduleCellCor[den,],
                 xLabels = names(cellprops),
                 yLabels = names(MEs)[den],
                 ySymbols = names(MEs)[den],
                 colorLabels = FALSE,
                 colors = blueWhiteRed(50),
                 textMatrix = textMatrix[den,],
                 setStdMargins = FALSE,
                 cex.text = 0.9,
                 zlim = c(-1,1),
                 main = paste("ME-celltype correlation"))
  par(parbackup)
  dev.off()
  figure<-figure+1
  
  # clustered version, with NA columns removed
  pdf(paste('fig_',q.i,'.',an.count,'.',figure,'_ME_cellprop_subset_correlation_clust.pdf',sep=""),height=pdf.options()$height*1,width=pdf.options()$width*1) 
  par(mar = c(2,2,2,5));par(cex.main=.5);
  heatmap(moduleCellCor[,!apply(moduleCellCor,2,is.na)[1,]],col=blueWhiteRed(50),main="Clustered module-cell relationships",mar=c(5,1),cexRow=.8,cexCol=.8)
  par(parbackup)
  dev.off()
  figure<-figure+1
  
  #edgewise module cell associations
  
  for (edge in 1:length(edges)){
    nSamples=length(edges[[edge]])
    
    #Correlate MEs and cellprops3
    datapoints<-nrow(cellprops3[which(!is.na(cellprops2$Th))%in%edges[[edge]],])
    
    moduleCellCor = cor(ME.edges.flow[[edge]],cellprops3[which(!is.na(cellprops2$Th))%in%edges[[edge]],], use = "p")
    retain<-!apply(moduleCellCor,2,is.na)[1,]
    moduleCellCor = moduleCellCor[,retain] #NA columns removed!
    moduleCellPvalue = corPvalueStudent(moduleCellCor, nSamples);
    moduleCellQvalue = matrix(p.adjust(moduleCellPvalue,method="BH"),nrow=nrow(moduleCellPvalue),dimnames=list(row.names(moduleCellPvalue),colnames(moduleCellPvalue)))
    #edges
    write.csv(moduleCellCor,file=file.path(dir.results,paste('table',q.i,'.',an.count,'.',i,'.',names(edges)[edge],".",names(questions)[[q.i]],'.',analysis,names(datalist)[i],"ModuleCell_withflow_Cor_edge",edge,".csv",sep="")))
    #edge annotation
    write.csv(moduleCellQvalue,file=file.path(dir.results,paste('table',q.i,'.',an.count,'.',i,names(edges)[edge],".",'.',names(questions)[[q.i]],'.',analysis,names(datalist)[i],"ModuleCell_withflow_Cor_qvalue_edge",edge,".csv",sep="")))
    
    #Display correlation matrix
    # Will display correlations and their BH-adjusted p-values
    textMatrix = paste(signif(moduleCellCor, 2), "\n(",
                       signif(moduleCellQvalue, 1), ")", sep = "");
    dim(textMatrix) = dim(moduleCellCor)
    
    pdf(paste('fig_',q.i,'.',an.count,'.',figure,'_ME_cellprop_correlation_withflow_edge',edge,'.pdf',sep=""),height=nrow(moduleCellCor)/2,width=ncol(moduleCellCor)*0.75) 
    par(mar = c(8,10,2,3));par(mfrow=c(1,1))
    # Display the correlation values within a heatmap plot
    labeledHeatmap(Matrix = moduleCellCor,
                   xLabels = colnames(cellprops3)[retain],
                   yLabels = names(ME.edges.flow[[edge]]),
                   ySymbols = names(ME.edges.flow[[edge]]),
                   colorLabels = FALSE,
                   colors = blueWhiteRed(50),
                   textMatrix = textMatrix,
                   setStdMargins = FALSE,
                   cex.text = 0.9,
                   zlim = c(-1,1),
                   main = paste("ME-celltype correlation",names(edges)[edge]))
    par(parbackup)
    figure<-figure+1
    dev.off()
    
    # clustered version, with NA columns removed
    
    par(mar = c(2,2,2,5));par(cex.main=.5);
    pdf(paste('fig_',q.i,'.',an.count,'.',figure,'_ME_cellprop_withflow_correlation_edge_clust',edge,'.pdf',sep=""),height=pdf.options()$height*1,width=pdf.options()$width*1) 
    heatmap(moduleCellCor[,!apply(moduleCellCor,2,is.na)[1,]],col=blueWhiteRed(50),main=paste("Clustered module-Cell_prop_flow relationships",names(edges)[edge]),mar=c(5,1),cexRow=.8,cexCol=.8)
    par(parbackup)
    dev.off()
    figure<-figure+1
  }
  
  ##test ratio (original proportions form deconvolution)
  
  #data prep
  cellprops3.cat<-cellprops3[order(rownames(cellprops3)),]
  combosub<-cellprops2_3[which(cellprops2_3$sample_name%in%rownames(cellprops3)),]
  combosub<-combosub[order(combosub$sample_name),]
  rownames(cellprops3.cat)<-combosub$comboname
  testsdata<-t(cellprops3.cat)
  
  #stats for comparisons (edge 5)
  tests<-data.frame()
  
  
  statlist<-list()
  ratiolist<-list()
  
  for (k in 1:nrow(testsdata)) {
    
    statlist[[k]]<-
      wilcox.test(testsdata[k,colnames(testsdata)==levs[3]|colnames(testsdata)==levs[4]],testsdata[k,colnames(testsdata)==levs[1]|colnames(testsdata)==levs[2]])$p.value
    
    ratiolist[[k]]<-
      median(na.omit(testsdata[k,colnames(testsdata)==levs[3]|colnames(testsdata)==levs[4]]))/median(na.omit(testsdata[k,colnames(testsdata)==levs[1]|colnames(testsdata)==levs[2]]))
    
  }
  tests<-rbind(tests,data.frame("type"="flat","names"=rownames(testsdata),"ratio"=as.vector(unlist(ratiolist)),"pval"=as.vector(unlist(statlist)),"man_B"=as.vector(unlist(statlist))*nrow(decmat.red),"bonf"=p.adjust(as.vector(unlist(statlist)),method="bonferroni"),"BH"=p.adjust(as.vector(unlist(statlist)),method="BH")))
  
  
  tests<-tests[with(tests,order(names)),]
  tests$edge<-rep(5,nrow(tests))
  cellprop_withflow_tests.flat<-tests
  write.csv(cellprop_withflow_tests.flat,file=file.path(dir.results,paste("cellprop_with_flow_Stats_",q.i,".",an.count,".",i,"by_",contrast.variable1,"_flat_wilcox.csv",sep="")))
  
  #stats for comparisons (edges 1-4)
  tests<-data.frame()
  
  
  for (the in list(c(3,1),c(4,2),c(2,1),c(4,3))){
    statlist<-list()
    ratiolist<-list()
    for (k in 1:nrow(testsdata)) {
      
      statlist[[k]]<-
        wilcox.test(testsdata[k,colnames(testsdata)==levs[the[[1]]]],testsdata[k,colnames(testsdata)==levs[the[[2]]]])$p.value
      
      ratiolist[[k]]<-
        median(na.omit(testsdata[k,colnames(testsdata)==levs[the[[1]]]]))/median(na.omit(testsdata[k,colnames(testsdata)==levs[the[[2]]]]))
      
    }
    tests<-rbind(tests,data.frame("type"=paste(as.character(the),collapse="."),"names"=rownames(testsdata),"ratio"=as.vector(unlist(ratiolist)),"pval"=as.vector(unlist(statlist)),"man_B"=as.vector(unlist(statlist))*nrow(testsdata),"bonf"=p.adjust(as.vector(unlist(statlist)),method="bonferroni"),"BH"=p.adjust(as.vector(unlist(statlist)),method="BH")))
  }
  
  tests<-tests[with(tests,order(names)),]
  tests$edge<-rep(1:4,nrow(tests)/4)
  cellprop_withflow_tests.edge<-tests
  write.csv(cellprop_withflow_tests.edge,file=file.path(dir.results,paste("cellprop_with_flow_Stats_",q.i,'.',an.count,'.',i,"by_",contrast.variable1,"_edges_wilcox.csv",sep="")))
  
  #combine and save
  cellprop_withflow_tests<-rbind(cellprop_withflow_tests.edge,cellprop_withflow_tests.flat)
  write.csv(cellprop_withflow_tests,file=file.path(dir.results,paste("cellprop_with_flow_Stats_",q.i,'.',an.count,'.',i,"by_",contrast.variable1,"_wilcox.csv",sep="")))
  
  
  ##end ratio
  
  
  
  ##test ratio (new proportions from deconvolution, where categories are combined to match flow categories) also known as network 2d
  
  #data prep
  cellprops3.cat<-cellprops3[order(rownames(cellprops3)),]
  combosub<-cellprops2_3[which(cellprops2_3$sample_name%in%rownames(cellprops3)),]
  combosub<-combosub[order(combosub$sample_name),]
  rownames(cellprops3.cat)<-combosub$comboname
  
  #new prop. matrix
  newCD4<-apply(cellprops3.cat[,1:2],1,sum)
  newCD8<-apply(cellprops3.cat[,3:4],1,sum)
  newB<-apply(cellprops3.cat[,5:10],1,sum)
  newNK<-apply(cellprops3.cat[,11:12],1,sum)
  newMono<-apply(cellprops3.cat[,13:14],1,sum)
  newDC<-apply(cellprops3.cat[,15:16],1,sum)
  newNeutro<-cellprops3.cat[,17]
  
  newprops<-cbind(newB,newDC,newMono,newNeutro,newNK,newCD4,newCD8)
  rownames(newprops)<-combosub$comboname
  testsdata<-t(newprops)
  
  #stats for comparisons (edge 5)
  tests<-data.frame()
  
  
  statlist<-list()
  ratiolist<-list()
  
  for (k in 1:nrow(testsdata)) {
    
    statlist[[k]]<-
      wilcox.test(testsdata[k,colnames(testsdata)==levs[3]|colnames(testsdata)==levs[4]],testsdata[k,colnames(testsdata)==levs[1]|colnames(testsdata)==levs[2]])$p.value
    
    ratiolist[[k]]<-
      median(na.omit(testsdata[k,colnames(testsdata)==levs[3]|colnames(testsdata)==levs[4]]))/median(na.omit(testsdata[k,colnames(testsdata)==levs[1]|colnames(testsdata)==levs[2]]))
    
  }
  tests<-rbind(tests,data.frame("type"="flat","names"=rownames(testsdata),"ratio"=as.vector(unlist(ratiolist)),"pval"=as.vector(unlist(statlist)),"man_B"=as.vector(unlist(statlist))*nrow(decmat.red),"bonf"=p.adjust(as.vector(unlist(statlist)),method="bonferroni"),"BH"=p.adjust(as.vector(unlist(statlist)),method="BH")))
  
  tests<-tests[with(tests,order(names)),]
  tests$edge<-rep(5,nrow(tests))
  cellprop_withflow_tests.flat.new<-tests
  write.csv(cellprop_withflow_tests.flat.new,file=file.path(dir.results,paste("cellprop_with_flow_Stats_NEW_",q.i,".",an.count,".",i,"by_",contrast.variable1,"_flat_wilcox.csv",sep="")))
  
  #stats for comparisons (edges 1-4)
  tests<-data.frame()
  
  
  for (the in list(c(3,1),c(4,2),c(2,1),c(4,3))){
    statlist<-list()
    ratiolist<-list()
    for (k in 1:nrow(testsdata)) {
      
      statlist[[k]]<-
        wilcox.test(testsdata[k,colnames(testsdata)==levs[the[[1]]]],testsdata[k,colnames(testsdata)==levs[the[[2]]]])$p.value
      
      ratiolist[[k]]<-
        median(na.omit(testsdata[k,colnames(testsdata)==levs[the[[1]]]]))/median(na.omit(testsdata[k,colnames(testsdata)==levs[the[[2]]]]))
      
    }
    tests<-rbind(tests,data.frame("type"=paste(as.character(the),collapse="."),"names"=rownames(testsdata),"ratio"=as.vector(unlist(ratiolist)),"pval"=as.vector(unlist(statlist)),"man_B"=as.vector(unlist(statlist))*nrow(testsdata),"bonf"=p.adjust(as.vector(unlist(statlist)),method="bonferroni"),"BH"=p.adjust(as.vector(unlist(statlist)),method="BH")))
  }
  
  tests<-tests[with(tests,order(names)),]
  tests$edge<-rep(1:4,nrow(tests)/4)
  cellprop_withflow_tests.edge.new<-tests
  write.csv(cellprop_withflow_tests.edge.new,file=file.path(dir.results,paste("cellprop_with_flow_Stats_NEW",q.i,'.',an.count,'.',i,"by_",contrast.variable1,"_edges_wilcox.csv",sep="")))
  
  #combine and save
  cellprop_withflow_tests.new<-rbind(cellprop_withflow_tests.edge.new,cellprop_withflow_tests.flat.new)
  write.csv(cellprop_withflow_tests.new,file=file.path(dir.results,paste("cellprop_with_flow_Stats_NEW",q.i,'.',an.count,'.',i,"by_",contrast.variable1,"_wilcox.csv",sep="")))
  
  
  ##end ratio
  ##end addiitonal test
  
  
  
}#end ifNETWORK 2c: Quantifying module cell-proportion associations (limited to samples with flow) ####
#12. GS/MM Gene relationship to trait and important modules: Gene Significance and Module Membership####

#We quantify associations of individual genes with our trait of interest (which depends on the question/ contrast) by defining Gene Significance GS 
#as (the absolute value of) the correlation between the gene and the trait. For each module, we also define a 
#quantitative measure of module membership MM as the correlation of the module eigengene and the gene expression 
#profile. This allows us to quantify the similarity of all genes on the array to every module.

#version 1: by contrast (other version would be by edge)

geneTraitSignificance<-list()
GSPvalue<-list()

for (contrastDirection in contrast.variable){
  
  #in this context 'class' refers to the contrast variable!
  class=as.data.frame(eval(parse(text=paste("datTraits$",contrastDirection))))
  names(class)=contrastDirection
  MEs0=bwMEs
  MEs = orderMEs(MEs0)
  # names (colors) of the modules (equivalent to using function signedKME())
  modNames = substring(names(MEs), 3)
  geneModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p"));
  MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));
  names(geneModuleMembership) = paste("MM", modNames, sep="");
  names(MMPvalue) = paste("p.MM", modNames, sep="");
  
  geneTraitSignificance[[contrastDirection]] = as.data.frame(cor(datExpr, class, use = "p"));
  GSPvalue[[contrastDirection]] = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance[[contrastDirection]]), nSamples));
  names(geneTraitSignificance[[contrastDirection]]) = paste("GS.", names(class), sep="");
  names(GSPvalue[[contrastDirection]]) = paste("p.GS.", names(class), sep="");
  
  #Intramodular analysis: identifying genes with high GS and MM
  
  #Using the GS and MM measures, we can identify genes that have a high significance for class as well as high module membership in interesting modules. We plot a scatterplot of Gene Significance vs. Module Membership in all modules:
  
  #add code to label the top 2 or 3 genes
  #then fish out the probes and make a heatmap consisting of the hub genes
  
  par(parbackup);
  toptwovector<-vector()
  tophubvector<-vector()
  toptwokMEvector<-vector()
  tophubkMEvector<-vector()
  
  #GS/MM plot
  pdf(paste('fig_',q.i,'.',an.count,'.',figure,'_GS_MM_plot_',contrastDirection,'.pdf',sep=""),height=pdf.options()$height*1,width=pdf.options()$width*1) 
  par(mfrow=c(3,3))
  for (mods in modNames){
    module = mods
    column = match(module, modNames);
    moduleGenes = bwmoduleColors==module;
    
    
    if(length(abs(geneModuleMembership[moduleGenes,column]))>1){#removed constraint for minimum module size which is 20
      
      #the actual GS/MM plot
      verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                         abs(geneTraitSignificance[[contrastDirection]][moduleGenes, 1]),
                         abline=T,abline.color = "red",
                         xlab = paste("MM", module, "module"),
                         ylab = paste("GS",contrastDirection),
                         main = paste(mods," MM vs. GS\n"),
                         cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)
      
      
      #extract top 1 and 2 probes to identify hub genes
      moddata<-geneModuleMembership[moduleGenes, column]
      names(moddata)<-row.names(geneModuleMembership[moduleGenes,])
      toptwo<-names(sort(moddata,decreasing=TRUE)[1:2])
      tophub<-names(sort(moddata,decreasing=TRUE)[1])
      toptwovector<-c(toptwovector,toptwo)
      tophubvector<-c(tophubvector,tophub)
      toptwokMEvector<-c(toptwokMEvector,sort(moddata,decreasing=TRUE)[1:2])
      tophubkMEvector<-c(tophubkMEvector,sort(moddata,decreasing=TRUE)[1])
    }
  }#end loop over modules
  par(parbackup)
  dev.off() 
  figure<-figure+1
  
  #info on the hub probes
  
  dicthub<-nuID2IlluminaID(toptwovector,species="Human")
  dicthub.entrez<-nuID2EntrezID(toptwovector,filterTh = NULL,lib.mapping='lumiHumanIDMapping', returnAllInfo = TRUE)
  dicthub3<-merge(dicthub,dicthub.entrez,by.x="nuID",by.y=0)
  dicthub4a<-cbind(as.character(sapply(modNames,function(x){rep(x,2)})),toptwokMEvector)
  dicthub4<-merge(dicthub4a,dicthub3,by.x=0,by.y="nuID")
  names(dicthub4)[c(1,2)]<-c("nuID","module")
  dicthub4<-dicthub4[order(dicthub4$module),]
  write.csv(dicthub4,file=file.path(dir.results,paste('table',q.i,'.',an.count,'.',i,'.',names(questions)[[q.i]],'.',analysis,names(datalist)[i],"hub_probes_for_modules_",contrastDirection,".csv",sep="")))
  
  #heatmap based on 2 hub genes (probes) 
  rsc<-matrix(cbind(as.character(dicthub4$module),as.character(dicthub4$module)),ncol=2)
  subsetModHubProbes<-is.element(rownames(exprs(data.q.fil)),toptwovector)
  pdf(paste('fig_',q.i,'.',an.count,'.',figure,'_heatmap_2_hub_probes.pdf',sep=""),height=pdf.options()$height*1,width=pdf.options()$width*1) 
  superHeatmap2_rsc(x=data.q.fil,y=subsetModHubProbes,phenomatrix=phenomatrix,scale="row",addTit="2 per module",labRow=dicthub4$ILMN_Gene,rsc)
  par(parbackup)
  dev.off()
  figure<-figure+1
  
  #heatmap based on 1 hub gene (probes)
  rsc<-matrix(cbind(as.character(dicthub4$module),as.character(dicthub4$module)),ncol=2)[seq(1,nrow(dicthub4),2),]
  subsetModHubProbes<-is.element(rownames(exprs(data.q.fil)),tophubvector)
  pdf(paste('fig_',q.i,'.',an.count,'.',figure,'_heatmap_1_hub_probe.pdf',sep=""),height=pdf.options()$height*1,width=pdf.options()$width*1) 
  superHeatmap2_rsc(x=data.q.fil,y=subsetModHubProbes,phenomatrix=phenomatrix,scale="row",addTit="1 per module",labRow=dicthub4$ILMN_Gene[seq(1,nrow(dicthub4),2)])
  par(parbackup)
  dev.off()
  figure<-figure+1
  
  rsc="none"
  
}#end GS/MM by contrast
#13. Summary output of network analysis results: annot and geneInfo####

#retrieve various annotations to attach to feature and module data
features<-featureNames(usedata[tv,])
features2<-nuID2probeID(features)
annotSYM<-as.character(getSYMBOL(features,"lumiHumanAll.db"))#need new annotSYM that matches targetID annotation. targetID is full of aliases!!!
annotSYM2<-featureData(usedata[features,])$TargetID#this is pretty rubbish, but usefule if no HGNC symbol exists
if(is.null(annotSYM2)){annotSYM2<-featureData(usedata[features,])$SYMBOL}
if(is.null(annotSYM2)){annotSYM2<-featureData(usedata[features,])$REF_ID}
if(is.null(annotSYM2)){annotSYM2<-featureData(usedata[features,])$ID_REF}
if(is.null(annotSYM2)){annotSYM2<-featureData(usedata[features,])$PROBE_ID}

testannot<-cbind(annotSYM,annotSYM2)
if(!is.null(annotSYM2)){
  discr<-testannot[which(testannot[,1]!=testannot[,2]),]
}

annotACC<-as.character(lookUp(features,"lumiHumanAll.db","ACCNUM"))
annotCHR<-as.character(lookUp(features,"lumiHumanAll.db","CHR"))
annotCHRLOC<-as.character(lookUp(features,"lumiHumanAll.db","CHRLOC"))#deprecated
annotENSEMBL<-as.character(lookUp(features,"lumiHumanAll.db","ENSEMBL"))
annotENTREZID<-as.character(lookUp(features,"lumiHumanAll.db","ENTREZID"))
annotGO<-as.character(lookUp(features,"lumiHumanAll.db","GO"))
annotREFSEQ<-as.character(lookUp(features,"lumiHumanAll.db","REFSEQ"))
annotUNIPROT<-as.character(lookUp(features,"lumiHumanAll.db","UNIPROT"))
annotPATH<-as.character(lookUp(features,"lumiHumanAll.db","PATH"))
annotOMIM<-as.character(lookUp(features,"lumiHumanAll.db","OMIM"))
annot<-data.frame(cbind(features,features2,annotSYM,annotSYM2,annotACC,annotCHR,annotCHRLOC,annotENSEMBL,annotENTREZID,annotGO,annotREFSEQ,annotUNIPROT,annotPATH,annotOMIM))

#save
write.csv(annot,file=file.path(dir.results,paste('table',q.i,'.',an.count,'.',i,'.',names(questions)[[q.i]],'.',analysis,names(datalist)[i],"module_annot.csv",sep="")))

#GeneInfo file for annotation of modules
probes = names(datExpr)
probes2annot = match(probes, annot$features)
print(paste("The following is the number or probes without annotation:",sum(is.na(probes2annot))))

# Create the starting data frame
geneInfo0 = data.frame(substanceBXH = probes,
                       probeID = annot$features2[probes2annot],
                       geneSymbol = annot$annotSYM[probes2annot],
                       HGNC_symbol = annot$annotSYM[probes2annot],
                       targetID = annot$annotSYM2[probes2annot],
                       EntrezID = annot$annotENTREZID[probes2annot],
                       moduleColor = bwmoduleColors,
                       geneTraitSignificance[[1]],
                       GSPvalue[[1]],
                       geneTraitSignificance[[2]],
                       GSPvalue[[2]],
                       OMIMid = annot$annotOMIM[probes2annot],
                       KEGGid = annot$annotPATH[probes2annot])
# Order modules by their significance for class
modOrder = order(-abs(cor(MEs, class, use = "p")));
# Add module membership information in the chosen order
for (mod in 1:ncol(geneModuleMembership)){
  oldNames = names(geneInfo0)
  geneInfo0 = data.frame(geneInfo0, geneModuleMembership[, modOrder[mod]],
                         MMPvalue[, modOrder[mod]]);
  names(geneInfo0) = c(oldNames, paste("MM.", modNames[modOrder[mod]], sep=""),
                       paste("p.MM.", modNames[modOrder[mod]], sep=""))
}
# Order the genes in the geneInfo variable first by module color, then by geneTraitSignificance for the first of the two contrast variables
geneOrder = order(geneInfo0$moduleColor, -abs(eval(parse(text=paste("geneInfo0$GS.",contrast.variable1,sep="")))));
geneInfo = geneInfo0[geneOrder, ]
write.csv(geneInfo,file=file.path(dir.results,paste('table',q.i,'.',an.count,'.',i,'.',names(questions)[[q.i]],'.',analysis,names(datalist)[i],"geneInfoClass2.csv",sep="")))

#Output genes of interest: GBPs and caspases (? expand to include other cell-death pathway genes)
casp.gbp<-geneInfo[which(geneInfo$geneSymbol%in%c("GBP1","GBP2","GBP3","GBP4","GBP5","GBP6","GBP7","CASP1","CASP2","CASP3","CASP4","CASP5","CASP6","CASP7","CASP8","CASP9")),]
write.csv(casp.gbp,file=file.path(dir.results,paste('table',q.i,'.',an.count,'.',i,'.',names(questions)[[q.i]],'.',analysis,names(datalist)[i],"geneInfoCASP_GBP.csv",sep="")))

modlin<-geneInfo[which(geneInfo$geneSymbol%in%c("IL32","TARP","IL2RB","PRF1","CD8A","RALGDS","GZMH","MCOLN2")),]
write.csv(modlin,file=file.path(dir.results,paste('table',q.i,'.',an.count,'.',i,'.',names(questions)[[q.i]],'.',analysis,names(datalist)[i],"geneInfoModlin8.csv",sep="")))
#14. Module/ signature relationships (sigenrich)####
modinfo<-geneInfo
modulecolors<-unique(as.character(modinfo$moduleColor))
signaturelist<-list(edge1.fac.TT.500,edge2.fac.TT.500,edge3.fac.TT.500,edge4.fac.TT.500,edge5.flat.TT.500)
#need for loop for each edge
for (edge in 1:length(edges)){
  #module signature enrichment plot (some variables duplicated at this point)
  signature<-signaturelist[[edge]]
  sigprobes<-rownames(signature)
  
  moduleprobes<-list()
  for (col in modulecolors){
    probes<-modinfo$substanceBXH[modinfo$moduleColor==col]
    moduleprobes[[col]]<-probes
  }
  
  allprobes<-as.character(modinfo$substanceBXH)
  allgenes<-as.character(modinfo$geneSymbol)
  enrichlist<-list()
  genesInMod<-list()
  enrichnormlist<-list()
  for (iter in 1:length(moduleprobes)){
    enrich<-sum(sigprobes%in%moduleprobes[[iter]])/length(sigprobes)
    #print(paste(names(moduleprobes)[iter],":",enrich))
    enrichlist[[names(moduleprobes)[iter]]]<-enrich
    
    enrichnorm=enrich*ncol(datExpr)/lapply(moduleprobes,length)[[iter]]
    enrichnormlist[[iter]]<-enrichnorm
    
    genes<-allgenes[allprobes%in%moduleprobes[[iter]]]
    genesInMod[[iter]]<-genes
  }
  names(genesInMod)<-modulecolors
  names(enrichnormlist)<-modulecolors
  sum(unlist(enrichlist))
  
  pdf(paste('fig_',q.i,'.',an.count,'.',figure,'_sigenrich_edge',edge,'.pdf',sep=""),height=pdf.options()$height*1,width=pdf.options()$width*1) 
  par(mfrow=c(3,1),mar=c(9,6,3,3))
  barplot(unlist(lapply(moduleprobes,length)),las=2,main="Module size (N probes)",col=modulecolors)
  barplot(unlist(enrichlist),las=2,main="Proportion of signature probes in module",col=modulecolors)
  barplot(unlist(enrichnormlist),las=2,main="Ratio of real to expected signature probe enrichment",col=modulecolors)
  abline(h=1,col="red")
  par(parbackup)
  dev.off()
  figure<-figure+1
  
  sigenrich<-data.frame(cbind(names(enrichlist),unlist(enrichlist)))
  colnames(sigenrich)<-c("module","% sig")
  write.csv(sigenrich,file=file.path(dir.results,paste('table',q.i,'.',an.count,'.',i,'.',names(questions)[[q.i]],'.',analysis,names(datalist)[i],"_node_module_sig_edge",edge,".csv",sep="")))
  
  sigenrichnorm<-data.frame(cbind(names(enrichnormlist),unlist(enrichnormlist)))
  colnames(sigenrichnorm)<-c("module","ratio_sig")
  write.csv(sigenrichnorm,file=file.path(dir.results,paste('table',q.i,'.',an.count,'.',i,'.',names(questions)[[q.i]],'.',analysis,names(datalist)[i],"_node_module_sig_norm_edge",edge,".csv",sep="")))
  
  #Move sigenrich properties to Neo4j
  if (neoinsert=="on"){
    
    nodes<-RNeo4j::nodes
    #network-specific cypher query
    graphdata=sigenrichnorm
    #query = paste("MERGE (wmod:wgcna {name:{wgcnamod},square:{squareg},edge:{edgeg},contrastvar:{contrastvarsg},contrast:{contrastg}}) ON CREATE SET wmod.sigenrich = ",modsig, " ON MATCH SET wmod.sigenrich = ",modsig,sep="")
    t = suppressMessages(newTransaction(graph))
    
    for (therow in 1:nrow(graphdata)) { #already in edgewise loop!
      the_wmod = as.character(graphdata[therow, ]$module)
      #modsig = as.numeric(as.character(graphdata[therow, ]$ratio_sig))        
      modsig = as.character(graphdata[therow, ]$ratio_sig)
      query = paste("MERGE (wmod:wgcna {name:{wgcnamod},square:{squareg},edge:{edgeg},contrastvar:{contrastvarsg},contrast:{contrastg}}) ON CREATE SET wmod.sigenrich = ",modsig, " ON MATCH SET wmod.sigenrich = ",modsig,sep="")
      
      suppressMessages(appendCypher(t, 
                                    query,
                                    wgcnamod = the_wmod,
                                    #modsig = the_modsig,
                                    #invariant node properties
                                    squareg = dataset.variables[[1]],
                                    edgeg = edge,
                                    contrastvarsg = contrastvars[[edge]],
                                    contrastg = paste(pathways[[edge]][[2]],"-",pathways[[edge]][[1]])                    
                                    
      ))
    }
    
    print("committing :)")
    suppressMessages(commit(t))
    
    #END NEO4J
    
  }#end neoinsert if
  
}#close edgewise sigenrich
#15. Output gene lists for use with online software and services.####
#superseded by neo4j

# Match probes in the data set to the probe IDs in the annotation file
probes = names(datExpr)
probes2annot = match(probes, annot$features)
# Get the corresponding Locus Link IDs
allLLIDs = as.numeric(annot$annotENTREZID[probes2annot]);
allLLIDs=allLLIDs[!is.na(allLLIDs)]
# $ Choose interesting modules
intModules = unique(bwmoduleColors)[modOrder]
#16. GO enrichment####
# As background in the enrichment analysis, we will use all probes in the analysis.
# fileName = file.path(dir.results,paste('ALL-LLID',q.i,'.',an.count,'.',i,'.',names(datalist)[[i]],"_entrez_ids-all.txt"));
# write.table(as.data.frame(allLLIDs), file = fileName,row.names = FALSE, col.names = FALSE)

# #GOenr = GOenrichmentAnalysis(bwmoduleColors,allLLIDs,organism = "human", nBestP = 20,includeOffspring=TRUE,verbose=100)
# GOenr = goenrich2(bwmoduleColors,allLLIDs,organism = "human", nBestP = 20,includeOffspring=FALSE,verbose=100)#the current version of the function is broken; I fixed it and renamed it goenrich2, which is sourced by setup script
# tab = GOenr$bestPTerms[[4]]$enrichment
# names(tab)
# write.table(tab,file=file.path(dir.results,paste('table',q.i,'.',an.count,'.',i,'.',names(datalist)[[i]],"_GOEnrichmentTable.csv",sep="")),sep=",",quote=TRUE,row.names=FALSE)
# 
# #on-screen
# keepCols = c(1, 2, 5, 6, 7, 12, 13);
# screenTab = tab[, keepCols];
# # Round the numeric columns to 2 decimal places:
# numCols = c(3, 4);
# screenTab[, numCols] = signif(apply(screenTab[, numCols], 2, as.numeric), 2)
# # Truncate the the term name to at most 150 characters
# screenTab[, 7] = substring(screenTab[, 7], 1, 100)
# # Shorten the column names:
# colnames(screenTab) = c("module", "size", "p-val", "Bonf", "nInTerm", "ont", "term name");
# rownames(screenTab) = NULL;
# # Set the width of R output. The reader should play with this number to obtain satisfactory output.
# options(width=40)
# # Finally, display the enrichment table:
# screenTab
#17. TOM plot####
#make TOMplot#only run this if lots of memory. don't include in RData

#problems:
#1. not complete network, but only one block at a time
#2. despite this, it is too big

# Calculate topological overlap anew: this could be done more efficiently by saving the TOM

#get hold of TOMs inRData files in ~

# calculated during module detection, but let us do it again here. #!No!
#dissTOM = 1-TOMsimilarityFromExpr(datExpr, power = 11);

# do this for eachj block, i.e.TOM?

##don't do TOM plots at this stage
load(file.path(dir.results,"TOM_blockwise-block.1.RData"))
TOM<-as.matrix(TOM)
dissTOM1=1-TOM
# Transform dissTOM with a power to make moderately strong connections more visible in the heatmap
plotTOM1 = dissTOM1^7;
# Set diagonal to NA for a nicer plot
diag(plotTOM1) = NA;
# Call the plot function
#sizeGrWindow(9,9)

# # #-- can't plot whole tom due to memory constraints...need to set ulimit -s 64. issues with C stack usage
# 
# plotDendroAndColors(bwnet$dendrograms[[1]], bwmoduleColors[bwnet$blockGenes[[1]]],
#                     "Module colors", main = "Gene dendrogram and module colors in block 1",
#                     dendroLabels = FALSE, hang = 0.03,
#                     addGuide = TRUE, guideHang = 0.05)
# 
# pdf(paste('fig_',q.i,'.',an.count,'.',figure,'_TOMPLOT.pdf',sep=""),height=pdf.options()$height*1,width=pdf.options()$width*1) 
# TOMplot(plotTOM1,bwnet$dendrograms[[1]], bwmoduleColors[bwnet$blockGenes[[1]]], main = "Network heatmap plot, all genes, block1")
# dev.off()
# 
# Cstack_info()

#Plot only some genes..
# #commented out from here
# nSelect = 1500
# # For reproducibility, we set the random seed
# set.seed(10);
# select = sample(dim(plotTOM1)[[1]], size = nSelect);
# selectTOM = dissTOM1[select, select];
# # There is no simple way of restricting a clustering tree to a subset of genes, so we must re-cluster.
# selectTree = flashClust(as.dist(selectTOM), method = "average")
# selectColors = bwmoduleColors[select];
# # Open a graphical window
# #sizeGrWindow(9,9)
# # Taking the dissimilarity to a power, say 10, makes the plot more informative by effectively changing
# # the color palette; setting the diagonal to NA also improves the clarity of the plot
# plotDiss = selectTOM^7;
# diag(plotDiss) = NA;
# TOMplot(plotDiss, selectTree, selectColors, main = "Network heatmap plot, selected genes")
# export.plot(file.prefix=paste('fig',q.i,'.',an.count,'.',i,'.',figure,names(questions)[[q.i]],'.',analysis,names(datalist)[i],"TOM_plot"),export.formats=export.formats.plots,height=height,width=(1 + sqrt(5))/2*height)
# par(parbackup)
# 
# figure<-figure+1
# 
# cmd1=cmdscale(as.dist(dissTOM1),2)
# sizeGrWindow(7, 6)
# par(mfrow=c(1,1))
# thecol<-bwmoduleColors[bwnet$blockGenes[[1]]]
# plot(cmd1, col=as.character(thecol), main="MDS plot",pch=4,xlab="Scaling Dimension 1", ylab="Scaling Dimension 2")

#rm(TOM)
#18. Visualize the network of eigengenes####

MEtraits<-list()
for (trait in 1:length(datTraits)){
  tmp<-as.data.frame(datTraits[[trait]])
  names(tmp)<-names(datTraits[trait])
  string=paste("ME network for ",names(datTraits[trait]),sep="")
  data=orderMEs(cbind(MEs,tmp))
  MEtraits[[trait]]=data
  names(MEtraits)[[trait]]=string
}

# Plot the relationships among the eigengenes and the trait
par(cex = 0.9)
for (m in 1:length(MEtraits)){
  pdf(paste('fig_',q.i,'.',an.count,'.',figure,'_ME_network_for',names(MEtraits)[[m]],'.pdf',sep=""),height=pdf.options()$height*1,width=pdf.options()$width*1)   
  plotEigengeneNetworks(MEtraits[[m]],names(MEtraits)[[m]], marDendro = c(0,12,2,7), marHeatmap = c(10,12,1,1), cex.lab = 0.8, xLabelsAngle = 90)
  dev.off() 
  figure<-figure+1
}
par(parbackup)
#19. Export to Cytoscape####

# Load topological overlap 
load(file.path(dir.results,"TOM_blockwise-block.1.RData"))
TOM<-as.matrix(TOM)

# Select modules
modules = intModules
probesM = features[bwnet$blockGenes[[1]]]

#main Cytoscape export loop
setwd(dir.results)
modNetList<-list()
for (cytmod in 1:length(modules)){
  
  # Select module probes and genes
  inModule = is.finite(match(bwmoduleColors[bwnet$blockGenes[[1]]], modules[cytmod]));#exchanged bwModuleColors for bwmoduleColors
  modProbes = probesM[inModule];
  modGenes = annot$annotSYM[match(modProbes, annot$features)];
  
  # Select the corresponding Topological Overlap
  modTOM = TOM[inModule, inModule];
  dimnames(modTOM) = list(modProbes, modProbes)
  
  # Export the network into edge and node list files Cytoscape can read
  moduleNetwork<-exportNetworkToCytoscape(modTOM,
                                          edgeFile = paste('CS',q.i,'.',an.count,'.',i,'.',names(datalist)[[i]],"CS-edges-", paste(bwmoduleColors[bwnet$blockGenes[[1]][inModule]][1], collapse="-"), ".txt", sep=""),
                                          nodeFile = paste('CS',q.i,'.',an.count,'.',i,'.',names(datalist)[[i]],"CS-nodes-", paste(bwmoduleColors[bwnet$blockGenes[[1]][inModule]][1], collapse="-"), ".txt", sep=""),
                                          weighted = TRUE,
                                          threshold = 0.0005,#??????????????????thisis where edges are lost
                                          nodeNames = modProbes,
                                          altNodeNames = modGenes,
                                          nodeAttr = bwmoduleColors[bwnet$blockGenes[[1]][inModule]]
  )
  modNetList[[cytmod]]<-moduleNetwork
  print(cytmod)
  print(paste("edgefile:",paste('CS',q.i,'.',an.count,'.',i,'.',names(datalist)[[i]],"CS-edges-", paste(bwmoduleColors[bwnet$blockGenes[[1]][inModule]][1], collapse="-"), ".txt", sep="")))
  print(paste("nodefile",paste('CS',q.i,'.',an.count,'.',i,'.',names(datalist)[[i]],"CS-nodes-", paste(bwmoduleColors[bwnet$blockGenes[[1]][inModule]][1], collapse="-"), ".txt", sep="")))
  print(bwmoduleColors[bwnet$blockGenes[[1]][inModule]][1])
  rm(modTOM)
  collectGarbage()
  gc()
}
names(modNetList)<-modules
#modNetList conatins all the individual module networks

setwd(dir.figures)

#Make a list of genes in each module for use later
modgeneslist<-list()  
for (iter in modules){
  inModule = is.finite(match(bwmoduleColors[bwnet$blockGenes[[1]]], iter));#exchanged bwModuleColors for bwmoduleColors
  modProbes = probesM[inModule];
  modGenes = annot$annotSYM[match(modProbes, annot$features)];
  modgeneslist[[iter]]<-as.character(modGenes)[!is.na(as.character(modGenes))];
}

#Make a list of entrezids in each module for use later
modentrezlist<-list() 
modNUIDlist<-list()
for (iter in modules){
  inModule = is.finite(match(bwmoduleColors[bwnet$blockGenes[[1]]], iter));#exchanged bwModuleColors for bwmoduleColors
  modProbes = probesM[inModule];
  modEntrez = annot$annotENTREZID[match(modProbes, annot$features)];
  modentrezlist[[iter]]<-as.character(modEntrez)[as.character(modEntrez)!="NA"];
  modNUIDlist[[iter]]<-as.character(modProbes)
}
#20. Plot module heatmap and eigengene expression####

#Not used but code kept. Edgewise ME barplots are saved instead (below)
datME=MEs

# for (module in intModules){
#   which.module=module
#   ME=datME[, paste("ME",which.module, sep="")]
#   par(mfrow=c(2,1), mar=c(0.3, 5.8, 10, 3))#these are a bit off
#   plotMat(t(scale(datExpr[,bwmoduleColors==which.module ]) ),nrgcols=50,rlabels=F,clabels=row.names(datExpr),rcols=which.module,ccols=as.character(classcolours1[[i]]),main=which.module, cex.main=1)
#   par(mar=c(4, 3.2, 0, 0.5))
#   barplot(ME, col=which.module, main="", cex.main=2,ylab="eigengene expression",xlab="array sample")
#   par(parbackup)
# }
#figure<-figure+1
setwd(dir.figures)
#subsets of MEs (not ME recalculated based on edge data subset): 4 individual edges
for (module in intModules){
  pdf(paste('fig_',q.i,'.',an.count,'.',figure,'_heatmap_barplot_edgeMEs_',module,'.pdf',sep=""),height=pdf.options()$height*1,width=pdf.options()$width*2.3)
  par(mfcol=c(2,4)) 
  for(edge in 1:4){
    edgeset=edges[[edge]]
    #(see above)
    which.module=module
    ME=MEs[edgeset, paste("ME",which.module, sep="")]
    par(mar=c(0.3, 7.4, 11, 2))
    #par(mfrow=c(2,1), mar=c(0.3, 5.5, 13, 2))
    plotMat(t(scale(datExpr[edges[[edge]],bwmoduleColors==which.module])),nrgcols=30,rlabels=F,clabels=row.names(datExpr)[edges[[edge]]],rcols=which.module,ccols=as.character(list(classcolours1[[i]],classcolours1[[i]],classcolours2[[i]],classcolours2[[i]])[[edge]])[edges[[edge]]],main=paste(names(edges)[edge],which.module), cex.main=2)
    par(mar=c(5, 5.7, 0, 0.7))
    barplot(ME, col=which.module, main="", cex.main=2,ylab="ME expression",cex.axis=1.5,cex.names=1.5,cex.lab=2)
  }#end edge loop, keep device open!
  dev.off()
  figure<-figure+1
}

#subsets of MEs (not ME recalculated based on edge data subset): edge 5
for (module in intModules){
  pdf(paste('fig_',q.i,'.',an.count,'.',figure,'_heatmap_barplot_MEs_',module,'.pdf',sep=""),height=pdf.options()$height*1,width=pdf.options()$width*1)
  par(mfcol=c(2,1)) 
  for(edge in 5:length(edges)){
    edgeset=edges[[edge]]
    #(see above)
    which.module=module
    ME=MEs[edgeset, paste("ME",which.module, sep="")]
    par(mar=c(0.3, 7.4, 11, 2))
    #par(mfrow=c(2,1), mar=c(0.3, 5.5, 13, 2))
    plotMat(t(scale(datExpr[edges[[edge]],bwmoduleColors==which.module])),nrgcols=30,rlabels=F,clabels=row.names(datExpr)[edges[[edge]]],rcols=which.module,ccols=as.character(list(classcolours1[[i]],classcolours1[[i]],classcolours2[[i]],classcolours2[[i]],classcolours1[[i]])[[edge]])[edges[[edge]]],main=paste(names(edges)[edge],which.module), cex.main=2)
    par(mar=c(5, 5.7, 0, 0.7))
    barplot(ME, col=which.module, main="", cex.main=2,ylab="ME expression",cex.axis=1.5,cex.names=1.5,cex.lab=2)
  }#end edge loop, keep device open!
  dev.off()
  figure<-figure+1
}

##--WITH matplots

for (module in intModules){
  hubgeneID<-as.character(dicthub4[which(dicthub4$module==module),"nuID"])
  pdf(paste('fig_',q.i,'.',an.count,'.',figure,'_heatmap_barplot_matplot_edgeMEs_',module,'.pdf',sep=""),height=pdf.options()$height*1,width=pdf.options()$width*2.3)
  par(mfcol=c(3,4)) 
  for(edge in 1:4){
    edgeset=edges[[edge]]
    #(see above)
    which.module=module
    ME=MEs[edgeset, paste("ME",which.module, sep="")]
    par(mar=c(0.3, 8.4, 11, 4)) 
    #par(mfrow=c(2,1), mar=c(0.3, 5.5, 13, 2))
    plotMat(t(scale(datExpr[edges[[edge]],bwmoduleColors==which.module])),nrgcols=30,rlabels=F,clabels=row.names(datExpr)[edges[[edge]]],rcols=which.module,ccols=as.character(list(classcolours1[[i]],classcolours1[[i]],classcolours2[[i]],classcolours2[[i]])[[edge]])[edges[[edge]]],main=paste(names(edges)[edge],which.module), cex.main=2)
    par(mar=c(5, 5.7, 0, 0.7))
    barplot(ME, col=which.module, main="", cex.main=2,ylab="ME expression",cex.axis=1.5,cex.names=1.5,cex.lab=2)
    matplot(datExpr[edges[[edge]],bwmoduleColors==which.module],type="l",col="gray80",lty=1,lwd=.1,ylab="intensity",cex.axis=1.5,cex.lab=2)
    lines(datExpr[edges[[edge]],bwmoduleColors==which.module][,hubgeneID[[1]]],col="red",lwd=3)
    lines(datExpr[edges[[edge]],bwmoduleColors==which.module][,hubgeneID[[2]]],col="blue",lwd=3)
  }#end edge loop, keep device open!
  dev.off()
  figure<-figure+1
}

for (module in intModules){
  hubgeneID<-as.character(dicthub4[which(dicthub4$module==module),"nuID"])
  pdf(paste('fig_',q.i,'.',an.count,'.',figure,'_heatmap_barplot_matplot_edgeMEs_',module,'.pdf',sep=""),height=pdf.options()$height*1,width=pdf.options()$width*1)
  par(mfcol=c(3,1)) 
  for(edge in 5:5){
    edgeset=edges[[edge]]
    #(see above)
    which.module=module
    ME=MEs[edgeset, paste("ME",which.module, sep="")]
    par(mar=c(0.3, 8.4, 11, 4))
    #par(mfrow=c(2,1), mar=c(0.3, 5.5, 13, 2))
    plotMat(t(scale(datExpr[edges[[edge]],bwmoduleColors==which.module])),nrgcols=30,rlabels=F,clabels=row.names(datExpr)[edges[[edge]]],rcols=which.module,ccols=as.character(list(classcolours1[[i]],classcolours1[[i]],classcolours2[[i]],classcolours2[[i]],classcolours1[[i]])[[edge]])[edges[[edge]]],main=paste(names(edges)[edge],which.module), cex.main=2)
    par(mar=c(5, 5.4, 0, 0.7))
    barplot(ME, col=which.module, main="", cex.main=2,ylab="ME expression",cex.axis=1.5,cex.names=1.5,cex.lab=2)
    matplot(datExpr[edges[[edge]],bwmoduleColors==which.module],type="l",col="gray80",lty=1,lwd=.1,ylab="intensity",cex.axis=1.5,cex.lab=2)
    lines(datExpr[edges[[edge]],bwmoduleColors==which.module][,hubgeneID[[1]]],col="red",lwd=3)
    lines(datExpr[edges[[edge]],bwmoduleColors==which.module][,hubgeneID[[2]]],col="blue",lwd=3)
  }#end edge loop, keep device open!
  dev.off()
  figure<-figure+1
}
##--##end
#21. MODULE NODE ANNOTATION calculation####
#NEW module node annotations: difference between median values of ME for each module, and the q-value for the difference

#for loop: contrast 1 and contrast 2
all.statlist.modules<-list()
all.difflist.modules<-list()
pathwaylist=eval(parse(text=paste("list(",contrast.variable1,"=c.q$pathway1,",contrast.variable2,"=c.q$pathway2)",sep="")))
classifierlist=list(classifier1,classifier2)

for (pathW in 1:length(pathwaylist)){
  statlist.modules<-list()
  difflist.modules<-list()
  for (mod in 1:length(intModules)) {
    which.module=intModules[mod]
    ME=datME[, paste("ME",which.module, sep="")]
    
    statlist.modules[[mod]]<-wilcox.test(
      ME[classifierlist[[pathW]]==pathwaylist[[pathW]][[2]]],ME[classifierlist[[pathW]]==pathwaylist[[pathW]][[1]]]
    )$p.value
    
    difflist.modules[[mod]]<-
      median(ME[classifierlist[[pathW]]==pathwaylist[[pathW]][[2]]])-
      median(ME[classifierlist[[pathW]]==pathwaylist[[pathW]][[1]]])
  }
  
  tests.modules<-data.frame("module"=intModules,"diffME"=as.vector(unlist(difflist.modules)),"pval"=as.vector(unlist(statlist.modules)),"man_B"=as.vector(unlist(statlist.modules))*length(intModules),"bonf"=p.adjust(as.vector(unlist(statlist.modules)),method="bonferroni"),"BH"=p.adjust(as.vector(unlist(statlist.modules)),method="BH"))
  write.csv(tests.modules,file=file.path(dir.results,paste("Module_stats_",names(pathwaylist)[[pathW]],".",q.i,'.',an.count,'.',i,"_",names(datalist)[[i]],"_","_wilcox.csv",sep="")))
  all.statlist.modules[[pathW]]<-statlist.modules
  all.difflist.modules[[pathW]]<-difflist.modules
}

#redone by edge! for loop: contrast 1 and contrast 2
all.statlist.modules<-list()
all.difflist.modules<-list()
all.proplist.modules<-list()
all.sumlist.modules<-list()
pathwaylist=eval(parse(text=paste("list(",contrast.variable1,"=c.q$pathway1,",contrast.variable2,"=c.q$pathway2)",sep="")))
classifierlist=list(classifier1,classifier2)
for (pathW in 1:length(pathwaylist)){
  for (edge in 1:length(edges)){
    edgeset<-edges[[edge]]
    pathway<-pathways[[edge]]
    classifier<-classifiers[[edge]]
    statlist.modules<-list()
    difflist.modules<-list()
    proplist.modules<-list()
    sumlist.modules<-list()
    for (mod in 1:length(intModules)) {
      which.module=intModules[mod]
      ME=datME[edgeset,][, paste("ME",which.module, sep="")]
      
      statlist.modules[[mod]]<-wilcox.test(
        ME[classifier==pathway[[2]]],ME[classifier==pathway[[1]]]
      )$p.value
      
      difflist.modules[[mod]]<-
        median(ME[classifier==pathway[[2]]])-
        median(ME[classifier==pathway[[1]]])
      
      proplist.modules[[mod]]<-
        (sum(ME[classifier==pathway[[2]]]>0)/length(ME[classifier==pathway[[2]]]))/
        (sum(ME[classifier==pathway[[1]]]>0)/length(ME[classifier==pathway[[1]]]))
      
      sumlist.modules[[mod]]<-
        (sum(ME[classifier==pathway[[2]]])/length(ME[classifier==pathway[[2]]]))/
        (sum(ME[classifier==pathway[[1]]])/length(ME[classifier==pathway[[1]]]))
    }
    
    tests.modules.edge<-data.frame("module"=intModules,"sumME"=as.vector(unlist(sumlist.modules)),"propUpRatio"=as.vector(unlist(proplist.modules)),"diffME"=as.vector(unlist(difflist.modules)),"pval"=as.vector(unlist(statlist.modules)),"man_B"=as.vector(unlist(statlist.modules))*length(intModules),"bonf"=p.adjust(as.vector(unlist(statlist.modules)),method="bonferroni"),"BH"=p.adjust(as.vector(unlist(statlist.modules)),method="BH"))
    write.csv(tests.modules.edge,file=file.path(dir.results,paste("Module_stats_",names(pathwaylist)[[pathW]],".",q.i,'.',an.count,'.',i,"_",names(datalist)[[i]],"_","_wilcox_edge",edge,".csv",sep="")))
    all.statlist.modules[[edge]]<-statlist.modules
    all.difflist.modules[[edge]]<-difflist.modules
    all.proplist.modules[[edge]]<-proplist.modules
    all.sumlist.modules[[edge]]<-sumlist.modules
  }
}#end pathW

#some diffME results
diffMEmatrix<-matrix(unlist(all.difflist.modules),nrow=length(intModules))
colnames(diffMEmatrix)<-paste("edge",1:5,sep="_")
rownames(diffMEmatrix)<-intModules

pheatmap(diffMEmatrix,filename=paste('fig_',q.i,'.',an.count,'.',figure,'_heatmap_edgeMEs_diffME.pdf',sep=""),height=pdf.options()$height*1,width=pdf.options()$width*1)
figure<-figure+1
#23. Functional enrichment analysis in modules####

#get EntrezIDs on a per-modules basis
moduleList<-list()
for (colour in as.character(unique(geneInfo$moduleColor))){
  mgr<-as.character(geneInfo$EntrezID[geneInfo$moduleColor==colour])
  mgr2<-mgr[mgr!="NA"]
  mgr3<-mgr2[!is.na(mgr2)]
  moduleList[[colour]]<-mgr3
}

#For the following, need to also generate gene lists for the results.

#pathway enrichment
res<-"fail"
try(res <- compareCluster(moduleList, fun = "enrichPathway"))
if (class(res)!="character"){
  resplot<-plot(res)
  resplot2<-resplot+theme(axis.text.x = element_text(angle = 270, hjust = 0))
  pdf(paste('fig_',q.i,'.',an.count,'.',figure,'_Module clusters Reactome.pdf',sep=""),height=pdf.options()$height*1.3,width=pdf.options()$width*1.3)
  print(resplot2)
  dev.off()
  figure<-figure+1}

#KEGG pathway enrichment
res2<-"fail"
res2 <- try(compareCluster(moduleList, fun = "enrichKEGG"))
if (class(res2)!="character" & class(res2)!="try-error"){
  if(nrow(summary(res2))>1){#added this as code fails with only one result
    resplot2<-plot(res2)
    resplot3<-resplot2+theme(axis.text.x = element_text(angle = 270, hjust = 0))
    pdf(paste('fig_',q.i,'.',an.count,'.',figure,'_Module clusters KEGG.pdf',sep=""),height=pdf.options()$height*1.3,width=pdf.options()$width*1.3)
    print(resplot3)
    dev.off()
    figure<-figure+1}
}

#Disease Ontology  enrichment
res3<-"fail"
try(res3 <- compareCluster(moduleList, fun = "enrichDO"))
if (class(res3)!="character"){
  if (nrow(summary(res3))>1){
    resplot4<-plot(res3)
    resplot5<-resplot4+theme(axis.text.x = element_text(angle = 270, hjust = 0))
    pdf(paste('fig_',q.i,'.',an.count,'.',figure,'_Module clusters DO.pdf',sep=""),height=pdf.options()$height*1.3,width=pdf.options()$width*1.3)
    print(resplot5)
    dev.off()
    figure<-figure+1}
}

#GO enrichment
res4<-"fail"
try(res4 <- compareCluster(moduleList,OrgDb='org.Hs.eg.db', fun = "enrichGO"))
if (class(res4)!="character"){
  if (nrow(summary(res4))>1){
    resplot6<-plot(res4)
    resplot7<-resplot6+theme(axis.text.x = element_text(angle = 270, hjust = 0))
    pdf(paste('fig_',q.i,'.',an.count,'.',figure,'_Module_enrichGO.pdf',sep=""),height=pdf.options()$height*1.3,width=pdf.options()$width*1.3)
    print(resplot7)
    dev.off()
    figure<-figure+1}
}
#24. NETWORK 3: module-reactome######
#Raw data
modreaclist<-list()
for (modIter in 1:length(moduleList)){
  theres<-c()
  try(
    theres <- enrichPathway(gene = moduleList[[modIter]], pvalueCutoff = 0.5,readable = T))
  modname<-names(moduleList)[[modIter]]
  try(write.csv(summary(theres),file=file.path(dir.results,paste("Module_Reactome_",names(moduleList)[modIter],"_",q.i,'.',an.count,'.',i,"_",names(datalist),".csv",sep=""))))
  
  try(modreaclist[[modIter]]<-cbind(rep(modname,nrow(summary(theres))),summary(theres)))
}

modreacdf<-data.frame()
for (member in 1:length(modreaclist)){
  print(nrow(modreaclist[[member]]))  
  if (nrow(modreaclist[[member]])>1&&class(modreaclist[[member]])!="NULL"){
    modreacdf<-rbind(modreacdf,modreaclist[[member]])
  }  
}
names(modreacdf)[1]<-"module"
try(write.csv(modreacdf,file=file.path(dir.results,paste("Module_Reactome_All","_",q.i,'.',an.count,'.',i,"_",names(datalist)[[i]],"_",".csv",sep=""))))

modreacdf.2<-data.frame()
for (member in 1:length(modreaclist)){
  if (nrow(modreaclist[[member]])>1&&class(modreaclist[[member]])!="NULL"){
    modreacdf.2<-rbind(modreacdf.2,modreaclist[[member]])
  }  
}
names(modreacdf.2)[1]<-"module"

#now select only the top pathway for each identical list of enriched genes
enriched.genes<-unique(modreacdf.2$geneID)
modreacdf.3<-data.frame()
for (genegroup in 1:length(enriched.genes)){
  to.select<-modreacdf.2[which(modreacdf.2$geneID==enriched.genes[genegroup]),]
  if (nrow(to.select)==1) {
    modreacdf.3<-rbind(modreacdf.3,to.select)
  }else{select.top<-to.select[which(to.select$pvalue==min(to.select$pvalue)),][1,]#select first if more than one with min(pvalue)
  modreacdf.3<-rbind(modreacdf.3,select.top)
  }
}
#Now Multiple Testing Adjustment by module
modreacdf.4<-data.frame()
mod.in.df<-unique(modreacdf.3$module)
for (modu in 1:length(mod.in.df)){
  #get module pathways
  to.adjust<-modreacdf.3[which(modreacdf.3$module==mod.in.df[modu]),]
  #get the pavals and adjust
  oldP<-to.adjust$pvalue
  newQ<-p.adjust(oldP,method="BH")
  to.adjust<-cbind(to.adjust,newQ)
  #select based on adjusted pvals
  adjusted<-to.adjust[which(to.adjust$newQ<0.05),]
  if(nrow(adjusted)==0 && nrow(to.adjust)>3){adjusted<-to.adjust[1:3,]}
  if(nrow(adjusted)>10){adjusted<-adjusted[1:10,]}
  #write to df
  modreacdf.4<-rbind(modreacdf.4,adjusted)
}
#25. NETWORK 4: module-cell_expr ####

#Blood atlases
try(bloodResults1<-userListEnrichment(geneR=geneInfo$geneSymbol,labelR=geneInfo$moduleColor,nameOut=file.path(dir.results,"BloodEnrichment1.csv"),omitCategories="grey",useBloodAtlases=TRUE,outputCorrectedPvalues=FALSE))
try(bloodResults2<-userListEnrichment(geneR=geneInfo$geneSymbol,labelR=geneInfo$moduleColor,nameOut=file.path(dir.results,"BloodEnrichment2.csv"),omitCategories="grey",useBloodAtlases=TRUE,outputCorrectedPvalues=TRUE))



if(neoinsert=="on"){
  #NewHaemAtlas only
  try(bloodResults3<-userListEnrichment(geneR=geneInfo$substanceBXH,labelR=geneInfo$moduleColor,fnIn=file.path(dir.results,"haemAtlas.csv"),catNmIn="HaemAtlas",nameOut=file.path(dir.results,"BloodEnrichment_HA1.csv"),omitCategories="grey",outputCorrectedPvalues=FALSE)) 
  try(bloodResults4<-userListEnrichment(geneR=geneInfo$substanceBXH,labelR=geneInfo$moduleColor,fnIn=file.path(dir.results,"haemAtlas.csv"),catNmIn="HaemAtlas",nameOut=file.path(dir.results,"BloodEnrichment_HA2.csv"),omitCategories="grey",outputCorrectedPvalues=TRUE)) 
  # gi<-read.csv("/Users/armindeffur/Dropbox/Latest_Run_Combo/version_2016-02-08_12_59_25/Q_2_berry.train/5_WGCNA_4k_tv/results/table2.5.1.Square2Berry.5_WGCNA_4k_tvBlood_TBgeneInfoClass2.csv")
  # res2<-userListEnrichment(geneR=as.character(gi$substanceBXH),labelR=as.character(gi$moduleColor),fnIn="ha4.csv",nameOut="newHA.csv",omitCategories = "grey",outputCorrectedPvalues = FALSE)
}else{#OldHaemAtlas only
  try(bloodResults3<-userListEnrichment(geneR=geneInfo$geneSymbol,labelR=geneInfo$moduleColor,fnIn=file.path(dir.results,"haemAtlas_old.csv"),catNmIn="HaemAtlas",nameOut=file.path(dir.results,"BloodEnrichment_HA1.csv"),omitCategories="grey",outputCorrectedPvalues=FALSE)) 
  try(bloodResults4<-userListEnrichment(geneR=geneInfo$geneSymbol,labelR=geneInfo$moduleColor,fnIn=file.path(dir.results,"haemAtlas_old.csv"),catNmIn="HaemAtlas",nameOut=file.path(dir.results,"BloodEnrichment_HA2.csv"),omitCategories="grey",outputCorrectedPvalues=TRUE)) 
}

#Abbas only
try(bloodResults5<-userListEnrichment(geneR=geneInfo$geneSymbol,labelR=geneInfo$moduleColor,fnIn=file.path(dir.results,"abbasvals.csv"),catNmIn="Abbas",nameOut=file.path(dir.results,"BloodEnrichment_abbas1.csv"),omitCategories="grey",outputCorrectedPvalues=FALSE)) 
try(bloodResults6<-userListEnrichment(geneR=geneInfo$geneSymbol,labelR=geneInfo$moduleColor,fnIn=file.path(dir.results,"abbasvals.csv"),catNmIn="Abbas",nameOut=file.path(dir.results,"BloodEnrichment_abbas2.csv"),omitCategories="grey",outputCorrectedPvalues=TRUE)) 
#26. NETWORK 5: module-PalazzoloWang ##################

if (length(moduleList)<60){
  
  try(PWResults<-userListEnrichment(geneR=geneInfo$geneSymbol,labelR=geneInfo$moduleColor,nameOut=file.path(dir.results,"PalazzoWangEnrichment1.csv"),omitCategories="grey",usePalazzoloWang=TRUE))
  try(PWResults<-userListEnrichment(geneR=geneInfo$geneSymbol,labelR=geneInfo$moduleColor,nameOut=file.path(dir.results,"PalazzoWangEnrichment2.csv"),omitCategories="grey",usePalazzoloWang=TRUE,outputCorrectedPvalues=FALSE))
  overlaps<-PWResults[["sigOverlaps"]]
  
  modules.PW<-unique(overlaps$InputCategories)
  
  PalWang.select<-data.frame()#this is quite elaborate now
  for (pwmod in 1:length(modules.PW)){
    modOverlap<-overlaps[which(overlaps$InputCategories==modules.PW[pwmod]),]
    #print(nrow(modOverlap))
    qvalues<-p.adjust(modOverlap$Pvalues,method="BH")
    modOverlap<-cbind(modOverlap,qvalues)
    modOverlap.select<-modOverlap[which(modOverlap$qvalues<0.01),]#select significant hits after BH correction
    if(nrow(modOverlap.select)>10){modOverlap.select<-modOverlap.select[1:10,]}#only take top 10
    if(nrow(modOverlap.select)<4){modOverlap.select<-modOverlap[1:min(4,nrow(modOverlap)),]}#only take top hit regardless of significance
    PalWang.select<-rbind(PalWang.select,modOverlap.select)
  }
  
}#end NETWORK 5 if statement
#27. NETWORK 6: module-ImmunePW ######

try(immuneResults1<-userListEnrichment(geneR=geneInfo$geneSymbol,labelR=geneInfo$moduleColor,nameOut=file.path(dir.results,"ImmunePathwayEnrichment1.csv"),omitCategories="grey",useImmunePathwayLists=TRUE,outputCorrectedPvalues=TRUE))
try(immuneResults2<-userListEnrichment(geneR=geneInfo$geneSymbol,labelR=geneInfo$moduleColor,nameOut=file.path(dir.results,"ImmunePathwayEnrichment2.csv"),omitCategories="grey",useImmunePathwayLists=TRUE,outputCorrectedPvalues=FALSE))

setwd(dir.figures)

#29. ROC curves for modules as classifiers####
par(mfrow=c(1,1),xpd=T)
list1<-vector()
list2<-vector()
aiclist1<-vector()
aiclist2<-vector()
estimatelist1<-vector()
estimatelist2<-vector()
siglist1<-vector()
siglist2<-vector()

pdf(paste('fig_',q.i,'.',an.count,'.',figure,'_modAUC1_ROC_curves.pdf',sep=""),height=pdf.options()$height*1.3,width=pdf.options()$width*1)
for (rociter in 1:length(modules)){
  predobj<-prediction(MEs[,order(colnames(MEs))][rociter],classifier1,label.ordering=c(pathway1[[1]][1],pathway1[[2]][1]))
  perf<-performance(predobj,"tpr","fpr")
  #logistic regression
  outcome<-abs(as.numeric(as.factor(classifier1))-2)
  mydata<-cbind(outcome=outcome,MEs[,order(colnames(MEs))][rociter])
  mylogit <- glm(outcome ~ mydata[,2], data = mydata, family = "binomial")
  MEaic<-mylogit$aic
  MEestimate<-coef(mylogit)[2]
  MEsig<-coef(summary(mylogit))[,4][2]
  #end logit
  par(bg="white",mai=c(1.2,1.5,1,1),mar=c(12.5,2.5,2.5,2.5),xpd=T)
  if (rociter==1){
    plot(perf,lwd=2,col=sort(modules)[rociter],main=contrast.variable[[1]])} else {
      plot(perf,lwd=2,col=sort(modules)[rociter],main=contrast.variable[[1]],add=T)}
  print(paste(modules[rociter],max(performance(predobj,"auc")@y.values[[1]],1-performance(predobj,"auc")@y.values[[1]])))
  
  #populate statistics lists
  list1[rociter]<-max(performance(predobj,"auc")@y.values[[1]],1-performance(predobj,"auc")@y.values[[1]])
  aiclist1[rociter]<-MEaic
  estimatelist1[rociter]<-MEestimate
  siglist1[rociter]<-MEsig
}
names(list1)<-sort(modules)
par(parbackup)
legend(0,-0.2,sort(modules),fill=sort(modules),ncol=4)
dev.off()
figure<-figure+1

pdf(paste('fig_',q.i,'.',an.count,'.',figure,'_modAUC2_ROC_curves.pdf',sep=""),height=pdf.options()$height*1.3,width=pdf.options()$width*1)
for (rociter in 1:length(modules)){
  predobj<-prediction(MEs[,order(colnames(MEs))][rociter],classifier2,label.ordering=c(pathway2[[1]][1],pathway2[[2]][1]))
  perf<-performance(predobj,"tpr","fpr")
  #logistic regression
  outcome<-abs(as.numeric(as.factor(classifier2))-2)
  mydata<-cbind(outcome=outcome,MEs[,order(colnames(MEs))][rociter])
  mylogit <- glm(outcome ~ mydata[,2], data = mydata, family = "binomial")
  MEaic<-mylogit$aic
  MEestimate<-coef(mylogit)[2]
  MEsig<-coef(summary(mylogit))[,4][2]
  #end logit
  par(bg="white",mai=c(1.2,1.5,1,1),mar=c(12.5,2.5,2.5,2.5),xpd=T)
  if (rociter==1){
    plot(perf,lwd=2,col=sort(modules)[rociter],main=contrast.variable[[2]])} else {
      plot(perf,lwd=2,col=sort(modules)[rociter],main=contrast.variable[[2]],add=T)}
  print(paste(modules[rociter],max(performance(predobj,"auc")@y.values[[1]],1-performance(predobj,"auc")@y.values[[1]])))
  
  #populate statistics lists
  list2[rociter]<-max(performance(predobj,"auc")@y.values[[1]],1-performance(predobj,"auc")@y.values[[1]])
  aiclist2[rociter]<-MEaic
  estimatelist2[rociter]<-MEestimate
  siglist2[rociter]<-MEsig
  
}
names(list1)<-sort(modules)
par(parbackup)
legend(0,-0.2,sort(modules),fill=sort(modules),ncol=4)
dev.off()
figure<-figure+1

modules.as.classifier<-data.frame(cbind(list1,estimatelist1,aiclist1,siglist1,list2,estimatelist2,aiclist2,siglist2))
colnames(modules.as.classifier)<-c("modAUC1","logit_est1","AIC1","P_1","modAUC2","logit_est2","AIC2","P_2")
write.csv(modules.as.classifier,file.path(dir.results,"modules_as_classifiers_AUC_glm.csv"))

for (edge in 1:5){
  
  if (neoinsert=="on"){
    
    nodes<-RNeo4j::nodes
    #network-specific cypher query
    graphdata=modules.as.classifier
    
    
    #query = paste("MERGE (wmod:wgcna {name:{wgcnamod},square:{squareg},edge:{edgeg},contrastvar:{contrastvarsg},contrast:{contrastg}}) ON MATCH SET wmod.modAUC1 = ",modAUC1,", wmod.modAUC2 = ",modAUC2, " ON CREATE SET wmod.modAUC1 = ",modAUC1, ",wmod.modAUC2 = ",modAUC2,sep="")
    t = suppressMessages(newTransaction(graph))
    
    for (therow in 1:nrow(graphdata)) { 
      wmod = rownames(graphdata)[therow]
      modAUC1 = graphdata[therow, ]$modAUC1
      modAUC2 = graphdata[therow, ]$modAUC2
      #new params:
      logit_est1 = graphdata[therow, ]$logit_est1
      logit_est2 = graphdata[therow, ]$logit_est2
      AIC1 = graphdata[therow, ]$AIC1
      AIC2 = graphdata[therow, ]$AIC2
      P_1 = graphdata[therow, ]$P_1
      P_2 = graphdata[therow, ]$P_2
      
      query = paste("MERGE (wmod:wgcna {name:{wgcnamod},square:{squareg},edge:{edgeg},contrastvar:{contrastvarsg},contrast:{contrastg}}) ON MATCH SET wmod.modAUC1 = ",modAUC1,", wmod.modAUC2 = ",modAUC2,", wmod.logit_est1 = ",logit_est1,", wmod.logit_est2 = ",logit_est2,", wmod.AIC1 = ",AIC1,", wmod.AIC2 = ",AIC2,", wmod.P_1 = ",P_1,", wmod.P_2 = ",P_2," ON CREATE SET wmod.modAUC1 = ",modAUC1, ", wmod.modAUC2 = ",modAUC2,", wmod.logit_est1 = ",logit_est1,", wmod.logit_est2 = ",logit_est2,", wmod.AIC1 = ",AIC1,", wmod.AIC2 = ",AIC2,", wmod.P_1 = ",P_1,", wmod.P_2 = ",P_2,sep="")
      print(query)
      suppressMessages(appendCypher(t, 
                                    query, 
                                    wgcnamod = wmod,
                                    #modAUC1 = graphdata[therow, ]$list1,
                                    #modAUC2 = graphdata[therow, ]$list2,
                                    #invariant node properties
                                    squareg = dataset.variables[[1]],
                                    edgeg = edge,
                                    contrastvarsg = contrastvars[[edge]],
                                    contrastg = paste(pathways[[edge]][[2]],"-",pathways[[edge]][[1]])
                                    
                                    
      ))
    }
    
    print("committing :)")
    suppressMessages(commit(t))
    
    #END NEO4J
    
  }#end neoinsert if
}#end edgewise loop
#30. Module networks: annotation and small networks####
data(BloodLists)
cellkinds<-unique(BloodLists[,2])
cellgenelist<-list()
cellprobelist<-list()
for (cellkind in 1:length(cellkinds)){
  cellgenes<-BloodLists[,1][BloodLists[,2]==cellkinds[cellkind]]
  cellgenes2<-unique(as.character(unlist(sapply(cellgenes,function(x){unlist(strsplit(x," /// "))}))))
  cellgenelist[[cellkinds[cellkind]]]<-cellgenes2
  cellprobevec<-c()
  for (cellgene in cellgenes2){
    cellprobes<-as.character(geneInfo[as.character(geneInfo$geneSymbol)[!is.na(as.character(geneInfo$geneSymbol))]==cellgene,"substanceBXH"])
    cellprobevec<-c(cellprobevec,cellprobes)
  }
  cellprobelist[[cellkinds[cellkind]]]<-cellprobevec
}

#make small_mod_networks for all modules and annotate them# does not work well as it removes key genes...
modulesNotGrey<-modules[modules!="grey"]
for (module in 1:length(modulesNotGrey)){
  colour<-modulesNotGrey[module]
  modnet<-read.table(file.path(dir.results,paste('CS',q.i,'.',an.count,'.',i,'.',names(datalist)[[i]],"CS-edges-", colour, ".txt", sep="")),header=TRUE)
  smallmodnet<-modnet[modnet$weight>2*median(as.numeric(modnet$weight)),]#network edge filter
  smallmodnetprobes<-unique(as.character(smallmodnet[,1]))
  write.table(smallmodnet,file.path(dir.results,paste("smallNet_",colour,".txt",sep="")),row.names=F)
  smallmodgenes<-as.character(unique(smallmodnet$fromAltName))
  write.table(smallmodgenes,file.path(dir.results,paste("smallNet_",colour,"_genes.txt")),row.names=F)
  probedic<-unique(smallmodnet[,c(1,5)])
  
  #map genes to Reactome pathways
  modnetdata<-modreacdf[modreacdf$module==modulesNotGrey[module],]#some modules are not enriched for reactome pw
  if(nrow(modnetdata)>1){
    frame<-c("pathway","probe")
    for (rowx in 1:nrow(modnetdata)){
      components<-list(modnetdata[rowx,"Description"],unlist(strsplit(modnetdata[rowx,"geneID"],"/")))
      for (genex in components[[2]]){
        pwprobes0<-as.character(probedic[as.character(probedic$fromAltName)==genex,1])
        pwprobes<-pwprobes0[!is.na(pwprobes0)]
        if(length(pwprobes)>0){
          for (probex in 1:length(pwprobes)){
            frame<-rbind(frame,c(components[[1]],pwprobes[probex]))
          }
        }
      }
    }
    
    if(class(frame)!="character"){#for the singular case of zero entries into frame...
      if(nrow(frame)>1){#for the singular case of one pathway and one probe...
        frame2<-data.frame(frame[1:nrow(frame),])
        frame3<-frame2[2:nrow(frame2),]
        colnames(frame3)<-frame[1,]
        write.csv(frame3,file.path(dir.results,paste("smallNet_",colour,"_pw_probes.txt")),row.names=F)
      }
    }
  }
  
}#end small_mod_net creation and annotation

#clean up
rm(TOM)
collectGarbage()
gc()

#6_Pathway analysis*****************************************####



analysis="6_Pathway"

print(paste("Analysis:",analysis))

an.count=6
figure<-1

resultlistPW<-list(
  "edge1.500"=edge1.fac.TT.500,
  "edge2.500"=edge2.fac.TT.500,
  "edge3.500"=edge3.fac.TT.500,
  "edge4.500"=edge4.fac.TT.500,
  "edge5.500"=edge5.flat.TT.500)
resultlistPW_colour<-list(
  "edge1.12K"=edge1.fac.TT.12K,
  "edge2.12K"=edge2.fac.TT.12K,
  "edge3.12K"=edge3.fac.TT.12K,
  "edge4.12K"=edge4.fac.TT.12K,
  "edge5.12K"=edge5.flat.TT.12K)
entrezlistPW<-entrez.resultlist[c(6,7,8,9,10)]
resultlistPW_colour<-entrez.resultlist[c(11,12,13,14,15)]
#0. Start edgewise for loop####
for (edge in 1:length(edges)){
  edgeset<-edges[[edge]]
  
  ## Define folder for storing results NB do this for each analysis
  dir.results <- file.path(dir.output.version,paste("Q_",q.i,"_",dataset.variables,sep=""),analysis,paste("edge",edge,sep=""), "results")
  ## Define folder for saving figures
  dir.figures <- file.path(dir.output.version,paste("Q_",q.i,"_",dataset.variables,sep=""),analysis,paste("edge",edge,sep=""), "figures")
  
  for (dir in c(dir.figures, dir.results)) {
    if (!file.exists(dir)) {
      dir.create(dir, recursive=T,showWarnings=T)
    }
  }
  
  for (dir in c(dir.figures, dir.results)) {
    if (!file.exists(dir)) {
      dir.create(dir, recursive=T,showWarnings=T)
    }
  }
  setwd(dir.figures)
  #1. Data (get from results above, used for colouring, not for pathway calculation)####
  r1<-resultlistPW[[edge]]
  r3<-r1
  r4<-resultlistPW_colour[[edge]]
  #2. Preparation####
  data(kegg.gs)
  data(kegg.sets.hs)
  data(kegg.gs.dise)
  data(carta.hs)
  data(go.gs)
  data(egSymb)
  
  #define empty pd list
  pd<-list()
  pd[[i]]<-pData(datalist[[i]][,edgeset])
  
  #subset
  ed0<-exprs(data.q[,edgeset])
  
  kegg.gs.sym<-lapply(kegg.gs,eg2sym) #pathways only
  kegg.sets.hs.sym<-lapply(kegg.sets.hs,eg2sym)# this contains both pathways and diseases
  kegg.gs.dise.sym<-lapply(kegg.gs.dise,eg2sym) #diseases only
  carta.hs.sym<-lapply(carta.hs,eg2sym)
  
  inkegg<-unique(unlist(kegg.gs.sym))
  inkegg.dise<-unique(unlist(kegg.gs.dise.sym))
  incarta<-unique(unlist(carta.hs.sym))
  
  inkegg.nuID<-targetID2nuID(inkegg,lib.mapping="lumiHumanIDMapping")[,7]
  inkegg.dise.nuID<-targetID2nuID(inkegg.dise,lib.mapping="lumiHumanIDMapping")[,7]
  incarta.nuID<-targetID2nuID(incarta,lib.mapping="lumiHumanIDMapping")[,7]
  
  inkeggbool<-is.element(rownames(ed0),inkegg.nuID)
  inkeggdisebool<-is.element(rownames(ed0),inkegg.dise.nuID)
  incartabool<-is.element(rownames(ed0),incarta.nuID)
  
  ed<-ed0[inkeggbool,]
  ed.dise<-ed0[inkeggdisebool,]
  ed.carta<-ed0[incartabool,]
  
  inboth<-nuID2targetID(row.names(ed),lib.mapping="lumiHumanIDMapping")
  inbothdise<-nuID2targetID(row.names(ed.dise),lib.mapping="lumiHumanIDMapping")
  inbothcarta<-nuID2targetID(row.names(ed.carta),lib.mapping="lumiHumanIDMapping")
  
  row.names(ed)<-inboth
  row.names(ed.dise)<-inbothdise
  row.names(ed.carta)<-inbothcarta
  
  #cases and controls: index vectors:
  #positions in contrast vector where value==ref and value==samp
  ref<-which(
    eval(parse(text=paste("datalist[[i]][,edgeset]$",contrastvars[[edge]],sep="")))
    ==pathways[[edge]][[1]])
  samp<-which(
    eval(parse(text=paste("datalist[[i]][,edgeset]$",contrastvars[[edge]],sep="")))
    ==pathways[[edge]][[2]])
  
  #the analysis is limited to genes that are found in KEGG, KEGG_disease, and BioCarta
  setwd(dir.figures)
  #KEGG
  kegg.p<-gage(ed,gsets=kegg.gs.sym,ref=ref,samp=samp,same.dir=T,compare="unpaired")
  write.csv(kegg.p,file=file.path(dir.results,paste(names(datalist)[[i]],"KEGG_Results.csv",sep="")))
  sigGS<-sigGeneSet(kegg.p,outname="KEGG",cexRow=.7,cexCol=1,mar=c(10,15),pdf.size = c(10,7))
  GSdata<-cbind(rownames(sigGS$stats),sprintf("%.3f",sigGS$stats[,1]))
  colnames(GSdata)<-c("Pathway","Result")
  write.csv(GSdata,file=file.path(dir.results,paste(names(datalist)[[i]],"KEGG_sig_Results.csv",sep="")))
  
  #KEGG disease
  kegg.p.dise<-gage(ed.dise,gsets=kegg.gs.dise.sym,ref=ref,samp=samp,same.dir=T,compare="unpaired")
  write.csv(kegg.p.dise,file=file.path(dir.results,paste(names(datalist)[[i]],"KEGG_DIS_Results.csv",sep="")))
  sigGS.dise<-sigGeneSet(kegg.p.dise,outname="KEGG_DIS",cexRow=.7,cexCol=1,mar=c(10,15),pdf.size = c(10,7))
  GSdata.dise<-cbind(rownames(sigGS.dise$stats),sprintf("%.3f",sigGS.dise$stats[,1]))
  colnames(GSdata.dise)<-c("Pathway","Result")
  write.csv(GSdata.dise,file=file.path(dir.results,paste(names(datalist)[[i]],"KEGG_DIS_sig_Results.csv",sep="")))
  
  #Biocarta
  carta.p<-gage(ed.carta,gsets=carta.hs.sym,ref=ref,samp=samp,same.dir=T,compare="unpaired")
  write.csv(carta.p,file=file.path(dir.results,paste(names(datalist)[[i]],"BioCarta_Results.csv",sep="")))
  sigGS.carta<-sigGeneSet(carta.p,outname="BioCarta",cexRow=.7,cexCol=1,mar=c(10,15),pdf.size = c(10,7))
  GSdata.carta<-cbind(rownames(sigGS.carta$stats),sprintf("%.3f",sigGS.carta$stats[,1]))
  colnames(GSdata.carta)<-c("Pathway","Result")
  write.csv(GSdata.carta,file=file.path(dir.results,paste(names(datalist)[[i]],"BioCarta_sig_Results.csv",sep="")))
  
  #xt<-xtable(GSdata,align=c("|p{1cm}","|p{9.5cm}|","p{1.5cm}|"))
  #print(xt,hline.after=-1:nrow(sigGS$stats),include.colnames=T,include.rownames=F,size="small",file="topPathways.tex")
  
  sigGS.KEGG<-sigGeneSet(kegg.p,outname=file.path(dir.results,"KEGG_hm"))
  sigGS.KEGG.dise<-sigGeneSet(kegg.p.dise,outname=file.path(dir.results,"KEGG_DIS_hm"))
  sigGS.carta<-sigGeneSet(carta.p,outname=file.path(dir.results,"BioCarta_hm"))
  
  genedata=r4$logFC
  names(genedata)<-r4$SYMBOL
  #names(genedata)<-r4$EntrezID
  if(is.null(names(genedata))){names(genedata)<-r4$Symbol}
  if(!is.null(r4$TargetID)){names(genedata)<-r4$TargetID}
  
  pw.KEGG<-as.character(sapply(row.names(sigGS.KEGG$stats),function(x){substr(x,1,8)}))
  pw.KEGG.dise<-as.character(sapply(row.names(sigGS.KEGG.dise$stats),function(x){substr(x,1,8)}))
  pw.carta<-as.character(sapply(row.names(sigGS.carta$stats),function(x){substr(x,1,8)}))
  
  data(gene.idtype.list)
  data("gene.idtype.bods")
  data(bods)
  
  setwd(dir.figures)
  
  names(pw.KEGG)<-row.names(sigGS.KEGG$stats)
  names(pw.KEGG.dise)<-row.names(sigGS.KEGG.dise$stats)
  
  pw.collection<-list(pw.KEGG,pw.KEGG.dise)
  names(pw.collection)<-c("KEGG","KEGG_DIS")
  
  #Redefine internal functions with conflicts
  getType<-KEGGgraph::getType
  nodes<-graph::nodes
  namecounter<-1
  
  #This crashes: 3- 8
  #3. Start loop over KEGG and KEGG_disease pathway collections####
  for(pwx in pw.collection){
    #ensure that all pathways are available vie KEGGREST api
    #pwx is unfiltered list
    #pw is filtered by pathway availability
    selvecKEGG<-logical()
    for (itera in 1:length(pwx)){
      a<-try(keggGet(pwx[itera]))
      if(class(a)=="try-error"){
        selvecKEGG[itera]<-FALSE
      }else{selvecKEGG[itera]<-TRUE}
    }
    pw<-pwx[which(selvecKEGG)]
    #4. Pathview####
    if (length(pw)>1){
        for (iterator in 1:length(pw)) {
        try(pathview(gene.data=genedata,pathway.id=as.character(pw[iterator]),species="hsa",gene.idtype="SYMBOL",kegg.dir=dir.kegg))
      }
      #5. KEGG_WGCNA enrichment heatmaps####Replace by Jaccard index??? todo
      mmat<-matrix(data=NA,nrow=length(modules),ncol=length(pw))
      for (ii in 1:length(pw)) {
        #overlap of pathways and modules (move into pathview loop when ready)
        pwg<-as.character(node.info(file.path(dir.kegg,paste(pw[ii],".xml",sep="")))$labels)
        mvec<-vector()
        for (j in 1:length(modgeneslist)){
          #m<-sum(modgeneslist[[j]]%in%pwg)/length(pwg);
          m<-sum(modgeneslist[[j]]%in%pwg)/length(modgeneslist[[j]]);
          #m<-sum(modgeneslist[[j]]%in%pwg);
          if (!is.na(m>0)){
            mvec<-c(mvec,m)}
          else mvec<-c(mvec,0)#each mvec is the enrichment for one pathway
        }
        mmat[,ii]<-mvec
      }
      rownames(mmat)<-modules
      colnames(mmat)<-pw
      colf=colorRampPalette(c("white","blue"))
      
      pdf(paste('fig_',q.i,'.',an.count,'.',figure,'_edge_',edge,'_pathwayEnrichMod_meth1.pdf',sep=""),height=pdf.options()$height*1.3,width=pdf.options()$width*1)
      heatmap_ad(mmat,scale="none",col=colf(40),main=paste("Pathway enrichment for modules",names(pw.collection)[[namecounter]],names(datalist)[[i]]),margins=c(7,7),cex.main=.8)
      dev.off()
      figure<-figure+1
      
      mmat2<-matrix(data=NA,nrow=length(modules),ncol=length(pw))
      for (iii in 1:length(pw)) {
        #overlap of pathways and modules (move into pathview loop when ready)
        pwg<-as.character(node.info(file.path(dir.kegg,paste(pw[iii],".xml",sep="")))$labels)
        mvec<-vector()
        for (j in 1:length(modgeneslist)){
          m<-sum(pwg%in%modgeneslist[[j]])/length(pwg);
          #m<-sum(pwg%in%modgeneslist[[j]])/length(modgeneslist[[j]]);
          #m<-sum(pwg%in%modgeneslist[[j]]);
          mvec<-c(mvec,m)#each mvec is the enrichment for one pathway
        }
        mmat2[,iii]<-mvec
      }
      rownames(mmat2)<-modules
      colnames(mmat2)<-pw
      colf=colorRampPalette(c("white","blue"))
      pdf(paste('fig_',q.i,'.',an.count,'.',figure,'_edge_',edge,'_pathwayEnrichMod_meth2.pdf',sep=""),height=pdf.options()$height*1.3,width=pdf.options()$width*1)
      heatmap_ad(mmat2,scale="none",col=colf(40),main=paste("Pathway enrichment for modules",names(pw.collection)[[namecounter]],names(datalist)[[i]]),margins=c(7,7),cex.main=.8)
      dev.off()
      figure<-figure+1
      #6. Pathenrich: normalised pathway enrichment for modules####
      if("grey"%in%modules)
      {mmat3<-matrix(data=NA,nrow=length(modules)-1,ncol=length(pw))#
      }else{mmat3<-matrix(data=NA,nrow=length(modules),ncol=length(pw))}#not grey}
      pwlength<-list()
      for (iv in 1:length(pw)) {
        #overlap of pathways and modules (move into pathview loop when ready)
        pwg<-as.character(node.info(file.path(dir.kegg,paste(pw[iv],".xml",sep="")))$labels)
        mvec<-vector()
        allprobes<-as.character(modinfo$substanceBXH)
        allgenes<-as.character(modinfo$geneSymbol)
        enrichlist<-list()
        genesInMod<-list()
        enrichnormlist<-list()
        pwlength[[iv]]<-length(pwg)
        for (jj in 1:length(modgeneslist)){
          #m<-sum(pwg%in%modgeneslist[[j]])/length(pwg);
          #m<-sum(pwg%in%modgeneslist[[j]])/length(modgeneslist[[j]]);
          #m<-sum(pwg%in%modgeneslist[[j]]);
          #mvec<-c(mvec,m)#each mvec is the enrichment for one pathway
          pwl<-length(pwg)
          enrich<-sum(pwg%in%modgeneslist[[jj]])/length(pwg)
          #print(paste(names(modgeneslist)[jj],":",enrich))
          enrichlist[[names(modgeneslist)[jj]]]<-enrich
          
          enrichnorm=enrich*12000/lapply(modgeneslist,length)[[jj]]
          enrichnormlist[[jj]]<-enrichnorm
          
          if(names(modgeneslist)[[jj]]!="grey"){
            mvec<-c(mvec,enrichnorm)}#each mvec is the enrichment for one pathway
          #genes<-allgenes[allprobes%in%moduleprobes[[iter]]]
          #genesInMod[[iter]]<-genes
          #names(genesInMod)<-modulecolors
          names(enrichnormlist)[[jj]]<-names(modgeneslist)[[jj]]
          sum(unlist(enrichlist))
          
          
        }
        mmat3[,iv]<-mvec
        pdf(paste('fig_',q.i,'.',an.count,'.',figure,'_edge_',edge,'module_pathway_enrichment_',pw[iv],'.pdf',sep=""),height=pdf.options()$height*1.3,width=pdf.options()$width*1)
        par(mfrow=c(3,1),mar=c(9,6,3,3))
        barplot(unlist(lapply(modgeneslist,length)),las=2,main="Module size (N genes)",col=names(modgeneslist))
        barplot(unlist(enrichlist),las=2,main=paste("Proportion of", names(pw)[iv] ,"genes in module"),col=names(modgeneslist))#fixed bug!
        barplot(unlist(enrichnormlist),las=2,main="Ratio of real to expected pathway gene enrichment",col=names(modgeneslist))
        abline(h=1,col="red")
        dev.off()
        par(parbackup)
        figure<-figure+1
        
      }
      #7. Heatmap based on pathenrich####
      truthvector<-lapply(modgeneslist,length)>1
      #rownames(mmat3)<-names(modgeneslist[truthvector])#here is a bug
      rownames(mmat3)<-modules[modules!="grey"]#here is a bug
      colnames(mmat3)<-pw
      colf=colorRampPalette(c("white","blue","red"))
      pdf(paste('fig_',q.i,'.',an.count,'.',figure,'_edge_',edge,'_pathwayEnrichMod_meth3.pdf',sep=""),height=pdf.options()$height*1.3,width=pdf.options()$width*1)
      heatmap_ad(mmat3,scale="none",col=colf(40),main=paste("Pathway enrichment for modules",names(pw.collection)[[namecounter]],names(datalist)[[i]]),margins=c(7,7),cex.main=.8)
      dev.off()
      figure<-figure+1
      
    } else {print("No significant pathways to show")}
    namecounter<-namecounter+1
    #8. Close pw.collection loop####
  }#close pw.collection loop
  #9. GO enrichment for top 500 de probes for the specified comparison (edge)####
  entrezIDs<-nuID2EntrezID(featureNames(datalist[[i]]),"lumiHumanIDMapping")
  #hei is haveEntrezID
  hei<-unique(as.character(entrezIDs[sapply(entrezIDs,function(x) x!="")]))
  
  genelistGO<-unique(entrezlistPW[[edge]]$EntrezID[!is.na(entrezlistPW[[edge]]$EntrezID)])
  
  p1.MF<-basicProfile(genelistGO,onto="MF",level=2,orgPackage="org.Hs.eg.db")
  p1.BP<-basicProfile(genelistGO,onto="BP",level=2,orgPackage="org.Hs.eg.db")
  p1.CC<-basicProfile(genelistGO,onto="CC",level=2,orgPackage="org.Hs.eg.db")
  
  printProfiles(p1.MF,percentage=T)
  printProfiles(p1.BP,percentage=T)
  printProfiles(p1.CC,percentage=T)
  
  pdf(paste('fig_',q.i,'.',an.count,'.',figure,'_edge_',edge,'_GO_MF_EnrichMod.pdf',sep=""),height=pdf.options()$height*1.3,width=pdf.options()$width*1)
  plotProfiles(p1.MF,aTitle=paste("DE genes",names(datalist)[[i]]),labelWidth=30)
  dev.off()
  figure<-figure+1
  
  pdf(paste('fig_',q.i,'.',an.count,'.',figure,'_edge_',edge,'_GO_BP_EnrichMod.pdf',sep=""),height=pdf.options()$height*1.3,width=pdf.options()$width*1)
  plotProfiles(p1.BP,aTitle=paste("DE genes",names(datalist)[[i]]),labelWidth=30)
  dev.off()
  figure<-figure+1
  
  pdf(paste('fig_',q.i,'.',an.count,'.',figure,'_edge_',edge,'_GO_CC_EnrichMod.pdf',sep=""),height=pdf.options()$height*1.3,width=pdf.options()$width*1)
  plotProfiles(p1.CC,aTitle=paste("DE genes",names(datalist)[[i]]),labelWidth=30)
  dev.off()
  figure<-figure+1
  #10. ReactomePA (max 1000 genes)####
  
  # make list of significantly DE genes and top 500 (max) DE genes for each cell subset
  genes.for.pathway.list=list(rownames(resultlistPW[[edge]]))
  names(genes.for.pathway.list)<-paste("Top500_edge",edge,sep="")
  
  setnames=names(genes.for.pathway.list)
  for (geneset in 1:length(genes.for.pathway.list)){
    
    
    rgenes<-genes.for.pathway.list[[geneset]]
    rgenesE<-data.frame(nuID2EntrezID(as.character(rgenes),filterTh = NULL,lib.mapping='lumiHumanIDMapping', returnAllInfo = TRUE))
    rgenesE2<-as.character(rgenesE$EntrezID)[as.character(rgenesE$EntrezID)!=""]
    # fc<-2^reslistTTsig[[i]]$logFC
    # 
    if (length(rgenesE2)>1000){
      rgenesE2<-rgenesE2[1:1000]#using all genes precipitates a bug in enrichPathway
      #   fc<-fc[1:1000]
    }
    rpa<-NA
    rpa <- try(enrichPathway(unique(rgenesE2), organism="human",pvalueCutoff = 0.4,readable = T))
    #head(summary(rpa))
    summary(rpa)
    
    #catch error conditions (NULL and NA content fpr rpa) and export result. This is subtle, and distinguishes between NULL and NA conditions
    #edit: added another error condition: rpa may succeed, but summary(rpa) is empty, I have no idea why
    
    succ<-TRUE
    if (length(rpa)==0){succ<-FALSE;print(paste(setnames[geneset],"rpa is of zero length"))} else if (is.na(rpa)){succ<-FALSE;print(paste(setnames[geneset],"rpa is NA"))} else if (class(rpa)=="try-error"){succ<-FALSE;print(paste(setnames[geneset],"rpa gives a try error"))} else if (nrow(summary(rpa))==0){succ<-FALSE;print(paste(setnames[geneset],"summary(rpa) is of zero length"))} else if (nrow(summary(rpa))==1){succ<-TRUE;print(paste(setnames[geneset],"nrow summary rpa is of length 1"))}#edited last condition
    
    if (succ) {
      write.csv(summary(rpa),file=file.path(dir.results,
                                            paste('table',q.i,'.',an.count,'.',i,"_",names(questions)[[q.i]],'_',analysis,names(datalist)[i],"Set_",names(genes.for.pathway.list)[geneset],"ReactomePA.csv",sep="")
      ))
      if(nrow(summary(rpa))>1){
        pdf(paste('fig_',q.i,'.',an.count,'.',figure,'_edge_',edge,'Set_',names(genes.for.pathway.list)[geneset],'reactomePA.pdf',sep=""),height=pdf.options()$height*1.3,width=pdf.options()$width*1.3)  
        theplot<-barplot(rpa,showCategory=8)
        print(theplot)
        dev.off()
        figure<-figure+1
      }#end if
      pdf(paste('fig_',q.i,'.',an.count,'.',figure,'_edge_',edge,'Set_',names(genes.for.pathway.list)[geneset],'cnetplot.pdf',sep=""),height=pdf.options()$height*1.3,width=pdf.options()$width*1)
      cnetplot(rpa, categorySize = "pvalue")#, foldChange = fc) # no FC info
      dev.off()
      figure<-figure+1
    } else {print(paste("no result for rpa: ",setnames[geneset]))}
    
  }#close geneset list
  #11. End edgewise pathway for loop####
}#end edgewise pathway for loop

#end pathway analysis

#7_CEGS modules********************************************####

analysis="7_CEGS"
an.count=7
figure<-1
#make relevant output folders
## Define folder for storing results NB do this for each analysis
dir.results <- file.path(dir.output.version,paste("Q_",q.i,"_",dataset.variables,sep=""),analysis, "results")
## Define folder for saving figures
dir.figures <- file.path(dir.output.version,paste("Q_",q.i,"_",dataset.variables,sep=""),analysis, "figures")

for (dir in c(dir.figures, dir.results)) {
  if (!file.exists(dir)) {
    dir.create(dir, recursive=T,showWarnings=T)
  }
}
setwd(dir.figures)

#make user list for later
data("HaemAtlas")
ha.2<-HaemAtlas[,c("SYMBOL","Type")]
colnames(ha.2)<-c("GENES","LISTS")
write.csv(ha.2,file.path(dir.results,"haemAtlas_old.csv"),row.names=F)
write.csv(ha.2,file.path(dir.results,"haemAtlas.csv"),row.names=F)#debug

if(neoinsert=="on"){
  ##new
  data(HaemAtlas)
  ha2<-HaemAtlas
  hanu<-HaemAtlas$nuID
  news<-nuID2targetID(hanu)
  ha2$news<-news
  ha2<-unique(ha2)
  
  # hanu2<-paste(hanu,collapse="','")
  # #need to add in platform specification
  # q<-paste("MATCH (p:PROBETYPE)-[r]-(s:SYMBOL) WHERE p.name IN ['",hanu2,"'] RETURN p.name AS nuID2, s.name AS SYMBOL2")
  # try(res<-cypher(graph,q))
  # ha3<-merge(ha2,res,by.x="nuID",by.y="nuID2",all.x=FALSE)
  # ha4<-ha3[,c("nuID","Type")]
  
  ha4<-ha2[,c("nuID","Type")]
  colnames(ha4)<-c("GENES","LISTS")
  ha4<-unique(ha4[order(ha4$LISTS),])
  write.csv(ha4,file.path(dir.results,"haemAtlas.csv"),row.names=F)
  ##\new
}

abbasvals<-cbind(as.character(marknames(Abbassym)),names(marknames(Abbassym)))
colnames(abbasvals)<-c("GENES","LISTS")
write.csv(abbasvals,file.path(dir.results,"abbasvals.csv"),row.names=F)
#1. start edgewise for loop####
diffExBaylor.edges<-list()
for (edge in 1:length(edges)){
  #2. Prepare DE analysis within each module####
  edgeset<-edges[[edge]]
  up.res.list<-list()
  down.res.list<-list()
  all.res.list<-list()
  alldata<-data.q[,edgeset]
  
  #don't re-use contrast and design matrices from DE
  
  TS<-factor(classifiers[[edge]])
  design<-model.matrix(~0+TS)
  colnames(design)<-levels(TS)
  
  cont.matrix<-makeContrasts(
    contrast.list=eval(parse(text=paste(pathways[[edge]][[2]],"-",pathways[[edge]][[1]],sep=""))),
    levels=design
  )
  #3. In each module (2008) do DE analysis and send result to module matrix####
  for (module in 1:length(modlist4)){
    #select module probes
    datasymbols<-probemap$SYMBOL
    if(is.null(datasymbols)){datasymbols<-fData(alldata)$TargetID}
    data<-alldata[datasymbols%in%modlist4[[module]],]
    #DE
    #Linear model: tb status
    fit.m<-lmFit(data,design)
    fit2.m<-contrasts.fit(fit.m,cont.matrix)
    fit2.m<-eBayes(fit2.m)
    res=topTable(fit2.m,number=Inf,adjust.method="none",p.value=0.05)
    up.probes<-res[res$logFC>0,]
    up.res<-nrow(up.probes)/nrow(exprs(data))
    down.probes<-res[res$logFC<0,]
    down.res<-nrow(down.probes)/nrow(exprs(data))
    all.res<-up.res-down.res
    up.res.list[[module]]<-up.res
    down.res.list[[module]]<-down.res
    all.res.list[[module]]<-all.res
    
    #within module correlation analysis
    #cor(data) or the transpose
  }
  #4. Generate module matrix####
  hmdata<-rep(NA,33)
  hmdata[c(1:8,12:20,22:31)]<-unlist(all.res.list)
  hmmatrix<-t(matrix(hmdata,ncol=3))
  
  pdf(paste('fig_',q.i,'.',an.count,'.',figure,'_edge_',edge,'CEGS_Matrix1.pdf',sep=""),height=pdf.options()$height*1,width=pdf.options()$width*1)
  labeledHeatmap(hmmatrix,xLabels=as.character(1:11),yLabels=c("M1","M2","M3"),colors=(blueWhiteRed(30)),zlim = c(-1,1),main=names(datalist)[i])
  dev.off()
  figure<-figure+1
  
  #correlate all module probes(----------)todo)
  #5. In each module (2013) do DE analysis and send result to module matrix####
  #make a scaffold to fit modules in
  for (module in 1:length(modlist4.2)){
    #select module probes
    data<-alldata[datasymbols%in%modlist4.2[[module]],]
    #DE
    #Linear model: tb status
    fit.m<-lmFit(data,design)
    fit2.m<-contrasts.fit(fit.m,cont.matrix)
    fit2.m<-eBayes(fit2.m)
    res=topTable(fit2.m,number=Inf,adjust.method="none",p.value=0.05)
    up.probes<-res[res$logFC>0,]
    up.res<-nrow(up.probes)/nrow(exprs(data))
    down.probes<-res[res$logFC<0,]
    down.res<-nrow(down.probes)/nrow(exprs(data))
    all.res<-up.res-down.res
    up.res.list[[module]]<-up.res
    down.res.list[[module]]<-down.res
    all.res.list[[module]]<-all.res
    
    #within module correlation analysis
    #cor(data) or the transpose
  }
  #6. Generate module matrix####
  hmdata<-rep(NA,260)
  hmdata<-unlist(all.res.list)
  dimensions<-strsplit(names(modlist4.2),".",fixed=T)
  coords.rows=unlist(lapply(dimensions,function(x){strsplit(x[1],"M",fixed=T)[[1]][2]}))
  coords.cols=unlist(lapply(dimensions,function(x){x[2]}))
  coords.all=cbind(as.numeric(coords.rows),as.numeric(coords.cols))
  hmmatrix.2013<-matrix(data=NA,nrow=9,ncol=max(coords.all[,2]))
  for (rows in 1:nrow(coords.all)){
    hmmatrix.2013[coords.all[rows,][1],coords.all[rows,][2]]<-hmdata[rows]
    
  }
  diffExBaylor<-cbind(as.character(names(modlist4.2)),unlist(all.res.list))
  diffExBaylor.edges[[edge]]<-diffExBaylor
  pdf(paste('fig_',q.i,'.',an.count,'.',figure,'_edge_',edge,'ChaussMatrix2.pdf',sep=""),height=pdf.options()$height*1,width=pdf.options()$width*1)
  labeledHeatmap(hmmatrix.2013,xLabels=1:ncol(hmmatrix.2013),yLabels=paste("M.",unique(coords.rows),sep=""),colors=(blueWhiteRed(30)),zlim = c(-1,1),main=names(datalist)[i])
  dev.off()
  figure<-figure+1
  
  cbind(names(modlist4.2),unlist(all.res.list))
}#end edgewise for loop###
#7. Baylor heatmaps####
val.list<-lapply(diffExBaylor.edges,function(x){as.numeric(x[,2])})
matBaylor<-matrix(as.numeric(unlist(val.list)),nrow=260)
colnames(matBaylor)<-paste("edge",1:5,sep="")
rownames(matBaylor)<-as.character(names(modlist4.2))

if(nrow(matBaylor[which(apply(matBaylor,1,sum)!=0),])>10){
  pheatmap(matBaylor[which(apply(matBaylor,1,sum)!=0),],main="non-zero rows all edges",filename=paste('fig_',q.i,'.',an.count,'.',figure,'_edge_',edge,'ChaussMatrix2_heatmap_nonzero.pdf',sep=""),height=pdf.options()$height*1,width=pdf.options()$width*1)
  figure<-figure+1
}

if(nrow(matBaylor[which(apply(matBaylor,1,function(x){max(abs(x))})>=0.2),])>10){
  pheatmap(matBaylor[which(apply(matBaylor,1,function(x){max(abs(x))})>=0.2),],main="at least 0.2 DE rows all edges",filename=paste('fig_',q.i,'.',an.count,'.',figure,'_edge_',edge,'ChaussMatrix2_heatmap_20pc.pdf',sep=""),height=pdf.options()$height*1,width=pdf.options()$width*1)
  figure<-figure+1
}
#8. NETWORK 8: Baylor module-reactome######

modNU<-list()
#Gene names for each module
for (mod in mod.names2){
  #data<-modgenes2[modgenes2[,1]==mod,7]
  dataBase<-illuminaHumanv2NUID[mappedkeys(illuminaHumanv2NUID)]
  data<-as.character(unlist(mget(as.character(modgenes2[modgenes2[,1]==mod,7]),dataBase)))
  modNU[[mod]]<-data
}

modEN<-lapply(lapply(lapply(modNU,nuID2EntrezID,lib.mapping="lumiHumanIDMapping"),as.character),function(x){x[which(x!="")]})

Baylor_modreaclist<-list()
for (modIter in 1:length(modEN)){
  theres<-c()
  try(
    theres <- enrichPathway(gene = modEN[[modIter]], pvalueCutoff = 0.5,readable = T))
  modname<-names(modEN)[[modIter]]
  #print(modname)
  #print(head(summary(theres)))
  try(write.csv(summary(theres),file=file.path(dir.results,paste("Baylor-Module_Reactome_",names(modEN)[modIter],"_",q.i,'.',an.count,'.',i,"_",names(datalist)[[i]],".csv",sep=""))))
  
  try(Baylor_modreaclist[[modIter]]<-cbind(rep(modname,nrow(summary(theres))),summary(theres)))
}

Baylor_modreacdf<-data.frame()
for (member in 1:length(Baylor_modreaclist)){
  #print(nrow(Baylor_modreaclist[[member]]))  
  if (nrow(Baylor_modreaclist[[member]])>1&&class(Baylor_modreaclist[[member]])!="NULL"){
    Baylor_modreacdf<-rbind(Baylor_modreacdf,Baylor_modreaclist[[member]])
  }  
}
names(Baylor_modreacdf)[1]<-"module"
try(write.csv(Baylor_modreacdf,file=file.path(dir.results,paste("Baylor-Module_Reactome_All","_",q.i,'.',an.count,'.',i,"_",names(datalist)[[i]],"_",".csv",sep=""))))

Baylor_modreacdf.2<-data.frame()
for (member in 1:length(Baylor_modreaclist)){
  print(nrow(Baylor_modreaclist[[member]]))  
  if (nrow(Baylor_modreaclist[[member]])>1&&class(Baylor_modreaclist[[member]])!="NULL"){
    Baylor_modreacdf.2<-rbind(Baylor_modreacdf.2,Baylor_modreaclist[[member]])
  }  
}
names(Baylor_modreacdf.2)[1]<-"module"

#now select only the top pathway for each identical list of enriched genes
Baylor_enriched.genes<-unique(Baylor_modreacdf.2$geneID)
Baylor_modreacdf.3<-data.frame()
for (genegroup in 1:length(Baylor_enriched.genes)){
  to.select<-Baylor_modreacdf.2[which(Baylor_modreacdf.2$geneID==Baylor_enriched.genes[genegroup]),]
  #print(nrow(to.select))
  if (nrow(to.select)==1) {
    Baylor_modreacdf.3<-rbind(Baylor_modreacdf.3,to.select)
  }else{
    select.top<-to.select[which(to.select$pvalue==min(to.select$pvalue)),][1,]#select first if more than one with min(pvalue)
    Baylor_modreacdf.3<-rbind(Baylor_modreacdf.3,select.top)
  }
  
}
#Now Multiple Testing Adjustment by module
Baylor_modreacdf.4<-data.frame()
Baylor_mod.in.df<-unique(Baylor_modreacdf.3$module)
for (modu in 1:length(Baylor_mod.in.df)){
  #get module pathways
  to.adjust<-Baylor_modreacdf.3[which(Baylor_modreacdf.3$module==Baylor_mod.in.df[modu]),]
  #get the pavals and adjust
  oldP<-to.adjust$pvalue
  newQ<-p.adjust(oldP,method="BH")
  to.adjust<-cbind(to.adjust,newQ)
  #select based on adjusted pvals
  adjusted<-to.adjust[which(to.adjust$newQ<0.05),]
  if(nrow(adjusted)==0 && nrow(to.adjust)>3){adjusted<-to.adjust[1:3,]}
  if(nrow(adjusted)>10){adjusted<-adjusted[1:10,]}
  #write to df
  Baylor_modreacdf.4<-rbind(Baylor_modreacdf.4,adjusted)
}
#9. NETWORK 9: module-cell_expr #########
#need to combine these somehow...
try(Baylor_bloodResults1<-userListEnrichment(geneR=modgenes2$Symbol,labelR=modgenes2$Module,nameOut=file.path(dir.results,"Baylor_BloodEnrichment1.csv"),omitCategories="grey",useBloodAtlases=TRUE,outputCorrectedPvalues=FALSE))
try(Baylor_bloodResults2<-userListEnrichment(geneR=modgenes2$Symbol,labelR=modgenes2$Module,nameOut=file.path(dir.results,"Baylor_BloodEnrichment2.csv"),omitCategories="grey",useBloodAtlases=TRUE,outputCorrectedPvalues=TRUE))

##

if(neoinsert=="on"){
  #NewHaemAtlas only
  convgenes<-IlluminaID2nuID(modgenes2$X)
  try(Baylor_bloodResults3<-userListEnrichment(geneR=convgenes[,"nuID"],labelR=modgenes2$Module,fnIn=file.path(dir.results,"haemAtlas.csv"),catNmIn="HaemAtlas",nameOut=file.path(dir.results,"Baylor_BloodEnrichment_HA1.csv"),omitCategories="grey",outputCorrectedPvalues=FALSE)) 
  try(Baylor_bloodResults4<-userListEnrichment(geneR=convgenes[,"nuID"],labelR=modgenes2$Module,fnIn=file.path(dir.results,"haemAtlas.csv"),catNmIn="HaemAtlas",nameOut=file.path(dir.results,"Baylor_BloodEnrichment_HA2.csv"),omitCategories="grey",outputCorrectedPvalues=TRUE))
  # gi<-read.csv("/Users/armindeffur/Dropbox/Latest_Run_Combo/version_2016-02-08_12_59_25/Q_2_berry.train/5_WGCNA_4k_tv/results/table2.5.1.Square2Berry.5_WGCNA_4k_tvBlood_TBgeneInfoClass2.csv")
  # res2<-userListEnrichment(geneR=as.character(gi$substanceBXH),labelR=as.character(gi$moduleColor),fnIn="ha4.csv",nameOut="newHA.csv",omitCategories = "grey",outputCorrectedPvalues = FALSE)
}else{#OldHaemAtlas only
  try(Baylor_bloodResults3<-userListEnrichment(geneR=modgenes2$Symbol,labelR=modgenes2$Module,fnIn=file.path(dir.results,"haemAtlas_old.csv"),catNmIn="HaemAtlas",nameOut=file.path(dir.results,"Baylor_BloodEnrichment_HA1.csv"),omitCategories="grey",outputCorrectedPvalues=FALSE)) 
  try(Baylor_bloodResults4<-userListEnrichment(geneR=modgenes2$Symbol,labelR=modgenes2$Module,fnIn=file.path(dir.results,"haemAtlas_old.csv"),catNmIn="HaemAtlas",nameOut=file.path(dir.results,"Baylor_BloodEnrichment_HA2.csv"),omitCategories="grey",outputCorrectedPvalues=TRUE))
}

##


try(Baylor_bloodResults5<-userListEnrichment(geneR=modgenes2$Symbol,labelR=modgenes2$Module,fnIn=file.path(dir.results,"abbasvals.csv"),catNmIn="Abbas",nameOut=file.path(dir.results,"Baylor_BloodEnrichment_abbas1.csv"),omitCategories="grey",outputCorrectedPvalues=FALSE)) 
try(Baylor_bloodResults6<-userListEnrichment(geneR=modgenes2$Symbol,labelR=modgenes2$Module,fnIn=file.path(dir.results,"abbasvals.csv"),catNmIn="Abbas",nameOut=file.path(dir.results,"Baylor_BloodEnrichment_abbas2.csv"),omitCategories="grey",outputCorrectedPvalues=TRUE)) 
#10. NETWORK 10: Baylor module-PalazzoloWang ##################

if (length(modEN)<300){#was 25
  
  try(Baylor_PWResults1<-userListEnrichment(geneR=modgenes2$Symbol,labelR=modgenes2$Module,nameOut=file.path(dir.results,"BaylorPalazzoWangEnrichment1.csv"),omitCategories="grey",usePalazzoloWang=TRUE))
  try(Baylor_PWResults2<-userListEnrichment(geneR=modgenes2$Symbol,labelR=modgenes2$Module,nameOut=file.path(dir.results,"Baylor_PalazzoWangEnrichment2.csv"),omitCategories="grey",usePalazzoloWang=TRUE,outputCorrectedPvalues=FALSE))
  overlaps<-Baylor_PWResults2[["sigOverlaps"]]
  nrow(overlaps)
  
  modules.PW<-unique(overlaps$InputCategories)
  
  Baylor_PalWang.select<-data.frame()#this is quite elaborate now
  for (pwmod in 1:length(modules.PW)){
    modOverlap<-overlaps[which(overlaps$InputCategories==modules.PW[pwmod]),]
    #print(nrow(modOverlap))
    qvalues<-p.adjust(modOverlap$Pvalues,method="BH")
    modOverlap<-cbind(modOverlap,qvalues)
    modOverlap.select<-modOverlap[which(modOverlap$qvalues<0.01),]#select significant hits after BH correction
    if(nrow(modOverlap.select)>10){modOverlap.select<-modOverlap.select[1:10,]}#only take top 10
    if(nrow(modOverlap.select)<4){modOverlap.select<-modOverlap[1:min(4,nrow(modOverlap)),]}#only take top hit regardless of significance
    Baylor_PalWang.select<-rbind(Baylor_PalWang.select,modOverlap.select)
  }
  
  #add BH column
}#end NETWORK 9 if statement
#11. NETWORK 11: Baylor module-ImmunePW ######

try(Baylor_immuneResults1<-userListEnrichment(geneR=modgenes2$Symbol,labelR=modgenes2$Module,nameOut=file.path(dir.results,"Baylor_ImmunePathwayEnrichment1.csv"),omitCategories="grey",useImmunePathwayLists=TRUE,outputCorrectedPvalues=TRUE))
try(Baylor_immuneResults2<-userListEnrichment(geneR=modgenes2$Symbol,labelR=modgenes2$Module,nameOut=file.path(dir.results,"Baylor_ImmunePathwayEnrichment2.csv"),omitCategories="grey",useImmunePathwayLists=TRUE,outputCorrectedPvalues=FALSE))
#add BH column
#12. NETWORK 12: WGCNA-module enrichment for Baylor modules####
#make a list and save as csv
BaylorEnrichmentList<-modgenes2[,c("Symbol","Module")]
names(BaylorEnrichmentList)<-c("GENES","LISTS")
write.csv(BaylorEnrichmentList,file.path(dir.results,"BaylorEnrichmentList.csv"),row.names=FALSE)

#do enrichment
try(WGCNA_BaylorResults1<-userListEnrichment(geneR=geneInfo$geneSymbol,labelR=geneInfo$moduleColor,fnIn=file.path(dir.results,"BaylorEnrichmentList.csv"),catNmIn="BaylorMod",nameOut=file.path(dir.results,"WGCNA_BaylorEnrichment1.csv"),omitCategories="grey",outputCorrectedPvalues=TRUE))
try(WGCNA_BaylorResults2<-userListEnrichment(geneR=geneInfo$geneSymbol,labelR=geneInfo$moduleColor,fnIn=file.path(dir.results,"BaylorEnrichmentList.csv"),catNmIn="BaylorMod",nameOut=file.path(dir.results,"WGCNA_BaylorEnrichment2.csv"),omitCategories="grey",outputCorrectedPvalues=FALSE))

setwd(dir.figures)
#end baylor



#8_Networks**************************************####

analysis="8_Networks"
an.count=8
figure<-1
#make relevant output folders
## Define folder for storing results NB do this for each analysis
dir.results <- file.path(dir.output.version,paste("Q_",q.i,"_",dataset.variables,sep=""),analysis, "results")
## Define folder for saving figures
dir.figures <- file.path(dir.output.version,paste("Q_",q.i,"_",dataset.variables,sep=""),analysis, "figures")

for (dir in c(dir.figures, dir.results)) {
  if (!file.exists(dir)) {
    dir.create(dir, recursive=T,showWarnings=T)
  }
}
setwd(dir.figures)

doNetworks<-TRUE
if(length(moduleList)>60){doNetworks<-FALSE}#skip this section for now if there are too many modules for sensible analysis. I need to find a way to bring down the number of modules for low N groups

#Start doNetworks if()####
if(doNetworks){
  if(cytoscapelink=="on"){
    # If you also want to choose a layout from R, a list  of available layouts can be accessed as follows:
    cy <- CytoscapeConnection(host="192.168.65.2",rpcPort=9000)
    hlp <-getLayoutNames(cy)
  }#end 
  #Start edgewise Cytoscape loop ####
  for (edge in 1:length(edges)){
    ## Define folder for storing results NB do this for each analysis
    dir.results <- file.path(dir.output.version,paste("Q_",q.i,"_",dataset.variables,sep=""),analysis,paste("edge",edge,sep=""), "results")
    ## Define folder for saving figures
    dir.figures <- file.path(dir.output.version,paste("Q_",q.i,"_",dataset.variables,sep=""),analysis,paste("edge",edge,sep=""), "figures")
    
    for (dir in c(dir.figures, dir.results)) {
      if (!file.exists(dir)) {
        dir.create(dir, recursive=T,showWarnings=T)
      }
    }
    setwd(dir.figures)
    
    #   #recalculate the correlation matrices for the relevant edge
    #   nSamples=length(edges[[edge]])
    #   #Correlate MEs and traits
    #   moduleTraitCor = cor(ME.edges[[edge]], datTraits[edges[[edge]],], use = "p");
    #   moduleTraitCor = moduleTraitCor[,!apply(moduleTraitCor,2,is.na)[1,]]
    #   moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);
    #   moduleTraitQvalue = matrix(p.adjust(moduleTraitPvalue,method="BH"),nrow=nrow(moduleTraitPvalue),dimnames=list(row.names(moduleTraitPvalue),colnames(moduleTraitPvalue)))
    #   #edges
    #   
    #   
    #   moduleCellCor = cor(ME.edges[[edge]], cellprops[edges[[edge]],], use = "p");
    #   moduleCellCor = moduleCellCor[,!apply(moduleCellCor,2,is.na)[1,]]
    #   moduleCellPvalue = corPvalueStudent(moduleCellCor, nSamples);
    #   moduleCellQvalue = matrix(p.adjust(moduleCellPvalue,method="BH"),nrow=nrow(moduleCellPvalue),dimnames=list(row.names(moduleCellPvalue),colnames(moduleCellPvalue)))
    
    print("helper 1")
    
    #NETWORK 1 helper####
    #recalculate edwise correlation data:
    nSamples=length(edges[[edge]])
    #Correlate MEs and traits
    moduleTraitCor = cor(ME.edges[[edge]], datTraits[edges[[edge]],], use = "p");
    moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);
    moduleTraitQvalue = matrix(p.adjust(moduleTraitPvalue,method="BH"),nrow=nrow(moduleTraitPvalue),dimnames=list(row.names(moduleTraitPvalue),colnames(moduleTraitPvalue)))
    
    q<-moduleTraitQvalue[,!apply(moduleTraitQvalue,2,is.na)[1,]]
    qm<-as.matrix(q,nrow=nrow(q),dimnames=dimnames(q))
    if(ncol(qm)==1){colnames(qm)<-colnames(moduleTraitQvalue)[which(!apply(moduleTraitQvalue,2,is.na)[1,])]}
    t<-moduleTraitCor[,!apply(moduleTraitCor,2,is.na)[1,]]
    q2<-as.matrix(q,nrow=nrow(q),dimnames=dimnames(q))
    if(ncol(q2)==1){colnames(q2)<-colnames(moduleTraitQvalue)[which(!apply(moduleTraitQvalue,2,is.na)[1,])]}
    t2<-as.matrix(t,nrow=nrow(t),dimnames=dimnames(t))
    if(ncol(t2)==1){colnames(t2)<-colnames(moduleTraitQvalue)[which(!apply(moduleTraitQvalue,2,is.na)[1,])]}
    q2[q2>0.05]<-1
    q2[q2<0.05]<-2
    q3<-q2-1
    t2[q3==0]<-0
    
    edgelist<-list()
    count<-0
    for (ix in 1:nrow(t2)){
      for (j in 1:ncol(t2)){
        if (t2[ix,j]!=0){
          count=count+1
          res<-c(row.names(t2)[ix],colnames(t2)[j],t2[ix,j],qm[ix,j])
          edgelist[[count]]<-res
        }
        
      }
    }
    
    if(length(edgelist)==0){
      q<-moduleTraitQvalue[,!apply(moduleTraitQvalue,2,is.na)[1,]]
      qm<-as.matrix(q,nrow=nrow(q),dimnames=dimnames(q))
      if(ncol(qm)==1){colnames(qm)<-colnames(moduleTraitQvalue)[which(!apply(moduleTraitQvalue,2,is.na)[1,])]}
      t<-moduleTraitCor[,!apply(moduleTraitCor,2,is.na)[1,]]
      q2<-as.matrix(q,nrow=nrow(q),dimnames=dimnames(q))
      if(ncol(q2)==1){colnames(q2)<-colnames(moduleTraitQvalue)[which(!apply(moduleTraitQvalue,2,is.na)[1,])]}
      t2<-as.matrix(t,nrow=nrow(t),dimnames=dimnames(t))
      if(ncol(t2)==1){colnames(t2)<-colnames(moduleTraitQvalue)[which(!apply(moduleTraitQvalue,2,is.na)[1,])]}
      q2[q2>0.15]<-1
      q2[q2<0.15]<-2
      q3<-q2-1
      t2[q3==0]<-0
      
      edgelist<-list()
      count<-0
      for (ix in 1:nrow(t2)){
        for (j in 1:ncol(t2)){
          if (t2[ix,j]!=0){
            count=count+1
            res<-c(row.names(t2)[ix],colnames(t2)[j],t2[ix,j],qm[ix,j])
            edgelist[[count]]<-res
          }
          
        }
      }
      
    }
    
    if(length(edgelist)==0){
      q<-moduleTraitQvalue[,!apply(moduleTraitQvalue,2,is.na)[1,]]
      qm<-as.matrix(q,nrow=nrow(q),dimnames=dimnames(q))
      if(ncol(qm)==1){colnames(qm)<-colnames(moduleTraitQvalue)[which(!apply(moduleTraitQvalue,2,is.na)[1,])]}
      t<-moduleTraitCor[,!apply(moduleTraitCor,2,is.na)[1,]]
      q2<-as.matrix(q,nrow=nrow(q),dimnames=dimnames(q))
      if(ncol(q2)==1){colnames(q2)<-colnames(moduleTraitQvalue)[which(!apply(moduleTraitQvalue,2,is.na)[1,])]}
      t2<-as.matrix(t,nrow=nrow(t),dimnames=dimnames(t))
      if(ncol(t2)==1){colnames(t2)<-colnames(moduleTraitQvalue)[which(!apply(moduleTraitQvalue,2,is.na)[1,])]}
      q2[q2>0.25]<-1
      q2[q2<0.25]<-2
      q3<-q2-1
      t2[q3==0]<-0
      
      edgelist<-list()
      count<-0
      for (ix in 1:nrow(t2)){
        for (j in 1:ncol(t2)){
          if (t2[ix,j]!=0){
            count=count+1
            res<-c(row.names(t2)[ix],colnames(t2)[j],t2[ix,j],qm[ix,j])
            edgelist[[count]]<-res
          }
          
        }
      }
      
    }
    
    if(length(edgelist)==0){
      q<-moduleTraitQvalue[,!apply(moduleTraitQvalue,2,is.na)[1,]]
      qm<-as.matrix(q,nrow=nrow(q),dimnames=dimnames(q))
      if(ncol(qm)==1){colnames(qm)<-colnames(moduleTraitQvalue)[which(!apply(moduleTraitQvalue,2,is.na)[1,])]}
      t<-moduleTraitCor[,!apply(moduleTraitCor,2,is.na)[1,]]
      q2<-as.matrix(q,nrow=nrow(q),dimnames=dimnames(q))
      if(ncol(q2)==1){colnames(q2)<-colnames(moduleTraitQvalue)[which(!apply(moduleTraitQvalue,2,is.na)[1,])]}
      t2<-as.matrix(t,nrow=nrow(t),dimnames=dimnames(t))
      if(ncol(t2)==1){colnames(t2)<-colnames(moduleTraitQvalue)[which(!apply(moduleTraitQvalue,2,is.na)[1,])]}
      q2[q2>1.1]<-1
      q2[q2<1.1]<-2
      q3<-q2-1
      t2[q3==0]<-0
      
      edgelist<-list()
      count<-0
      for (ix in 1:nrow(t2)){
        for (j in 1:ncol(t2)){
          if (t2[ix,j]!=0){
            count=count+1
            res<-c(row.names(t2)[ix],colnames(t2)[j],t2[ix,j],qm[ix,j])
            edgelist[[count]]<-res
          }
          
        }
      }
      
    }
    
    edgemm<-t(matrix(unlist(edgelist),nrow=4))
    colnames(edgemm)<-c("module","pheno","weight","qvalue")
    newmod<-apply(edgemm,1,function(x){x[1]<-substr(x[1],3,nchar(x[1]))})
    edgemm[,1]<-newmod
    
    write.csv(edgemm,file=file.path(dir.results,paste('table',q.i,'.',an.count,'.',i,'.',names(questions)[[q.i]],'.',analysis,names(datalist)[i],"module_pheno.csv",sep="")))
    
    phenos<-unique(edgemm[,2])
    phenos2<-cbind(phenos,rep("pheno",length(phenos)))
    colnames(phenos2)<-c("pheno","kind")
    write.csv(phenos2,file=file.path(dir.results,paste('table',q.i,'.',an.count,'.',i,'.',names(questions)[[q.i]],'.',analysis,names(datalist)[i],"node_pheno.csv",sep="")))
    
    #node annotation
    modinfo<-geneInfo
    modules<-as.character(unique(modinfo$moduleColor))
    node.modules<-data.frame(cbind(modules),rep("module",length(modules)))
    colnames(node.modules)<-c("module","kind")
    write.csv(node.modules,file=file.path(dir.results,paste('table',q.i,'.',an.count,'.',i,'.',names(questions)[[q.i]],'.',analysis,names(datalist)[i],"node_module.csv",sep="")))
    #END NETWORK 1 helper
    
    print("igraph 1")
    #NETWORK 1 igraph####
    
    # Data format: dataframe with 3 variables; variables 1 & 2 correspond to interactions; variable 3 corresponds to the weight of interaction
    dataSet0 <- data.frame(edgemm)
    dataSet0<-merge(dataSet0,diffMEmatrix[,edge], by.x="module",by.y="row.names",all.x=T,all.y=F)
    dataSet0<-transform(dataSet0,weight=as.numeric(as.character(weight)))
    colnames(dataSet0)<-c("module","pheno","weight","qvalue","diffME")
    dataSet<-cbind(dataSet0,abs(dataSet0$weight)^2)
    colnames(dataSet)<-c("module","pheno","weight","qvalue","diffME","R_sq")
    # Create a graph. Use simplify to ensure that there are no duplicated edges or self loops
    gD <- graph.data.frame(dataSet[,c(1,2,3,4,6)], directed=FALSE)
    #not all modules have associations with phenotype variables, hence:
    modulesInSet<-levels(as.factor(dataSet$module))
    phenosInSet<-levels(as.factor(dataSet$pheno))
    
    # Print number of nodes and edges
    vcount(gD)
    ecount(gD)
    
    V(gD)[modulesInSet]$kind="module"
    V(gD)[phenosInSet]$kind="pheno"
    
    # Check the attributes
    summary(gD)
    
    gD <- set.edge.attribute(gD, "weight", index = E(gD), value = 0.0)
    gD <- set.edge.attribute(gD, "R_sq", index = E(gD), value = 0.0)
    gD <- set.edge.attribute(gD, "qvalue", index = E(gD), value = 0.0)
    
    #gD <- set.edge.attribute(gD, "similarity", index = E(gD), value = 0)
    
    # The order of interactions in gD is not the same as it is in dataSet or as it is in the edge list,
    # and for that reason these values cannot be assigned directly
    
    E(gD)[as.character(dataSet$pheno) %--% as.character(dataSet$module)]$R_sq <- as.numeric(as.character(dataSet$R_sq))
    E(gD)[as.character(dataSet$pheno) %--% as.character(dataSet$module)]$weight <- as.numeric(as.character(dataSet$weight))
    E(gD)[as.character(dataSet$pheno) %--% as.character(dataSet$module)]$qvalue <- as.numeric(as.character(dataSet$qvalue))
    
    # Check the attributes
    summary(gD)
    
    #END IGRAPH
    #NETWORK 1 neo4j####
    if (neoinsert=="on"){
      #redefine nodes here so that RNeo4j works
      nodes<-RNeo4j::nodes
      
      #get the dataframe to use for relationships and nodes
      graphdata<-dataSet
      #network-specific cypher query #NB this does not return all modules for the edge, only those that have a significant pheno association
      # query = "
      # MERGE (wmod:wgcna {name:{wgcnamod},square:{squareg},edge:{edgeg},contrastvar:{contrastvarsg},contrast:{contrastg}})
      # MERGE (phen:pheno {name:{pheno},square:{squareg},edge:{edgeg},contrastvar:{contrastvarsg},contrast:{contrastg}})
      # CREATE (wmod)-[:correlates {weight:{weightg},qvalue:{qval},Rsq:{rsqg},square:{squareg},edge:{edgeg},contrastvar:{contrastvarsg},contrast:{contrastg}}]->(phen)
      # "
      t = suppressMessages(newTransaction(graph))
      
      for (therow in 1:nrow(graphdata)) {
        wmod = as.character(graphdata[therow, ]$module)
        phen = as.character(graphdata[therow, ]$pheno)
        diffMEg = as.numeric(as.character(graphdata[therow,]$diffME))
        query = paste("MERGE (wmod:wgcna {name:{wgcnamod},square:{squareg},edge:{edgeg},contrastvar:{contrastvarsg},contrast:{contrastg}}) ON CREATE SET wmod.diffME = ",diffMEg," ON MATCH SET wmod.diffME = ",diffMEg," MERGE (phen:pheno {name:{pheno},square:{squareg},edge:{edgeg},contrastvar:{contrastvarsg},contrast:{contrastg}})  CREATE (wmod)-[:correlates {weight:{weightg},qvalue:{qval},Rsq:{rsqg},square:{squareg},edge:{edgeg},contrastvar:{contrastvarsg},contrast:{contrastg}}]->(phen)",sep="")
        
        suppressMessages(appendCypher(t, 
                                      query, 
                                      wgcnamod = wmod, 
                                      pheno = phen,
                                      squareg = dataset.variables[[1]],
                                      edgeg = edge,
                                      contrastvarsg = contrastvars[[edge]],
                                      contrastg = paste(pathways[[edge]][[2]],"-",pathways[[edge]][[1]]),
                                      qval = as.numeric(as.character(graphdata[therow,]$qvalue)),
                                      weightg = as.numeric(as.character(graphdata[therow,]$weight)),
                                      rsqg = as.numeric(as.character(graphdata[therow,]$R_sq))#,
                                      #diffMEg = as.numeric(as.character(graphdata[therow,]$diffME))
        ))
      }
      
      suppressMessages(commit(t))
      
      #END NEO4J
    }#end neoinsert if
    #NETWORK 1 graphNEL####
    
    gD.cyt <- igraph.to.graphNEL(gD)
    
    # We have to create attributes for graphNEL
    # We'll keep the same name, so the values are passed from igraph node attributes
    
    # gD.cyt <- initNodeAttribute(gD.cyt, 'label', 'char', "x")
    gD.cyt <- initNodeAttribute(gD.cyt, 'kind', 'char', "y")
    # gD.cyt <- initNodeAttribute(gD.cyt, 'degree', 'numeric', 0) 
    # gD.cyt <- initNodeAttribute(gD.cyt, 'betweenness', 'numeric', 0) 
    gD.cyt <- initNodeAttribute(gD.cyt, 'diffME', 'numeric', 0.0) 
    # edge attributes
    gD.cyt <- initEdgeAttribute (gD.cyt, "qvalue", 'numeric', 0)
    gD.cyt <- initEdgeAttribute (gD.cyt, "weight", 'numeric', 0)
    #gD.cyt <- initEdgeAttribute (gD.cyt, "similarity", 'numeric', 0)
    gD.cyt <- initEdgeAttribute (gD.cyt, "R_sq", 'numeric', 0.0)
    gD.cyt <- initEdgeAttribute (gD.cyt, "edgeType", "char", "undefined")
    #NETWORK 1 RCytoscape####
    if(cytoscapelink=="on"){
      #open Cytoscape connection
      
      # Now we can create a new graph window in cytoscape
      # Be sure that CytoscapeRPC plugin is activated
      gDCW <- new.CytoscapeWindow("NETWORK 1", graph = gD.cyt, overwriteWindow = TRUE,host="192.168.65.2",rpcPort=9000)
      getAllEdgeAttributes(gDCW)
      # We can display graph, with defaults color/size scheme
      displayGraph(gDCW)
      
      # # If you also want to choose a layout from R, a list  of available layouts can be accessed as follows:
      # cy <- CytoscapeConnection(host="192.168.65.2",rpcPort=9000)
      # hlp <-getLayoutNames(cy)
      # See properties for the given layout
      # getLayoutPropertyNames(cy, hlp[6])
      # Apply values to some of the properties and plot the layout
      setLayoutProperties (gDCW, hlp[6], list (edge_attribute = 'R_sq'))
      #another change: use fruchterman-rheingold [19]
      layoutNetwork(gDCW, hlp[6])
      
      # Now, we can define our own default color/size scheme
      setDefaultBackgroundColor(gDCW, '#FFFFFF')
      setDefaultEdgeColor(gDCW, '#CDC9C9')
      setDefaultEdgeLineWidth(gDCW, 1)
      setDefaultNodeBorderColor(gDCW, '#000000')
      setDefaultNodeBorderWidth(gDCW, 3)
      setDefaultNodeShape(gDCW, 'ellipse')
      setDefaultNodeColor(gDCW, '#87CEFA')
      #setDefaultNodeSize(gDCW, 60)
      #setDefaultNodeFontSize(gDCW, 20)
      #setDefaultNodeLabelColor(gDCW, '#000000')
      
      # And we can replot it 
      redraw(gDCW)       
      
      # Rules for node colors, node sizes, and edge colors
      
      #Rules
      ##Send diffME to Cytoscape
      new.nodes<-intModules
      new.diffME<-unlist(all.difflist.modules[[edge]])
      setNodeAttributesDirect(gDCW,"diffME","floating",new.nodes,new.diffME)
      
      #Colour modules by up or down
      control.points <- c (-1.0, 0.0, 1.0)
      node.colors <- rgb(t(col2rgb(c("black","blue","white","red","green"))),maxColorValue = 255)
      setNodeColorRule(gDCW, node.attribute.name='diffME', control.points,node.colors, mode='interpolate')
      redraw (gDCW)
      
      #pheno nodes are red
      phenocol<-rgb(t(col2rgb(c("red"))),maxColorValue = 255)
      setNodeColorDirect (gDCW,phenos,phenocol)
      
      #invariant node shapes
      data.values <- c ("module","pheno")
      node.shapes<-c("round_rect","ellipse")
      setNodeShapeRule(gDCW,node.attribute.name="kind",data.values,node.shapes,default.shape="ellipse")
      
      #invariant node dimensions
      lockNodeDimensions(gDCW,FALSE)
      node.width <- c(40,20) 
      #setNodeSizeRule(gDCW, 'kind',data.values,node.sizes,mode="lookup")
      setNodeWidthDirect(gDCW, modules,60)
      setNodeWidthDirect(gDCW, phenos,50)
      newlabelcol<-rgb(t(col2rgb(c("lightblue"))),maxColorValue = 255)
      setNodeLabelColorDirect(gDCW, phenos, newlabelcol)
      
      control.points <- c (-1.0, 0.0, 1.0)
      edge.colors <- rgb(t(col2rgb(c("black","blue","white","red","green"))),maxColorValue = 255)
      setEdgeColorRule(gDCW, edge.attribute.name='weight', control.points,edge.colors, mode='interpolate')
      redraw (gDCW)
    }#end optional cytoscape
    
    #NETWORK 2 helper####
    print("helper 2")
    nSamples=length(edges[[edge]])
    #Correlate MEs and cellprops
    moduleCellCor = cor(ME.edges[[edge]], cellprops[edges[[edge]],], use = "p");
    moduleCellPvalue = corPvalueStudent(moduleCellCor, nSamples);
    moduleCellQvalue = matrix(p.adjust(moduleCellPvalue,method="BH"),nrow=nrow(moduleCellPvalue),dimnames=list(row.names(moduleCellPvalue),colnames(moduleCellPvalue)))
    
    r<-moduleCellQvalue[,!apply(moduleCellQvalue,2,is.na)[1,]] 
    s<-moduleCellCor[,!apply(moduleCellCor,2,is.na)[1,]]
    r2<-as.matrix(r,nrow=nrow(r),dimnames=dimnames(r))
    s2<-as.matrix(s,nrow=nrow(s),dimnames=dimnames(s))
    r2[r2>0.05]<-1
    r2[r2<0.05]<-2
    r3<-r2-1
    s2[r3==0]<-0
    
    edgelist2<-list()
    count<-0
    for (ix in 1:nrow(s2)){
      for (j in 1:ncol(s2)){
        if (s2[ix,j]!=0){
          count=count+1
          res<-c(row.names(s2)[ix],colnames(s2)[j],s2[ix,j],r[ix,j])
          edgelist2[[count]]<-res
        }
        
      }
    }
    
    if(length(edgelist2)==0){
      r<-moduleCellQvalue[,!apply(moduleCellQvalue,2,is.na)[1,]] 
      s<-moduleCellCor[,!apply(moduleCellCor,2,is.na)[1,]]
      r2<-as.matrix(r,nrow=nrow(r),dimnames=dimnames(r))
      s2<-as.matrix(s,nrow=nrow(s),dimnames=dimnames(s))
      r2[r2>0.25]<-1
      r2[r2<0.25]<-2
      r3<-r2-1
      s2[r3==0]<-0
      edgelist2<-list()
      count<-0
      for (ix in 1:nrow(s2)){
        for (j in 1:ncol(s2)){
          if (s2[ix,j]!=0){
            count=count+1
            res<-c(row.names(s2)[ix],colnames(s2)[j],s2[ix,j],r[ix,j])
            edgelist2[[count]]<-res
          }
          
        }
      }  
      
    }
    
    ###TMP fix for bad datasets: no multiple testing
    
    if(length(edgelist2)==0){
      print(paste("no multiple testing for edge",edge))
      r<-moduleCellPvalue[,!apply(moduleCellPvalue,2,is.na)[1,]] 
      s<-moduleCellCor[,!apply(moduleCellCor,2,is.na)[1,]]
      r2<-as.matrix(r,nrow=nrow(r),dimnames=dimnames(r))
      s2<-as.matrix(s,nrow=nrow(s),dimnames=dimnames(s))
      r2[r2>0.25]<-1
      r2[r2<0.25]<-2
      r3<-r2-1
      s2[r3==0]<-0
      edgelist2<-list()
      count<-0
      for (ix in 1:nrow(s2)){
        for (j in 1:ncol(s2)){
          if (s2[ix,j]!=0){
            count=count+1
            res<-c(row.names(s2)[ix],colnames(s2)[j],s2[ix,j],r[ix,j])
            edgelist2[[count]]<-res
          }
          
        }
      }  
      
    }
    ##end fix
    
    edgemm2<-t(matrix(unlist(edgelist2),nrow=4))
    colnames(edgemm2)<-c("module","cell","weight","qvalue")
    newmod2<-apply(edgemm2,1,function(x){x[1]<-substr(x[1],3,nchar(x[1]))})
    edgemm2[,1]<-newmod2
    write.csv(edgemm2,file=file.path(dir.results,paste('table',q.i,'.',an.count,'.',i,'.',names(questions)[[q.i]],'.',analysis,names(datalist)[i],"module_cellprop.csv",sep="")))
    
    cells<-unique(edgemm2[,2])
    cells2<-cbind(cells,rep("cell_prop",length(cells)))
    colnames(cells2)<-c("cell_prop","kind")
    write.csv(cells2,file=file.path(dir.results,paste('table',q.i,'.',an.count,'.',i,'.',names(questions)[[q.i]],'.',analysis,names(datalist)[i],"node_cellprop.csv",sep="")))
    
    #END NETWORK 2 helper
    
    print("igraph 2")
    #NETWORK 2 igraph####
    
    #BEGIN igraph 2
    
    # Read a data set. 
    # Data format: dataframe with 3 variables; variables 1 & 2 correspond to interactions; variable 3 corresponds to the weight of interaction
    dataSet0 <- data.frame(edgemm2)
    dataSet0<-merge(dataSet0,diffMEmatrix[,edge], by.x="module",by.y="row.names",all.x=T,all.y=F)
    dataSet0<-transform(dataSet0,weight=as.numeric(as.character(weight)))
    colnames(dataSet0)<-c("module","cell","weight","qvalue","diffME")
    dataSet<-cbind(dataSet0,abs(dataSet0$weight)^2)
    colnames(dataSet)<-c("module","cell","weight","qvalue","diffME","R_sq")
    modulesInSet<-levels(as.factor(dataSet$module))
    cellsInSet<-levels(as.factor(dataSet$cell))
    # Create a graph. Use simplify to ensure that there are no duplicated edges or self loops
    gD2 <- graph.data.frame(dataSet[,c(1,2,3,4,6)], directed=FALSE)
    
    # Print number of nodes and edges
    vcount(gD2)
    ecount(gD2)
    # 
    # # Calculate some node properties and node similarities that will be used to illustrate 
    # # different plotting abilities
    # 
    # # Calculate degree for all nodes
    # degAll <- degree(gD2, v = V(gD2), mode = "all")
    # 
    # # Calculate betweenness for all nodes
    # betAll <- betweenness(gD2, v = V(gD2), directed = FALSE,weights=dataSet$R_sq) / (((vcount(gD2) - 1) * (vcount(gD2)-2)) / 2)
    # betAll.norm <- (betAll - min(betAll))/(max(betAll) - min(betAll))
    # # rm(betAll)
    # 
    # #Calculate Dice similarities between all pairs of nodes
    # #dsAll <- similarity.dice(gD, vids = V(gD), mode = "all")
    
    # Add new node/edge attributes based on the calculated node properties/similarities
    
    # gD2 <- set.vertex.attribute(gD2, "degree", index = V(gD2), value = degAll)
    # gD2 <- set.vertex.attribute(gD2, "betweenness", index = V(gD2), value = betAll.norm)
    V(gD2)[modulesInSet]$kind="module"
    V(gD2)[cellsInSet]$kind="cell_prop"
    
    list.vertex.attributes(gD2)
    
    # Check the attributes
    summary(gD2)
    
    #F1 <- function(x) {data.frame(similarity = dsAll[which(V(gD)$name == as.character(x$module)), which(V(gD)$name == as.character(x$pheno))])}
    #dataSet.ext <- ddply(.data=dataSet, .variables=c("module", "pheno", "weight"), .fun=function(x) data.frame(F1(x)))
    dataSet.ext<-dataSet
    
    gD2 <- set.edge.attribute(gD2, "weight", index = E(gD2), value = 0.0)
    gD2 <- set.edge.attribute(gD2, "R_sq", index = E(gD2), value = 0.0)
    gD2 <- set.edge.attribute(gD2, "qvalue", index = E(gD2), value = 0.0)
    #gD <- set.edge.attribute(gD, "similarity", index = E(gD), value = 0)
    
    # The order of interactions in gD is not the same as it is in dataSet or as it is in the edge list,
    # and for that reason these values cannot be assigned directly
    
    E(gD2)[as.character(dataSet.ext$module) %--% as.character(dataSet.ext$cell)]$R_sq <- dataSet.ext$R_sq
    E(gD2)[as.character(dataSet.ext$module) %--% as.character(dataSet.ext$cell)]$weight <- dataSet$weight
    E(gD2)[as.character(dataSet.ext$module) %--% as.character(dataSet.ext$cell)]$qvalue <- as.numeric(as.character(dataSet$qvalue))
    #E(gD)[as.character(dataSet.ext$module) %--% as.character(dataSet.ext$pheno)]$similarity <- as.numeric(dataSet.ext$similarity)
    
    # Check the attributes
    summary(gD2)
    
    #END IGRAPH
    #NETWORK 2 neo4j####
    if (neoinsert=="on"){
      #get the dataframe to use for relationships and nodes
      graphdata=dataSet
      t = suppressMessages(newTransaction(graph))
      
      #remove rows where ratio or p-value are NA or infinite; these are at best hard to interpret, but more likely meaningless results. The full results are still in the exported csv
      
      cellproptests<-cellproptests[!is.na(cellproptests$ratio),]
      cellproptests<-cellproptests[!is.na(cellproptests$pval),]
      cellproptests<-cellproptests[cellproptests$ratio!=Inf,]
      
      for (therow in 1:nrow(graphdata)) {
        wmod = as.character(graphdata[therow, ]$module)
        cp = as.character(graphdata[therow, ]$cell)
        
        #data for some cell types may have been removed above! Some of these reults will be empty.
        ratioX = as.numeric(as.character(cellproptests[which(cellproptests$names==cp&cellproptests$edge==edge),]$ratio))
        diffPX = as.numeric(as.character(cellproptests[which(cellproptests$names==cp&cellproptests$edge==edge),]$pval))
        diffQX = as.numeric(as.character(cellproptests[which(cellproptests$names==cp&cellproptests$edge==edge),]$BH))
        
        #code below only sets ratio, diffP and diffQ properties on cellprop nodes if the values exist. Correlation with cellprop is still included as result, even if ratio is absent or meaningless.
        if (length(ratioX)==0|length(diffPX)==0|length(diffQX)==0){
          diffMEg = as.numeric(as.character(graphdata[therow,]$diffME))
          query = paste("MERGE (wmod:wgcna {name:{wgcnamod},square:{squareg},edge:{edgeg},contrastvar:{contrastvarsg},contrast:{contrastg}}) ON CREATE SET wmod.diffME = ",diffMEg," ON MATCH SET wmod.diffME = ",diffMEg," MERGE (cp:cellprop {name:{cellprop},square:{squareg},edge:{edgeg},contrastvar:{contrastvarsg},contrast:{contrastg}}) CREATE (wmod)-[:correlates {weight:{weightg},qvalue:{qval},Rsq:{rsqg},square:{squareg},edge:{edgeg},contrastvar:{contrastvarsg},contrast:{contrastg}}]->(cp)",sep="")
          #print(paste("case1",ratioX,diffPX,diffQX))
          suppressMessages(appendCypher(t, 
                                        query, 
                                        #nodes
                                        wgcnamod = wmod, 
                                        cellprop = cp,
                                        #invariant node properties
                                        squareg = dataset.variables[[1]],
                                        edgeg = edge,
                                        contrastvarsg = contrastvars[[edge]],
                                        contrastg = paste(pathways[[edge]][[2]],"-",pathways[[edge]][[1]]),
                                        qval = as.numeric(as.character(graphdata[therow,]$qvalue)),
                                        weightg = as.numeric(as.character(graphdata[therow,]$weight)),
                                        rsqg = as.numeric(as.character(graphdata[therow,]$R_sq))#,
                                        #                ratiog = ratioX,
                                        #                diffPg = diffPX,
                                        #                diffQg = diffQX
          ))
        }else{
          diffMEg = as.numeric(as.character(graphdata[therow,]$diffME))
          query = paste("MERGE (wmod:wgcna {name:{wgcnamod},square:{squareg},edge:{edgeg},contrastvar:{contrastvarsg},contrast:{contrastg}}) ON CREATE SET wmod.diffME = ",diffMEg," ON MATCH SET wmod.diffME = ",diffMEg," MERGE (cp:cellprop {name:{cellprop},square:{squareg},edge:{edgeg},contrastvar:{contrastvarsg},contrast:{contrastg},ratio:{ratiog},diffP:{diffPg},diffQ:{diffQg}}) CREATE (wmod)-[:correlates {weight:{weightg},qvalue:{qval},Rsq:{rsqg},square:{squareg},edge:{edgeg},contrastvar:{contrastvarsg},contrast:{contrastg}}]->(cp)",sep="")
          #print(paste("case2",ratioX,diffPX,diffQX))
          suppressMessages(appendCypher(t, 
                                        query, 
                                        #nodes
                                        wgcnamod = wmod, 
                                        cellprop = cp,
                                        #invariant node properties
                                        squareg = dataset.variables[[1]],
                                        edgeg = edge,
                                        contrastvarsg = contrastvars[[edge]],
                                        contrastg = paste(pathways[[edge]][[2]],"-",pathways[[edge]][[1]]),
                                        qval = as.numeric(as.character(graphdata[therow,]$qvalue)),
                                        weightg = as.numeric(as.character(graphdata[therow,]$weight)),
                                        rsqg = as.numeric(as.character(graphdata[therow,]$R_sq)),
                                        ratiog = ratioX,
                                        diffPg = diffPX,
                                        diffQg = diffQX
          ))
        }
      }
      suppressMessages(commit(t))
      
      #END NEO4J
    }#end neoinsert if
    #NETWORK 2 graphNEL####
    
    gD2.cyt <- igraph.to.graphNEL(gD2)
    
    # We have to create attributes for graphNEL
    # We'll keep the same name, so the values are passed from igraph
    # node attributes
    # gD2.cyt <- initNodeAttribute(gD2.cyt, 'label', 'char', "x")
    gD2.cyt <- initNodeAttribute(gD2.cyt, 'kind', 'char', "y")
    # gD2.cyt <- initNodeAttribute(gD2.cyt, 'degree', 'numeric', 0) 
    # gD2.cyt <- initNodeAttribute(gD2.cyt, 'betweenness', 'numeric', 0) 
    gD2.cyt <- initNodeAttribute (gD2.cyt, "diffME", 'numeric', 0.0)
    # edge attributes
    gD2.cyt <- initEdgeAttribute (gD2.cyt, "weight", 'numeric', 0)
    gD2.cyt <- initEdgeAttribute (gD2.cyt, "qvalue", 'numeric', 0)
    #gD2.cyt <- initEdgeAttribute (gD2.cyt, "weight2", 'numeric', 0)
    #gD.cyt <- initEdgeAttribute (gD.cyt, "similarity", 'numeric', 0)
    gD2.cyt <- initEdgeAttribute (gD2.cyt, "R_sq", 'numeric', 0.0)
    gD2.cyt <- initEdgeAttribute (gD2.cyt, "edgeType", "char", "undefined")
    #NETWORK 2 RCytoscape####
    if(cytoscapelink=="on"){
      #open Cytoscape connection
      
      # Now we can create a new graph window in cytoscape
      # Be sure that CytoscapeRPC plugin is activated
      gDCW2 <- new.CytoscapeWindow("NETWORK 2", graph = gD2.cyt, overwriteWindow = TRUE,host="192.168.65.2",rpcPort=9000)
      
      ##Send diffME to Cytoscape
      displayGraph(gDCW2)
      new.nodes<-intModules[as.character(tests.modules$module)%in%as.character(unique(dataSet$module))]
      new.diffME<-unlist(all.difflist.modules[[edge]])[as.character(tests.modules$module)%in%as.character(unique(dataSet$module))]
      
      setNodeAttributesDirect(gDCW2,"diffME","numeric",new.nodes,new.diffME)
      # getNodeAttribute(gDCW2,"salmon","diffME")
      getAllNodeAttributes(gDCW2)
      getAllNodes(gDCW2)
      
      #Colour modules by up or down
      control.points <- c (-1.0, 0.0, 1.0)
      node.colors <- rgb(t(col2rgb(c("black","blue","white","red","green"))),maxColorValue = 255)
      setNodeColorRule(gDCW2, node.attribute.name='diffME', control.points,node.colors, mode='interpolate')
      redraw (gDCW2)
      
      # # If you also want to choose a layout from R, a list  of available layouts can be accessed as follows:
      # cy <- CytoscapeConnection(host="192.168.65.2",rpcPort=9000)
      # hlp <-getLayoutNames(cy)
      
      setLayoutProperties (gDCW2, hlp[6], list (edge_attribute = 'R_sq'))
      layoutNetwork(gDCW2, hlp[6])
      
      # Finally, we can define rules for node colors, node sizes, and edge colors
      #pheno nodes are red
      cellcol<-rgb(t(col2rgb(c("green"))),maxColorValue = 255)
      setNodeColorDirect(gDCW2,cells,cellcol)
      
      lockNodeDimensions(gDCW2,FALSE)
      setNodeWidthDirect(gDCW2, modules,60)
      setNodeWidthDirect(gDCW2, cells,50)
      newlabelcol<-rgb(t(col2rgb(c("black"))),maxColorValue = 255)
      setNodeLabelColorDirect(gDCW2, cells, newlabelcol)
      
      control.points <- c (-1.0, 0.0, 1.0)
      edge.colors <- rgb(t(col2rgb(c("black","blue","white","red","green"))),maxColorValue = 255)
      setEdgeColorRule(gDCW2, edge.attribute.name='weight', control.points,edge.colors, mode='interpolate')
      redraw (gDCW2)
    }#end optional cytoscape
    
    #NETWORK 2b conditional####
    if (dataset.variables=='berry.test'){
      #NETWORK 2b prep####
      #calc datapoints for the edge
      narow<-!apply(cellprops2m[edges[[edge]],],2,is.na)[,1]
      datapoints<-nrow(cellprops2m[edges[[edge]],][narow,])
      #Only add correlations where we are reasonably certain, i.e. at least 5 datapoints for the correlation
      if(datapoints>5){
        #NETWORK 2b helper####
        nSamples=length(edges[[edge]])
        
        #Correlate MEs and cellprops2m (this is now the same as in WGCNA section)
        narow<-!apply(cellprops2m[edges[[edge]],],2,is.na)[,1]
        datapoints<-nrow(cellprops2m[edges[[edge]],][narow,])
        moduleCellCor = cor(ME.edges[[edge]], cellprops2m[edges[[edge]],], use = "p");
        retain<-!apply(moduleCellCor,2,is.na)[1,]
        moduleCellCor = moduleCellCor[,retain] #NA columns removed!
        moduleCellPvalue = corPvalueStudent(moduleCellCor, nSamples);
        moduleCellQvalue = matrix(p.adjust(moduleCellPvalue,method="BH"),nrow=nrow(moduleCellPvalue),dimnames=list(row.names(moduleCellPvalue),colnames(moduleCellPvalue)))
        #edges
        write.csv(moduleCellCor,file=file.path(dir.results,paste('table',q.i,'.',an.count,'.',i,'.',names(edges)[edge],".",names(questions)[[q.i]],'.',analysis,names(datalist)[i],"ModuleFlowCellCor_edge",edge,".csv",sep="")))
        
        r<-moduleCellQvalue[,!apply(moduleCellQvalue,2,is.na)[1,]] 
        s<-moduleCellCor[,!apply(moduleCellCor,2,is.na)[1,]]
        r2<-as.matrix(r,nrow=nrow(r),dimnames=dimnames(r))
        s2<-as.matrix(s,nrow=nrow(s),dimnames=dimnames(s))
        r2[r2>0.05]<-1
        r2[r2<0.05]<-2
        r3<-r2-1
        s2[r3==0]<-0
        
        edgelist2<-list()
        count<-0
        for (ix in 1:nrow(s2)){
          for (j in 1:ncol(s2)){
            if (s2[ix,j]!=0){
              count=count+1
              res<-c(row.names(s2)[ix],colnames(s2)[j],s2[ix,j],r[ix,j])
              edgelist2[[count]]<-res
            }
            
          }
        }
        
        if(length(edgelist2)==0){
          r<-moduleCellQvalue[,!apply(moduleCellQvalue,2,is.na)[1,]] 
          s<-moduleCellCor[,!apply(moduleCellCor,2,is.na)[1,]]
          r2<-as.matrix(r,nrow=nrow(r),dimnames=dimnames(r))
          s2<-as.matrix(s,nrow=nrow(s),dimnames=dimnames(s))
          r2[r2>0.25]<-1
          r2[r2<0.25]<-2
          r3<-r2-1
          s2[r3==0]<-0
          edgelist2<-list()
          count<-0
          for (ix in 1:nrow(s2)){
            for (j in 1:ncol(s2)){
              if (s2[ix,j]!=0){
                count=count+1
                res<-c(row.names(s2)[ix],colnames(s2)[j],s2[ix,j],r[ix,j])
                edgelist2[[count]]<-res
              }
              
            }
          }  
          
        }
        
        edgemm2<-t(matrix(unlist(edgelist2),nrow=4))
        colnames(edgemm2)<-c("module","flowcell","weight","qvalue")
        newmod2<-apply(edgemm2,1,function(x){x[1]<-substr(x[1],3,nchar(x[1]))})
        edgemm2[,1]<-newmod2
        write.csv(edgemm2,file=file.path(dir.results,paste('table',q.i,'.',an.count,'.',i,'.',names(questions)[[q.i]],'.',analysis,names(datalist)[i],"module_flow_cellprop.csv",sep="")))
        
        cells<-unique(edgemm2[,2])
        cells2<-cbind(cells,rep("flow_cell_prop",length(cells)))
        colnames(cells2)<-c("flow_cell_prop","kind")
        write.csv(cells2,file=file.path(dir.results,paste('table',q.i,'.',an.count,'.',i,'.',names(questions)[[q.i]],'.',analysis,names(datalist)[i],"node_flow_cellprop.csv",sep="")))
        
        #END NETWORK 2 helper
        
        print("igraph 2")
        #NETWORK 2b igraph####
        
        #BEGIN igraph 2
        
        # Read a data set. 
        # Data format: dataframe with 3 variables; variables 1 & 2 correspond to interactions; variable 3 corresponds to the weight of interaction
        dataSet0 <- data.frame(edgemm2)
        dataSet0<-merge(dataSet0,diffMEmatrix[,edge], by.x="module",by.y="row.names",all.x=T,all.y=F)
        dataSet0<-transform(dataSet0,weight=as.numeric(as.character(weight)))
        colnames(dataSet0)<-c("module","flowcell","weight","qvalue","diffME")
        dataSet<-cbind(dataSet0,abs(dataSet0$weight)^2)
        colnames(dataSet)<-c("module","flowcell","weight","qvalue","diffME","R_sq")
        modulesInSet<-levels(as.factor(dataSet$module))
        cellsInSet<-levels(as.factor(dataSet$flowcell))
        # Create a graph. Use simplify to ensure that there are no duplicated edges or self loops
        gD2b <- graph.data.frame(dataSet[,c(1,2,3,4,6)], directed=FALSE)
        
        # Print number of nodes and edges
        vcount(gD2b)
        ecount(gD2b)
        # 
        # # Calculate some node properties and node similarities that will be used to illustrate 
        # # different plotting abilities
        # 
        # # Calculate degree for all nodes
        # degAll <- degree(gD2, v = V(gD2), mode = "all")
        # 
        # # Calculate betweenness for all nodes
        # betAll <- betweenness(gD2, v = V(gD2), directed = FALSE,weights=dataSet$R_sq) / (((vcount(gD2) - 1) * (vcount(gD2)-2)) / 2)
        # betAll.norm <- (betAll - min(betAll))/(max(betAll) - min(betAll))
        # # rm(betAll)
        # 
        # #Calculate Dice similarities between all pairs of nodes
        # #dsAll <- similarity.dice(gD, vids = V(gD), mode = "all")
        
        # Add new node/edge attributes based on the calculated node properties/similarities
        
        # gD2 <- set.vertex.attribute(gD2, "degree", index = V(gD2), value = degAll)
        # gD2 <- set.vertex.attribute(gD2, "betweenness", index = V(gD2), value = betAll.norm)
        V(gD2b)[modulesInSet]$kind="module"
        V(gD2b)[cellsInSet]$kind="flowcell"
        
        list.vertex.attributes(gD2b)
        
        # Check the attributes
        summary(gD2b)
        
        #F1 <- function(x) {data.frame(similarity = dsAll[which(V(gD)$name == as.character(x$module)), which(V(gD)$name == as.character(x$pheno))])}
        #dataSet.ext <- ddply(.data=dataSet, .variables=c("module", "pheno", "weight"), .fun=function(x) data.frame(F1(x)))
        dataSet.ext<-dataSet
        
        gD2b <- set.edge.attribute(gD2b, "weight", index = E(gD2b), value = 0.0)
        gD2b <- set.edge.attribute(gD2b, "R_sq", index = E(gD2b), value = 0.0)
        gD2b <- set.edge.attribute(gD2b, "qvalue", index = E(gD2b), value = 0.0)
        #gD <- set.edge.attribute(gD, "similarity", index = E(gD), value = 0)
        
        # The order of interactions in gD is not the same as it is in dataSet or as it is in the edge list,
        # and for that reason these values cannot be assigned directly
        
        E(gD2b)[as.character(dataSet.ext$module) %--% as.character(dataSet.ext$flowcell)]$R_sq <- dataSet.ext$R_sq
        E(gD2b)[as.character(dataSet.ext$module) %--% as.character(dataSet.ext$flowcell)]$weight <- dataSet$weight
        E(gD2b)[as.character(dataSet.ext$module) %--% as.character(dataSet.ext$flowcell)]$qvalue <- as.numeric(as.character(dataSet$qvalue))
        
        # Check the attributes
        summary(gD2)
        
        #END IGRAPH
        #NETWORK 2b neo4j####
        if (neoinsert=="on"){
          #get the dataframe to use for relationships and nodes
          graphdata=dataSet
          t = suppressMessages(newTransaction(graph))
          
          #remove rows where ratio or p-value are NA or infinite; these are at best hard to interpret, but more likely meaningless results. The full results are still in the exported csv
          flowcellproptests<-flowcellproptests[!is.na(flowcellproptests$ratio),]
          flowcellproptests<-flowcellproptests[!is.na(flowcellproptests$pval),]
          flowcellproptests<-flowcellproptests[flowcellproptests$ratio!=Inf,]
          
          for (therow in 1:nrow(graphdata)) {
            wmod = as.character(graphdata[therow, ]$module)
            cp = as.character(graphdata[therow, ]$flowcell)
            
            #data for some cell types may have been removed above! Some of these reults will be empty.
            ratioX = as.numeric(as.character(flowcellproptests[which(flowcellproptests$names==cp&flowcellproptests$edge==edge),]$ratio))
            diffPX = as.numeric(as.character(flowcellproptests[which(flowcellproptests$names==cp&flowcellproptests$edge==edge),]$pval))
            diffQX = as.numeric(as.character(flowcellproptests[which(flowcellproptests$names==cp&flowcellproptests$edge==edge),]$BH))
            
            #code below only sets ratio, diffP and diffQ properties on flowcellprop nodes if the values exist. Correlation with flowcellprop is still included as result, even if ratio is absent or meaningless.
            if (length(ratioX)==0|length(diffPX)==0|length(diffQX)==0){
              diffMEg = as.numeric(as.character(graphdata[therow,]$diffME))
              query = paste("MERGE (wmod:wgcna {name:{wgcnamod},square:{squareg},edge:{edgeg},contrastvar:{contrastvarsg},contrast:{contrastg}}) ON CREATE SET wmod.diffME = ",diffMEg," ON MATCH SET wmod.diffME = ",diffMEg," MERGE (cp:flowcellprop {name:{flowcellprop},square:{squareg},edge:{edgeg},contrastvar:{contrastvarsg},contrast:{contrastg}}) CREATE (wmod)-[:correlates {weight:{weightg},qvalue:{qval},Rsq:{rsqg},square:{squareg},edge:{edgeg},contrastvar:{contrastvarsg},contrast:{contrastg}}]->(cp)",sep="")
              #print(paste("case1",ratioX,diffPX,diffQX))
              suppressMessages(appendCypher(t, 
                                            query, 
                                            #nodes
                                            wgcnamod = wmod, 
                                            flowcellprop = cp,
                                            #invariant node properties
                                            squareg = dataset.variables[[1]],
                                            edgeg = edge,
                                            contrastvarsg = contrastvars[[edge]],
                                            contrastg = paste(pathways[[edge]][[2]],"-",pathways[[edge]][[1]]),
                                            qval = as.numeric(as.character(graphdata[therow,]$qvalue)),
                                            weightg = as.numeric(as.character(graphdata[therow,]$weight)),
                                            rsqg = as.numeric(as.character(graphdata[therow,]$R_sq))#,
                                            #                ratiog = ratioX,
                                            #                diffPg = diffPX,
                                            #                diffQg = diffQX
              ))
            }else{
              diffMEg = as.numeric(as.character(graphdata[therow,]$diffME))
              query = paste("MERGE (wmod:wgcna {name:{wgcnamod},square:{squareg},edge:{edgeg},contrastvar:{contrastvarsg},contrast:{contrastg}}) ON CREATE SET wmod.diffME = ",diffMEg," ON MATCH SET wmod.diffME = ",diffMEg," MERGE (cp:flowcellprop {name:{flowcellprop},square:{squareg},edge:{edgeg},contrastvar:{contrastvarsg},contrast:{contrastg},ratio:{ratiog},diffP:{diffPg},diffQ:{diffQg}}) CREATE (wmod)-[:correlates {weight:{weightg},qvalue:{qval},Rsq:{rsqg},square:{squareg},edge:{edgeg},contrastvar:{contrastvarsg},contrast:{contrastg}}]->(cp)",sep="")
              #print(paste("case2",ratioX,diffPX,diffQX))
              suppressMessages(appendCypher(t, 
                                            query, 
                                            #nodes
                                            wgcnamod = wmod, 
                                            flowcellprop = cp,
                                            #invariant node properties
                                            squareg = dataset.variables[[1]],
                                            edgeg = edge,
                                            contrastvarsg = contrastvars[[edge]],
                                            contrastg = paste(pathways[[edge]][[2]],"-",pathways[[edge]][[1]]),
                                            qval = as.numeric(as.character(graphdata[therow,]$qvalue)),
                                            weightg = as.numeric(as.character(graphdata[therow,]$weight)),
                                            rsqg = as.numeric(as.character(graphdata[therow,]$R_sq)),
                                            ratiog = ratioX,
                                            diffPg = diffPX,
                                            diffQg = diffQX
              ))
            }
          }
          suppressMessages(commit(t))
          
          #END NEO4J
        }#end neoinsert if
        #NETWORK 2b graphNEL####
        
        gD2b.cyt <- igraph.to.graphNEL(gD2b)
        
        # We have to create attributes for graphNEL
        # We'll keep the same name, so the values are passed from igraph
        # node attributes
        # gD2.cyt <- initNodeAttribute(gD2.cyt, 'label', 'char', "x")
        gD2b.cyt <- initNodeAttribute(gD2b.cyt, 'kind', 'char', "y")
        # gD2.cyt <- initNodeAttribute(gD2.cyt, 'degree', 'numeric', 0) 
        # gD2.cyt <- initNodeAttribute(gD2.cyt, 'betweenness', 'numeric', 0) 
        gD2b.cyt <- initNodeAttribute (gD2b.cyt, "diffME", 'numeric', 0.0)
        # edge attributes
        gD2b.cyt <- initEdgeAttribute (gD2b.cyt, "weight", 'numeric', 0)
        gD2b.cyt <- initEdgeAttribute (gD2b.cyt, "qvalue", 'numeric', 0)
        #gD2.cyt <- initEdgeAttribute (gD2.cyt, "weight2", 'numeric', 0)
        #gD.cyt <- initEdgeAttribute (gD.cyt, "similarity", 'numeric', 0)
        gD2b.cyt <- initEdgeAttribute (gD2b.cyt, "R_sq", 'numeric', 0.0)
        gD2b.cyt <- initEdgeAttribute (gD2b.cyt, "edgeType", "char", "undefined")
        #NETWORK 2b RCytoscape####
        if(cytoscapelink == "on"){
          #open Cytoscape connection
          
          # Now we can create a new graph window in cytoscape
          # Be sure that CytoscapeRPC plugin is activated
          gDCW2b <- new.CytoscapeWindow("NETWORK 2b", graph = gD2b.cyt, overwriteWindow = TRUE,host="192.168.65.2",rpcPort=9000)
          
          ##Send diffME to Cytoscape
          displayGraph(gDCW2b)
          new.nodes<-intModules[as.character(tests.modules$module)%in%as.character(unique(dataSet$module))]
          new.diffME<-unlist(all.difflist.modules[[edge]])[as.character(tests.modules$module)%in%as.character(unique(dataSet$module))]
          
          setNodeAttributesDirect(gDCW2b,"diffME","numeric",new.nodes,new.diffME)
          # getNodeAttribute(gDCW2,"salmon","diffME")
          getAllNodeAttributes(gDCW2b)
          getAllNodes(gDCW2b)
          
          #Colour modules by up or down
          control.points <- c (-1.0, 0.0, 1.0)
          node.colors <- rgb(t(col2rgb(c("black","blue","white","red","green"))),maxColorValue = 255)
          setNodeColorRule(gDCW2b, node.attribute.name='diffME', control.points,node.colors, mode='interpolate')
          redraw (gDCW2b)
          
          # # If you also want to choose a layout from R, a list  of available layouts can be accessed as follows:
          # cy <- CytoscapeConnection(host="192.168.65.2",rpcPort=9000)
          # hlp <-getLayoutNames(cy)
          
          setLayoutProperties (gDCW2b, hlp[6], list (edge_attribute = 'R_sq'))
          layoutNetwork(gDCW2b, hlp[6])
          
          # Finally, we can define rules for node colors, node sizes, and edge colors
          #pheno nodes are red
          cellcol<-rgb(t(col2rgb(c("green"))),maxColorValue = 255)
          setNodeColorDirect(gDCW2b,cells,cellcol)
          
          lockNodeDimensions(gDCW2b,FALSE)
          setNodeWidthDirect(gDCW2b, modules,60)
          setNodeWidthDirect(gDCW2b, cells,50)
          newlabelcol<-rgb(t(col2rgb(c("black"))),maxColorValue = 255)
          setNodeLabelColorDirect(gDCW2b, cells, newlabelcol)
          
          control.points <- c (-1.0, 0.0, 1.0)
          edge.colors <- rgb(t(col2rgb(c("black","blue","white","red","green"))),maxColorValue = 255)
          setEdgeColorRule(gDCW2b, edge.attribute.name='weight', control.points,edge.colors, mode='interpolate')
          redraw (gDCW2b)
        }#end optional cytoscape}
        #End optional network####
      }#end datapoints in flowcellprops if
    }#end NETWORK 2b conditional
    
    
    #NETWORK 2c conditional####
    if (dataset.variables=='berry.test'){
      #NETWORK 2c prep####
      #NETWORK 2c helper####
      # Define numbers of genes and samples
      nGenes = ncol(datExpr);
      nSamples = nrow(datExpr);
      # Recalculate edgewise MEs with color labels and select those with flow info only
      MEs0=bwMEs
      MEs = orderMEs(MEs0)
      ME.edges.flow=list()
      for (tempedge in 1:length(edges)){
        edgeset<-edges[[tempedge]]
        flowset<-which(!is.na(cellprops2$Th))
        edgeflow<-intersect(edgeset,flowset)
        MEs0.edge = MEs[edgeflow,]
        MEs.edge = orderMEs(MEs0.edge)
        ME.edges.flow[[names(edges)[tempedge]]]<-MEs.edge
      }
      cellprops3<-t(as.matrix(coef(decdat1)))[which(!is.na(cellprops2$Th)),]
      
      nSamples=length(edges[[edge]])
      
      #Correlate MEs and cellprops3
      datapoints<-nrow(cellprops3[which(!is.na(cellprops2$Th))%in%edges[[edge]],])
      
      moduleCellCor = cor(ME.edges.flow[[edge]],cellprops3[which(!is.na(cellprops2$Th))%in%edges[[edge]],], use = "p")
      retain<-!apply(moduleCellCor,2,is.na)[1,]
      moduleCellCor = moduleCellCor[,retain] #NA columns removed!
      moduleCellPvalue = corPvalueStudent(moduleCellCor, nSamples);
      moduleCellQvalue = matrix(p.adjust(moduleCellPvalue,method="BH"),nrow=nrow(moduleCellPvalue),dimnames=list(row.names(moduleCellPvalue),colnames(moduleCellPvalue)))
      #edges
      write.csv(moduleCellCor,file=file.path(dir.results,paste('table',q.i,'.',an.count,'.',i,'.',names(edges)[edge],".",names(questions)[[q.i]],'.',analysis,names(datalist)[i],"ModuleCell_withflow_Cor_edge",edge,".csv",sep="")))
      #edge annotation
      write.csv(moduleCellQvalue,file=file.path(dir.results,paste('table',q.i,'.',an.count,'.',i,names(edges)[edge],".",'.',names(questions)[[q.i]],'.',analysis,names(datalist)[i],"ModuleCell_withflow_Cor_qvalue_edge",edge,".csv",sep="")))
      
      r<-moduleCellQvalue[,!apply(moduleCellQvalue,2,is.na)[1,]] 
      s<-moduleCellCor[,!apply(moduleCellCor,2,is.na)[1,]]
      r2<-as.matrix(r,nrow=nrow(r),dimnames=dimnames(r))
      s2<-as.matrix(s,nrow=nrow(s),dimnames=dimnames(s))
      r2[r2>0.05]<-1
      r2[r2<0.05]<-2
      r3<-r2-1
      s2[r3==0]<-0
      
      edgelist2<-list()
      count<-0
      for (ix in 1:nrow(s2)){
        for (j in 1:ncol(s2)){
          if (s2[ix,j]!=0){
            count=count+1
            res<-c(row.names(s2)[ix],colnames(s2)[j],s2[ix,j],r[ix,j])
            edgelist2[[count]]<-res
          }
          
        }
      }
      
      if(length(edgelist2)==0){
        r<-moduleCellQvalue[,!apply(moduleCellQvalue,2,is.na)[1,]] 
        s<-moduleCellCor[,!apply(moduleCellCor,2,is.na)[1,]]
        r2<-as.matrix(r,nrow=nrow(r),dimnames=dimnames(r))
        s2<-as.matrix(s,nrow=nrow(s),dimnames=dimnames(s))
        r2[r2>0.25]<-1
        r2[r2<0.25]<-2
        r3<-r2-1
        s2[r3==0]<-0
        edgelist2<-list()
        count<-0
        for (ix in 1:nrow(s2)){
          for (j in 1:ncol(s2)){
            if (s2[ix,j]!=0){
              count=count+1
              res<-c(row.names(s2)[ix],colnames(s2)[j],s2[ix,j],r[ix,j])
              edgelist2[[count]]<-res
            }
            
          }
        }  
        
      }
      
      edgemm2<-t(matrix(unlist(edgelist2),nrow=4))
      colnames(edgemm2)<-c("module","cellpropF","weight","qvalue")
      newmod2<-apply(edgemm2,1,function(x){x[1]<-substr(x[1],3,nchar(x[1]))})
      edgemm2[,1]<-newmod2
      write.csv(edgemm2,file=file.path(dir.results,paste('table',q.i,'.',an.count,'.',i,'.',names(questions)[[q.i]],'.',analysis,names(datalist)[i],"module_flow_cellprop.csv",sep="")))
      
      cells<-unique(edgemm2[,2])
      cells2<-cbind(cells,rep("cellpropF",length(cells)))
      colnames(cells2)<-c("cellpropF","kind")
      write.csv(cells2,file=file.path(dir.results,paste('table',q.i,'.',an.count,'.',i,'.',names(questions)[[q.i]],'.',analysis,names(datalist)[i],"node_flow_cellprop.csv",sep="")))
      
      #END NETWORK 2 helper
      
      print("igraph 2")
      #NETWORK 2c igraph####
      
      #BEGIN igraph 2
      
      # Read a data set. 
      # Data format: dataframe with 3 variables; variables 1 & 2 correspond to interactions; variable 3 corresponds to the weight of interaction
      dataSet0 <- data.frame(edgemm2)
      dataSet0<-merge(dataSet0,diffMEmatrix[,edge], by.x="module",by.y="row.names",all.x=T,all.y=F)
      dataSet0<-transform(dataSet0,weight=as.numeric(as.character(weight)))
      colnames(dataSet0)<-c("module","cellpropF","weight","qvalue","diffME")
      dataSet<-cbind(dataSet0,abs(dataSet0$weight)^2)
      colnames(dataSet)<-c("module","cellpropF","weight","qvalue","diffME","R_sq")
      modulesInSet<-levels(as.factor(dataSet$module))
      cellsInSet<-levels(as.factor(dataSet$cellpropF))
      # Create a graph. Use simplify to ensure that there are no duplicated edges or self loops
      gD2c <- graph.data.frame(dataSet[,c(1,2,3,4,6)], directed=FALSE)
      
      # Print number of nodes and edges
      vcount(gD2c)
      ecount(gD2c)
      # 
      # # Calculate some node properties and node similarities that will be used to illustrate 
      # # different plotting abilities
      # 
      # # Calculate degree for all nodes
      # degAll <- degree(gD2, v = V(gD2), mode = "all")
      # 
      # # Calculate betweenness for all nodes
      # betAll <- betweenness(gD2, v = V(gD2), directed = FALSE,weights=dataSet$R_sq) / (((vcount(gD2) - 1) * (vcount(gD2)-2)) / 2)
      # betAll.norm <- (betAll - min(betAll))/(max(betAll) - min(betAll))
      # # rm(betAll)
      # 
      # #Calculate Dice similarities between all pairs of nodes
      # #dsAll <- similarity.dice(gD, vids = V(gD), mode = "all")
      
      # Add new node/edge attributes based on the calculated node properties/similarities
      
      # gD2 <- set.vertex.attribute(gD2, "degree", index = V(gD2), value = degAll)
      # gD2 <- set.vertex.attribute(gD2, "betweenness", index = V(gD2), value = betAll.norm)
      V(gD2c)[modulesInSet]$kind="module"
      V(gD2c)[cellsInSet]$kind="cellpropF"
      
      list.vertex.attributes(gD2c)
      
      # Check the attributes
      summary(gD2c)
      
      #F1 <- function(x) {data.frame(similarity = dsAll[which(V(gD)$name == as.character(x$module)), which(V(gD)$name == as.character(x$pheno))])}
      #dataSet.ext <- ddply(.data=dataSet, .variables=c("module", "pheno", "weight"), .fun=function(x) data.frame(F1(x)))
      dataSet.ext<-dataSet
      
      gD2c <- set.edge.attribute(gD2c, "weight", index = E(gD2c), value = 0.0)
      gD2c <- set.edge.attribute(gD2c, "R_sq", index = E(gD2c), value = 0.0)
      gD2c <- set.edge.attribute(gD2c, "qvalue", index = E(gD2c), value = 0.0)
      #gD <- set.edge.attribute(gD, "similarity", index = E(gD), value = 0)
      
      # The order of interactions in gD is not the same as it is in dataSet or as it is in the edge list,
      # and for that reason these values cannot be assigned directly
      
      E(gD2c)[as.character(dataSet.ext$module) %--% as.character(dataSet.ext$cellpropF)]$R_sq <- dataSet.ext$R_sq
      E(gD2c)[as.character(dataSet.ext$module) %--% as.character(dataSet.ext$cellpropF)]$weight <- dataSet$weight
      E(gD2c)[as.character(dataSet.ext$module) %--% as.character(dataSet.ext$cellpropF)]$qvalue <- as.numeric(as.character(dataSet$qvalue))
      
      # Check the attributes
      summary(gD2c)
      
      #END IGRAPH
      #NETWORK 2c neo4j####
      if (neoinsert=="on"){
        #get the dataframe to use for relationships and nodes
        graphdata=dataSet
        t = suppressMessages(newTransaction(graph))
        
        #remove rows where ratio or p-value are NA or infinite; these are at best hard to interpret, but more likely meaningless results. The full results are still in the exported csv
        cellprop_withflow_tests<-cellprop_withflow_tests[!is.na(cellprop_withflow_tests$ratio),]
        cellprop_withflow_tests<-cellprop_withflow_tests[!is.na(cellprop_withflow_tests$pval),]
        cellprop_withflow_tests<-cellprop_withflow_tests[cellprop_withflow_tests$ratio!=Inf,]
        
        for (therow in 1:nrow(graphdata)) {
          wmod = as.character(graphdata[therow, ]$module)
          cp = as.character(graphdata[therow, ]$cellpropF)
          
          #data for some cell types may have been removed above! Some of these reults will be empty.
          ratioX = as.numeric(as.character(cellprop_withflow_tests[which(cellprop_withflow_tests$names==cp&cellprop_withflow_tests$edge==edge),]$ratio))
          diffPX = as.numeric(as.character(cellprop_withflow_tests[which(cellprop_withflow_tests$names==cp&cellprop_withflow_tests$edge==edge),]$pval))
          diffQX = as.numeric(as.character(cellprop_withflow_tests[which(cellprop_withflow_tests$names==cp&cellprop_withflow_tests$edge==edge),]$BH))
          
          #code below only sets ratio, diffP and diffQ properties on cellpropFprop nodes if the values exist. Correlation with cellpropFprop is still included as result, even if ratio is absent or meaningless.
          if (length(ratioX)==0|length(diffPX)==0|length(diffQX)==0){
            diffMEg = as.numeric(as.character(graphdata[therow,]$diffME))
            query = paste("MERGE (wmod:wgcna {name:{wgcnamod},square:{squareg},edge:{edgeg},contrastvar:{contrastvarsg},contrast:{contrastg}}) ON CREATE SET wmod.diffME = ",diffMEg," ON MATCH SET wmod.diffME = ",diffMEg," MERGE (cp:cellpropFprop {name:{cellpropFprop},square:{squareg},edge:{edgeg},contrastvar:{contrastvarsg},contrast:{contrastg}}) CREATE (wmod)-[:correlates {weight:{weightg},qvalue:{qval},Rsq:{rsqg},square:{squareg},edge:{edgeg},contrastvar:{contrastvarsg},contrast:{contrastg}}]->(cp)",sep="")
            #print(paste("case1",ratioX,diffPX,diffQX))
            suppressMessages(appendCypher(t, 
                                          query, 
                                          #nodes
                                          wgcnamod = wmod, 
                                          cellpropFprop = cp,
                                          #invariant node properties
                                          squareg = dataset.variables[[1]],
                                          edgeg = edge,
                                          contrastvarsg = contrastvars[[edge]],
                                          contrastg = paste(pathways[[edge]][[2]],"-",pathways[[edge]][[1]]),
                                          qval = as.numeric(as.character(graphdata[therow,]$qvalue)),
                                          weightg = as.numeric(as.character(graphdata[therow,]$weight)),
                                          rsqg = as.numeric(as.character(graphdata[therow,]$R_sq))#,
                                          #                ratiog = ratioX,
                                          #                diffPg = diffPX,
                                          #                diffQg = diffQX
            ))
          }else{
            diffMEg = as.numeric(as.character(graphdata[therow,]$diffME))
            query = paste("MERGE (wmod:wgcna {name:{wgcnamod},square:{squareg},edge:{edgeg},contrastvar:{contrastvarsg},contrast:{contrastg}}) ON CREATE SET wmod.diffME = ",diffMEg," ON MATCH SET wmod.diffME = ",diffMEg," MERGE (cp:cellpropFprop {name:{cellpropFprop},square:{squareg},edge:{edgeg},contrastvar:{contrastvarsg},contrast:{contrastg},ratio:{ratiog},diffP:{diffPg},diffQ:{diffQg}}) CREATE (wmod)-[:correlates {weight:{weightg},qvalue:{qval},Rsq:{rsqg},square:{squareg},edge:{edgeg},contrastvar:{contrastvarsg},contrast:{contrastg}}]->(cp)",sep="")
            #print(paste("case2",ratioX,diffPX,diffQX))
            suppressMessages(appendCypher(t, 
                                          query, 
                                          #nodes
                                          wgcnamod = wmod, 
                                          cellpropFprop = cp,
                                          #invariant node properties
                                          squareg = dataset.variables[[1]],
                                          edgeg = edge,
                                          contrastvarsg = contrastvars[[edge]],
                                          contrastg = paste(pathways[[edge]][[2]],"-",pathways[[edge]][[1]]),
                                          qval = as.numeric(as.character(graphdata[therow,]$qvalue)),
                                          weightg = as.numeric(as.character(graphdata[therow,]$weight)),
                                          rsqg = as.numeric(as.character(graphdata[therow,]$R_sq)),
                                          ratiog = ratioX,
                                          diffPg = diffPX,
                                          diffQg = diffQX
            ))
          }
        }
        suppressMessages(commit(t))
        
        #END NEO4J
      }#end neoinsert if
      #NETWORK 2c graphNEL####
      
      gD2c.cyt <- igraph.to.graphNEL(gD2c)
      
      # We have to create attributes for graphNEL
      # We'll keep the same name, so the values are passed from igraph
      # node attributes
      # gD2.cyt <- initNodeAttribute(gD2.cyt, 'label', 'char', "x")
      gD2c.cyt <- initNodeAttribute(gD2c.cyt, 'kind', 'char', "y")
      # gD2.cyt <- initNodeAttribute(gD2.cyt, 'degree', 'numeric', 0) 
      # gD2.cyt <- initNodeAttribute(gD2.cyt, 'betweenness', 'numeric', 0) 
      gD2c.cyt <- initNodeAttribute (gD2c.cyt, "diffME", 'numeric', 0.0)
      # edge attributes
      gD2c.cyt <- initEdgeAttribute (gD2c.cyt, "weight", 'numeric', 0)
      gD2c.cyt <- initEdgeAttribute (gD2c.cyt, "qvalue", 'numeric', 0)
      #gD2.cyt <- initEdgeAttribute (gD2.cyt, "weight2", 'numeric', 0)
      #gD.cyt <- initEdgeAttribute (gD.cyt, "similarity", 'numeric', 0)
      gD2c.cyt <- initEdgeAttribute (gD2c.cyt, "R_sq", 'numeric', 0.0)
      gD2c.cyt <- initEdgeAttribute (gD2c.cyt, "edgeType", "char", "undefined")
      #NETWORK 2c RCytoscape####
      if(cytoscapelink == "on"){
        #open Cytoscape connection
        
        # Now we can create a new graph window in cytoscape
        # Be sure that CytoscapeRPC plugin is activated
        gDCW2c <- new.CytoscapeWindow("NETWORK 2c", graph = gD2c.cyt, overwriteWindow = TRUE,host="192.168.65.2",rpcPort=9000)
        
        ##Send diffME to Cytoscape
        displayGraph(gDCW2c)
        new.nodes<-intModules[as.character(tests.modules$module)%in%as.character(unique(dataSet$module))]
        new.diffME<-unlist(all.difflist.modules[[edge]])[as.character(tests.modules$module)%in%as.character(unique(dataSet$module))]
        
        setNodeAttributesDirect(gDCW2c,"diffME","numeric",new.nodes,new.diffME)
        # getNodeAttribute(gDCW2,"salmon","diffME")
        getAllNodeAttributes(gDCW2c)
        getAllNodes(gDCW2c)
        
        #Colour modules by up or down
        control.points <- c (-1.0, 0.0, 1.0)
        node.colors <- rgb(t(col2rgb(c("black","blue","white","red","green"))),maxColorValue = 255)
        setNodeColorRule(gDCW2c, node.attribute.name='diffME', control.points,node.colors, mode='interpolate')
        redraw (gDCW2c)
        
        # # If you also want to choose a layout from R, a list  of available layouts can be accessed as follows:
        # cy <- CytoscapeConnection(host="192.168.65.2",rpcPort=9000)
        # hlp <-getLayoutNames(cy)
        
        setLayoutProperties (gDCW2c, hlp[6], list (edge_attribute = 'R_sq'))
        layoutNetwork(gDCW2c, hlp[6])
        
        # Finally, we can define rules for node colors, node sizes, and edge colors
        #pheno nodes are red
        cellcol<-rgb(t(col2rgb(c("green"))),maxColorValue = 255)
        setNodeColorDirect (gDCW2c,cells,cellcol)
        
        lockNodeDimensions(gDCW2c,FALSE)
        setNodeWidthDirect(gDCW2c, modules,60)
        setNodeWidthDirect(gDCW2c, cells,50)
        newlabelcol<-rgb(t(col2rgb(c("black"))),maxColorValue = 255)
        setNodeLabelColorDirect(gDCW2c, cells, newlabelcol)
        
        control.points <- c (-1.0, 0.0, 1.0)
        edge.colors <- rgb(t(col2rgb(c("black","blue","white","red","green"))),maxColorValue = 255)
        setEdgeColorRule(gDCW2c, edge.attribute.name='weight', control.points,edge.colors, mode='interpolate')
        redraw (gDCW2c)
        #End optional network####
      }#end optional cytoscape
    }#end NETWORK 2c conditional
    
    #NETWORK 3 helper####
    print("helper 3")
    reactomeinfo<-modreacdf
    reactome<-as.character(unique(reactomeinfo$Description))
    node.reactome<-data.frame(cbind(reactome,rep("reactomePW",length(reactome))))
    colnames(node.reactome)<-c("reactomePW","kind")
    write.csv(node.reactome,file.path(dir.results,"node_reactomePW.csv"),row.names=FALSE)
    write.csv(reactomeinfo,file=file.path(dir.results,"all_reactomePW.csv"))
    
    print("igraph 3")
    #NETWORK 3 igraph####
    
    #pathway stringency selection
    pdf(paste('fig_',q.i,'.',an.count,'.',figure,'_reactome_p_q_vals.pdf',sep=""),height=pdf.options()$height*1,width=pdf.options()$width*1)
    plot(sort(modreacdf$pvalue),pch=18,col="green",main="P and q values for reactome pathways",ylim=c(0,0.3))
    points(sort(modreacdf$qvalue),col="blue",pch=20)
    points(sort(modreacdf$p.adjust),col="orange",pch=20)
    legend("topleft",c("P","P.adj","q"),fill=c("green","orange","blue"))
    abline(h=0.05,col="red")
    dev.off()
    figure<-figure+1
    
    #Previous code
    # dataSet0 <- modreacdf[modreacdf$module!="grey",c("module","Description","qvalue")]
    # dataSet0$Description<-gsub("~","approx.",gsub(")","*",gsub("(","*",dataSet0$Description,fixed=TRUE),fixed=TRUE),fixed=TRUE)
    # dataSet0<-transform(dataSet0,qvalue=as.numeric(as.character(qvalue)))
    # colnames(dataSet0)<-c("module","ID","qvalue")#we will use qvalue as weight
    # #dataSet0[is.na(dataSet0$qvalue),]<-1
    # dataSet1<-dataSet0[!is.na(dataSet0$qvalue),]#why are qvalues NA???
    # dataSet<-dataSet1[dataSet1$qvalue<0.05,]
    # modulesInSet<-unique(as.character(dataSet$module))
    # reacIDInSet<-unique(as.character(dataSet$ID))
    # # Create a graph. Use simplify to ensure that there are no duplicated edges or self loops
    # gD3 <- simplify(graph.data.frame(dataSet, directed=FALSE))
    
    dataSet0 <- modreacdf.4[modreacdf.4$module!="grey",c("module","Description","newQ")]
    dataSet0$Description<-gsub(",",";;",gsub("~","approx.",gsub(")","*",gsub("(","*",dataSet0$Description,fixed=TRUE),fixed=TRUE),fixed=TRUE),fixed=TRUE)
    dataSet0<-merge(dataSet0,diffMEmatrix[,edge], by.x="module",by.y="row.names",all.x=T,all.y=F)
    colnames(dataSet0)<-c("module","Description","qvalue","diffME")
    dataSet0<-transform(dataSet0,qvalue=as.numeric(as.character(dataSet0$qvalue)))
    
    #NB I deleted "weight" this may break things
    
    colnames(dataSet0)<-c("module","ID","qvalue","diffME")
    #dataSet0[is.na(dataSet0$qvalue),]<-1
    dataSet1<-dataSet0[!is.na(dataSet0$qvalue),]#why are qvalues NA???
    #dataSet<-dataSet1[dataSet1$qvalue<0.05,]#not required here as values are already mtc
    dataSet<-dataSet1
    modulesInSet<-unique(as.character(dataSet$module))
    reacIDInSet<-unique(as.character(dataSet$ID))
    # Create a graph. Use simplify to ensure that there are no duplicated edges or self loops
    gD3 <- graph.data.frame(dataSet[,1:3], directed=FALSE)
    
    # Print number of nodes and edges
    vcount(gD3)
    ecount(gD3)
    
    # Calculate some node properties and node similarities that will be used to illustrate 
    # different plotting abilities
    
    # Calculate degree for all nodes
    #degAll <- degree(gD3, v = V(gD3), mode = "all")
    
    #no betweenness calcualtions: igraph is very crashy for some reason
    # # Calculate betweenness for all nodes (?igraph bug with crash for zero weights)
    # betAll <- betweenness(gD3, v = V(gD3), directed = FALSE,weights=NULL) / (((vcount(gD3) - 1) * (vcount(gD3)-2)) / 2)
    # betAll.norm <- (betAll - min(betAll))/(max(betAll) - min(betAll))
    
    # #Calculate Dice similarities between all pairs of nodes
    # dsAll <- similarity.dice(gD, vids = V(gD), mode = "all")
    
    # Add new node/edge attributes based on the calculated node properties/similarities
    
    modinreac<-modules[modules%in%as.character(dataSet$module)]
    reacID<-unique(as.character(dataSet$ID))
    #gD3 <- set.vertex.attribute(gD3, "degree", index = V(gD3), value = degAll)
    #gD3 <- set.vertex.attribute(gD3, "betweenness", index = V(gD3), value = betAll.norm)
    V(gD3)[modulesInSet]$kind="module"
    #V(gD3)[modinreac]$kind="module"
    V(gD3)[reacIDInSet]$kind="reactomeID"#edit this
    #V(gD3)[reacID]$kind="reactomeID"#edit this
    
    # Check the attributes
    summary(gD3)
    
    #F1 <- function(x) {data.frame(similarity = dsAll[which(V(gD)$name == as.character(x$module)), which(V(gD)$name == as.character(x$pheno))])}
    #dataSet.ext <- ddply(.data=dataSet, .variables=c("module", "pheno", "weight"), .fun=function(x) data.frame(F1(x)))
    dataSet.ext<-dataSet
    
    #gD3 <- set.edge.attribute(gD3, "weight", index = E(gD), value = 0.0)
    gD3 <- set.edge.attribute(gD3, "qvalue", index = E(gD3), value = 0.0)
    #gD <- set.edge.attribute(gD, "similarity", index = E(gD), value = 0)
    
    # The order of interactions in gD is not the same as it is in dataSet or as it is in the edge list,
    # and for that reason these values cannot be assigned directly
    
    #E(gD)[as.character(dataSet.ext$module) %--% as.character(dataSet.ext$pheno)]$R_sq <- dataSet.ext$R_sq
    E(gD3)[as.character(dataSet.ext$module) %--% as.character(dataSet.ext$ID)]$qvalue <- dataSet$qvalue#crash!
    #E(gD)[as.character(dataSet.ext$module) %--% as.character(dataSet.ext$pheno)]$similarity <- as.numeric(dataSet.ext$similarity)
    
    # Check the attributes
    summary(gD3)
    plot(gD3)
    #END IGRAPH
    #NETWORK 3 neo4j####
    if (neoinsert=="on"){
      #get the dataframe to use for relationships and nodes
      graphdata=dataSet
      #network-specific cypher query
      # query = "
      # MERGE (wmod:wgcna {name:{wgcnamod},square:{squareg},edge:{edgeg},contrastvar:{contrastvarsg},contrast:{contrastg},diffME:{diffMEg}})
      # MERGE (reac:reactomePW {name:{reacPW}})
      # CREATE (wmod)-[:enriched {qvalue:{qval},square:{squareg},edge:{edgeg},contrastvar:{contrastvarsg},contrast:{contrastg}}]->(reac)
      # "
      t = suppressMessages(newTransaction(graph))
      
      for (therow in 1:nrow(graphdata)) {
        wmod = as.character(graphdata[therow, ]$module)
        reac = as.character(graphdata[therow, ]$ID)
        diffMEg = as.numeric(as.character(graphdata[therow,]$diffME))
        query = paste("MERGE (wmod:wgcna {name:{wgcnamod},square:{squareg},edge:{edgeg},contrastvar:{contrastvarsg},contrast:{contrastg}}) ON CREATE SET wmod.diffME = ",diffMEg," ON MATCH SET wmod.diffME = ",diffMEg," MERGE (reac:reactomePW {name:{reacPW}}) CREATE (wmod)-[:enriched {qvalue:{qval},square:{squareg},edge:{edgeg},contrastvar:{contrastvarsg},contrast:{contrastg}}]->(reac)",sep="")
        suppressMessages(appendCypher(t, 
                                      query, 
                                      #nodes
                                      wgcnamod = wmod, 
                                      reacPW = reac,
                                      #invariant node properties
                                      squareg = dataset.variables[[1]],
                                      edgeg = edge,
                                      contrastvarsg = contrastvars[[edge]],
                                      contrastg = paste(pathways[[edge]][[2]],"-",pathways[[edge]][[1]]),
                                      qval = as.numeric(as.character(graphdata[therow,]$qvalue))#,
                                      #diffMEg = as.numeric(as.character(graphdata[therow,]$diffME))
        ))
      }
      
      suppressMessages(commit(t))
      
      #END NEO4J
    }#end neoinsert if
    #NETWORK 3 graphNEL####
    
    gD3.cyt <- igraph.to.graphNEL(gD3)
    
    # We have to create attributes for graphNEL
    # We'll keep the same name, so the values are passed from igraph
    # node attributes
    #gD3.cyt <- initNodeAttribute(gD3.cyt, 'label', 'char', "x")
    gD3.cyt <- initNodeAttribute(gD3.cyt, 'kind', 'char', "y")
    #gD3.cyt <- initNodeAttribute(gD3.cyt, 'degree', 'numeric', 0) 
    #gD3.cyt <- initNodeAttribute(gD3.cyt, 'betweenness', 'numeric', 0) 
    # edge attributes
    gD3.cyt <- initEdgeAttribute (gD3.cyt, "weight", 'numeric', 0)
    gD3.cyt <- initEdgeAttribute (gD3.cyt, "qvalue", 'numeric', 0)
    #gD3.cyt <- initEdgeAttribute (gD3.cyt, "R_sq", 'numeric', 0.0)
    gD3.cyt <- initEdgeAttribute (gD3.cyt, "edgeType", "char", "undefined")
    #NETWORK 3 RCytoscape####
    if(cytoscapelink == "on"){
      #open Cytoscape connection
      
      # Now we can create a new graph window in cytoscape
      # Be sure that CytoscapeRPC plugin is activated
      gDCW3 <- new.CytoscapeWindow("NETWORK 3", graph = gD3.cyt, overwriteWindow = TRUE,host="192.168.65.2",rpcPort=9000)
      getAllEdgeAttributes(gDCW3)
      # We can display graph, with defaults color/size scheme
      displayGraph(gDCW3)
      
      # # If you also want to choose a layout from R, a list  of available layouts can be accessed as follows:
      # cy <- CytoscapeConnection(host="192.168.65.2",rpcPort=9000)
      # hlp <-getLayoutNames(cy)
      # We'll select the layour number 18 - "fruchterman-rheingold" layout 
      # See properties for the given layout
      # getLayoutPropertyNames(cy, hlp[6])
      # Apply values to some of the properties and plot the layout
      #setLayoutProperties (gDCW, hlp[6], list (edge_attribute = 'similarity', iterations = 1000))
      setLayoutProperties (gDCW3, hlp[6], list (edge_attribute = 'qvalue'))
      layoutNetwork(gDCW3, hlp[6])
      
      # # Now, we can define our own default color/size scheme
      # setDefaultBackgroundColor(gDCW3, '#FFFFFF')
      # setDefaultEdgeColor(gDCW3, '#CDC9C9')
      # setDefaultEdgeLineWidth(gDCW3, 1)
      # setDefaultNodeBorderColor(gDCW3, '#000000')
      # setDefaultNodeBorderWidth(gDCW3, 3)
      # setDefaultNodeShape(gDCW3, 'ellipse')
      # setDefaultNodeColor(gDCW3, '#87CEFA')
      # #setDefaultNodeSize(gDCW, 60)
      # #setDefaultNodeFontSize(gDCW, 20)
      # #setDefaultNodeLabelColor(gDCW, '#000000')
      
      # And we can replot it 
      #redraw(gDCW3)       
      
      # Rules for node colors, node sizes, and edge colors
      
      #Rules
      # ##Send diffME to Cytoscape
      # new.nodes<-as.character(tests.modules$module)
      # new.diffME<-tests.modules$diffME
      # setNodeAttributesDirect(gDCW,"diffME","floating",new.nodes,new.diffME)
      
      # #Colour modules by up or down
      # control.points <- c (-1.0, 0.0, 1.0)
      # node.colors <- rgb(t(col2rgb(c("black","blue","white","red","green"))),maxColorValue = 255)
      # setNodeColorRule(gDCW, node.attribute.name='diffME', control.points,node.colors, mode='interpolate')
      # redraw (gDCW)
      
      #reactome nodes are lightyellow
      reactocol<-rgb(t(col2rgb(c("lightyellow"))),maxColorValue = 255)
      setNodeColorDirect (gDCW3,reacID,reactocol)
      redraw (gDCW3)
      
      # #invariant node shapes
      # data.values <- c ("module","pheno")
      # node.shapes<-c("round_rect","ellipse")
      # setNodeShapeRule(gDCW,node.attribute.name="kind",data.values,node.shapes,default.shape="ellipse")
      
      #invariant node dimensions
      lockNodeDimensions(gDCW3,FALSE)
      #node.width <- c(40,20) 
      #setNodeSizeRule(gDCW, 'kind',data.values,node.sizes,mode="lookup")
      #setNodeWidthDirect(gDCW3, modules,60)
      setNodeWidthDirect(gDCW3, reacID,100)
      #newlabelcol<-rgb(t(col2rgb(c("lightblue"))),maxColorValue = 255)
      #setNodeLabelColorDirect(gDCW, phenos, newlabelcol)
      redraw (gDCW3)
      
      # control.points <- c (-1.0, 0.0, 1.0)
      # edge.colors <- rgb(t(col2rgb(c("black","blue","white","red","green"))),maxColorValue = 255)
      # setEdgeColorRule(gDCW, edge.attribute.name='weight', control.points,edge.colors, mode='interpolate')
      # redraw (gDCW)
    }#end optional cytoscape
    
    #NETWORK 4####
    
    print("helper 4")
    
    #I should iterate over all three result sets (BloodLists, HaemAtlas and Abbas)
    #Also, BH adjustment should be used, therefore use unadjusted method in userListEnrichment code
    #Start iterative cellEx####
    if (nrow(bloodResults5$sigOverlaps)>0){
      cellinfolist<-list(bloodResults1,bloodResults3,bloodResults5)
    }else{cellinfolist<-list(bloodResults1,bloodResults3)}
    
    
    for (methlist in 1:length(cellinfolist)){
      #NETWORK 4 helper####
      cellinfo<-cellinfolist[[methlist]]$sigOverlaps[,c(1,2,4)]#Now uncorrected!
      cellinfo.sig<-cellinfolist[[methlist]]$pValues
      newCellNames<-sapply(as.character(cellinfo.sig$UserDefinedCategories),function(x){as.character(strsplit(x,"_")[[1]][1])})
      names(newCellNames)<-NULL
      newCellNames2<-paste(newCellNames,"_xp_",methlist,sep="")
      
      cellinfo.sig$UserDefinedCategories<-newCellNames2
      cell.expr<-as.character(unique(cellinfo.sig$UserDefinedCategories))
      node.cell.expr<-data.frame(cbind(cell.expr,rep("cell_expr",length(cell.expr))))
      colnames(node.cell.expr)<-c("cell_expr","kind")
      write.csv(node.cell.expr,file.path(dir.results,"node_cell_expr.csv"))
      
      print("igraph 4")
      #NETWORK 4 igraph####
      
      # Data format: dataframe with 3 variables; variables 1 & 2 correspond to interactions; variable 3 corresponds to the weight of interaction
      dataSet0 <- cellinfo.sig
      q.values<-p.adjust(cellinfo.sig$Pvalues,method="BH")
      
      #colnames(dataSet0)<-c("module","Description","qvalue","diffME")
      #dataSet0<-transform(dataSet0,qvalue=as.numeric(as.character(dataSet0$qvalue)))
      
      dataSet02<-cbind(cellinfo.sig,q.values)
      dataSet03<-transform(dataSet02,q.values=as.numeric(as.character(q.values)))
      dataSet04<-transform(dataSet03,CorrectedPvalues=as.numeric(as.character(CorrectedPvalues)))
      dataSet05<-dataSet04[dataSet04$q.values<0.05,c(1,2,7)]
      if(nrow(dataSet05)<2){dataSet05<-dataSet04[dataSet04$q.values<1,c(1,2,7)]}
      if(nrow(dataSet05)<2){dataSet05<-dataSet04[dataSet04$Pvalues<1,c(1,2,7)]}#much less stringent
      colnames(dataSet05)<-c("module","cell","qvalue")#we will use qvalue as weight
      dataSet05<-merge(dataSet05,diffMEmatrix[,edge], by.x="module",by.y="row.names",all.x=T,all.y=F)
      dataSet<-dataSet05
      colnames(dataSet)<-c("module","cell","qvalue","diffME")
      dataSet[is.na(dataSet)]<-1
      cell.expr2<-as.character(unique(dataSet$cell))
      modulesInSet<-levels(as.factor(dataSet$module))
      
      # Use simplify to ensure that there are no duplicated edges or self loops
      gD4 <- graph.data.frame(dataSet[,1:3], directed=FALSE)
      
      # Print number of nodes and edges
      vcount(gD4)
      ecount(gD4)
      
      # Add new node/edge attributes based on the calculated node properties/similarities
      
      modincellexp<-modules[modules%in%as.character(dataSet$module)]
      cellexpID<-unique(as.character(dataSet$cell))
      #gD4 <- set.vertex.attribute(gD4, "degree", index = V(gD4), value = degAll)
      #gD4 <- set.vertex.attribute(gD4, "betweenness", index = V(gD4), value = betAll.norm)
      V(gD4)[modincellexp]$kind="module"
      #V(gD4)[modulesInSet]$kind="module"
      V(gD4)[cellexpID]$kind="cellexpID"#edit this
      
      # Check the attributes
      summary(gD4)
      
      #F1 <- function(x) {data.frame(similarity = dsAll[which(V(gD)$name == as.character(x$module)), which(V(gD)$name == as.character(x$pheno))])}
      #dataSet.ext <- ddply(.data=dataSet, .variables=c("module", "pheno", "weight"), .fun=function(x) data.frame(F1(x)))
      dataSet.ext<-dataSet
      
      #gD4 <- set.edge.attribute(gD4, "weight", index = E(gD), value = 0.0)
      gD4 <- set.edge.attribute(gD4, "qvalue", index = E(gD4), value = 0.0)
      #gD <- set.edge.attribute(gD, "similarity", index = E(gD), value = 0)
      
      # The order of interactions in gD is not the same as it is in dataSet or as it is in the edge list,
      # and for that reason these values cannot be assigned directly
      
      #E(gD)[as.character(dataSet.ext$module) %--% as.character(dataSet.ext$pheno)]$R_sq <- dataSet.ext$R_sq
      E(gD4)[as.character(dataSet.ext$module) %--% as.character(dataSet.ext$cell)]$qvalue <- dataSet$qvalue#crash!
      #E(gD)[as.character(dataSet.ext$module) %--% as.character(dataSet.ext$pheno)]$similarity <- as.numeric(dataSet.ext$similarity)
      
      # Check the attributes
      summary(gD4)
      
      #END IGRAPH
      #NETWORK 4 neo4j####
      if (neoinsert=="on"){
        #get the dataframe to use for relationships and nodes
        graphdata=dataSet
        #network-specific cypher query
        # query = "
        # MERGE (wmod:wgcna {name:{wgcnamod},square:{squareg},edge:{edgeg},contrastvar:{contrastvarsg},contrast:{contrastg},diffME:{diffMEg}})
        # MERGE (cellx:cellEx {name:{cellExpr}})
        # CREATE (wmod)-[:enriched {qvalue:{qval},square:{squareg},edge:{edgeg},contrastvar:{contrastvarsg},contrast:{contrastg}}]->(cellx)
        # "
        t = suppressMessages(newTransaction(graph))
        
        for (therow in 1:nrow(graphdata)) {
          wmod = as.character(graphdata[therow, ]$module)
          cellx = as.character(graphdata[therow, ]$cell)
          diffMEg = as.numeric(as.character(graphdata[therow,]$diffME))
          query = paste("MERGE (wmod:wgcna {name:{wgcnamod},square:{squareg},edge:{edgeg},contrastvar:{contrastvarsg},contrast:{contrastg}}) ON CREATE SET wmod.diffME = ",diffMEg," ON MATCH SET wmod.diffME = ",diffMEg," MERGE (cellx:cellEx {name:{cellExpr}}) CREATE (wmod)-[:enriched {qvalue:{qval},square:{squareg},edge:{edgeg},contrastvar:{contrastvarsg},contrast:{contrastg}}]->(cellx)",sep="")
          suppressMessages(appendCypher(t, 
                                        query, 
                                        #nodes
                                        wgcnamod = wmod, #name
                                        cellExpr = cellx, #name of instance 
                                        #invariant node properties
                                        squareg = dataset.variables[[1]],
                                        edgeg = edge,
                                        contrastvarsg = contrastvars[[edge]],
                                        contrastg = paste(pathways[[edge]][[2]],"-",pathways[[edge]][[1]]),
                                        qval = as.numeric(as.character(graphdata[therow,]$qvalue))#,
                                        #diffMEg = as.numeric(as.character(graphdata[therow,]$diffME))
          ))
        }
        
        suppressMessages(commit(t))
        
        #END NEO4J
      }#end neoinsert if
      #NETWORK 4 graphNEL####
      
      gD4.cyt <- igraph.to.graphNEL(gD4)
      
      # We have to create attributes for graphNEL
      # We'll keep the same name, so the values are passed from igraph
      # node attributes
      #gD4.cyt <- initNodeAttribute(gD4.cyt, 'label', 'char', "x")
      gD4.cyt <- initNodeAttribute(gD4.cyt, 'kind', 'char', "y")
      #gD4.cyt <- initNodeAttribute(gD4.cyt, 'degree', 'numeric', 0) 
      #gD4.cyt <- initNodeAttribute(gD4.cyt, 'betweenness', 'numeric', 0) 
      # edge attributes
      gD4.cyt <- initEdgeAttribute (gD4.cyt, "weight", 'numeric', 0)
      gD4.cyt <- initEdgeAttribute (gD4.cyt, "qvalue", 'numeric', 0)
      #gD4.cyt <- initEdgeAttribute (gD4.cyt, "R_sq", 'numeric', 0.0)
      gD4.cyt <- initEdgeAttribute (gD4.cyt, "edgeType", "char", "undefined")
      #NETWORK 4 RCytoscape####
      if(cytoscapelink == "on"){
        #open Cytoscape connection
        
        # Now we can create a new graph window in cytoscape
        # Be sure that CytoscapeRPC plugin is activated
        gDCW4 <- new.CytoscapeWindow("NETWORK 4", graph = gD4.cyt, overwriteWindow = TRUE,host="192.168.65.2",rpcPort=9000)
        getAllEdgeAttributes(gDCW4)
        # We can display graph, with defaults color/size scheme
        displayGraph(gDCW4)
        
        # # If you also want to choose a layout from R, a list  of available layouts can be accessed as follows:
        # cy <- CytoscapeConnection(host="192.168.65.2",rpcPort=9000)
        # hlp <-getLayoutNames(cy)
        # We'll select the layour number 18 - "fruchterman-rheingold" layout 
        # See properties for the given layout
        # getLayoutPropertyNames(cy, hlp[6])
        # Apply values to some of the properties and plot the layout
        #setLayoutProperties (gDCW, hlp[6], list (edge_attribute = 'similarity', iterations = 1000))
        setLayoutProperties (gDCW4, hlp[6], list (edge_attribute = 'qvalue'))
        layoutNetwork(gDCW4, hlp[6])
        
        # # Now, we can define our own default color/size scheme
        # setDefaultBackgroundColor(gDCW4, '#FFFFFF')
        # setDefaultEdgeColor(gDCW4, '#CDC9C9')
        # setDefaultEdgeLineWidth(gDCW4, 1)
        # setDefaultNodeBorderColor(gDCW4, '#000000')
        # setDefaultNodeBorderWidth(gDCW4, 3)
        # setDefaultNodeShape(gDCW4, 'ellipse')
        # setDefaultNodeColor(gDCW4, '#87CEFA')
        # #setDefaultNodeSize(gDCW, 60)
        # #setDefaultNodeFontSize(gDCW, 20)
        # #setDefaultNodeLabelColor(gDCW, '#000000')
        # 
        # # And we can replot it 
        # redraw(gDCW4)       
        
        # Rules for node colors, node sizes, and edge colors
        
        #Rules
        # ##Send diffME to Cytoscape
        # new.nodes<-as.character(tests.modules$module)
        # new.diffME<-tests.modules$diffME
        # setNodeAttributesDirect(gDCW,"diffME","floating",new.nodes,new.diffME)
        
        # #Colour modules by up or down
        # control.points <- c (-1.0, 0.0, 1.0)
        # node.colors <- rgb(t(col2rgb(c("black","blue","white","red","green"))),maxColorValue = 255)
        # setNodeColorRule(gDCW, node.attribute.name='diffME', control.points,node.colors, mode='interpolate')
        # redraw (gDCW)
        
        #cell_exp nodes are light
        cellexpcol<-rgb(t(col2rgb(c("azure"))),maxColorValue = 255)
        setNodeColorDirect (gDCW4,cell.expr2,cellexpcol)
        redraw (gDCW4)
        
        # #invariant node shapes
        # data.values <- c ("module","pheno")
        # node.shapes<-c("round_rect","ellipse")
        # setNodeShapeRule(gDCW,node.attribute.name="kind",data.values,node.shapes,default.shape="ellipse")
        
        #invariant node dimensions
        lockNodeDimensions(gDCW4,FALSE)
        #node.width <- c(40,20) 
        #setNodeSizeRule(gDCW, 'kind',data.values,node.sizes,mode="lookup")
        #setNodeWidthDirect(gDCW4, modules,60)
        setNodeWidthDirect(gDCW4, reacID,100)
        #newlabelcol<-rgb(t(col2rgb(c("lightblue"))),maxColorValue = 255)
        #setNodeLabelColorDirect(gDCW, phenos, newlabelcol)
        redraw (gDCW4)
        
        # control.points <- c (-1.0, 0.0, 1.0)
        # edge.colors <- rgb(t(col2rgb(c("black","blue","white","red","green"))),maxColorValue = 255)
        # setEdgeColorRule(gDCW, edge.attribute.name='weight', control.points,edge.colors, mode='interpolate')
        # redraw (gDCW)
      }#end optional cytoscape
      #End iterative cellEx####
    }#end iterative cellEx
    
    print("helper 5")
    
    #NETWORK 5(it may not exist if moduleList >60)####
    #Start NETWORK 5 if clause####
    if (length(moduleList)<60){
      #NETWORK 5 helper####
      palWanginfo<-PWResults$sigOverlaps[,c(1,2,4)]#Bonferroni-corrected only
      palWanginfo.sig<-PWResults$pValues#Changed to new version
      palWanginfo.sig<-PalWang.select
      
      newPWNames<-sapply(as.character(palWanginfo.sig$UserDefinedCategories),function(x){as.character(strsplit(x,"__")[[1]][1])})
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
      
      palWanginfo.sig$UserDefinedCategories_2<-newPWNames2
      palWanginfo.sig$UserDefinedCategories_long<-newPWNames3
      write.csv(palWanginfo.sig,file.path(dir.results,"PalWangPW_selected.csv"))
      
      
      palWang<-as.character(unique(palWanginfo.sig$UserDefinedCategories_2))
      #palWang2<-sapply(palWang,function(f){gsub(","," ",f)})
      #palWang3<-gsub(","," ",palWang)
      node.palWang<-data.frame(cbind(palWang,rep("palWangPW",length(palWang))))
      colnames(node.palWang)<-c("palWangPW","kind")
      write.csv(node.palWang,file.path(dir.results,"node_palWangPW.csv"))
      
      print("igraph 5")
      #NETWORK 5 igraph####
      
      # Data format: dataframe with 3 variables; variables 1 & 2 correspond to interactions; variable 3 corresponds to the weight of interaction
      dataSet0 <- palWanginfo.sig
      dataSet0$UserDefinedCategories<-gsub(",",";;",gsub(")","*",gsub("(","*",dataSet0$UserDefinedCategories,fixed=TRUE),fixed=TRUE),fixed=TRUE)
      dataSet0$UserDefinedCategories_2<-gsub(",",";;",gsub(")","*",gsub("(","*",dataSet0$UserDefinedCategories_2,fixed=TRUE),fixed=TRUE),fixed=TRUE)
      dataSet0$UserDefinedCategories_long<-gsub(",",";;",gsub("~","approx.",gsub(")","*",gsub("(","*",dataSet0$UserDefinedCategories_long,fixed=TRUE),fixed=TRUE),fixed=TRUE),fixed=TRUE)
      #below is old code, not required with better filtering more proximally
      # q.values<-p.adjust(palWanginfo.sig$Pvalues,method="BH")
      # dataSet02<-cbind(dataSet0,q.values)
      # dataSet03<-transform(dataSet02,q.values=as.numeric(as.character(q.values)))
      # dataSet04<-transform(dataSet03,CorrectedPvalues=as.numeric(as.character(CorrectedPvalues)))
      # dataSet05<-dataSet04[dataSet04$q.values<0.03,c(1,2,7)]#more stringent!
      dataSet05<-dataSet0[,c(1,7,5)]
      colnames(dataSet05)<-c("module","PalWang","qvalue")#we will use qvalue as weight
      
      dataSet05<-merge(dataSet05,diffMEmatrix[,edge], by.x="module",by.y="row.names",all.x=T,all.y=F)
      dataSet<-dataSet05
      colnames(dataSet)<-c("module","PalWang","qvalue","diffME")
      modulesInSet<-unique(as.character(dataSet$module))
      PalWang2<-as.character(unique(dataSet$PalWang))
      
      dataSet[is.na(dataSet)]<-1
      # Create a graph. Use simplify to ensure that there are no duplicated edges or self loops
      gD5 <- graph.data.frame(dataSet[,1:3], directed=FALSE)
      
      # Print number of nodes and edges
      vcount(gD5)
      ecount(gD5)
      
      # Calculate some node properties and node similarities that will be used to illustrate 
      # different plotting abilities
      
      # Calculate degree for all nodes
      #degAll <- degree(gD5, v = V(gD5), mode = "all")
      
      #no betweenness calcualtions: igraph is very crashy for some reason
      # # Calculate betweenness for all nodes (?igraph bug with crash for zero weights)
      # betAll <- betweenness(gD5, v = V(gD5), directed = FALSE,weights=NULL) / (((vcount(gD5) - 1) * (vcount(gD5)-2)) / 2)
      # betAll.norm <- (betAll - min(betAll))/(max(betAll) - min(betAll))
      
      # #Calculate Dice similarities between all pairs of nodes
      # dsAll <- similarity.dice(gD, vids = V(gD), mode = "all")
      
      # Add new node/edge attributes based on the calculated node properties/similarities
      
      modinPW<-modules[modules%in%as.character(dataSet$module)]
      PW.ID<-unique(as.character(dataSet$PalWang))
      #gD5 <- set.vertex.attribute(gD5, "degree", index = V(gD5), value = degAll)
      #gD5 <- set.vertex.attribute(gD5, "betweenness", index = V(gD5), value = betAll.norm)
      V(gD5)[modulesInSet]$kind="module"
      V(gD5)[PalWang2]$kind="PalWang"#edit this
      
      # Check the attributes
      summary(gD5)
      
      #F1 <- function(x) {data.frame(similarity = dsAll[which(V(gD)$name == as.character(x$module)), which(V(gD)$name == as.character(x$pheno))])}
      #dataSet.ext <- ddply(.data=dataSet, .variables=c("module", "pheno", "weight"), .fun=function(x) data.frame(F1(x)))
      dataSet.ext<-dataSet
      
      #gD5 <- set.edge.attribute(gD5, "weight", index = E(gD), value = 0.0)
      gD5 <- set.edge.attribute(gD5, "qvalue", index = E(gD5), value = 0.0)
      #gD <- set.edge.attribute(gD, "similarity", index = E(gD), value = 0)
      
      # The order of interactions in gD is not the same as it is in dataSet or as it is in the edge list,
      # and for that reason these values cannot be assigned directly
      
      #E(gD)[as.character(dataSet.ext$module) %--% as.character(dataSet.ext$pheno)]$R_sq <- dataSet.ext$R_sq
      E(gD5)[as.character(dataSet.ext$module) %--% as.character(dataSet.ext$PalWang)]$qvalue <- dataSet$qvalue#crash!
      #E(gD)[as.character(dataSet.ext$module) %--% as.character(dataSet.ext$pheno)]$similarity <- as.numeric(dataSet.ext$similarity)
      
      # Check the attributes
      summary(gD5)
      gD5s<-simplify(gD5)
      summary(gD5s)
      #END IGRAPH
      #NETWORK 5 neo4j####
      if (neoinsert=="on"){
        #get the dataframe to use for relationships and nodes
        graphdata=dataSet
        #network-specific cypher query
        # query = "
        # MERGE (wmod:wgcna {name:{wgcnamod},square:{squareg},edge:{edgeg},contrastvar:{contrastvarsg},contrast:{contrastg},diffME:{diffMEg}})
        # MERGE (pwm:PalWangPW {name:{pwinstance}})
        # CREATE (wmod)-[:enriched {qvalue:{qval},square:{squareg},edge:{edgeg},contrastvar:{contrastvarsg},contrast:{contrastg}}]->(pwm)
        # "
        t = suppressMessages(newTransaction(graph))
        
        for (therow in 1:nrow(graphdata)) {
          wmod = as.character(graphdata[therow, ]$module)
          pwm = as.character(graphdata[therow, ]$PalWang)
          diffMEg = as.numeric(as.character(graphdata[therow,]$diffME))
          query = paste("MERGE (wmod:wgcna {name:{wgcnamod},square:{squareg},edge:{edgeg},contrastvar:{contrastvarsg},contrast:{contrastg}}) ON CREATE SET wmod.diffME = ",diffMEg," ON MATCH SET wmod.diffME = ",diffMEg," MERGE (pwm:PalWangPW {name:{pwinstance}}) CREATE (wmod)-[:enriched {qvalue:{qval},square:{squareg},edge:{edgeg},contrastvar:{contrastvarsg},contrast:{contrastg}}]->(pwm)",sep="")
          suppressMessages(appendCypher(t, 
                                        query, 
                                        #nodes
                                        wgcnamod = wmod, #name
                                        pwinstance = pwm, #name of instance 
                                        #invariant node properties
                                        squareg = dataset.variables[[1]],
                                        edgeg = edge,
                                        contrastvarsg = contrastvars[[edge]],
                                        contrastg = paste(pathways[[edge]][[2]],"-",pathways[[edge]][[1]]),
                                        qval = as.numeric(as.character(graphdata[therow,]$qvalue))#,
                                        #diffMEg = as.numeric(as.character(graphdata[therow,]$diffME))
          ))
        }
        
        suppressMessages(commit(t))
        
        #END NEO4J
      }#end neoinsert if
      #NETWORK 5 graphNEL####
      
      gD5.cyt <- igraph.to.graphNEL(gD5)
      
      # We have to create attributes for graphNEL
      # We'll keep the same name, so the values are passed from igraph
      # node attributes
      #gD5.cyt <- initNodeAttribute(gD5.cyt, 'label', 'char', "x")
      gD5.cyt <- initNodeAttribute(gD5.cyt, 'kind', 'char', "y")
      gD5.cyt <- initNodeAttribute(gD5.cyt, 'diffME', 'numeric', 0) 
      #gD5.cyt <- initNodeAttribute(gD5.cyt, 'betweenness', 'numeric', 0) 
      # edge attributes
      gD5.cyt <- initEdgeAttribute (gD5.cyt, "weight", 'numeric', 0)
      gD5.cyt <- initEdgeAttribute (gD5.cyt, "qvalue", 'numeric', 0)
      #gD5.cyt <- initEdgeAttribute (gD5.cyt, "R_sq", 'numeric', 0.0)
      gD5.cyt <- initEdgeAttribute (gD5.cyt, "edgeType", "char", "undefined")
      #NETWORK 5 RCytoscape####
      if(cytoscapelink == "on"){
        #open Cytoscape connection
        
        # Now we can create a new graph window in cytoscape
        # Be sure that CytoscapeRPC plugin is activated
        gDCW5 <- new.CytoscapeWindow("NETWORK 5", graph = gD5.cyt, overwriteWindow = TRUE,host="192.168.65.2",rpcPort=9000)
        getAllEdgeAttributes(gDCW5)
        # We can display graph, with defaults color/size scheme
        displayGraph(gDCW5)
        
        # # If you also want to choose a layout from R, a list  of available layouts can be accessed as follows:
        # cy <- CytoscapeConnection(host="192.168.65.2",rpcPort=9000)
        # hlp <-getLayoutNames(cy)
        # We'll select the layour number 18 - "fruchterman-rheingold" layout 
        # See properties for the given layout
        # getLayoutPropertyNames(cy, hlp[6])
        # Apply values to some of the properties and plot the layout
        #setLayoutProperties (gDCW, hlp[6], list (edge_attribute = 'similarity', iterations = 1000))
        setLayoutProperties (gDCW5, hlp[6], list (edge_attribute = 'qvalue'))
        layoutNetwork(gDCW5, hlp[6])
        
        # # Now, we can define our own default color/size scheme
        # setDefaultBackgroundColor(gDCW5, '#FFFFFF')
        # setDefaultEdgeColor(gDCW5, '#CDC9C9')
        # setDefaultEdgeLineWidth(gDCW5, 1)
        # setDefaultNodeBorderColor(gDCW5, '#000000')
        # setDefaultNodeBorderWidth(gDCW5, 3)
        # setDefaultNodeShape(gDCW5, 'ellipse')
        # setDefaultNodeColor(gDCW5, '#87CEFA')
        # #setDefaultNodeSize(gDCW, 60)
        # #setDefaultNodeFontSize(gDCW, 20)
        # #setDefaultNodeLabelColor(gDCW, '#000000')
        # 
        # # And we can replot it 
        # redraw(gDCW5)       
        
        # Rules for node colors, node sizes, and edge colors
        
        #Rules
        # ##Send diffME to Cytoscape
        # new.nodes<-as.character(tests.modules$module)
        # new.diffME<-tests.modules$diffME
        # setNodeAttributesDirect(gDCW,"diffME","floating",new.nodes,new.diffME)
        
        # #Colour modules by up or down
        # control.points <- c (-1.0, 0.0, 1.0)
        # node.colors <- rgb(t(col2rgb(c("black","blue","white","red","green"))),maxColorValue = 255)
        # setNodeColorRule(gDCW, node.attribute.name='diffME', control.points,node.colors, mode='interpolate')
        # redraw (gDCW)
        
        # #invariant node shapes
        # data.values <- c ("module","pheno")
        # node.shapes<-c("round_rect","ellipse")
        # setNodeShapeRule(gDCW,node.attribute.name="kind",data.values,node.shapes,default.shape="ellipse")
        
        #invariant node dimensions
        lockNodeDimensions(gDCW5,FALSE)
        #node.width <- c(40,20) 
        #setNodeSizeRule(gDCW, 'kind',data.values,node.sizes,mode="lookup")
        #setNodeWidthDirect(gDCW5, modules,60)
        setNodeWidthDirect(gDCW5, PW.ID,100)
        #newlabelcol<-rgb(t(col2rgb(c("lightblue"))),maxColorValue = 255)
        #setNodeLabelColorDirect(gDCW, phenos, newlabelcol)
        redraw (gDCW5)
        
        # control.points <- c (-1.0, 0.0, 1.0)
        # edge.colors <- rgb(t(col2rgb(c("black","blue","white","red","green"))),maxColorValue = 255)
        # setEdgeColorRule(gDCW, edge.attribute.name='weight', control.points,edge.colors, mode='interpolate')
        # redraw (gDCW)
      }#end optional cytoscape
      #end NETWORK 5 if statement##########
    }#end NETWORK 5 if statement#####
    
    #NETWORK 6####
    
    print("helper 6")
    #NETWORK 6 helper####
    immuneinfo<-immuneResults1$sigOverlaps[,c(1,2,4)]
    immuneinfo.sig<-immuneResults1$pValues
    immune<-as.character(unique(immuneinfo$UserDefinedCategories))
    node.immune<-data.frame(cbind(immune,rep("immunePW",length(immune))))
    colnames(node.immune)<-c("immunePW","kind")
    write.csv(node.immune,file.path(dir.results,"node_immunePW.csv"))
    
    print("igraph 6")
    #NETWORK 6 igraph####
    
    # Data format: dataframe with 3 variables; variables 1 & 2 correspond to interactions; variable 3 corresponds to the weight of interaction
    dataSet0 <- immuneinfo.sig
    dataSet0$UserDefinedCategories<-gsub(",",";;",gsub(")","*",gsub("(","*",dataSet0$UserDefinedCategories,fixed=TRUE),fixed=TRUE),fixed=TRUE)
    q.values<-p.adjust(immuneinfo.sig$Pvalues,method="BH")
    dataSet02<-cbind(dataSet0,q.values)
    dataSet03<-transform(dataSet02,q.values=as.numeric(as.character(q.values)))
    dataSet04<-transform(dataSet03,CorrectedPvalues=as.numeric(as.character(CorrectedPvalues)))
    dataSet05<-dataSet04[dataSet04$q.values<0.05,c(1,2,7)]
    
    if(nrow(dataSet05)<3){dataSet05<-dataSet04[dataSet04$Pvalues<0.025,c(1,2,7)]}#loosen restrictions
    network6<-TRUE
    if(nrow(dataSet05)<3){network6<-FALSE}
    
    if(network6){####
      dataSet05<-merge(dataSet05,diffMEmatrix[,edge], by.x="InputCategories",by.y="row.names",all.x=T,all.y=F)
      colnames(dataSet05)<-c("module","immunePW","qvalue","diffME")#we will use qvalue as weight
      dataSet<-dataSet05
      colnames(dataSet)<-c("module","immunePW","qvalue","diffME")
      immune2<-as.character(unique(dataSet$immunePW))
      modulesInSet<-unique(as.character(dataSet$module))
      dataSet[is.na(dataSet)]<-1
      # Create a graph. Use simplify to ensure that there are no duplicated edges or self loops
      gD6 <- graph.data.frame(dataSet[,1:3], directed=FALSE)
      
      # Print number of nodes and edges
      vcount(gD6)
      ecount(gD6)
      
      # Calculate some node properties and node similarities that will be used to illustrate 
      # different plotting abilities
      
      # Calculate degree for all nodes
      #degAll <- degree(gD6, v = V(gD6), mode = "all")
      
      #no betweenness calcualtions: igraph is very crashy for some reason
      # # Calculate betweenness for all nodes (?igraph bug with crash for zero weights)
      # betAll <- betweenness(gD6, v = V(gD6), directed = FALSE,weights=NULL) / (((vcount(gD6) - 1) * (vcount(gD6)-2)) / 2)
      # betAll.norm <- (betAll - min(betAll))/(max(betAll) - min(betAll))
      
      # #Calculate Dice similarities between all pairs of nodes
      # dsAll <- similarity.dice(gD, vids = V(gD), mode = "all")
      
      # Add new node/edge attributes based on the calculated node properties/similarities
      
      modinimmune<-modules[modules%in%as.character(dataSet$module)]
      immune.ID<-unique(as.character(dataSet$immunePW))
      #gD6 <- set.vertex.attribute(gD6, "degree", index = V(gD6), value = degAll)
      #gD6 <- set.vertex.attribute(gD6, "betweenness", index = V(gD6), value = betAll.norm)
      V(gD6)[modulesInSet]$kind="module"
      V(gD6)[immune.ID]$kind="immunePW"#edit this
      
      # Check the attributes
      summary(gD6)
      
      #F1 <- function(x) {data.frame(similarity = dsAll[which(V(gD)$name == as.character(x$module)), which(V(gD)$name == as.character(x$pheno))])}
      #dataSet.ext <- ddply(.data=dataSet, .variables=c("module", "pheno", "weight"), .fun=function(x) data.frame(F1(x)))
      dataSet.ext<-dataSet
      
      #gD6 <- set.edge.attribute(gD6, "weight", index = E(gD), value = 0.0)
      gD6 <- set.edge.attribute(gD6, "qvalue", index = E(gD6), value = 0.0)
      #gD <- set.edge.attribute(gD, "similarity", index = E(gD), value = 0)
      
      # The order of interactions in gD is not the same as it is in dataSet or as it is in the edge list,
      # and for that reason these values cannot be assigned directly
      
      #E(gD)[as.character(dataSet.ext$module) %--% as.character(dataSet.ext$pheno)]$R_sq <- dataSet.ext$R_sq
      E(gD6)[as.character(dataSet.ext$module) %--% as.character(dataSet.ext$immunePW)]$qvalue <- dataSet$qvalue#crash!
      #E(gD)[as.character(dataSet.ext$module) %--% as.character(dataSet.ext$pheno)]$similarity <- as.numeric(dataSet.ext$similarity)
      
      # Check the attributes
      summary(gD6)
      gD6s<-simplify(gD6)
      summary(gD6s)
      
      # #igraph analysis code
      # wc <- walktrap.community(gD)
      # modularity(wc)
      # membership(wc)
      # plot(wc, gD)
      
      #END IGRAPH
      #NETWORK 6 neo4j####
      if (neoinsert=="on"){
        #get the dataframe to use for relationships and nodes
        graphdata=dataSet
        #network-specific cypher query
        # query = "
        # MERGE (wmod:wgcna {name:{wgcnamod},square:{squareg},edge:{edgeg},contrastvar:{contrastvarsg},contrast:{contrastg},diffME:{diffMEg}})
        # MERGE (ipw:ImmunePW {name:{ipwinstance}})
        # CREATE (wmod)-[:enriched {qvalue:{qval},square:{squareg},edge:{edgeg},contrastvar:{contrastvarsg},contrast:{contrastg}}]->(ipw)
        # "
        t = suppressMessages(newTransaction(graph))
        
        for (therow in 1:nrow(graphdata)) {
          wmod = as.character(graphdata[therow, ]$module)
          ipw = as.character(graphdata[therow, ]$immunePW)
          diffMEg = as.numeric(as.character(graphdata[therow,]$diffME))
          query = paste("MERGE (wmod:wgcna {name:{wgcnamod},square:{squareg},edge:{edgeg},contrastvar:{contrastvarsg},contrast:{contrastg}}) ON CREATE SET wmod.diffME = ",diffMEg," ON MATCH SET wmod.diffME = ",diffMEg," MERGE (ipw:ImmunePW {name:{ipwinstance}}) CREATE (wmod)-[:enriched {qvalue:{qval},square:{squareg},edge:{edgeg},contrastvar:{contrastvarsg},contrast:{contrastg}}]->(ipw)",sep="")
          suppressMessages(appendCypher(t, 
                                        query, 
                                        #nodes
                                        wgcnamod = wmod, #name
                                        ipwinstance = ipw, #name of instance 
                                        #invariant node properties
                                        squareg = dataset.variables[[1]],
                                        edgeg = edge,
                                        contrastvarsg = contrastvars[[edge]],
                                        contrastg = paste(pathways[[edge]][[2]],"-",pathways[[edge]][[1]]),
                                        qval = as.numeric(as.character(graphdata[therow,]$qvalue))#,
                                        #diffMEg = as.numeric(as.character(graphdata[therow,]$diffME))
          ))
        }
        
        suppressMessages(commit(t))
        
        #END NEO4J
      }#end neoinsert if
      #NETWORK 6 graphNEL####
      
      gD6.cyt <- igraph.to.graphNEL(gD6)
      
      # We have to create attributes for graphNEL
      # We'll keep the same name, so the values are passed from igraph
      # node attributes
      #gD6.cyt <- initNodeAttribute(gD6.cyt, 'label', 'char', "x")
      gD6.cyt <- initNodeAttribute(gD6.cyt, 'kind', 'char', "y")
      gD6.cyt <- initNodeAttribute(gD6.cyt, 'diffME', 'numeric', 0) 
      #gD6.cyt <- initNodeAttribute(gD6.cyt, 'betweenness', 'numeric', 0) 
      # edge attributes
      gD6.cyt <- initEdgeAttribute (gD6.cyt, "weight", 'numeric', 0)
      gD6.cyt <- initEdgeAttribute (gD6.cyt, "qvalue", 'numeric', 0)
      #gD6.cyt <- initEdgeAttribute (gD6.cyt, "R_sq", 'numeric', 0.0)
      gD6.cyt <- initEdgeAttribute (gD6.cyt, "edgeType", "char", "undefined")
      #NETWORK 6 RCytoscape####
      if(cytoscapelink == "on"){
        #open Cytoscape connection
        
        # Now we can create a new graph window in cytoscape
        # Be sure that CytoscapeRPC plugin is activated
        gDCW6 <- new.CytoscapeWindow("NETWORK 6", graph = gD6.cyt, overwriteWindow = TRUE,host="192.168.65.2",rpcPort=9000)
        getAllEdgeAttributes(gDCW6)
        # We can display graph, with defaults color/size scheme
        displayGraph(gDCW6)
        
        # # If you also want to choose a layout from R, a list  of available layouts can be accessed as follows:
        # cy <- CytoscapeConnection(host="192.168.65.2",rpcPort=9000)
        # hlp <-getLayoutNames(cy)
        # We'll select the layour number 18 - "fruchterman-rheingold" layout 
        # See properties for the given layout
        # getLayoutPropertyNames(cy, hlp[6])
        # Apply values to some of the properties and plot the layout
        #setLayoutProperties (gDCW, hlp[6], list (edge_attribute = 'similarity', iterations = 1000))
        setLayoutProperties (gDCW6, hlp[6], list (edge_attribute = 'qvalue'))
        layoutNetwork(gDCW6, hlp[6])
        
        # # Now, we can define our own default color/size scheme
        # setDefaultBackgroundColor(gDCW6, '#FFFFFF')
        # setDefaultEdgeColor(gDCW6, '#CDC9C9')
        # setDefaultEdgeLineWidth(gDCW6, 1)
        # setDefaultNodeBorderColor(gDCW6, '#000000')
        # setDefaultNodeBorderWidth(gDCW6, 3)
        # setDefaultNodeShape(gDCW6, 'ellipse')
        # setDefaultNodeColor(gDCW6, '#87CEFA')
        # #setDefaultNodeSize(gDCW, 60)
        # #setDefaultNodeFontSize(gDCW, 20)
        # #setDefaultNodeLabelColor(gDCW, '#000000')
        # 
        # # And we can replot it 
        # redraw(gDCW6)       
        
        # Rules for node colors, node sizes, and edge colors
        
        #Rules
        # ##Send diffME to Cytoscape
        # new.nodes<-as.character(tests.modules$module)
        # new.diffME<-tests.modules$diffME
        # setNodeAttributesDirect(gDCW,"diffME","floating",new.nodes,new.diffME)
        
        # #Colour modules by up or down
        # control.points <- c (-1.0, 0.0, 1.0)
        # node.colors <- rgb(t(col2rgb(c("black","blue","white","red","green"))),maxColorValue = 255)
        # setNodeColorRule(gDCW, node.attribute.name='diffME', control.points,node.colors, mode='interpolate')
        # redraw (gDCW)
        
        # #invariant node shapes
        # data.values <- c ("module","pheno")
        # node.shapes<-c("round_rect","ellipse")
        # setNodeShapeRule(gDCW,node.attribute.name="kind",data.values,node.shapes,default.shape="ellipse")
        
        #invariant node dimensions
        lockNodeDimensions(gDCW6,FALSE)
        #node.width <- c(40,20) 
        #setNodeSizeRule(gDCW, 'kind',data.values,node.sizes,mode="lookup")
        #setNodeWidthDirect(gDCW6, modules,60)
        setNodeWidthDirect(gDCW6, immune.ID,100)
        #newlabelcol<-rgb(t(col2rgb(c("lightblue"))),maxColorValue = 255)
        #setNodeLabelColorDirect(gDCW, phenos, newlabelcol)
        redraw (gDCW6)
        
        # control.points <- c (-1.0, 0.0, 1.0)
        # edge.colors <- rgb(t(col2rgb(c("black","blue","white","red","green"))),maxColorValue = 255)
        # setEdgeColorRule(gDCW, edge.attribute.name='weight', control.points,edge.colors, mode='interpolate')
        # redraw (gDCW)
      }#end optional cytoscape
    }#end if(network6) #network 6 code is evaluated only if it exists
    
    
    #NETWORK 8####
    print("helper 8")
    #NETWORK 8 helper####
    reactomeinfo<-Baylor_modreacdf
    reactome<-as.character(unique(reactomeinfo$Description))
    node.reactome<-data.frame(cbind(reactome,rep("reactomePW",length(reactome))))
    colnames(node.reactome)<-c("reactomePW","kind")
    write.csv(node.reactome,file.path(dir.results,"Baylor_node_reactomePW.csv"),row.names=FALSE)
    write.csv(reactomeinfo,file=file.path(dir.results,"Baylor_all_reactomePW.csv"))
    
    
    baylor_modules<-as.character(names(modlist4.2))
    node.modules<-data.frame(cbind(baylor_modules),rep("module",length(baylor_modules)))
    colnames(node.modules)<-c("Baylormod","kind")
    write.csv(node.modules,file=file.path(dir.results,paste('table',q.i,'.',an.count,'.',i,'.',names(questions)[[q.i]],'.',analysis,names(datalist)[i],"baylor_node_module.csv",sep="")))
    
    print("igraph 8")
    #NETWORK 8 igraph####
    
    dataSet0 <- Baylor_modreacdf.4[Baylor_modreacdf.4$module!="grey",c("module","Description","newQ")]
    dataSet0$Description<-gsub(",",";;",gsub("~","approx.",gsub(")","*",gsub("(","*",dataSet0$Description,fixed=TRUE),fixed=TRUE),fixed=TRUE),fixed=TRUE)
    dataSet0<-transform(dataSet0,newQ=as.numeric(as.character(newQ)))
    colnames(dataSet0)<-c("Baylormod","ID","qvalue")#we will use qvalue as weight
    #dataSet0[is.na(dataSet0$qvalue),]<-1
    dataSet1<-dataSet0[!is.na(dataSet0$qvalue),]#why are qvalues NA???
    #dataSet<-dataSet1[dataSet1$qvalue<0.05,]#not required here as values are already mtc
    
    dataSet1<-merge(dataSet1,matBaylor[,edge], by.x="Baylormod",by.y="row.names",all.x=T,all.y=F)
    colnames(dataSet1)<-c("Baylormod","ID","qvalue","diffEX")
    dataSet<-dataSet1
    
    modulesInSet<-unique(as.character(dataSet$module))
    reacIDInSet<-unique(as.character(dataSet$ID))
    # Create a graph. Use simplify to ensure that there are no duplicated edges or self loops
    gD8 <- graph.data.frame(dataSet[,1:3], directed=FALSE)
    
    # Print number of nodes and edges
    vcount(gD8)
    ecount(gD8)
    
    # Add new node/edge attributes based on the calculated node properties/similarities
    modinreac<-unique(dataSet$module)
    reacID<-unique(as.character(dataSet$ID))
    V(gD8)[modulesInSet]$kind="Baylormod"
    V(gD8)[reacIDInSet]$kind="reactomeID"#edit this
    
    
    # Check the attributes
    summary(gD8)
    
    #F1 <- function(x) {data.frame(similarity = dsAll[which(V(gD)$name == as.character(x$module)), which(V(gD)$name == as.character(x$pheno))])}
    #dataSet.ext <- ddply(.data=dataSet, .variables=c("module", "pheno", "weight"), .fun=function(x) data.frame(F1(x)))
    dataSet.ext<-dataSet
    gD8 <- set.edge.attribute(gD8, "qvalue", index = E(gD8), value = 0.0)
    
    # The order of interactions in gD is not the same as it is in dataSet or as it is in the edge list,
    # and for that reason these values cannot be assigned directly
    
    #E(gD)[as.character(dataSet.ext$module) %--% as.character(dataSet.ext$pheno)]$R_sq <- dataSet.ext$R_sq
    E(gD8)[as.character(dataSet.ext$Baylormod) %--% as.character(dataSet.ext$ID)]$qvalue <- dataSet$qvalue#crash!
    #E(gD)[as.character(dataSet.ext$module) %--% as.character(dataSet.ext$pheno)]$similarity <- as.numeric(dataSet.ext$similarity)
    
    # Check the attributes
    summary(gD8)
    plot(gD8)
    #END IGRAPH
    #NETWORK 8 neo4j####
    if (neoinsert=="on"){
      #get the dataframe to use for relationships and nodes
      graphdata=dataSet
      #network-specific cypher query
      query = "
      MERGE (bmod:baylor {name:{baylormod},square:{squareg},edge:{edgeg},contrastvar:{contrastvarsg},contrast:{contrastg},diffEX:{diffEXg}})
      MERGE (reac:reactomePW {name:{reacPW}})
      CREATE (bmod)-[:enriched {qvalue:{qval},square:{squareg},edge:{edgeg},contrastvar:{contrastvarsg},contrast:{contrastg}}]->(reac)
      "
      t = suppressMessages(newTransaction(graph))
      
      for (therow in 1:nrow(graphdata)) {
        bmod = as.character(graphdata[therow, ]$Baylormod)
        reac = as.character(graphdata[therow, ]$ID)
        
        suppressMessages(appendCypher(t, 
                                      query, 
                                      #nodes
                                      baylormod = bmod, 
                                      reacPW = reac,
                                      #invariant node properties
                                      squareg = dataset.variables[[1]],
                                      edgeg = edge,
                                      contrastvarsg = contrastvars[[edge]],
                                      contrastg = paste(pathways[[edge]][[2]],"-",pathways[[edge]][[1]]),
                                      qval = as.numeric(as.character(graphdata[therow,]$qvalue)),
                                      diffEXg = as.numeric(as.character(graphdata[therow,]$diffEX))
        ))
      }
      
      suppressMessages(commit(t))
      
      #END NEO4J
    }#end neoinsert if
    #NETWORK 8 graphNEL####
    
    gD8.cyt<-igraph.to.graphNEL(gD8)
    
    # We have to create attributes for graphNEL
    # We'll keep the same name, so the values are passed from igraph
    # node attributes
    
    gD8.cyt <- initNodeAttribute(gD8.cyt, 'kind', 'char', "y")
    gD8.cyt <- initNodeAttribute(gD8.cyt, 'diffEX', 'numeric', 0) 
    #gD3.cyt <- initNodeAttribute(gD3.cyt, 'betweenness', 'numeric', 0)
    #gD7.cyt <- initNodeAttribute(gD7.cyt, 'diffME', 'numeric', 0.0) 
    # edge attributes
    gD8.cyt <- initEdgeAttribute (gD8.cyt, "weight", 'numeric', 0)
    gD8.cyt <- initEdgeAttribute (gD8.cyt, "qvalue", 'numeric', 0)
    #gD3.cyt <- initEdgeAttribute (gD3.cyt, "R_sq", 'numeric', 0.0)
    gD8.cyt <- initEdgeAttribute (gD8.cyt, "edgeType", "char", "undefined")
    #NETWORK 8 RCytoscape####
    if(cytoscapelink == "on"){
      #open Cytoscape connection
      
      # Now we can create a new graph window in cytoscape
      # Be sure that CytoscapeRPC plugin is activated
      gDCW8 <- new.CytoscapeWindow("NETWORK 8", graph = gD8.cyt, overwriteWindow = TRUE,host="192.168.65.2",rpcPort=9000)
      getAllEdgeAttributes(gDCW8)
      # We can display graph, with defaults color/size scheme
      displayGraph(gDCW8)
      
      # # If you also want to choose a layout from R, a list  of available layouts can be accessed as follows:
      # cy <- CytoscapeConnection(host="192.168.65.2",rpcPort=9000)
      # hlp <-getLayoutNames(cy)
      # We'll select the layour number 18 - "fruchterman-rheingold" layout 
      # See properties for the given layout
      # getLayoutPropertyNames(cy, hlp[6])
      # Apply values to some of the properties and plot the layout
      #setLayoutProperties (gDCW, hlp[6], list (edge_attribute = 'similarity', iterations = 1000))
      setLayoutProperties (gDCW8, hlp[6], list (edge_attribute = 'qvalue'))
      layoutNetwork(gDCW8, hlp[6])
      # 
      # # Now, we can define our own default color/size scheme
      # setDefaultBackgroundColor(gDCW7, '#FFFFFF')
      # setDefaultEdgeColor(gDCW7, '#CDC9C9')
      # setDefaultEdgeLineWidth(gDCW7, 1)
      # setDefaultNodeBorderColor(gDCW7, '#000000')
      # setDefaultNodeBorderWidth(gDCW7, 3)
      # setDefaultNodeShape(gDCW7, 'ellipse')
      # setDefaultNodeColor(gDCW7, '#87CEFA')
      # #setDefaultNodeSize(gDCW, 60)
      # #setDefaultNodeFontSize(gDCW, 20)
      # #setDefaultNodeLabelColor(gDCW, '#000000')
      
      # And we can replot it 
      redraw(gDCW8)       
      
      # Rules for node colors, node sizes, and edge colors
      # 
      # #Rules
      ##Send diffME to Cytoscape
      new.nodes<-unlist(diffExBaylor.edges[[edge]])[,1]
      new.diffEX<-unlist(diffExBaylor.edges[[edge]])[,2]
      setNodeAttributesDirect(gDCW8,"diffEX","floating",new.nodes,new.diffEX)
      
      #Colour modules by up or down
      control.points <- c (-1.0, 0.0, 1.0)
      node.colors <- rgb(t(col2rgb(c("black","blue","white","red","green"))),maxColorValue = 255)
      setNodeColorRule(gDCW8, node.attribute.name='diffEX', control.points,node.colors, mode='interpolate')
      redraw (gDCW8)
      
      #reactome nodes are lightyellow
      reactocol<-rgb(t(col2rgb(c("lightyellow"))),maxColorValue = 255)
      setNodeColorDirect (gDCW8,reacID,reactocol)
      redraw (gDCW8)
      
      # #invariant node shapes
      # data.values <- c ("module","pheno")
      # node.shapes<-c("round_rect","ellipse")
      # setNodeShapeRule(gDCW,node.attribute.name="kind",data.values,node.shapes,default.shape="ellipse")
      
      #invariant node dimensions
      lockNodeDimensions(gDCW8,FALSE)
      #node.width <- c(40,20) 
      #setNodeSizeRule(gDCW, 'kind',data.values,node.sizes,mode="lookup")
      #setNodeWidthDirect(gDCW3, modules,60)
      setNodeWidthDirect(gDCW8, reacID,100)
      #newlabelcol<-rgb(t(col2rgb(c("lightblue"))),maxColorValue = 255)
      #setNodeLabelColorDirect(gDCW, phenos, newlabelcol)
      redraw (gDCW8)
      
      # control.points <- c (-1.0, 0.0, 1.0)
      # edge.colors <- rgb(t(col2rgb(c("black","blue","white","red","green"))),maxColorValue = 255)
      # setEdgeColorRule(gDCW, edge.attribute.name='weight', control.points,edge.colors, mode='interpolate')
      # redraw (gDCW)
    }#end optional cytoscape
    
    
    #NETWORK 9####
    
    cellinfolist=list(Baylor_bloodResults1,Baylor_bloodResults3,Baylor_bloodResults5)
    #cellinfolist2=list(bloodResults2,bloodResults4,bloodResults6)
    
    print("helper 9")
    for (methlist in 1:length(cellinfolist)){
      #NETWORK 9 helper####
      cellinfo<-cellinfolist[[methlist]]$sigOverlaps[,c(1,2,4)]#not Bonferroni-corrected only
      cellinfo.sig<-cellinfolist[[methlist]]$pValues
      newCellNames<-sapply(as.character(cellinfo.sig$UserDefinedCategories),function(x){as.character(strsplit(x,"_")[[1]][1])})
      names(newCellNames)<-NULL
      newCellNames2<-paste(newCellNames,"_xp_",methlist,sep="")
      cellinfo.sig$UserDefinedCategories<-newCellNames2
      # cellinfo<-read.csv("/Users/armindeffur/Desktop/tmp/Analysis_printing/Version 6_TBPC_compartment/edges/mod_cell_exp.csv")
      cell.expr<-as.character(unique(cellinfo.sig$UserDefinedCategories))
      node.cell.expr<-data.frame(cbind(cell.expr,rep("cell_expr",length(cell.expr))))
      colnames(node.cell.expr)<-c("cell_expr","kind")
      write.csv(node.cell.expr,file.path(dir.results,"Baylor_node_cell_expr.csv"))
      #NETWORK 9 igraph####
      
      # Data format: dataframe with 3 variables; variables 1 & 2 correspond to interactions; variable 3 corresponds to the weight of interaction
      dataSet0 <- cellinfo.sig
      q.values<-p.adjust(cellinfo.sig$Pvalues,method="BH")
      dataSet02<-cbind(cellinfo.sig,q.values)
      dataSet03<-transform(dataSet02,q.values=as.numeric(as.character(q.values)))
      dataSet04<-transform(dataSet03,CorrectedPvalues=as.numeric(as.character(CorrectedPvalues)))
      dataSet05<-dataSet04[dataSet04$q.values<0.05,c(1,2,7)]
      colnames(dataSet05)<-c("Baylormod","cell","qvalue")#we will use qvalue as weight
      
      dataSet05<-merge(dataSet05,matBaylor[,edge], by.x="Baylormod",by.y="row.names",all.x=T,all.y=F)
      colnames(dataSet05)<-c("Baylormod","cell","qvalue","diffEX")
      dataSet<-dataSet05
      dataSet[is.na(dataSet)]<-1
      cell.expr2<-as.character(unique(dataSet$cell))
      modulesInSet<-unique(as.character(dataSet$Baylormod))
      
      # Create a graph. Use simplify to ensure that there are no duplicated edges or self loops
      #gD9<- simplify(graph.data.frame(dataSet, directed=FALSE))
      gD9<- graph.data.frame(dataSet[,1:3], directed=FALSE)#allow duplicate edges. see if this breaks things.
      
      # Print number of nodes and edges
      vcount(gD9)
      ecount(gD9)
      
      # Add new node/edge attributes based on the calculated node properties/similarities
      
      modincellexp<-as.character(unique(dataSet$Baylormod))
      cellexpID<-unique(as.character(dataSet$cell))
      V(gD9)[modincellexp]$kind="Baylormod"
      V(gD9)[cellexpID]$kind="cellexpID"#edit this
      
      # Check the attributes
      summary(gD9)
      
      #F1 <- function(x) {data.frame(similarity = dsAll[which(V(gD)$name == as.character(x$module)), which(V(gD)$name == as.character(x$pheno))])}
      #dataSet.ext <- ddply(.data=dataSet, .variables=c("module", "pheno", "weight"), .fun=function(x) data.frame(F1(x)))
      dataSet.ext<-dataSet
      
      #gD4 <- set.edge.attribute(gD4, "weight", index = E(gD), value = 0.0)
      gD9 <- set.edge.attribute(gD9, "qvalue", index = E(gD9), value = 0.0)
      #gD <- set.edge.attribute(gD, "similarity", index = E(gD), value = 0)
      
      # The order of interactions in gD is not the same as it is in dataSet or as it is in the edge list,
      # and for that reason these values cannot be assigned directly
      
      #E(gD)[as.character(dataSet.ext$module) %--% as.character(dataSet.ext$pheno)]$R_sq <- dataSet.ext$R_sq
      E(gD9)[as.character(dataSet.ext$cell) %--% as.character(dataSet.ext$Baylormod)]$qvalue <- dataSet$qvalue#crash!
      #E(gD)[as.character(dataSet.ext$module) %--% as.character(dataSet.ext$pheno)]$similarity <- as.numeric(dataSet.ext$similarity)
      
      # Check the attributes
      summary(gD9)
      #NETWORK 9 neo4j####
      if (neoinsert=="on"){
        #get the dataframe to use for relationships and nodes
        graphdata=dataSet
        #network-specific cypher query
        query = "
        MERGE (bmod:baylor {name:{baylormod},square:{squareg},edge:{edgeg},contrastvar:{contrastvarsg},contrast:{contrastg},diffEX:{diffEXg}})
        MERGE (cellx:cellEx {name:{cellExpr}})
        CREATE (bmod)-[:enriched {qvalue:{qval},square:{squareg},edge:{edgeg},contrastvar:{contrastvarsg},contrast:{contrastg}}]->(cellx)
        "
        t = suppressMessages(newTransaction(graph))
        
        for (therow in 1:nrow(graphdata)) {
          bmod = as.character(graphdata[therow, ]$Baylormod)
          cellx = as.character(graphdata[therow, ]$cell)
          
          suppressMessages(appendCypher(t, 
                                        query, 
                                        #nodes
                                        baylormod = bmod, 
                                        cellExpr = cellx, #name of instance 
                                        #invariant node properties
                                        squareg = dataset.variables[[1]],
                                        edgeg = edge,
                                        contrastvarsg = contrastvars[[edge]],
                                        contrastg = paste(pathways[[edge]][[2]],"-",pathways[[edge]][[1]]),
                                        qval = as.numeric(as.character(graphdata[therow,]$qvalue)),
                                        diffEXg = as.numeric(as.character(graphdata[therow,]$diffEX))
          ))
        }
        
        suppressMessages(commit(t))
        
        #END NEO4J
      }#end neoinsert 
      #NETWORK 9 graphNEL####
      
      gD9.cyt <- igraph.to.graphNEL(gD9)
      
      # We have to create attributes for graphNEL
      # We'll keep the same name, so the values are passed from igraph
      # node attributes
      #gD9.cyt <- initNodeAttribute(gD9.cyt, 'label', 'char', "x")
      gD9.cyt <- initNodeAttribute(gD9.cyt, 'kind', 'char', "y")
      gD9.cyt <- initNodeAttribute(gD9.cyt, 'diffEX', 'numeric', 0) 
      #gD4.cyt <- initNodeAttribute(gD4.cyt, 'betweenness', 'numeric', 0) 
      # edge attributes
      gD9.cyt <- initEdgeAttribute (gD9.cyt, "weight", 'numeric', 0)
      gD9.cyt <- initEdgeAttribute (gD9.cyt, "qvalue", 'numeric', 0)
      #gD4.cyt <- initEdgeAttribute (gD4.cyt, "R_sq", 'numeric', 0.0)
      gD9.cyt <- initEdgeAttribute (gD9.cyt, "edgeType", "char", "undefined")
      #NETWORK 9 RCytoscape####
      if(cytoscapelink == "on"){
        #open Cytoscape connection
        
        # Now we can create a new graph window in cytoscape
        # Be sure that CytoscapeRPC plugin is activated
        gDCW9 <- new.CytoscapeWindow("NETWORK 9", graph = gD9.cyt, overwriteWindow = TRUE,host="192.168.65.2",rpcPort=9000)
        getAllEdgeAttributes(gDCW9)
        # We can display graph, with defaults color/size scheme
        displayGraph(gDCW9)
        
        # # If you also want to choose a layout from R, a list  of available layouts can be accessed as follows:
        # cy <- CytoscapeConnection(host="192.168.65.2",rpcPort=9000)
        # hlp <-getLayoutNames(cy)
        # We'll select the layour number 18 - "fruchterman-rheingold" layout 
        # See properties for the given layout
        # getLayoutPropertyNames(cy, hlp[6])
        # Apply values to some of the properties and plot the layout
        #setLayoutProperties (gDCW, hlp[6], list (edge_attribute = 'similarity', iterations = 1000))
        setLayoutProperties (gDCW9, hlp[6], list (edge_attribute = 'qvalue'))
        layoutNetwork(gDCW9, hlp[6])
        
        #Rules
        ##Send diffME to Cytoscape
        new.nodes<-unlist(diffExBaylor.edges[[edge]])[,1]
        new.diffEX<-unlist(diffExBaylor.edges[[edge]])[,2]
        setNodeAttributesDirect(gDCW9,"diffEX","floating",new.nodes,new.diffEX)
        
        #Colour modules by up or down
        control.points <- c (-1.0, 0.0, 1.0)
        node.colors <- rgb(t(col2rgb(c("black","blue","white","red","green"))),maxColorValue = 255)
        setNodeColorRule(gDCW9, node.attribute.name='diffEX', control.points,node.colors, mode='interpolate')
        redraw (gDCW9)
        
        #cell_exp nodes are light
        cellexpcol<-rgb(t(col2rgb(c("azure"))),maxColorValue = 255)
        setNodeColorDirect (gDCW9,cell.expr2,cellexpcol)
        redraw (gDCW9)
        
        # #invariant node shapes
        # data.values <- c ("module","pheno")
        # node.shapes<-c("round_rect","ellipse")
        # setNodeShapeRule(gDCW,node.attribute.name="kind",data.values,node.shapes,default.shape="ellipse")
        
        #invariant node dimensions
        lockNodeDimensions(gDCW9,FALSE)
        #node.width <- c(40,20) 
        #setNodeSizeRule(gDCW, 'kind',data.values,node.sizes,mode="lookup")
        #setNodeWidthDirect(gDCW4, modules,60)
        setNodeWidthDirect(gDCW9, reacID,100)
        #newlabelcol<-rgb(t(col2rgb(c("lightblue"))),maxColorValue = 255)
        #setNodeLabelColorDirect(gDCW, phenos, newlabelcol)
        redraw (gDCW9)
      }#end optional cytoscape
    }#end iterative cellEX for Baylor
    
    #NETWORK 10 helper#### 
    
    palWanginfo<-Baylor_PWResults2$sigOverlaps[,c(1,2,4)]#not Bonferroni-corrected only
    #palWanginfo.sig<-Baylor_PWResults2$pValues#Changed to new version
    palWanginfo.sig<-Baylor_PalWang.select
    
    newPWNames<-sapply(as.character(palWanginfo.sig$UserDefinedCategories),function(x){as.character(strsplit(x,"__")[[1]][1])})
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
    
    palWanginfo.sig$UserDefinedCategories_2<-newPWNames2
    palWanginfo.sig$UserDefinedCategories_long<-newPWNames3
    write.csv(palWanginfo.sig,file.path(dir.results,"Baylor_PalWangPW_selected.csv"))
    
    
    palWang<-as.character(unique(palWanginfo.sig$UserDefinedCategories_2))
    #palWang2<-sapply(palWang,function(f){gsub(","," ",f)})
    #palWang3<-gsub(","," ",palWang)
    node.palWang<-data.frame(cbind(palWang,rep("palWangPW",length(palWang))))
    colnames(node.palWang)<-c("palWangPW","kind")
    write.csv(node.palWang,file.path(dir.results,"Baylor_node_palWangPW.csv"))
    #NETWORK 10 igraph####
    
    # Data format: dataframe with 3 variables; variables 1 & 2 correspond to interactions; variable 3 corresponds to the weight of interaction
    dataSet0 <- palWanginfo.sig
    dataSet0$UserDefinedCategories<-gsub(",",";;",gsub(")","*",gsub("(","*",dataSet0$UserDefinedCategories,fixed=TRUE),fixed=TRUE),fixed=TRUE)
    dataSet0$UserDefinedCategories_2<-gsub(",",";;",gsub(")","*",gsub("(","*",dataSet0$UserDefinedCategories_2,fixed=TRUE),fixed=TRUE),fixed=TRUE)
    dataSet0$UserDefinedCategories_long<-gsub(",",";;",gsub(")","*",gsub("(","*",dataSet0$UserDefinedCategories_long,fixed=TRUE),fixed=TRUE),fixed=TRUE)
    #below is old code, not required with better filtering more proximally
    # q.values<-p.adjust(palWanginfo.sig$Pvalues,method="BH")
    # dataSet02<-cbind(dataSet0,q.values)
    # dataSet03<-transform(dataSet02,q.values=as.numeric(as.character(q.values)))
    # dataSet04<-transform(dataSet03,CorrectedPvalues=as.numeric(as.character(CorrectedPvalues)))
    # dataSet05<-dataSet04[dataSet04$q.values<0.03,c(1,2,7)]#more stringent!
    dataSet05<-dataSet0[,c(1,7,5)]
    colnames(dataSet05)<-c("Baylormod","PalWang","qvalue")#we will use qvalue as weight
    dataSet05<-merge(dataSet05,matBaylor[,edge], by.x="Baylormod",by.y="row.names",all.x=T,all.y=F)
    colnames(dataSet05)<-c("Baylormod","PalWang","qvalue","diffEX")
    
    dataSet<-dataSet05
    modulesInSet<-unique(as.character(dataSet$Baylormod))
    
    colnames(dataSet)<-c("Baylormod","PalWang","qvalue","diffEX")
    PalWang2<-as.character(unique(dataSet$PalWang))
    
    dataSet[is.na(dataSet)]<-1
    # Create a graph. Use simplify to ensure that there are no duplicated edges or self loops
    gD10 <- graph.data.frame(dataSet[,1:3], directed=FALSE)
    
    # Print number of nodes and edges
    vcount(gD10)
    ecount(gD10)
    
    # Add new node/edge attributes based on the calculated node properties/similarities
    
    modinPW<-modules[modules%in%as.character(dataSet$Baylormod)]
    PW.ID<-unique(as.character(dataSet$PalWang))
    V(gD10)[modulesInSet]$kind="Baylormod"
    V(gD10)[PalWang2]$kind="PalWang"#edit this
    
    # Check the attributes
    summary(gD10)
    
    #F1 <- function(x) {data.frame(similarity = dsAll[which(V(gD)$name == as.character(x$module)), which(V(gD)$name == as.character(x$pheno))])}
    #dataSet.ext <- ddply(.data=dataSet, .variables=c("module", "pheno", "weight"), .fun=function(x) data.frame(F1(x)))
    dataSet.ext<-dataSet
    
    #gD5 <- set.edge.attribute(gD5, "weight", index = E(gD), value = 0.0)
    gD10 <- set.edge.attribute(gD10, "qvalue", index = E(gD10), value = 0.0)
    #gD <- set.edge.attribute(gD, "similarity", index = E(gD), value = 0)
    
    # The order of interactions in gD is not the same as it is in dataSet or as it is in the edge list,
    # and for that reason these values cannot be assigned directly
    
    #E(gD)[as.character(dataSet.ext$module) %--% as.character(dataSet.ext$pheno)]$R_sq <- dataSet.ext$R_sq
    E(gD10)[as.character(dataSet.ext$Baylormod) %--% as.character(dataSet.ext$PalWang)]$qvalue <- dataSet$qvalue#crash!
    #E(gD)[as.character(dataSet.ext$module) %--% as.character(dataSet.ext$pheno)]$similarity <- as.numeric(dataSet.ext$similarity)
    
    # Check the attributes
    summary(gD10)
    gD5s<-simplify(gD10)
    summary(gD10)
    #END IGRAPH
    #NETWORK 10 neo4j####
    if (neoinsert=="on"){
      #get the dataframe to use for relationships and nodes
      graphdata=dataSet
      #network-specific cypher query
      query = "
      MERGE (bmod:baylor {name:{baylormod},square:{squareg},edge:{edgeg},contrastvar:{contrastvarsg},contrast:{contrastg},diffEX:{diffEXg}})
      MERGE (pwm:PalWangPW {name:{pwinstance}})
      CREATE (bmod)-[:enriched {qvalue:{qval},square:{squareg},edge:{edgeg},contrastvar:{contrastvarsg},contrast:{contrastg}}]->(pwm)
      "
      t = suppressMessages(newTransaction(graph))
      
      for (therow in 1:nrow(graphdata)) {
        bmod = as.character(graphdata[therow, ]$Baylormod)
        pwm = as.character(graphdata[therow, ]$PalWang)
        
        suppressMessages(appendCypher(t, 
                                      query, 
                                      #nodes
                                      baylormod = bmod, 
                                      pwinstance = pwm, #name of instance 
                                      #invariant node properties
                                      squareg = dataset.variables[[1]],
                                      edgeg = edge,
                                      contrastvarsg = contrastvars[[edge]],
                                      contrastg = paste(pathways[[edge]][[2]],"-",pathways[[edge]][[1]]),
                                      qval = as.numeric(as.character(graphdata[therow,]$qvalue)),
                                      diffEXg = as.numeric(as.character(graphdata[therow,]$diffEX))
        ))
      }
      
      suppressMessages(commit(t))
      
      #END NEO4J
    }#end neoinsert if
    #NETWORK 10 graphNEL####
    
    gD10.cyt <- igraph.to.graphNEL(gD10)
    
    # We have to create attributes for graphNEL
    # We'll keep the same name, so the values are passed from igraph
    # node attributes
    #gD10.cyt <- initNodeAttribute(gD10.cyt, 'label', 'char', "x")
    gD10.cyt <- initNodeAttribute(gD10.cyt, 'kind', 'char', "y")
    gD10.cyt <- initNodeAttribute(gD10.cyt, 'diffEX', 'numeric', 0) 
    #gD5.cyt <- initNodeAttribute(gD5.cyt, 'betweenness', 'numeric', 0) 
    # edge attributes
    gD10.cyt <- initEdgeAttribute (gD10.cyt, "weight", 'numeric', 0)
    gD10.cyt <- initEdgeAttribute (gD10.cyt, "qvalue", 'numeric', 0)
    #gD5.cyt <- initEdgeAttribute (gD5.cyt, "R_sq", 'numeric', 0.0)
    gD10.cyt <- initEdgeAttribute (gD10.cyt, "edgeType", "char", "undefined")
    #NETWORK 10 RCytoscape####
    if(cytoscapelink == "on"){
      # Now we can create a new graph window in cytoscape
      # Be sure that CytoscapeRPC plugin is activated
      gDCW10 <- new.CytoscapeWindow("NETWORK 10", graph = gD10.cyt, overwriteWindow = TRUE,host="192.168.65.2",rpcPort=9000)
      getAllEdgeAttributes(gDCW10)
      # We can display graph, with defaults color/size scheme
      displayGraph(gDCW10)
      
      # # If you also want to choose a layout from R, a list  of available layouts can be accessed as follows:
      # cy <- CytoscapeConnection(host="192.168.65.2",rpcPort=9000)
      # hlp <-getLayoutNames(cy)
      # We'll select the layour number 18 - "fruchterman-rheingold" layout 
      # See properties for the given layout
      # getLayoutPropertyNames(cy, hlp[6])
      # Apply values to some of the properties and plot the layout
      #setLayoutProperties (gDCW, hlp[6], list (edge_attribute = 'similarity', iterations = 1000))
      setLayoutProperties (gDCW10, hlp[6], list (edge_attribute = 'qvalue'))
      layoutNetwork(gDCW10, hlp[6])
      
      # # Now, we can define our own default color/size scheme
      # setDefaultBackgroundColor(gDCW5, '#FFFFFF')
      # setDefaultEdgeColor(gDCW5, '#CDC9C9')
      # setDefaultEdgeLineWidth(gDCW5, 1)
      # setDefaultNodeBorderColor(gDCW5, '#000000')
      # setDefaultNodeBorderWidth(gDCW5, 3)
      # setDefaultNodeShape(gDCW5, 'ellipse')
      # setDefaultNodeColor(gDCW5, '#87CEFA')
      # #setDefaultNodeSize(gDCW, 60)
      # #setDefaultNodeFontSize(gDCW, 20)
      # #setDefaultNodeLabelColor(gDCW, '#000000')
      # 
      # # And we can replot it 
      # redraw(gDCW5)       
      
      # Rules for node colors, node sizes, and edge colors
      
      #Rules
      ##Send diffEX to Cytoscape
      new.nodes<-unlist(diffExBaylor.edges[[edge]])[,1]
      new.diffEX<-unlist(diffExBaylor.edges[[edge]])[,2]
      setNodeAttributesDirect(gDCW10,"diffEX","floating",new.nodes,new.diffEX)
      
      #Colour modules by up or down
      control.points <- c (-1.0, 0.0, 1.0)
      node.colors <- rgb(t(col2rgb(c("black","blue","white","red","green"))),maxColorValue = 255)
      setNodeColorRule(gDCW10, node.attribute.name='diffEX', control.points,node.colors, mode='interpolate')
      redraw (gDCW10)
      
      
      
      # #invariant node shapes
      # data.values <- c ("module","pheno")
      # node.shapes<-c("round_rect","ellipse")
      # setNodeShapeRule(gDCW,node.attribute.name="kind",data.values,node.shapes,default.shape="ellipse")
      
      #invariant node dimensions
      lockNodeDimensions(gDCW10,FALSE)
      #node.width <- c(40,20) 
      #setNodeSizeRule(gDCW, 'kind',data.values,node.sizes,mode="lookup")
      #setNodeWidthDirect(gDCW5, modules,60)
      setNodeWidthDirect(gDCW10, PW.ID,100)
      #newlabelcol<-rgb(t(col2rgb(c("lightblue"))),maxColorValue = 255)
      #setNodeLabelColorDirect(gDCW, phenos, newlabelcol)
      redraw (gDCW10)
    }#end optional cytoscape
    
    #NETWORK 11 helper####
    
    immuneinfo<-Baylor_immuneResults2$sigOverlaps[,c(1,2,4)]
    immuneinfo.sig<-Baylor_immuneResults2$pValues
    immune<-as.character(unique(immuneinfo$UserDefinedCategories))
    node.immune<-data.frame(cbind(immune,rep("immunePW",length(immune))))
    colnames(node.immune)<-c("immunePW","kind")
    write.csv(node.immune,file.path(dir.results,"Baylor_node_immunePW.csv"))
    
    print("igraph 11")
    #NETWORK 11 igraph####
    
    # Data format: dataframe with 3 variables; variables 1 & 2 correspond to interactions; variable 3 corresponds to the weight of interaction
    dataSet0 <- immuneinfo.sig
    dataSet0$UserDefinedCategories<-gsub(",",";;",gsub(")","*",gsub("(","*",dataSet0$UserDefinedCategories,fixed=TRUE),fixed=TRUE),fixed=TRUE)
    q.values<-p.adjust(immuneinfo.sig$Pvalues,method="BH")
    dataSet02<-cbind(dataSet0,q.values)
    dataSet03<-transform(dataSet02,q.values=as.numeric(as.character(q.values)))
    dataSet04<-transform(dataSet03,CorrectedPvalues=as.numeric(as.character(CorrectedPvalues)))
    dataSet05<-dataSet04[dataSet04$q.values<0.05,c(1,2,7)]
    if(nrow(dataSet05)<3){dataSet05<-dataSet04[dataSet04$Pvaluess<0.05,c(1,2,7)]}#loosen restrictions i.e.ignore multiple testing
    
    dataSet05<-merge(dataSet05,matBaylor[,edge], by.x="InputCategories",by.y="row.names",all.x=T,all.y=F)
    dataSet<-dataSet05
    
    colnames(dataSet05)<-c("Baylormod","immunePW","qvalue","diffEX")#we will use qvalue as weight
    dataSet<-dataSet05
    colnames(dataSet)<-c("Baylormod","immunePW","qvalue","diffEX")
    immune2<-as.character(unique(dataSet$immunePW))
    modulesInSet<-unique(as.character(dataSet$Baylormod))
    dataSet[is.na(dataSet)]<-1
    # Create a graph. Use simplify to ensure that there are no duplicated edges or self loops
    gD11 <- graph.data.frame(dataSet[,1:3], directed=FALSE)
    
    # Print number of nodes and edges
    vcount(gD11)
    ecount(gD11)
    
    # Calculate some node properties and node similarities that will be used to illustrate 
    # different plotting abilities
    
    # Calculate degree for all nodes
    #degAll <- degree(gD6, v = V(gD6), mode = "all")
    
    #no betweenness calcualtions: igraph is very crashy for some reason
    # # Calculate betweenness for all nodes (?igraph bug with crash for zero weights)
    # betAll <- betweenness(gD6, v = V(gD6), directed = FALSE,weights=NULL) / (((vcount(gD6) - 1) * (vcount(gD6)-2)) / 2)
    # betAll.norm <- (betAll - min(betAll))/(max(betAll) - min(betAll))
    
    # #Calculate Dice similarities between all pairs of nodes
    # dsAll <- similarity.dice(gD, vids = V(gD), mode = "all")
    
    # Add new node/edge attributes based on the calculated node properties/similarities
    
    modinimmune<-modules[modules%in%as.character(dataSet$Baylormod)]
    immune.ID<-unique(as.character(dataSet$immunePW))
    #gD6 <- set.vertex.attribute(gD6, "degree", index = V(gD6), value = degAll)
    #gD6 <- set.vertex.attribute(gD6, "betweenness", index = V(gD6), value = betAll.norm)
    V(gD11)[modulesInSet]$kind="Baylormod"
    V(gD11)[immune.ID]$kind="immunePW"#edit this
    
    # Check the attributes
    summary(gD11)
    
    #F1 <- function(x) {data.frame(similarity = dsAll[which(V(gD)$name == as.character(x$module)), which(V(gD)$name == as.character(x$pheno))])}
    #dataSet.ext <- ddply(.data=dataSet, .variables=c("module", "pheno", "weight"), .fun=function(x) data.frame(F1(x)))
    dataSet.ext<-dataSet
    
    #gD6 <- set.edge.attribute(gD6, "weight", index = E(gD), value = 0.0)
    gD11 <- set.edge.attribute(gD11, "qvalue", index = E(gD11), value = 0.0)
    #gD <- set.edge.attribute(gD, "similarity", index = E(gD), value = 0)
    
    # The order of interactions in gD is not the same as it is in dataSet or as it is in the edge list,
    # and for that reason these values cannot be assigned directly
    
    #E(gD)[as.character(dataSet.ext$module) %--% as.character(dataSet.ext$pheno)]$R_sq <- dataSet.ext$R_sq
    E(gD11)[as.character(dataSet.ext$Baylormod) %--% as.character(dataSet.ext$immunePW)]$qvalue <- dataSet$qvalue#crash!
    #E(gD)[as.character(dataSet.ext$module) %--% as.character(dataSet.ext$pheno)]$similarity <- as.numeric(dataSet.ext$similarity)
    
    # Check the attributes
    summary(gD11)
    gD11s<-simplify(gD11)
    summary(gD11s)
    
    # #igraph analysis code
    # wc <- walktrap.community(gD)
    # modularity(wc)
    # membership(wc)
    # plot(wc, gD)
    
    #END IGRAPH
    #NETWORK 11 neo4j####
    if (neoinsert=="on"){
      #get the dataframe to use for relationships and nodes
      graphdata=dataSet
      #network-specific cypher query
      query = "
      MERGE (bmod:baylor {name:{baylormod},square:{squareg},edge:{edgeg},contrastvar:{contrastvarsg},contrast:{contrastg},diffEX:{diffEXg}})
      MERGE (ipw:ImmunePW {name:{ipwinstance}})
      CREATE (bmod)-[:enriched {qvalue:{qval},square:{squareg},edge:{edgeg},contrastvar:{contrastvarsg},contrast:{contrastg}}]->(ipw)
      "
      t = suppressMessages(newTransaction(graph))
      
      for (therow in 1:nrow(graphdata)) {
        bmod = as.character(graphdata[therow, ]$Baylormod)
        ipw = as.character(graphdata[therow, ]$immunePW)
        
        suppressMessages(appendCypher(t, 
                                      query, 
                                      #nodes
                                      baylormod = bmod, 
                                      ipwinstance = ipw, #name of instance 
                                      #invariant node properties
                                      squareg = dataset.variables[[1]],
                                      edgeg = edge,
                                      contrastvarsg = contrastvars[[edge]],
                                      contrastg = paste(pathways[[edge]][[2]],"-",pathways[[edge]][[1]]),
                                      qval = as.numeric(as.character(graphdata[therow,]$qvalue)),
                                      diffEXg = as.numeric(as.character(graphdata[therow,]$diffEX))
        ))
      }
      
      suppressMessages(commit(t))
      
      #END NEO4J
    }#end neoinsert if
    #NETWORK 11 graphNEL####
    
    gD11.cyt <- igraph.to.graphNEL(gD11)
    
    # We have to create attributes for graphNEL
    # We'll keep the same name, so the values are passed from igraph
    # node attributes
    #gD11.cyt <- initNodeAttribute(gD11.cyt, 'label', 'char', "x")
    gD11.cyt <- initNodeAttribute(gD11.cyt, 'kind', 'char', "y")
    gD11.cyt <- initNodeAttribute(gD11.cyt, 'diffEX', 'numeric', 0) 
    #gD6.cyt <- initNodeAttribute(gD6.cyt, 'betweenness', 'numeric', 0) 
    # edge attributes
    gD11.cyt <- initEdgeAttribute (gD11.cyt, "weight", 'numeric', 0)
    gD11.cyt <- initEdgeAttribute (gD11.cyt, "qvalue", 'numeric', 0)
    #gD6.cyt <- initEdgeAttribute (gD6.cyt, "R_sq", 'numeric', 0.0)
    gD11.cyt <- initEdgeAttribute (gD11.cyt, "edgeType", "char", "undefined")
    #NETWORK 11 RCytoscape####
    if(cytoscapelink == "on"){
      #open Cytoscape connection
      
      # Now we can create a new graph window in cytoscape
      # Be sure that CytoscapeRPC plugin is activated
      gDCW11 <- new.CytoscapeWindow("NETWORK 11", graph = gD11.cyt, overwriteWindow = TRUE,host="192.168.65.2",rpcPort=9000)
      getAllEdgeAttributes(gDCW6)
      # We can display graph, with defaults color/size scheme
      displayGraph(gDCW11)
      
      # # If you also want to choose a layout from R, a list  of available layouts can be accessed as follows:
      # cy <- CytoscapeConnection(host="192.168.65.2",rpcPort=9000)
      # hlp <-getLayoutNames(cy)
      # We'll select the layour number 18 - "fruchterman-rheingold" layout 
      # See properties for the given layout
      # getLayoutPropertyNames(cy, hlp[6])
      # Apply values to some of the properties and plot the layout
      #setLayoutProperties (gDCW, hlp[6], list (edge_attribute = 'similarity', iterations = 1000))
      setLayoutProperties (gDCW11, hlp[6], list (edge_attribute = 'qvalue'))
      layoutNetwork(gDCW11, hlp[6])
      
      # # Now, we can define our own default color/size scheme
      # setDefaultBackgroundColor(gDCW6, '#FFFFFF')
      # setDefaultEdgeColor(gDCW6, '#CDC9C9')
      # setDefaultEdgeLineWidth(gDCW6, 1)
      # setDefaultNodeBorderColor(gDCW6, '#000000')
      # setDefaultNodeBorderWidth(gDCW6, 3)
      # setDefaultNodeShape(gDCW6, 'ellipse')
      # setDefaultNodeColor(gDCW6, '#87CEFA')
      # #setDefaultNodeSize(gDCW, 60)
      # #setDefaultNodeFontSize(gDCW, 20)
      # #setDefaultNodeLabelColor(gDCW, '#000000')
      # 
      # # And we can replot it 
      # redraw(gDCW6)       
      
      # Rules for node colors, node sizes, and edge colors
      
      #Rules
      ##Send diffME to Cytoscape
      new.nodes<-unlist(diffExBaylor.edges[[edge]])[,1]
      new.diffEX<-unlist(diffExBaylor.edges[[edge]])[,2]
      setNodeAttributesDirect(gDCW11,"diffEX","floating",new.nodes,new.diffEX)
      
      #Colour modules by up or down
      control.points <- c (-1.0, 0.0, 1.0)
      node.colors <- rgb(t(col2rgb(c("black","blue","white","red","green"))),maxColorValue = 255)
      setNodeColorRule(gDCW11, node.attribute.name='diffEX', control.points,node.colors, mode='interpolate')
      redraw (gDCW11)
      
      
      
      # #invariant node shapes
      # data.values <- c ("module","pheno")
      # node.shapes<-c("round_rect","ellipse")
      # setNodeShapeRule(gDCW,node.attribute.name="kind",data.values,node.shapes,default.shape="ellipse")
      
      #invariant node dimensions
      lockNodeDimensions(gDCW11,FALSE)
      #node.width <- c(40,20) 
      #setNodeSizeRule(gDCW, 'kind',data.values,node.sizes,mode="lookup")
      #setNodeWidthDirect(gDCW6, modules,60)
      setNodeWidthDirect(gDCW11, immune.ID,100)
      #newlabelcol<-rgb(t(col2rgb(c("lightblue"))),maxColorValue = 255)
      #setNodeLabelColorDirect(gDCW, phenos, newlabelcol)
      redraw (gDCW11)
      
      # control.points <- c (-1.0, 0.0, 1.0)
      # edge.colors <- rgb(t(col2rgb(c("black","blue","white","red","green"))),maxColorValue = 255)
      # setEdgeColorRule(gDCW, edge.attribute.name='weight', control.points,edge.colors, mode='interpolate')
      # redraw (gDCW)
    }#end optional cytoscape
    
    #NETWORK 12 helper####
    
    modEnrich<-WGCNA_BaylorResults2$sigOverlaps[,c(1,2,4)]
    modEnrich.sig<-WGCNA_BaylorResults2$pValues
    modsenriched<-as.character(unique(modEnrich$UserDefinedCategories))
    
    newmodenrichNames_a<-sapply(as.character(modEnrich$UserDefinedCategories),function(x){as.character(strsplit(x,"__")[[1]][1])})
    names(newmodenrichNames_a)<-NULL
    modEnrich$UserDefinedCategories<-newmodenrichNames_a
    
    newmodenrichNames<-sapply(as.character(modEnrich.sig$UserDefinedCategories),function(x){as.character(strsplit(x,"__")[[1]][1])})
    names(newmodenrichNames)<-NULL
    
    modEnrich.sig$UserDefinedCategories_2<-newmodenrichNames
    
    write.csv(modEnrich.sig,file.path(dir.results,"modEnrich_selected.csv"))
    
    node.modenrich<-data.frame(cbind(modsenriched,rep("BaylorMod",length(modsenriched))))
    colnames(node.modenrich)<-c("BaylorMod","kind")
    write.csv(node.modenrich,file.path(dir.results,"WGCNA_Baylor_node_enrich.csv"))
    
    print("igraph 12")
    #NETWORK 12 igraph####
    dataSet0<-modEnrich.sig
    q.values<-p.adjust(modEnrich.sig$Pvalues,method="BH")
    dataSet02<-cbind(dataSet0,q.values)
    dataSet03<-transform(dataSet02,q.values=as.numeric(as.character(q.values)))
    dataSet04<-transform(dataSet03,CorrectedPvalues=as.numeric(as.character(CorrectedPvalues)))
    dataSet05<-dataSet04[dataSet04$q.values<0.05,c(1,7,8)]
    
    dataSet05<-merge(dataSet05,diffMEmatrix[,edge], by.x="InputCategories",by.y="row.names",all.x=T,all.y=F)
    dataSet05<-merge(dataSet05,matBaylor[,edge], by.x="UserDefinedCategories_2",by.y="row.names",all.x=T,all.y=F)
    
    dataSet<-dataSet05
    colnames(dataSet)<-c("Baylormod","module","qvalue","diffME","diffEX")
    Baylormod2<-as.character(unique(dataSet$Baylormod))
    modulesInSet<-unique(as.character(dataSet$module))
    dataSet[is.na(dataSet)]<-1
    # Create a graph. Use simplify to ensure that there are no duplicated edges or self loops
    gD12 <- graph.data.frame(dataSet[,1:3], directed=FALSE)
    
    # Print number of nodes and edges
    vcount(gD12)
    ecount(gD12)
    
    
    # Add new node/edge attributes based on the calculated node properties/similarities
    
    modinenrich<-modules[modules%in%as.character(dataSet$module)]
    enrich.ID<-unique(as.character(dataSet$Baylormod))
    #gD6 <- set.vertex.attribute(gD6, "degree", index = V(gD6), value = degAll)
    #gD6 <- set.vertex.attribute(gD6, "betweenness", index = V(gD6), value = betAll.norm)
    V(gD12)[modulesInSet]$kind="module"
    V(gD12)[enrich.ID]$kind="Baylormod"#edit this
    
    # Check the attributes
    summary(gD12)
    
    dataSet.ext<-dataSet
    
    gD12 <- set.edge.attribute(gD12, "qvalue", index = E(gD12), value = 0.0)
    
    # The order of interactions in gD is not the same as it is in dataSet or as it is in the edge list,
    # and for that reason these values cannot be assigned directly
    
    E(gD12)[as.character(dataSet.ext$module) %--% as.character(dataSet.ext$Baylormod)]$qvalue <- dataSet$qvalue#crash!
    
    # Check the attributes
    summary(gD12)
    gD12s<-simplify(gD12)
    summary(gD12s)
    
    
    #END IGRAPH
    #NETWORK 12 neo4j####
    if (neoinsert=="on"){
      #get the dataframe to use for relationships and nodes
      graphdata=dataSet
      #network-specific cypher query
      # query = "
      # MERGE (wmod:wgcna {name:{wgcnamod},square:{squareg},edge:{edgeg},contrastvar:{contrastvarsg},contrast:{contrastg},diffME:{diffMEg}})
      # MERGE (bmod:baylor {name:{baylormod},square:{squareg},edge:{edgeg},contrastvar:{contrastvarsg},contrast:{contrastg},diffEX:{diffEXg}})
      # CREATE (wmod)-[:enriched {qvalue:{qval},square:{squareg},edge:{edgeg},contrastvar:{contrastvarsg},contrast:{contrastg}}]->(bmod)
      # "
      t = suppressMessages(newTransaction(graph))
      
      for (therow in 1:nrow(graphdata)) {
        bmod = as.character(graphdata[therow, ]$Baylormod)
        wmod = as.character(graphdata[therow, ]$module)
        diffMEg = as.numeric(as.character(graphdata[therow,]$diffME))
        query = paste("MERGE (wmod:wgcna {name:{wgcnamod},square:{squareg},edge:{edgeg},contrastvar:{contrastvarsg},contrast:{contrastg}}) ON CREATE SET wmod.diffME = ",diffMEg," ON MATCH SET wmod.diffME = ",diffMEg," MERGE (bmod:baylor {name:{baylormod},square:{squareg},edge:{edgeg},contrastvar:{contrastvarsg},contrast:{contrastg},diffEX:{diffEXg}}) CREATE (wmod)-[:enriched {qvalue:{qval},square:{squareg},edge:{edgeg},contrastvar:{contrastvarsg},contrast:{contrastg}}]->(bmod)",sep="")
        suppressMessages(appendCypher(t, 
                                      query, 
                                      #nodes
                                      baylormod = bmod, 
                                      wgcnamod = wmod, #name
                                      #invariant node properties
                                      squareg = dataset.variables[[1]],
                                      edgeg = edge,
                                      contrastvarsg = contrastvars[[edge]],
                                      contrastg = paste(pathways[[edge]][[2]],"-",pathways[[edge]][[1]]),
                                      qval = as.numeric(as.character(graphdata[therow,]$qvalue)),
                                      diffEXg = as.numeric(as.character(graphdata[therow,]$diffEX))#,
                                      #diffMEg = as.numeric(as.character(graphdata[therow,]$diffME))
        ))
      }
      
      suppressMessages(commit(t))
      
      #END NEO4J
    }#end neoinsert if
    #NETWORK 12 graphNEL####
    
    gD12.cyt <- igraph.to.graphNEL(gD12)
    
    # We have to create attributes for graphNEL
    # We'll keep the same name, so the values are passed from igraph
    # node attributes
    #gD12.cyt <- initNodeAttribute(gD12.cyt, 'label', 'char', "x")
    gD12.cyt <- initNodeAttribute(gD12.cyt, 'kind', 'char', "y")
    gD12.cyt <- initNodeAttribute(gD12.cyt, 'diffME', 'numeric', 0) 
    gD12.cyt <- initNodeAttribute(gD12.cyt, 'diffEX', 'numeric', 0) 
    #gD6.cyt <- initNodeAttribute(gD6.cyt, 'betweenness', 'numeric', 0) 
    # edge attributes
    gD12.cyt <- initEdgeAttribute (gD12.cyt, "weight", 'numeric', 0)
    gD12.cyt <- initEdgeAttribute (gD12.cyt, "qvalue", 'numeric', 0)
    #gD6.cyt <- initEdgeAttribute (gD6.cyt, "R_sq", 'numeric', 0.0)
    gD12.cyt <- initEdgeAttribute (gD12.cyt, "edgeType", "char", "undefined")
    #NETWORK 12 RCytoscape####
    if(cytoscapelink == "on"){
      #open Cytoscape connection
      
      # Now we can create a new graph window in cytoscape
      # Be sure that CytoscapeRPC plugin is activated
      gDCW12 <- new.CytoscapeWindow("NETWORK 12", graph = gD12.cyt, overwriteWindow = TRUE,host="192.168.65.2",rpcPort=9000)
      getAllEdgeAttributes(gDCW12)
      # We can display graph, with defaults color/size scheme
      displayGraph(gDCW12)
      
      # # If you also want to choose a layout from R, a list  of available layouts can be accessed as follows:
      # cy <- CytoscapeConnection(host="192.168.65.2",rpcPort=9000)
      # hlp <-getLayoutNames(cy)
      # We'll select the layour number 18 - "fruchterman-rheingold" layout 
      # See properties for the given layout
      # getLayoutPropertyNames(cy, hlp[6])
      # Apply values to some of the properties and plot the layout
      #setLayoutProperties (gDCW, hlp[6], list (edge_attribute = 'similarity', iterations = 1000))
      setLayoutProperties (gDCW12, hlp[6], list (edge_attribute = 'qvalue'))
      layoutNetwork(gDCW12, hlp[6])
      
      # Rules for node colors, node sizes, and edge colors
      
      #Rules
      
      #Send diffME to Cytoscape
      new.nodes<-as.character(tests.modules$module)
      new.diffME<-tests.modules$diffME
      setNodeAttributesDirect(gDCW,"diffME","floating",new.nodes,new.diffME)
      
      ##Send diffEX to Cytoscape
      new.nodes<-unlist(diffExBaylor.edges[[edge]])[,1]
      new.diffEX<-unlist(diffExBaylor.edges[[edge]])[,2]
      setNodeAttributesDirect(gDCW12,"diffEX","floating",new.nodes,new.diffEX)
      
      #Colour modules by up or down
      control.points <- c (-1.0, 0.0, 1.0)
      node.colors <- rgb(t(col2rgb(c("black","blue","white","red","green"))),maxColorValue = 255)
      setNodeColorRule(gDCW12, node.attribute.name='diffME', control.points,node.colors, mode='interpolate')
      setNodeColorRule(gDCW12, node.attribute.name='diffEX', control.points,node.colors, mode='interpolate')
      redraw (gDCW12)
      
      
      
      # #invariant node shapes
      # data.values <- c ("module","pheno")
      # node.shapes<-c("round_rect","ellipse")
      # setNodeShapeRule(gDCW,node.attribute.name="kind",data.values,node.shapes,default.shape="ellipse")
      
      #invariant node dimensions
      lockNodeDimensions(gDCW12,FALSE)
      #node.width <- c(40,20) 
      #setNodeSizeRule(gDCW, 'kind',data.values,node.sizes,mode="lookup")
      #setNodeWidthDirect(gDCW6, modules,60)
      setNodeWidthDirect(gDCW12, enrich.ID,100)
      #newlabelcol<-rgb(t(col2rgb(c("lightblue"))),maxColorValue = 255)
      #setNodeLabelColorDirect(gDCW, phenos, newlabelcol)
      redraw (gDCW12)
      
      # control.points <- c (-1.0, 0.0, 1.0)
      # edge.colors <- rgb(t(col2rgb(c("black","blue","white","red","green"))),maxColorValue = 255)
      # setEdgeColorRule(gDCW, edge.attribute.name='weight', control.points,edge.colors, mode='interpolate')
      # redraw (gDCW)
      
    }#end optional cytoscape
    
    #Merge Networks#####
    if(cytoscapelink == "on"){
      ##
      if (length(moduleList)<60&&network6){
        graphlist<-list(gD2.cyt,gD3.cyt,gD4.cyt,gD5.cyt,gD6.cyt,gD7.cyt,gD8.cyt,gD9.cyt,gD10.cyt,gD11.cyt,gD12.cyt,gD13.cyt)
        windowlist<-list(gDCW2,gDCW3,gDCW4,gDCW5,gDCW6,gDCW7,gDCW8,gDCW9,gDCW10,gDCW11,gDCW12,gDCW13)}
      
      if (length(moduleList)<60&&!network6){
        graphlist<-list(gD2.cyt,gD3.cyt,gD4.cyt,gD5.cyt,gD7.cyt,gD8.cyt,gD9.cyt,gD10.cyt,gD11.cyt,gD12.cyt,gD13.cyt)
        windowlist<-list(gDCW2,gDCW3,gDCW4,gDCW5,gDCW7,gDCW8,gDCW9,gDCW10,gDCW11,gDCW12,gDCW13)}
      
      if (!length(moduleList)<60&&network6){
        graphlist<-list(gD2.cyt,gD3.cyt,gD4.cyt,gD6.cyt,gD7.cyt,gD8.cyt,gD9.cyt,gD10.cyt,gD11.cyt,gD12.cyt,gD13.cyt)
        windowlist<-list(gDCW2,gDCW3,gDCW4,gDCW6,gDCW7,gDCW8,gDCW9,gDCW10,gDCW11,gDCW12,gDCW13)}
      
      if (!length(moduleList)<60&&!network6){
        graphlist<-list(gD2.cyt,gD3.cyt,gD4.cyt,gD7.cyt,gD8.cyt,gD9.cyt,gD10.cyt,gD11.cyt,gD12.cyt,gD13.cyt)
        windowlist<-list(gDCW2,gDCW3,gDCW4,gDCW7,gDCW8,gDCW9,gDCW10,gDCW11,gDCW12,gDCW13)}
      
      ##
      
      print("merging networks")
      
      window.name<-"MERGE_1"
      gDCWM <- new.CytoscapeWindow(window.name, graph=gD.cyt,overwriteWindow = TRUE,host="192.168.65.2",rpcPort=9000)
      displayGraph (gDCWM)
      
      
      
      for(graphiter in 1:length(graphlist)){
        print(paste("merge",graphiter))
        if(numEdges(graphlist[[graphiter]])>1){
          addGraphToGraph(gDCWM, graphlist[[graphiter]])
          all.edges<-getAllEdges(gDCWM)
          #add undirected edge attribute to all edges
          setEdgeAttributesDirect(gDCWM,"interaction","char",all.edges,rep("undirected",length(all.edges)))
          new.edges<-getAllEdges(windowlist[[graphiter]])
          new.edges2<-rev(gsub("unspecified","NA",new.edges))
          
          new.weight<-getAllEdgeAttributes(windowlist[[graphiter]])$weight
          setEdgeAttributesDirect(gDCWM,"weight","floating",new.edges2,new.weight)#wrong order now fixed
          
          new.R_sq<-getAllEdgeAttributes(windowlist[[graphiter]])$R_sq
          setEdgeAttributesDirect(gDCWM,"R_sq","floating",new.edges2,new.R_sq)#wrong order now fixed
          
          new.qvalue<-getAllEdgeAttributes(windowlist[[graphiter]])$qvalue
          setEdgeAttributesDirect(gDCWM,"qvalue","floating",new.edges2,new.qvalue)#wrong order now fixed
          
          redraw(gDCWM)
          layoutNetwork(gDCWM)
          setEdgeAttributes(gDCWM,"weight")  
          control.points <- c (-1.0,-0.75, 0.0, 0.75,1.0)
          edge.colors <- rgb(t(col2rgb(c("black","blue","lightblue","white","pink","red","green"))),maxColorValue = 255)
          setEdgeColorRule(gDCWM, edge.attribute.name='weight', control.points,edge.colors, mode='interpolate')
          redraw (gDCWM)
          setLayoutProperties (gDCWM, hlp[6], list (edge_attribute = 'weight'))
          layoutNetwork(gDCWM)
        }else{print(paste("graph",graphiter,"has less than 2 edges and can't be added"))}
      }
      
      #resend diffMEs as merging loses this info in some cases
      new.nodes_WGCNA<-intModules
      new.diffME_WGCNA<-unlist(all.difflist.modules[[edge]])
      setNodeAttributesDirect(gDCWM,"diffME","floating",new.nodes_WGCNA,new.diffME_WGCNA)
      
      new.nodes_Baylor<-unlist(diffExBaylor.edges[[edge]])[,1]
      new.diffME_Baylor<-unlist(diffExBaylor.edges[[edge]])[,2]
      setNodeAttributesDirect(gDCWM,"diffEX","floating",new.nodes_Baylor,new.diffME_Baylor)
      
      #Colour modules by up or down
      control.points <- c (-1.0, 0.0, 1.0)
      node.colors <- rgb(t(col2rgb(c("black","blue","white","red","green"))),maxColorValue = 255)
      setNodeColorRule(gDCWM, node.attribute.name='diffME', control.points,node.colors, mode='interpolate')
      redraw (gDCWM)
      layoutNetwork(gDCWM)
      
      # all.edges<-getAllEdges(gDCWM)
      # setEdgeAttributesDirect(gDCWM,"interaction","char",all.edges,rep("undirected",length(all.edges)))
      # current.edges<-getAllEdges(gDCWM)
      # 
      # new.edges2<-rev(gsub("NA","unspecified",new.edges))
      # setEdgeAttributesDirect(gDCWM,"ID","char",all.edges,rep("undirected",length(all.edges)))
      
      #this may work, possibly remove. But nice to have... (requires mounting of volume, implemented in docker-compose file)
      #the issue here is that saveNetwork() writes to the filesystem on which cytoscape is running. In this case it is a different container! The solution is to mount the main ANIMA directory in that container as well, and automagically the file will appear on the host filesystem as required. Specifically, ANIMA is mounted in the root of the cytoscape container, which is where session files are written to.
      
      try(saveNetwork(gDCWM,file.path(gsub("/home/rstudio/","",dir.results,fixed=TRUE),"session.gml"),format='gml'))
      
      gDCWM.copy<-existing.CytoscapeWindow("MERGE_1",host="192.168.65.2",rpcPort=9000,copy.graph.from.cytoscape.to.R=TRUE)
      #gM.cyt<-getGraphFromCyWindow(gDCWM,"MERGE_1")
      
      edgeAttr<-getAllEdgeAttributes(gDCWM.copy)
      nodeAttr<-getAllNodeAttributes(gDCWM.copy)
      
      edgeAttr2<-edgeAttr[,2:7]
      nodeAttr2<-nodeAttr[,c(2,4,6,7,8)]
      
      write.csv(edgeAttr,file=file.path(dir.results,paste('table',q.i,'.',an.count,'.',i,'.',names(questions)[[q.i]],'.',analysis,names(datalist)[i],"edgeAttr_edge",edge,".csv",sep="")))
      write.csv(nodeAttr2,file=file.path(dir.results,paste('table',q.i,'.',an.count,'.',i,'.',names(questions)[[q.i]],'.',analysis,names(datalist)[i],"nodeAttr_edge",edge,".csv",sep="")))
      
      deleteAllWindows(cy)
      rm(list=ls()[grep("gDC.*",ls())])
      gc()
      
    }#end optional cytoscape
    
    #End edgewise Cytoscape####
  }#end edgewise cytoscape
  #end doNetworks if loop####
  #end doNetworks if loop####
  }#end doNetworks if loop####

#NETWORK 14: Map module probes to Neo4j####
if (neoinsert=="on"){
  nodes<-RNeo4j::nodes
  #network-specific cypher query
  
  write.csv(geneInfo,file.path(dir.figures,"tmp.csv"))
  if(class(datalist[[i]])=="ExpressionSet"){
    for (edge in 1:5){
      query = paste("LOAD CSV WITH HEADERS FROM 'file://",file.path(dir.figures,"tmp.csv"),"' AS csvLine 
                    MATCH (wmod:wgcna {name:csvLine.moduleColor,square:'",dataset.variables[[1]],"',edge:toInt(",edge,")})
                    MATCH (probe:PROBE {name:csvLine.substanceBXH,square:'",dataset.variables[[1]],"',edge:toInt(",edge,")})
                    CREATE (probe)-[:mapsTo]->(wmod)              
                    MERGE (gene:SYMBOL {name:csvLine.HGNC_symbol}) 
                    MERGE (probe)-[:mapsTo]->(gene)",sep="")
      cypher(graph,query)
      queryclean<-"MATCH (s:SYMBOL {name:'NA'}) OPTIONAL MATCH (s)-[r]-() DELETE s, r"
      cypher(graph,queryclean)
      
    }
  }else if(class(datalist[[i]])=="LumiBatch"){
    for (edge in 1:5){
      query = paste("LOAD CSV WITH HEADERS FROM 'file://",file.path(dir.figures,"tmp.csv"),"' AS csvLine 
                    MATCH (wmod:wgcna {name:csvLine.moduleColor,square:'",dataset.variables[[1]],"',edge:toInt(",edge,")})
                    MATCH (probe:PROBE {name:csvLine.substanceBXH,square:'",dataset.variables[[1]],"',edge:toInt(",edge,")})
                    CREATE (probe)-[:mapsTo]->(wmod)              
                    MERGE (gene:SYMBOL {name:csvLine.targetID}) 
                    MERGE (probe)-[:mapsTo]->(gene)",sep="")
      cypher(graph,query)
      queryclean<-"MATCH (s:SYMBOL {name:'NA'}) OPTIONAL MATCH (s)-[r]-() DELETE s, r"
      cypher(graph,queryclean)
      
    }  
  }
  
  
}#end neoinsert if

#old
#debug create (probe)- changed to merge
#it has now been changed back to create for a new run

# query = "
# MATCH (wmod:wgcna {name:{wgcnamod},square:"",edge:{edgeg}})
# MATCH (probe:PROBE {name:{nuID},square:{squareg},edge:{edgeg}})
# MERGE (gene:SYMBOL {name:{targetID}}) 
# CREATE (probe)-[:mapsTo]->(wmod)
# CREATE (probe)-[:mapsTo]->(gene)
# "
# blocks<-split(1:nrow(geneInfo),ceiling(seq_along(1:nrow(geneInfo))/1000))
# for (edge in 1:5){#FIXED!
#   for (block in blocks){
#     t = suppressMessages(newTransaction(graph))
#     graphdata=geneInfo[block,]
#     for (therow in 1:nrow(graphdata)) { #this can get extremely slow
#       probe = as.character(graphdata[therow, ]$substanceBXH)
#       wmod = as.character(graphdata[therow, ]$moduleColor)
#       if (is.na(graphdata[therow, ]$geneSymbol)){
#         gene = as.character(graphdata[therow, ]$targetID)#use secondary name derived from chip annotation rather than illumina.db lookup. this will be non HGNC in many cases
#       } else {gene = as.character(graphdata[therow, ]$geneSymbol)}#modify# this causes problems as it is not the same as targetID}
#       
#       
#       suppressMessages(appendCypher(t, 
#                    query, 
#                    wgcnamod = wmod,
#                    nuID = probe, 
#                    targetID = gene,
#                    #invariant node properties
#                    squareg = dataset.variables[[1]],
#                    edgeg = edge               
#       ))
#     }
#     
#     #print("committing :)")
#     print(paste("Edge:",edge,"Batch:", ceiling(therow/1000),"of",ceiling(nrow(graphdata)/1000), "committed."))
#     suppressMessages(commit(t))
#   }#end blocks
#   #END NEO4J
# }#end NEO4J edges


#NETWORK 15: NEO4J CODE TO WRITE INTERPROBE NETWORKS####
#debug: used merge
#back to create
if (neoinsert=="on"){
  nodes<-RNeo4j::nodes
  #network-specific cypher query
  for(cytmod in modNetList){
    graphdata_big<-cytmod$edgeData[order(-cytmod$edgeData$weight),]
    print(paste("graphdata_big has",nrow(graphdata_big),"rows"))
    graphdata_all<-graphdata_big[1:min(max(c(nrow(cytmod$nodeData)*3,ceiling(nrow(graphdata_big)/10))),2000),]#this is still experimental and requires optimisation
    print(nrow(cytmod$nodeData)*3)
    print(ceiling(nrow(graphdata_big)/10))
    print(paste("graphdata_all has",nrow(graphdata_all),"rows"))
    
    write.csv(graphdata_all,file.path(dir.figures,"tmp.csv"))
    for (edge in 1:5){
      query = paste("LOAD CSV WITH HEADERS FROM 'file://",file.path(dir.figures,"tmp.csv"),"' AS csvLine 
                    MATCH (probe1:PROBE {name:csvLine.fromNode,square:'",dataset.variables[[1]],"',edge:toInt(",edge,")})
                    MATCH (probe2:PROBE {name:csvLine.toNode,square:'",dataset.variables[[1]],"',edge:toInt(",edge,")})
                    CREATE (probe1)-[:TOMCOR {TOMweight:toFloat(csvLine.weight),square:'",dataset.variables[[1]],"',edge:toInt(",edge,")}]->(probe2)",sep="")
      cypher(graph,query)
    }#end NEO4J edges
    
    
    
    # query = "
    # MATCH (probe1:PROBE {name:{probe1g},square:{squareg},edge:{edgeg}})
    # MATCH (probe2:PROBE {name:{probe2g},square:{squareg},edge:{edgeg}})
    # CREATE (probe1)-[:TOMCOR {TOMweight:{TOMweightg},square:{squareg},edge:{edgeg}}]->(probe2)
    # "
    # t = suppressMessages(newTransaction(graph))
    # blocks<-split(1:nrow(graphdata_all),ceiling(seq_along(1:nrow(graphdata_all))/1000))
    # for (edge in 1:5){
    #   
    # 
    #   for (block in blocks){
    #     t = suppressMessages(newTransaction(graph))
    #     graphdata=graphdata_all[block,]
    #     for (therow in 1:nrow(graphdata)) { #this can get extremely slow
    #       probe1 = as.character(graphdata[therow, ]$fromNode)
    #       probe2 = as.character(graphdata[therow, ]$toNode)
    #       TOMweight = as.numeric(as.character(graphdata[therow,]$weight))
    #       suppressMessages(appendCypher(t, 
    #                                     query, 
    #                                     probe1g = probe1,
    #                                     probe2g = probe2,
    #                                     TOMweightg = TOMweight,
    #                                     #invariant node properties
    #                                     squareg = dataset.variables[[1]],
    #                                     edgeg = edge)
    #       
    #       )
    #     }
    #     
    #     #print("committing :)")
    #     print(paste("Edge:",edge,"Batch:", ceiling(therow/1000),"of",ceiling(nrow(graphdata)/1000), "committed."))
    #     suppressMessages(commit(t))
    #   }#end blocks
    #   #END NEO4J
    # }#end NEO4J edges
  }#end modNetList loop
  }#end neoinsert if

#9_Baylor and WGCNA module correspondence*********************************####
analysis="9_Baylor_WGCNA_Modules"
an.count=9
figure<-1
dir.results <- file.path(dir.output.version,paste("Q_",q.i,"_",dataset.variables,sep=""),analysis, "results")
## Define folder for saving figures
dir.figures <- file.path(dir.output.version,paste("Q_",q.i,"_",dataset.variables,sep=""),analysis, "figures")

for (dir in c(dir.figures, dir.results)) {
  if (!file.exists(dir)) {
    dir.create(dir, recursive=T,showWarnings=T)
  }
}
setwd(dir.figures)

sum(unique(unlist(modlist4))%in%unique(unlist(modgeneslist)))/length(unique(unlist(modlist4)))
sum(unique(unlist(modentrez4))%in%unique(unlist(modentrezlist)))/length(unique(unlist(modentrez4)))

capturematrix<-matrix(data=NA,nrow=length(modlist4),ncol=length(modgeneslist))
dimnames(capturematrix)<-list(names(modlist4),names(modgeneslist))
for (baylormodule in 1:length(modlist4)){
  for(wgcnamodule in 1:length(modgeneslist)){
    capturematrix[baylormodule,wgcnamodule]<-sum(modlist4[[baylormodule]]%in%modgeneslist[[wgcnamodule]])/length(modlist4[[baylormodule]])
  }
}
hmcolf<-colorRampPalette(c("white","black"))

pdf(paste('fig_',q.i,'.',an.count,'.',figure,'_Baylor_WGCNA.pdf',sep=""),height=pdf.options()$height*1,width=pdf.options()$width*1)
heatmap(capturematrix,scale="none",col=hmcolf(50))
dev.off()
figure<-figure+1

capturematrix2<-matrix(data=NA,nrow=length(modentrez4),ncol=length(modentrezlist))
dimnames(capturematrix2)<-list(names(modentrez4),names(modentrezlist))
for (baylormodule in 1:length(modentrez4)){
  theBaylorPresenceProportion<-sum(modentrez4[[baylormodule]]%in%unique(unlist(modentrezlist)))/length(modentrez4[[baylormodule]])
  print(paste(names(modentrez4)[[baylormodule]],theBaylorPresenceProportion))
  for(wgcnamodule in 1:length(modentrezlist)){
    
    thecmres<-sum(modentrez4[[baylormodule]]%in%modentrezlist[[wgcnamodule]])*theBaylorPresenceProportion
    capturematrix2[baylormodule,wgcnamodule]<-thecmres
  }
}
pdf(paste('fig_',q.i,'.',an.count,'.',figure,'_Baylor_WGCNA2.pdf',sep=""),height=pdf.options()$height*1,width=pdf.options()$width*1)
heatmap(capturematrix2,scale="none",col=hmcolf(50),main="BaylorPresenceNormalised")
dev.off()

##
#Jaccard index
capturematrix3<-matrix(data=NA,nrow=length(modentrez4),ncol=length(modentrezlist))
dimnames(capturematrix3)<-list(names(modentrez4),names(modentrezlist))
for (baylormodule in 1:length(modentrez4)){
  #normalize for gene presence in the Baylor module list
  theBaylorPresenceProportion<-sum(modentrez4[[baylormodule]]%in%unique(unlist(modentrezlist)))/length(modentrez4[[baylormodule]])
  print(paste(names(modentrez4)[[baylormodule]],theBaylorPresenceProportion))
  for(wgcnamodule in 1:length(modentrezlist)){
    jaccard.index<-sum(modentrez4[[baylormodule]]%in%modentrezlist[[wgcnamodule]])/length(union(modentrez4[[baylormodule]],modentrezlist[[wgcnamodule]]))/theBaylorPresenceProportion
    capturematrix3[baylormodule,wgcnamodule]<-jaccard.index
  }
}
pdf(paste('fig_',q.i,'.',an.count,'.',figure,'_Baylor_WGCNA3_Jaccard.pdf',sep=""),height=pdf.options()$height*1,width=pdf.options()$width*1)
heatmap(capturematrix3,scale="none",col=hmcolf(50),main="JaccardIndex_normalised")
dev.off()

##

#10_Baylor2013 and WGCNA module correspondence*********************************####
analysis="10_Baylor2013_WGCNA_Modules"
an.count=9
figure<-1
dir.results <- file.path(dir.output.version,paste("Q_",q.i,"_",dataset.variables,sep=""),analysis, "results")
## Define folder for saving figures
dir.figures <- file.path(dir.output.version,paste("Q_",q.i,"_",dataset.variables,sep=""),analysis, "figures")

for (dir in c(dir.figures, dir.results)) {
  if (!file.exists(dir)) {
    dir.create(dir, recursive=T,showWarnings=T)
  }
}
setwd(dir.figures)

sum(unique(unlist(modlist4.2))%in%unique(unlist(modgeneslist)))/length(unique(unlist(modlist4.2)))#gene level modlist 4.2 is (new) baylor modgeneslist is wgcna
sum(unique(unlist(modNU))%in%unique(unlist(modNUIDlist)))/length(unique(unlist(modNUIDlist)))#nuID level modNU is (new) Baylor modNUIDlist is wgcna

capturematrix<-matrix(data=NA,nrow=length(modlist4.2),ncol=length(modgeneslist))
dimnames(capturematrix)<-list(names(modlist4.2),names(modgeneslist))
for (baylormodule in 1:length(modlist4.2)){
  for(wgcnamodule in 1:length(modgeneslist)){
    capturematrix[baylormodule,wgcnamodule]<-sum(modlist4.2[[baylormodule]]%in%modgeneslist[[wgcnamodule]])/length(modlist4.2[[baylormodule]])
  }
}
hmcolf<-colorRampPalette(c("white","black"))
pdf(paste('fig_',q.i,'.',an.count,'.',figure,'_Baylor_2013_WGCNA.pdf',sep=""),height=pdf.options()$height*1,width=pdf.options()$width*1)
heatmap(capturematrix,scale="none",col=hmcolf(50),main="gene level")
dev.off()
figure<-figure+1

#enrichnorm=enrich*12000/lapply(moduleprobes,length)[[iter]]
#enrichnormlist[[iter]]<-enrichnorm

capturematrix2<-matrix(data=NA,nrow=length(modNU),ncol=length(modNUIDlist))
dimnames(capturematrix2)<-list(names(modNU),names(modNUIDlist))
for (baylormodule in 1:length(modNU)){
  theBaylorPresenceProportion<-sum(modNU[[baylormodule]]%in%unique(unlist(modNUIDlist)))/length(modNU[[baylormodule]])
  print(paste(names(modNU)[[baylormodule]],theBaylorPresenceProportion))
  for(wgcnamodule in 1:length(modNUIDlist)){
    
    thecmres<-sum(modNU[[baylormodule]]%in%modNUIDlist[[wgcnamodule]])*theBaylorPresenceProportion
    capturematrix2[baylormodule,wgcnamodule]<-thecmres
  }
}
pdf(paste('fig_',q.i,'.',an.count,'.',figure,'_Baylor_2013_WGCNA2.pdf',sep=""),height=pdf.options()$height*1,width=pdf.options()$width*1)
heatmap(capturematrix2,scale="none",col=hmcolf(50),main="BaylorPresenceNormalised")
dev.off()##
#Jaccard index
capturematrix3<-matrix(data=NA,nrow=length(modNU),ncol=length(modNUIDlist))
dimnames(capturematrix3)<-list(names(modNU),names(modNUIDlist))
for (baylormodule in 1:length(modNU)){
  #normalize for gene presence in the Baylor module list
  theBaylorPresenceProportion<-sum(modNU[[baylormodule]]%in%unique(unlist(modNUIDlist)))/length(modNU[[baylormodule]])
  print(paste(names(modNU)[[baylormodule]],theBaylorPresenceProportion))
  for(wgcnamodule in 1:length(modNUIDlist)){
    jaccard.index<-sum(modNU[[baylormodule]]%in%modNUIDlist[[wgcnamodule]])/length(union(modNU[[baylormodule]],modNUIDlist[[wgcnamodule]]))#/theBaylorPresenceProportion
    capturematrix3[baylormodule,wgcnamodule]<-jaccard.index
  }
}
pdf(paste('fig_',q.i,'.',an.count,'.',figure,'_Baylor_2013_WGCNA3.pdf',sep=""),height=pdf.options()$height*1,width=pdf.options()$width*1)
heatmap(capturematrix3,scale="none",col=hmcolf(50),main="JaccardIndex_not_normalised")
dev.off()
##

#Exit inner loop####
#}#close inner loop dataset (i)

#end of questions
#end of all analysis for the question, now comes housekeeping in preparation for the next question


#testobject
#testobject<-matrix(data=1,nrow=10000,ncol=10000)

#all variables up to now
allvar<-ls()

#only the variables created in the for loop
tmpvar<-allvar[!(allvar%in%mastervar2)]

#remove the loop-specific variables to reclaim memory and perform garbage collection


print("memory_used before purge inside loop")
print(mem_used())
print(system("free"))

rm(list=tmpvar)

collectGarbage()
gc()

print("memory_used after purge inside loop")
print(mem_used())
print(system("free"))
#mail
try(sendmail("armin.deffur@icloud.com", subject=paste("Notification about question",q.i),message=paste("Calculations for question:",q.i," finished!",sep=""), password="rmail"))
#Exit outer loop####