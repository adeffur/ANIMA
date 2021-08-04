library(neo2R)
library(ggpubr)
source(file.path("/home/rstudio/source_data/questions.R"),local=TRUE)
source("/home/rstudio/source_data/setlist.R",local=TRUE)
graph<-startGraph(graphstring,database="animadb",username="neo4j",password="anima")

# 1. select the square

square<-squarephenofun()

square<-"blood.PCF.defPC"
square<-"blood.PCF.probPC"
square<-"nonLTBI.TBPC"

# 2. get the groups 

query<-paste("MATCH (p:personPheno {square:'",square,"'})-[r]-(vst:varsubtype)-[r2]-(vt:vartype) RETURN DISTINCT vst.name AS vst, vt. name as vt",sep="")

groups<-cypher(graph,query)
vst.list<-unique(groups$vst)
vt.list<-unique(groups$vt)

#some reactive ui input depends on square

# 3. select the group (can be vt or vst)

#decide by vt or by vst

# 4. get the data for the group and square
#select fill variable
fillvar<-fillvarfun()
# fillvar<-"Categories"
# fillvar<-"class1"
# fillvar<-"class2"

bygroup<-bygroupfun()
# bygroup<-"vst"
# bygroup<-"vt"

if(bygroup=="vst"){
  #by vst
  group<-phenovstfun()
  
  # group<-"Cytokine"
  # group<-"neutrophil_prot"
  # group<-"MMP/inh"
  # group<-"NETs"
  # group<-"autoantibodies"
  # group<-"ECM"
  # group<-"TB_response"
  query<-paste("MATCH (p:personPheno {square:'",square,"'})-[r]-(vst:varsubtype {name:'",group,"'}) RETURN DISTINCT p.personName as person, p.class1 as class1, p.class2 as class2, p.name as pheno, p.value as value",sep="")
  
}else if(bygroup=="vt"){
  #by vt
  group<-phenovtfun()
  
  # group<-"Pathology"
  # group<-"Clinical"
  query<-paste("MATCH (p:personPheno {square:'",square,"'})-[r]-(vst:varsubtype)-[r2]-(vt:vartype {name:'",group,"'}) RETURN DISTINCT p.personName as person, p.class1 as class1, p.class2 as class2, p.name as pheno, p.value as value")
  
}



#query
res<-cypher(graph,query)

# 5. determine number of classes (2 or 4)
Categories<-interaction(res$class1,res$class2)
res$Categories<-Categories
res$numeric<-as.numeric(res$value)
res$character<-as.character(res$value)
data<-res

# 6. determine if categorical or numeric, if both split into 2 datasets

# 7. generate ggplot on grid; user may specify grid?; also a separate plot type for categorical and numeric, mixing two types on facet does not make sense

#need to identify from data:
# categorical vs numeric
# paired or unpaired
# number of classes that can be compared



complevs<-eval(parse(text=paste("levels(as.factor(data$",fillvar,"))",sep="")))
if (length(complevs)==4){
  my_comparisons <- list( c(complevs[1], complevs[2]),c(complevs[1], complevs[3]),c(complevs[2], complevs[4]),c(complevs[3], complevs[4]) )
}else{
  my_comparisons <- list( c(complevs[1], complevs[2]) )
}

plotfig<-ggplot(data, aes(x=eval(parse(text=fillvar)), y=numeric, fill=eval(parse(text=fillvar)))) + 
  geom_boxplot(aes(alpha=0.3 )) +
  scale_fill_manual(values=c("lightsteelblue", "mistyrose", "moccasin","lightgreen")) +
  geom_dotplot(mapping=aes(x=eval(parse(text=fillvar)), y=numeric, color=eval(parse(text=fillvar)), fill = eval(parse(text=fillvar))),binaxis='y', stackdir='center',stackratio=0.7,method='dotdensity',inherit.aes = FALSE) +
  scale_color_manual(values=c("darkblue", "darkred", "orange","darkgreen")) +
  #geom_jitter() +
  #stat_summary(fun.y=median, geom="point", shape=18, color="red") +
  theme_minimal() +
  theme(axis.title.y = element_text(size = rel(2.5), angle = 90, color = "slategrey",margin = margin(t = 0, r = 50, b = 50, l = 50))) +
  theme(axis.title.x = element_text(size = rel(2.5), color = "slategrey",margin = margin(t = 50, r = 50, b = 0, l = 50))) + 
  theme(axis.text = element_text(size = rel(0.8))) +
  theme(legend.text = element_text(size = rel(1.0))) +
  theme(legend.title = element_text(size = rel(1.0))) +
  theme(legend.text=element_text(size=rel(1.0))) +
  theme(legend.key.size = unit(1,"cm")) +
  labs(y = "phenofun2()") +
  ggtitle(paste("Cytokines","by class")) +
  theme(plot.title = element_text(size=rel(2.5),hjust = 0.5,color="maroon")) +
  theme(strip.text.x = element_text(size = 12, colour = "blue", angle = 0)) +
  stat_compare_means(comparisons = my_comparisons,method="wilcox.test",label.y.npc="center") + # Add pairwise comparisons p-value
  #stat_compare_means(label.y = 50) +     # Add global p-value
  scale_y_continuous(expand=expand_scale(mult=c(0,0.2)))
  
  

potfig<-plotfig + facet_wrap(vars(pheno),scales="free",ncol=3)

print(plotfig)

########fix igraphplotter

query.base<-"MATCH (c:CELL)-[r0]-(x1)-[r1]-(n:wgcna {square:'blood.PCF.probTB',edge:5,name:'black'})-[r2]-(x2) WHERE (x1:cellEx OR x1:cellprop) AND (x2:reactomePW OR x2:ImmunePW OR x2:pheno OR x2:PalWangPW)"
query.base<-"MATCH (c:CELL)-[r0]-(x1)-[r1]-(n:wgcna {square:'blood.PCF.probTB',edge:5,name:'blue'})-[r2]-(x2) WHERE (x1:cellEx OR x1:cellprop) AND (x2:reactomePW OR x2:ImmunePW OR x2:pheno OR x2:PalWangPW)"
nodelist<-c("c","x1","n","x2")
edgetrips<-list(c("c","r0","x1"),c("x1","r1","n"),c("n","r2","x2"))
igr<-igraph_plotter(query.base,nodelist,edgetrips,rimpar="diffME",plot=TRUE,csv=FALSE,prefix="/home",filename = "igraphWGCNAannot",plotd3=FALSE)



