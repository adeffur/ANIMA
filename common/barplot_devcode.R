# 1. select the square
print("phenogrid start::###########")

square<-squarephenofun()
#square<-"blood.PCF.defPC"
print("square:")
print(square)

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
  query<-paste("MATCH (ps:person)-[r0]-(p:personPheno {square:'",square,"'})-[r]-(vst:varsubtype {name:'",group,"'}) RETURN DISTINCT p.personName as person, ps.class1 as class1, ps.class2 as class2, p.name as pheno, p.value as value",sep="")
  
}else if(bygroup=="vt"){
  #by vt
  group<-phenovtfun()
  
  # group<-"Pathology"
  # group<-"Clinical"
  query<-paste("MATCH (ps:person)-[r0]-(p:personPheno {square:'",square,"'})-[r]-(vst:varsubtype)-[r2]-(vt:vartype {name:'",group,"'}) RETURN DISTINCT p.personName as person, ps.class1 as class1, ps.class2 as class2, p.name as pheno, p.value as value",sep="")
  
}



#query
res<-cypher(graph,query)
res$value<-as.numeric(res$value)
Categories<-interaction(res$class1,res$class2)
res$Categories<-Categories

df.wide <- pivot_wider(res, names_from = pheno, values_from = value) 
# 5. determine number of classes (2 or 4)
# Categories<-interaction(res$class1,res$class2)
# res$Categories<-Categories
# res$numeric<-as.numeric(res$value)
# res$character<-as.character(res$value)
# data<-res
# df.wide <- pivot_wider(data[,1:4], names_from = var, values_from = value) 
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

myData<-aggregate(df.wide$IL.8,by=list(cats=df.wide$Categories),FUN=function(x){c(median=median(x,na.rm=T),iqr=IQR(x,na.rm=T),n=length(x))})
myData<-do.call(data.frame,myData)
colnames(myData)<-c("Category","median","iqr","n")

dodge <- position_dodge(width = 0.9)
limits <- aes(ymax = median + iqr,
              ymin = median - iqr)
p <- ggplot(data = myData, aes(x = Category, y = median, fill = Category))
p + geom_bar(stat = "identity", position = dodge) +
  geom_errorbar(limits, position = dodge, width = 0.25) +
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(),
        axis.title.x=element_blank())
# p <- ggplot(data = myData, aes(x = names, y = mean, fill = names))
# p + geom_bar(stat = "identity", position = dodge) +
#    geom_errorbar(limits, position = dodge, width = 0.25) +

  #plotfig<-ggplot(data, aes(x=eval(parse(text=fillvar)), y=numeric, fill=eval(parse(text=fillvar)))) + 
  plotfig<-ggplot(myData, aes(x=eval(parse(text=fillvar)), y=median, fill=eval(parse(text=fillvar)))) + 
  #geom_boxplot(aes(alpha=0.3 ),width=.25) +
  geom_bar(stat = "identity") +
  #geom_col(stat = "median") +
  #geom_errorbar(limits,position=position_dodge(0.9),width = 0.25) +
  #geom_boxplot(width=.25) +
  #scale_fill_manual(values=c("lightsteelblue", "mistyrose", "moccasin","lightgreen")) +
  #geom_dotplot(mapping=aes(x=eval(parse(text=fillvar)), y=numeric, color=eval(parse(text=fillvar)), fill = eval(parse(text=fillvar))),binaxis='y', stackdir='center',stackratio=0.7,method='dotdensity',inherit.aes = FALSE) +
  #scale_color_manual(values=c("darkblue", "darkred", "orange","darkgreen")) +
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
  labs(y = "value") +
  labs(x = "groupings") +
  ggtitle(paste(group,"values by",fillvar)) +
  theme(plot.title = element_text(size=rel(2.5),hjust = 0.5,color="maroon")) +
  theme(strip.text.x = element_text(size = 12, colour = "blue", angle = 0)) +
  stat_compare_means(comparisons = my_comparisons,method="wilcox.test",label.y.npc="center") + # Add pairwise comparisons p-value
  #stat_compare_means(label.y = 50) +     # Add global p-value
  scale_y_continuous(expand=expand_scale(mult=c(0,0.2)))



plotfig<-plotfig + facet_wrap(vars(pheno),scales="free",ncol=3)

print(plotfig)