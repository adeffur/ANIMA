# 0. Initialize
library(dplyr)
source("~/scripts/ANIMA_setup.R")
source(file.path(dir.data_root,"questions.R")) #this is specific to the project; each project will get its own set of questions tied to specific data.
source("~/source_data/setlist.R")
graph<-startGraph(graphstring)#using docker for Mac!

table1num<-data.frame("vartype"=character(),"varsubtype"=character(),"var"=character(),"All n"=integer(),"All median"=numeric(),"All IQR"=numeric(),
                       "grp1 n"=integer(),"grp1 median"=numeric(),"grp1 IQR"=numeric(),"grp2 n"=numeric(),"grp2 median"=numeric(),"grp2 IQR"=numeric(),"testname"=character(),"P-value"=numeric())

table1cat<-data.frame("vartype"=character(),"varsubtype"=character(),"var"=character(),"All n"=integer(),"grp1"=character(),"grp2"=character(),"testname"=character(),"P-value"=numeric())

count<-0
count2<-0
# 1. Get vartypes for square
query<-"MATCH (pp:personPheno {square:'blood.PCF.defPC'})-[r]-(vst:varsubtype)-[r2]-(vt:vartype) RETURN DISTINCT(vt.name)"
res<-cypher(graph,query)
vartypes<-res[,1]

# 2. For loop over vartypes
for (vt in vartypes){
# 3. Get varsubtypes for vartype

  query<-paste("MATCH (pp:personPheno {square:'blood.PCF.defPC'})-[r]-(vst:varsubtype)-[r2]-(vt:vartype) WHERE vt.name = '",vt,"' RETURN DISTINCT(vst.name)",sep="")
  res2<-cypher(graph,query)
  varsubtypes<-res2[,1]
  
# 4. For loop over varsubtypes
  for (vst in varsubtypes){
# 5. Get variables for varsubtype
    query<-paste("MATCH (pp:personPheno {square:'blood.PCF.defPC'})-[r]-(vst:varsubtype)-[r2]-(vt:vartype) WHERE vt.name = '",vt,"' AND vst.name = '",vst,"' RETURN DISTINCT(pp.name)",sep="")
    res3<-cypher(graph,query)
    vars<-res3[,1]
# 6. For loop over variables
    for (var in vars){
      #count<-count+1
# 7. Get data for the variable
      query<-paste("MATCH (pp:personPheno {square:'blood.PCF.defPC'})-[r]-(vst:varsubtype)-[r2]-(vt:vartype) WHERE vt.name = '",vt,"' AND vst.name = '",vst,"' AND pp.name = '",var,"' RETURN 
                   pp.personName AS sampleID, pp.class1 AS class1, pp.class2 AS class2, pp.value AS value",sep="")
      data<-cypher(graph,query)
     
# 8. Check if the class1 subsets are identical. This is the case for matched samples (blood/fluid where the variable isn't different between blood and fluid, but not otherwise) select only data for one class1 type
      c1<-sort(data$value[which(data$class1==unique(data$class1)[1])])
      c2<-sort(data$value[which(data$class1==unique(data$class1)[2])])
      needtosplit<="no"
      if (length(c1)==length(c2)){
        if (sum(c1==c2)==length(c1)){
          data<-data[which(data$class1==sort(unique(data$class1))[1]),]
          print("debug1")
          needtosplit<-"no"
        }else{
          data<-data
          print("debug2")
          needtosplit<-"yes"
        }
      }#end lengthcheck
      #changeindex<-which(c1!=c2)
# 9. Define variable class (numeric or character)
    datavec<-as.numeric(data$value)
    if (sum(is.na(datavec))==length(datavec)) {
      data$value<-as.character(data$value)
    }else{
      data$value<-as.numeric(data$value)
    }
# 10. if numeric
    if (class(data$value)=="numeric" & needtosplit == "no"){
      
      count<-count+1#remove NA
      data<-data[which(!is.na(data$value)),]
      
      data$class2<-as.factor(data$class2)
      alln<-nrow(data)
      allmed<-median(data$value)
      alliqr<-IQR(data$value)
      
      
      sum<-data %>%
        group_by(class2) %>% 
        summarise(n=n(),median = median(value),IQR = IQR(value))
      
      class2names<-as.character(sum$class2)
      
      statname="kruskal wallis"
      stats<-kruskal.test(value~class2,data=data)
      
      table1num[count,]<-c(vt,vst,var,alln,allmed,alliqr,sum[1,2],sum[1,3],sum[1,4],sum[2,2],sum[2,3],sum[2,4],stats$method,stats$p.value)
      
      
    }else if (class(data$value)=="numeric" & needtosplit == "yes"){
      #split data in 2 
      splitclasses<-sort(unique(data$class1))
      
      #PartA
      count<-count+1#remove NA
      data<-data[which(!is.na(data$value)),]
      dataS1<-data[which(data$class1==splitclasses[1]),]
      
      
      dataS1$class2<-as.factor(dataS1$class2)
      alln<-nrow(dataS1)
      allmed<-median(dataS1$value)
      alliqr<-IQR(dataS1$value)
      
      
      sum<-dataS1 %>%
        group_by(class2) %>% 
        summarise(n=n(),median = median(value),IQR = IQR(value))
      
      class2names<-as.character(sum$class2)
      
      statname="kruskal wallis"
      stats<-"fail"
      
      try(stats<-kruskal.test(value~class2,data=dataS1))
      if(stats!="fail"){
        table1num[count,]<-c(vt,vst,paste(var,splitclasses[1],sep="."),alln,allmed,alliqr,sum[1,2],sum[1,3],sum[1,4],sum[2,2],sum[2,3],sum[2,4],stats$method,stats$p.value)
      }else{
        table1num[count,]<-c(vt,vst,paste(var,splitclasses[1],sep="."),alln,allmed,alliqr,sum[1,2],sum[1,3],sum[1,4],sum[2,2],sum[2,3],sum[2,4],"NA",NA)
      }
      
      
      
      
      
      #PartB
      count<-count+1#remove NA
      data<-data[which(!is.na(data$value)),]
      dataS2<-data[which(data$class1==splitclasses[2]),]
      
      
      dataS2$class2<-as.factor(dataS2$class2)
      alln<-nrow(dataS2)
      allmed<-median(dataS2$value)
      alliqr<-IQR(dataS2$value)
      
      
      sum<-dataS2 %>%
        group_by(class2) %>% 
        summarise(n=n(),median = median(value),IQR = IQR(value))
      
      class2names<-as.character(sum$class2)
      
      statname="kruskal wallis"
      stats<-"fail"
      
      try(stats<-kruskal.test(value~class2,data=dataS2))
      if(stats!="fail"){
        table1num[count,]<-c(vt,vst,paste(var,splitclasses[2],sep="."),alln,allmed,alliqr,sum[1,2],sum[1,3],sum[1,4],sum[2,2],sum[2,3],sum[2,4],stats$method,stats$p.value)
      }else{
        table1num[count,]<-c(vt,vst,paste(var,splitclasses[2],sep="."),alln,allmed,alliqr,sum[1,2],sum[1,3],sum[1,4],sum[2,2],sum[2,3],sum[2,4],"NA",NA)
      }
      
      print("hey!")
    }else if (class(data$value)=="character" & needtosplit == "no"){
      
      #remove NA
      data<-data[which(data$value!=""),]
      data<-data[which(data$value!="n.a"),]
      
      if (nrow(data)>0){
        data$class2<-as.factor(data$class2)
        alln<-nrow(data)
        data$value<-as.factor(data$value)
        ft<-ftable(class2~value,data=data)
        if(nrow(ft)==2){
          ft[,1]<-ft[,1]/colSums(ft)[1]
          ft[,2]<-ft[,2]/colSums(ft)[2]
          statname="Chi-squared test"
          stats<-chisq.test(ft)
          
          ftdf<-as.data.frame((ft))
          class2names<-as.character(levels(ftdf$class2))
          ftdf1<-ftdf[which(ftdf$class2==class2names[1]),]
          ftdf2<-ftdf[which(ftdf$class2==class2names[2]),]
          
          grp1<-paste(paste(ftdf1$value,ftdf1$Freq,sep=":"),collapse="__")
          grp2<-paste(paste(ftdf2$value,ftdf2$Freq,sep=":"),collapse="__")
          
          count2<-count2+1
          table1cat[count2,]<-c(vt,vst,var,alln,grp1,grp2,stats$method,stats$p.value)
          }#endif
      
      }#endif
      
      
    }else if (class(data$value)=="character" & needtosplit == "yes"){
      print("need to split cat!")
      
    }else{NULL}

  
  
    }#end vars loop
  }#end varsubtypes loop  
}#end vartypes loop
colnames(table1num)<-gsub("grp1",class2names[1],colnames(table1num))
colnames(table1num)<-gsub("grp2",class2names[2],colnames(table1num))
colnames(table1cat)<-gsub("grp1",class2names[1],colnames(table1cat))
colnames(table1cat)<-gsub("grp2",class2names[2],colnames(table1cat))
