install.packages("table1")
library(table1)
#library(dplyr)
# 0. Initialize



source("~/scripts/ANIMA_setup.R")
source(file.path(dir.data_root,"questions.R")) #this is specific to the project; each project will get its own set of questions tied to specific data.
source("~/source_data/setlist.R")
graph<-startGraph(graphstring,database="animadb",username="neo4j",password="anima")

#function
pvalue <- function(x, ...) {
  # Construct vectors of data y, and groups (strata) g
  y <- unlist(x)
  g <- factor(rep(1:length(x), times=sapply(x, length)))
  if (is.numeric(y)) {
    # For numeric variables, perform a non-parametric test
    p <- kruskal.test(y ~ g)$p.value
  } else {
    # For categorical variables, perform a chi-squared test of independence
    p <- chisq.test(table(y, g))$p.value
  }
  # Format the p-value, using an HTML entity for the less-than sign.
  # The initial empty string places the output on the line below the variable label.
  #c("", sub("<", "&lt;", format.pval(p, digits=3, eps=0.001)))
  c("", sub("<", "<", format.pval(p, digits=3, eps=0.001)))
}

html_to_pdf <- function(html_file, pdf_file) {
  cmd <- sprintf("pandoc --pdf-engine=pdflatex %s -t latex -o %s", html_file, pdf_file)
  "pandoc reports/7/report.html -t latex -o reports/7/report.pdf"
  system(cmd)
}

#Clinical characteristics####
#get data in dataframe

#blood.pcf.defpc####
square<-"blood.PCF.defPC"
query<-paste("MATCH (pr:person)-[r0]-(p:personPheno {square:'",square,"'})-[r]-(v1:varsubtype)-[r2]-(v2:vartype) WHERE (v2.name IN ['Clinical','Pathology','Imaging','Demographics'] AND p.personName =~ '.*?blood') RETURN p.personName AS name,p.name AS var,p.value AS value, v1.name AS vst, v2.name AS VT",sep="")
query<-paste("MATCH (pr:person)-[r0]-(p:personPheno {square:'",square,"'})-[r]-(v1:varsubtype)-[r2]-(v2:vartype) WHERE (v2.name IN ['Clinical','Pathology','Imaging','Demographics'] AND p.personName =~ '.*?blood') RETURN p.personName AS name,pr.class2 as HIV,p.name AS var,p.value AS value, v1.name AS vst, v2.name AS VT",sep="")
res<-cypher(graph,query)
df.wide <- pivot_wider(res[,1:4], names_from = var, values_from = value) 
print(colnames(df.wide))
df.wide[, c(5,8,9,10,11,12,13:16)] <- sapply(df.wide[, c(5,8,9,10,11,12,13:16)], as.numeric)
#labels<-c("HIV","Ethnicity","Age","Sex","Weight loss","Night sweats","Temperature","Heart rate","Systolic BP","Diastolic BP","Haemoglobin","WBC count","Platelet count","CD4 count")

df.wide$ETHNICITY<-factor(df.wide$ETHNICITY,labels=c("Black","Coloured","Indian"))
df.wide$SEX<-factor(df.wide$SEX,labels=c("Female","Male"))
df.wide$WtLoss<-factor(df.wide$WtLoss,levels=c(0,1),labels=c("No","Yes"))
df.wide$NtSweat<-factor(df.wide$NtSweat,levels=c(0,1),labels=c("No","Yes"))
df.wide$TB.class<-factor(df.wide$TB.class,levels=c("definite","probable"),labels<-c("Definite TB-PC","Probable TB-PC"))
df.wide$HIV<-factor(df.wide$HIV,levels=c("HIVneg","HIVpos"),labels<-c("HIV uninfected","HIV infected"))

table1::label(df.wide$ETHNICITY)<-"Ethnicity"
table1::label(df.wide$SEX)<-"Sex"
table1::label(df.wide$WtLoss)<-"Weight loss"
table1::label(df.wide$NtSweat)<-"Night sweats"
table1::label(df.wide$TB.class)<-"TB class"
table1::label(df.wide$HIV)<-"HIV status"

table1::label(df.wide$AGE_CALC)<-"Age"
table1::label(df.wide$Temp)<-"Body temperature"
table1::label(df.wide$Pulse)<-"Pulse rate"
table1::label(df.wide$BPsystolic)<-"Systolic blood pressure"
table1::label(df.wide$BPdiastolic)<-"Diastolic blood pressure"
table1::label(df.wide$Hb)<-"Haemoglobin concentration"
table1::label(df.wide$WCC)<-"White blood cell count"
table1::label(df.wide$PLT)<-"Platelet count"
table1::label(df.wide$CD4)<-"CD4 count"


table1(~ ETHNICITY + AGE_CALC + SEX + WtLoss + NtSweat + Temp + Pulse + BPsystolic + BPdiastolic + Hb + WCC + PLT +CD4 | TB.class + HIV , data=df.wide, overall=F, extra.col=list(`P-value`=pvalue),topclass="Rtable1-zebra")
table1(~ ETHNICITY + AGE_CALC + SEX + WtLoss + NtSweat + Temp + Pulse | TB.class + HIV , data=df.wide, overall=F, extra.col=list(`P-value`=pvalue),topclass="Rtable1-zebra")
table1(~ BPsystolic + BPdiastolic + Hb + WCC + PLT +CD4 | TB.class + HIV , data=df.wide, overall=F, extra.col=list(`P-value`=pvalue),topclass="Rtable1-zebra")
table1(~ ETHNICITY + SEX + AGE_CALC + CD4 | TB.class + HIV , data=df.wide, overall=F, extra.col=list(`P-value`=pvalue),topclass="Rtable1-zebra")

#html_to_pdf("/home/rstudio/output/table1.html","/home/rstudio/output/table1.pdf")


#blood.pcf.probpc####
square<-"blood.PCF.probPC"
query<-paste("MATCH (pr:person)-[r0]-(p:personPheno {square:'",square,"'})-[r]-(v1:varsubtype)-[r2]-(v2:vartype) WHERE (v2.name IN ['Clinical','Pathology','Imaging','Demographics'] AND p.personName =~ '.*?blood') RETURN p.personName AS name,p.name AS var,p.value AS value, v1.name AS vst, v2.name AS VT",sep="")
query<-paste("MATCH (pr:person)-[r0]-(p:personPheno {square:'",square,"'})-[r]-(v1:varsubtype)-[r2]-(v2:vartype) WHERE (v2.name IN ['Clinical','Pathology','Imaging','Demographics'] AND p.personName =~ '.*?blood') RETURN p.personName AS name,pr.class2 as HIV,p.name AS var,p.value AS value, v1.name AS vst, v2.name AS VT",sep="")
res<-cypher(graph,query)
df.wide <- pivot_wider(res[,1:4], names_from = var, values_from = value) 
print(colnames(df.wide))

df.wide$ETHNICITY<-factor(df.wide$ETHNICITY,levels=c("black","coloured"),labels=c("Black","Coloured"))
df.wide$SEX<-factor(df.wide$SEX,labels=c("Female","Male"))
df.wide$WtLoss<-factor(df.wide$WtLoss,levels=c(0,1),labels=c("No","Yes"))
df.wide$NtSweat<-factor(df.wide$NtSweat,levels=c(0,1),labels=c("No","Yes"))
df.wide$TB.class<-factor(df.wide$TB.class,levels=c("definite","probable"),labels<-c("Definite TB-PC","Probable TB-PC"))
df.wide$HIV<-factor(df.wide$HIV,levels=c("HIVneg","HIVpos"),labels<-c("HIV uninfected","HIV infected"))

df.wide[, c(5,8,9,10,11,12,13:16)] <- sapply(df.wide[, c(5,8,9,10,11,12,13:16)], as.numeric)

table1::label(df.wide$ETHNICITY)<-"Ethnicity"
table1::label(df.wide$SEX)<-"Sex"
table1::label(df.wide$WtLoss)<-"Weight loss"
table1::label(df.wide$NtSweat)<-"Night sweats"
table1::label(df.wide$TB.class)<-"TB class"
table1::label(df.wide$HIV)<-"HIV status"

table1::label(df.wide$AGE_CALC)<-"Age"
table1::label(df.wide$Temp)<-"Body temperature"
table1::label(df.wide$Pulse)<-"Pulse rate"
table1::label(df.wide$BPsystolic)<-"Systolic blood pressure"
table1::label(df.wide$BPdiastolic)<-"Diastolic blood pressure"
table1::label(df.wide$Hb)<-"Haemoglobin concentration"
table1::label(df.wide$WCC)<-"White blood cell count"
table1::label(df.wide$PLT)<-"Platelet count"
table1::label(df.wide$CD4)<-"CD4 count"



#labels<-c("HIV","Ethnicity","Age","Sex","Weight loss","Night sweats","Temperature","Heart rate","Systolic BP","Diastolic BP","Haemoglobin","WBC count","Platelet count","CD4 count")
table1(~ ETHNICITY + AGE_CALC + SEX + WtLoss + NtSweat + Temp + Pulse + BPsystolic + BPdiastolic + Hb + WCC + PLT +CD4 | TB.class + HIV , data=df.wide, overall=F, extra.col=list(`P-value`=pvalue),topclass="Rtable1-zebra")
#html_to_pdf("/home/rstudio/output/table1.html","/home/rstudio/output/table1.pdf")



#prob.def.blood####
square<-"prob.def.blood"
query<-paste("MATCH (pr:person)-[r0]-(p:personPheno {square:'",square,"'})-[r]-(v1:varsubtype)-[r2]-(v2:vartype) WHERE (v2.name IN ['Clinical','Pathology','Imaging','Demographics'] AND p.personName =~ '.*?blood') RETURN p.personName AS name,pr.class2 as HIV,p.name AS var,p.value AS value, v1.name AS vst, v2.name AS VT",sep="")
res<-cypher(graph,query)
df.wide <- pivot_wider(res[,1:4], names_from = var, values_from = value) 
print(colnames(df.wide))
df.wide[, c(5,8,9,10,11,12,13:16)] <- sapply(df.wide[, c(5,8,9,10,11,12,13:16)], as.numeric)

df.wide$ETHNICITY<-factor(df.wide$ETHNICITY,levels=c("black","coloured","indian"),labels=c("Black","Coloured","Indian"))
df.wide$SEX<-factor(df.wide$SEX,labels=c("Female","Male"))
df.wide$WtLoss<-factor(df.wide$WtLoss,levels=c(0,1),labels=c("No","Yes"))
df.wide$NtSweat<-factor(df.wide$NtSweat,levels=c(0,1),labels=c("No","Yes"))
df.wide$TB.class<-factor(df.wide$TB.class,levels=c("definite","probable"),labels<-c("Definite TB-PC","Probable TB-PC"))
df.wide$HIV<-factor(df.wide$HIV,levels=c("HIVneg","HIVpos"),labels<-c("HIV uninfected","HIV infected"))

table1::label(df.wide$ETHNICITY)<-"Ethnicity"
table1::label(df.wide$SEX)<-"Sex"
table1::label(df.wide$WtLoss)<-"Weight loss"
table1::label(df.wide$NtSweat)<-"Night sweats"
table1::label(df.wide$TB.class)<-"TB class"
table1::label(df.wide$HIV)<-"HIV status"

table1::label(df.wide$AGE_CALC)<-"Age"
table1::label(df.wide$Temp)<-"Body temperature"
table1::label(df.wide$Pulse)<-"Pulse rate"
table1::label(df.wide$BPsystolic)<-"Systolic blood pressure"
table1::label(df.wide$BPdiastolic)<-"Diastolic blood pressure"
table1::label(df.wide$Hb)<-"Haemoglobin concentration"
table1::label(df.wide$WCC)<-"White blood cell count"
table1::label(df.wide$PLT)<-"Platelet count"
table1::label(df.wide$CD4)<-"CD4 count"


table1(~ ETHNICITY + AGE_CALC + SEX + WtLoss + NtSweat + Temp + Pulse + BPsystolic + BPdiastolic + Hb + WCC + PLT +CD4 | TB.class + HIV , data=df.wide, overall=F, extra.col=list(`P-value`=pvalue),topclass="Rtable1-zebra")
#html_to_pdf("/home/rstudio/output/table1.html","/home/rstudio/output/table1.pdf")



#nonLTBI.PTB####
square<-"nonLTBI.PTB"
query<-paste("MATCH (pr:person)-[r0]-(p:personPheno {square:'",square,"'})-[r]-(v1:varsubtype)-[r2]-(v2:vartype) WHERE (v2.name IN ['Clinical','Pathology','Imaging','Demographics']) RETURN p.personName AS name,p.name AS var,p.value AS value, v1.name AS vst, v2.name AS VT",sep="")
query<-paste("MATCH (pr:person)-[r0]-(p:personPheno {square:'",square,"'})-[r]-(v1:varsubtype)-[r2]-(v2:vartype) WHERE (v2.name IN ['Clinical','Pathology','Imaging','Demographics']) RETURN p.personName AS name,pr.class2 as HIV,p.name AS var,p.value AS value, v1.name AS vst, v2.name AS VT",sep="")
res<-cypher(graph,query)
df.wide <- pivot_wider(res[,1:4], names_from = var, values_from = value) 
print(colnames(df.wide))

df.wide$ETHNICITY<-factor(df.wide$ETHNICITY,levels=c("black"),labels=c("Black"))
df.wide$SEX<-factor(df.wide$SEX,labels=c("Female","Male"))
df.wide$WtLoss<-factor(df.wide$WtLoss,levels=c(0,1),labels=c("No","Yes"))
df.wide$NtSweat<-factor(df.wide$NtSweat,levels=c(0,1),labels=c("No","Yes"))
df.wide$TB.class<-factor(df.wide$TB.class,levels=c("definite","probable"),labels<-c("Definite TB-PC","Probable TB-PC"))
df.wide$HIV<-factor(df.wide$HIV,levels=c("HIVneg","HIVpos"),labels<-c("HIV uninfected","HIV infected"))

table1::label(df.wide$ETHNICITY)<-"Ethnicity"
table1::label(df.wide$SEX)<-"Sex"
table1::label(df.wide$WtLoss)<-"Weight loss"
table1::label(df.wide$NtSweat)<-"Night sweats"
table1::label(df.wide$TB.class)<-"TB class"
table1::label(df.wide$HIV)<-"HIV status"

table1::label(df.wide$AGE_CALC)<-"Age"
table1::label(df.wide$Temp)<-"Body temperature"
table1::label(df.wide$Pulse)<-"Pulse rate"
table1::label(df.wide$BPsystolic)<-"Systolic blood pressure"
table1::label(df.wide$BPdiastolic)<-"Diastolic blood pressure"
table1::label(df.wide$Hb)<-"Haemoglobin concentration"
table1::label(df.wide$WCC)<-"White blood cell count"
table1::label(df.wide$PLT)<-"Platelet count"
table1::label(df.wide$CD4)<-"CD4 count"

df.wide[, c(5,8,9,10,11,12,13:16)] <- sapply(df.wide[, c(5,8,9,10,11,12,13:16)], as.numeric)
#labels<-c("HIV","Ethnicity","Age","Sex","Weight loss","Night sweats","Temperature","Heart rate","Systolic BP","Diastolic BP","Haemoglobin","WBC count","Platelet count","CD4 count")
table1(~ ETHNICITY + AGE_CALC + SEX + WtLoss + NtSweat + Temp + Pulse + BPsystolic + BPdiastolic + Hb + WCC + PLT +CD4 | HIV , data=df.wide, overall=F,topclass="Rtable1-zebra")
#html_to_pdf("/home/rstudio/output/table1.html","/home/rstudio/output/table1.pdf")



#Cytokines
square<-"blood.PCF.defPC"
query<-paste("MATCH (pr:person)-[r0]-(p:personPheno {square:'",square,"'})-[r]-(v1:varsubtype) WHERE (v1.name IN ['Cytokine']) RETURN p.personName AS name,pr.class2 as HIV,p.name AS var,p.value AS value",sep="")
res<-cypher(graph,query)
df.wide <- pivot_wider(res[,1:4], names_from = var, values_from = value)
df.wide<-df.wide %>%
  separate(name, c("person", "compartment"), " ")


print(colnames(df.wide))
df.wide[, c(4:18)] <- sapply(df.wide[, c(4:18)], as.numeric)
#labels<-c("HIV","Ethnicity","Age","Sex","Weight loss","Night sweats","Temperature","Heart rate","Systolic BP","Diastolic BP","Haemoglobin","WBC count","Platelet count","CD4 count")
table1(~ IL.1.beta + IL.6 + IL.18 + TNF.a + IFN.gamma.y + IP.10 + IL.13 + IL.17  + IL.22 + IL.23.A + TGF.beta.1 + IL.10 + IL.8 + IL.12p70 | HIV + compartment , data=df.wide, overall=F, extra.col=list(`P-value`=pvalue),topclass="Rtable1-zebra")

