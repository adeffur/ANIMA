#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#
#install packages until docker image updated:
 #install.packages("shinycssloaders")
 # install.packages("ggpubr")
 # devtools::install_github("ricardo-bion/ggradar",dependencies=TRUE)
 # devtools::install_github("rstudio/gt",dependencies=TRUE)
 #install.packages("shinythemes")
 #install.packages("neo2R")
#remotes::install_github("rstudio/shinymeta")
#install.packages("ggnewscale")

#install.packages("huxtable")

#install.packages("reactlog")
#library(reactlog)
#reactlog_enable()
#File paths
tabledir="/home/rstudio/output/project/tabular"
cytodir="/home/rstudio/output/project/cytoscape"
venndir="/home/rstudio/output/project/vennDiagrams"
figdir="/home/rstudio/output/project/linePlots"
igraphdir="/home/rstudio/output/project/igraphs"
datadir="/home/rstudio/output/RData"
#cellcordir="/home/rstudio/output/project/cellcor"
modcordir="/home/rstudio/output/project/modcor"

options(stringsAsFactors = FALSE)

for (dir in c(tabledir, cytodir,venndir,figdir,igraphdir,datadir)) {
  if (!file.exists(dir)) {
    dir.create(dir, recursive=T,showWarnings=T)
  }
} 
#library(shinythemes)
library(shiny)
library(shinymeta)
library(ggplot2)
library(scales)
library(networkD3)
library(htmlwidgets)
library(igraph)
library(WGCNA)
#library(raster)

library(plyr)

library(xtable)

library(DT)
library(data.table)
library(pheatmap)
library(ReactomePA)
library(illuminaHumanv3.db)
library(illuminaHumanv4.db)

library("lumi", lib.loc="/usr/local/lib/R/site-library")
library(Heatplus)
# library(wordcloud)
# library(tm)
# library(grImport)

#install/packages("shinycssloaders")
library(shinycssloaders)
library(ggrepel)
library(ggpubr)

library(ggradar)
library(gridExtra)
library(cowplot)
library(dplyr)
library(tidyr)
library(huxtable)

library(table1)
library(kableExtra)

#library(textreadr)
scripts<-system(paste("ls","/home/rstudio/scripts"),intern=TRUE)
datafiles<-system(paste("ls",datadir),intern=TRUE)
excl<-grep("86|393",datafiles)
datafiles<-datafiles[!1:length(datafiles)%in%excl]
source(file.path("/home/rstudio/source_data/questions.R"),local=TRUE)
source("/home/rstudio/source_data/setlist.R",local=TRUE)
#library(RNeo4j)
library(neo2R)
#graph<-startGraph(graphstring)#!!!!! note the address!!

print("HERE WE ARE")
graph<-startGraph(graphstring,database="animadb",username="neo4j",password="anima")
print("connected!!!")
shadowtext <- function(x, y=NULL, labels, col='white', bg='black', 
                       theta= seq(0, 2*pi, length.out=50), r=0.07, ... ) {
  
  xy <- xy.coords(x,y)
  xo <- r*strwidth('A')
  yo <- r*strheight('A')
  
  # draw background text with small shift in x and y in background colour
  for (i in theta) {
    text( xy$x + cos(i)*xo, xy$y + sin(i)*yo, labels, col=bg, ... )
  }
  # draw actual text in exact xy position in foreground colour
  text(xy$x, xy$y, labels, col=col, ... )
}

doIntracor<-TRUE

#File paths2
#Projects root (scripts and output)
dir.root <- '/home/rstudio'
dir.common <- file.path(dir.root,"common")
#Output folders
dir.home<-file.path(dir.root,"ANIMA")
dir.main <- file.path(dir.root, 'output')

## Define folder for saving RData files
dir.rdata <- file.path(dir.main,'RData')
#print(paste("Data will be saved to", dir.rdata))

## Define folder for re-using KEGG png and xml files
dir.kegg <- file.path(dir.main,'KEGG')
#print(paste("KEGG Data will be saved to", dir.kegg))

#setwd("dockerhost")
for (dir in c(dir.main, dir.rdata,dir.kegg)) {
  if (!file.exists(dir)) {
    dir.create(dir, recursive=T,showWarnings=T)
  }
}

#Functions

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


#R-scripts root
dir.R.files <- file.path(dir.common, "03_Functions")
## R utilities folder
dir.util <- file.path(dir.R.files, 'util')

#Data root
dir.data_root<-file.path(dir.root,'source_data')

#annotation and module data
dir.annot<-file.path(dir.root,"source_data")

#Individual subject distributions based on 2 modules

#buildlisting<-system(paste("ls","/home/rstudio/output/build/"),intern=TRUE)

#source(file.path(dir.R.files,"igraph_plotter_newer.R"),local=TRUE)
source(file.path(dir.R.files,"igraph_plotter_4.R"),local=TRUE)
source(file.path(dir.R.files,"moduleMeta.R"),local=TRUE)
source(file.path("/home/rstudio/source_data/questions.R"),local=TRUE)
source(file.path(dir.R.files,"probe_boxplot4.R"),local=TRUE)
source(file.path(dir.R.files,"probe_boxplot5.R"),local=TRUE)
source(file.path(dir.R.files,"probe_boxplot6.R"),local=TRUE)
source(file.path(dir.R.files,"mwat.R"),local=TRUE)

##Required for G1

#new
source("/home/rstudio/source_data/setlist.R",local=TRUE)
#end new

edgelist<-1:5

query<-"MATCH (c:CELL) RETURN DISTINCT c.name AS cellname"
cellnamelist<-cypher(graph,query)
cellnames<-cellnamelist$cellname

##End required for G1

#Required for boxplot
#dataBXP<-cypher(graph,"MATCH (n:wgcna) RETURN DISTINCT n.square AS sets")
#datafilesBXP<-gsub("train\\.LTBI","train_LTBI",gsub("berry\\.","berry_",paste(datadir,"/",dataBXP$sets,".RData",sep="")))
datafilesBXP<-paste(datadir,"/",unlist2(setlistnameslist),".RData",sep="")
print("begin loading datasets")
for (file in datafilesBXP){
  load(file)
  print(paste("loaded",file))
  }
print("finished loading datasets")



#End required for boxplot

###Required for MM
data<-cypher(graph,"MATCH (n:wgcna) RETURN DISTINCT n.square AS sets")
#all<-1:258
all<-1:260
#all<-1:900
print("debug start")
print(data$sets)
neutromods<-cypher(graph,"MATCH (b:baylor {edge:'5'})-[r]-(c:cellEx)-[r2]-(c2:CELL) WHERE c2.name = 'neutrophil' RETURN DISTINCT b.name")$b.name
monomods<-cypher(graph,"MATCH (b:baylor {edge:'5'})-[r]-(c:cellEx)-[r2]-(c2:CELL) WHERE c2.name = 'monocyte' RETURN DISTINCT b.name")$b.name
bcellmods<-cypher(graph,"MATCH (b:baylor {edge:'5'})-[r]-(c:cellEx)-[r2]-(c2:CELL) WHERE c2.name = 'B-cell' RETURN DISTINCT b.name")$b.name
cd4mods<-cypher(graph,"MATCH (b:baylor {edge:'5'})-[r]-(c:cellEx)-[r2]-(c2:CELL) WHERE c2.name = 'CD4+ T-cell' RETURN DISTINCT b.name")$b.name
#cd8mods<-cypher(graph,"MATCH (b:baylor {edge:'5'})-[r]-(c:cellEx)-[r2]-(c2:CELL) WHERE c2.name = 'CD8+ T-cell' RETURN DISTINCT b.name")$b.name #(only one result!)
pltmods<-cypher(graph,"MATCH (b:baylor {edge:'5'})-[r]-(c:cellEx)-[r2]-(c2:CELL) WHERE c2.name = 'platelet' RETURN DISTINCT b.name")$b.name
nkmods<-cypher(graph,"MATCH (b:baylor {edge:'5'})-[r]-(c:cellEx)-[r2]-(c2:CELL) WHERE c2.name = 'NK-cell' RETURN DISTINCT b.name")$b.name
lymmods<-cypher(graph,"MATCH (b:baylor {edge:'5'})-[r]-(c:cellEx)-[r2]-(c2:CELL) WHERE c2.name = 'lymphocyte' RETURN DISTINCT b.name")$b.name
rbcmods<-cypher(graph,"MATCH (b:baylor {edge:'5'})-[r]-(c:cellEx)-[r2]-(c2:CELL) WHERE c2.name = 'RBC' RETURN DISTINCT b.name")$b.name

#functions and pathways
ifnmods<-cypher(graph,"MATCH (b:baylor {edge:'5'})-[r]-(pw) WHERE pw.name =~ '(?i).*interferon.*|.*IFN.*|.*ISG15.*' RETURN DISTINCT b.name")$b.name
lysosomemods<-cypher(graph,"MATCH (b:baylor {edge:'5'})-[r]-(pw) WHERE (pw:reactomePW OR pw:ImmunePW OR pw:PalWangPW) AND pw.name =~ '(?i).*lysosom[a,e].*' RETURN DISTINCT b.name")$b.name
phagocytmods<-cypher(graph,"MATCH (b:baylor {edge:'5'})-[r]-(pw) WHERE (pw:reactomePW OR pw:ImmunePW OR pw:PalWangPW) AND pw.name =~ '(?i).*phagocyt[o,i].*' RETURN DISTINCT b.name")$b.name
antprocmods<-cypher(graph,"MATCH (b:baylor {edge:'5'})-[r]-(pw) WHERE (pw:reactomePW OR pw:ImmunePW OR pw:PalWangPW) AND pw.name =~ '(?i).*antigen.*proc.*pres.*' RETURN DISTINCT b.name")$b.name
glycolmods<-cypher(graph,"MATCH (b:baylor {edge:'5'})-[r]-(pw) WHERE (pw:reactomePW OR pw:ImmunePW OR pw:PalWangPW) AND pw.name =~ '(?i).*glycoly[s,t].*' RETURN DISTINCT b.name")$b.name
oxphosmods<-cypher(graph,"MATCH (b:baylor {edge:'5'})-[r]-(pw) WHERE (pw:reactomePW OR pw:ImmunePW OR pw:PalWangPW) AND pw.name =~ '(?i).*oxidat.*phos.*' RETURN DISTINCT b.name")$b.name
ribosomemods<-cypher(graph,"MATCH (b:baylor {edge:'5'})-[r]-(pw) WHERE (pw:reactomePW OR pw:ImmunePW OR pw:PalWangPW) AND pw.name =~ '(?i).*ribosom[a,e].*' RETURN DISTINCT b.name")$b.name
lipidmods<-cypher(graph,"MATCH (b:baylor {edge:'5'})-[r]-(pw) WHERE (pw:reactomePW OR pw:ImmunePW OR pw:PalWangPW) AND pw.name =~ '(?i).*lipid.*|.*cholesterol.*' RETURN DISTINCT b.name")$b.name
prolifmods<-cypher(graph,"MATCH (b:baylor {edge:'5'})-[r]-(pw) WHERE (pw:reactomePW OR pw:ImmunePW OR pw:PalWangPW) AND pw.name =~ '(?i).*cell.*cycle*|.*mitosis*|.*checkpoint.*|.*cyclin.*' RETURN DISTINCT b.name")$b.name

allmods<-cypher(graph,"MATCH (b:baylor) RETURN DISTINCT b.name as all")
allmods<-allmods$all

subset.list<-list("All"=all,"Neutrophils"=neutromods,"Monocytes"=monomods,"CD4 T cells"=cd4mods,"B cells"=bcellmods,"NK cells"=nkmods,"Lymphocytes"=lymmods,
                  "Red blood cells"=rbcmods,"Platelets"=pltmods,"Interferon signaling"=ifnmods,"Lysosome"=lysosomemods,
                  "Phagocytosis"=phagocytmods,"Antigen proc. pres."=antprocmods,"Glycolysis"=glycolmods,"Oxidative phosphorylation"=oxphosmods,"Ribosome"=ribosomemods,"Lipids"=lipidmods,"Cell proliferation"=prolifmods)


subset.list2<-list("All"=allmods,"Neutrophils"=neutromods,"Monocytes"=monomods,"CD4 T cells"=cd4mods,"B cells"=bcellmods,"NK cells"=nkmods,"Lymphocytes"=lymmods,
                  "Red blood cells"=rbcmods,"Platelets"=pltmods,"Interferon signaling"=ifnmods,"Lysosome"=lysosomemods,
                  "Phagocytosis"=phagocytmods,"Antigen proc. pres."=antprocmods,"Glycolysis"=glycolmods,"Oxidative phosphorylation"=oxphosmods,"Ribosome"=ribosomemods,"Lipids"=lipidmods,"Cell proliferation"=prolifmods)

sdf<-matrix(nrow=(length(subset.list2)-1),ncol=length(allmods),data=0)

rown<-names(subset.list2)[2:(length(subset.list2))]
print(rown)
print(dim(sdf))

row.names(sdf)<-rown
colnames(sdf)<-allmods

#dimnames(sdf)<-list(rown,allmods)
print(sdf)
for (sub in 2:length(subset.list)){
  name<-names(subset.list)[sub]
  sdf[name,]<-jitter(as.numeric(allmods%in%subset.list[[sub]]),factor=1,amount=0.0001)
  sdf[name,]<-as.character(allmods%in%subset.list[[sub]])
}

print(sdf)
sdf<-t(sdf)
sdf<-as.data.frame(sdf,row.names = allmods)
#datafiles<-paste(data$sets,".R",sep="")
datafiles<-paste(unlist2(setlistnameslist),".R",sep="")

#new
data$sets<-sort(data$sets)
#//

datasetmatrix<-data.frame(matrix(nrow=length(data$sets),ncol=5,data=1),row.names=data$sets)
colnames(datasetmatrix)<-c("edge1","edge2","edge3","edge4","edge5")
print("datasetmatrix")
print(datasetmatrix)


for (set in data$sets){
  line=character()
  for (edge in 1:5){
  edgesChar<-cypher(graph,paste("MATCH (n:wgcna {square:'",set,"',edge:'",edge,"'}) RETURN DISTINCT n.contrast AS comparison",sep=""))
  line<-c(line,edgesChar$comparison)
  }#edge
  
  groupline<-character()
  for (group in c(3,1)){
    builder<-line[group]
    groupA<-unlist(strsplit(builder," - "))[2]
    groupB<-unlist(strsplit(builder," - "))[1]
    groupAB<-c(groupA,groupB)
    groupline<-c(groupline,groupAB)
    }#groups
  groupline<-c(groupline,"All")
  
  for (edge in 1:5){
    datasetmatrix[set,edge]<-paste(line[edge],groupline[edge],sep="\n")
  }#edge
}#set



##end required for MM

#Required for G2
#setlislist stuff in G1 dependencies, it is newer
data<-cypher(graph,"MATCH (n:wgcna) RETURN DISTINCT n.square AS sets")
#datafiles<-paste(data$sets,".R",sep="")
datafiles<-paste(unlist2(setlistnameslist),".R",sep="")
#load(file=file.path(cellcordir,"csl.RData"))#precomputed dependency

#datasetmatrix<-data.frame(matrix(nrow=length(data$sets),ncol=5,data=1),row.names=data$sets)
#print(rownames(datasetmatrix))
#colnames(datasetmatrix)<-c("edge1","edge2","edge3","edge4","edge5")


allcellsDF<-cypher(graph,"MATCH (c1:CELL)-[r]-(c:cellEx) RETURN DISTINCT c1.name as cellgroup, c.name as cell")
allcells<-apply(allcellsDF,1,function(x){paste(x[1],x[2],sep=".")})
print(sort(allcells))
#End required for G2

##Required for Machines
#Individual inflammasomes####
annot<-read.csv(file.path("/home/rstudio/source_data","inflammasomes.csv"))
annot2<-unique(annot[,2:3])
annot2<-annot2[which(annot2$Symbol!=""),]
generegexlist<-list(
  "NLRP1"="NLRP1|PYCARD|CASP1|GSDMD|IL1|IL18|SYK|MAPK8|CARD16|CARD17|CARD18|CASP12",
  "NLRP3"="NLRP3|NEK7|PYCARD|CASP1|GSDMD|IL1|IL18|BRCC3|NOS2|POP1|POP2|SYK|MAPK8|CARD16|CARD17|CARD18|CASP12",
  "NLRC4"="NAIP|NLRC4|PYCARD|CASP1|GSDMD|IL1|IL18|PRKC|SYK|MAPK8|CARD16|CARD17|CARD18|CASP12",
  "AIM2"="AIM2|GBP2|GBP5|PYCARD|CASP1|GSDMD|IL1|IL18|TP53BP1|POP3|SYK|MAPK8|CARD16|CARD17|CARD18|CASP12",
  "Pyrin"="MEFV|PYCARD|CASP1|GSDMD|IL1|IL18|SYK|MAPK8|CARD16|CARD17|CARD18|CASP12",
  "nonCan"="CASP4|CASP5|GBP1|GBP2|GBP3|GBP4|GBP5|GBP6|GBP7|GSDMD",
  "NLRP6"="NLRP6|PYCARD|CASP1|GSDMD|IL1|IL18|SYK|MAPK8|CARD16|CARD17|CARD18|CASP12",
  "IFI16"="IFI16|PYCARD|CASP1|GSDMD|IL1|IL18|SYK|MAPK8|CARD16|CARD17|CARD18|CASP12",
  "NLRP7"="NLRP7|PYCARD|CASP1|GSDMD|IL1|IL18|SYK|MAPK8|CARD16|CARD17|CARD18|CASP12",
  "NLRP12"="NLRP12|PYCARD|CASP1|GSDMD|IL1|IL18|SYK|MAPK8|CARD16|CARD17|CARD18|CASP12",
  "RIG-I"="DDX58|PYCARD|CASP1|GSDMD|IL1|IL18|SYK|MAPK8|CARD16|CARD17|CARD18|CASP12"
)
squares<-cypher(graph,"MATCH (n:wgcna) RETURN DISTINCT n.square as squares")$squares
labels<-squares
labelsvec<-paste(squares,collapse="|")
print("machines ok")
#Required for cellcore3
#source(file.path(dir.R.files,"cellCore3.R"),local=TRUE)
#source(file.path(dir.R.files,"cellCore2.R"),local=TRUE)
source(file.path(dir.R.files,"cellCore4.R"),local=TRUE)

query<-"MATCH (c:cellEx) RETURN DISTINCT c.name AS cellname2"
cellnamelist2<-cypher(graph,query)
cellnames2<-cellnamelist2$cellname2

query<-"MATCH (c:cellprop) RETURN DISTINCT c.name AS cellname3"
cellnamelist3<-cypher(graph,query)
cellnames3<-cellnamelist3$cellname3

print("cellcore3 ok")

#load all module data
load("/home/rstudio/output/build/allmodules.RData")   
  
# Define UI ####
ui<-fluidPage(
  #shinythemes::themeSelector(),
              titlePanel(paste("ANIMA:",project)),
              sidebarLayout(
                sidebarPanel(
                  tabsetPanel(
                    tabPanel("Phenotype",
                             tags$head(tags$style("#Phenotype{height:90vh !important;}")),
                             tags$head(tags$style("#phenoPlot{height:45vh !important;}")),
                             #tags$head(tags$style("#Phenotype{height:90vh !important;}")),
                             h3(paste('Phenotype:',project),style="color:#FF0000"),
                             fluidRow(column(4,uiOutput("choosesquarePheno"),offset=.33),column(3,uiOutput("chooseVST"),offset=.33),column(4,uiOutput("chooseVT"),offset=.33)),
                             fluidRow(column(4,uiOutput("choosePheno"),offset=.33),column(4,selectInput("fillvarselect","choose fillvar",c("class1","class2","Categories")),offset=.33),column(4,selectInput("bygroupselect","choose var. group",c("vst","vt"),multiple = FALSE),offset=.33)),
                             fluidRow(column(4,selectInput("phenoclasses","choose classes",c("class1","class2","both")),offset=.33),column(4,uiOutput("chooseTable1Data"),offset=.33)),
                             fluidRow(column(4,sliderInput("phenoplotheight","phenoplotheight",500,4500,500,50),offset=.33),column(3,downloadButton("tablenumcsv", "TableNUM")
),column(3,downloadButton("tablecatcsv", "TableCAT")
))
                    ),
                    tabPanel("Gene",
                      tags$head(tags$style("#boxplot{height:100vh !important; overflow-y: scroll}")),
                      h3(paste('Individual probe intensities:',project),style="color:#FF0000"),
                      fluidRow(column(8,uiOutput("choosesquareBXP"),offset=.33),
                               column(4,selectInput("orderP","Boxplot order",list("Group contrast 1"="c(1,3,2,4)","Group contrast 2"="c(1,2,3,4)"),selected=1,multiple=FALSE))),
                      fluidRow(
                        column(8,textInput("regex","Write regex",value="STAT.*")),
                        column(4,sliderInput("plotheight","plot height",1000,6000,1000,500),offset=0.33)
                      ),
                      fluidRow(
                        column(4,sliderInput("ymin","y min",0,20,6,1)),
                        column(4,sliderInput("ymax","y max",0,20,14,1),offset=1)),
                      fluidRow(checkboxInput('returncsv2', 'output csv of differential expression?', FALSE))
                               ),
                    tabPanel("Inflammasomes",
                             tags$head(tags$style("#inflammasomes{height:90vh !important;}")),
                             h3(paste('Inflammasomes:',project),style="color:#FF0000"),
                             fluidRow(column(4,uiOutput("choosesquareMACHINE"),offset=.33)),
                             fluidRow(
                               column(8,selectInput("inflammasome","Select inflammasome",names(generegexlist),selected = names(generegexlist)[1],multiple=FALSE)),
                               column(4,sliderInput("plotheightMachine","plotheightMachine",500,1500,700,150))
                             )
                    ),
                    tabPanel("Signatures",
                             tags$head(tags$style("#signatures{height:90vh !important;}")),
                             tags$head(tags$style("#moduleSankey{height:85vh !important;}")),
                             h3(paste('Signatures:',project),style="color:#FF0000"),
                             fluidRow(
                               #column(4,uiOutput("choosesquareSIG"),offset=.33),
                               
                                 column(4,selectInput("studySIG","Choose study",studies,selected=studies[1],multiple=FALSE)),
                                 column(4,uiOutput("selectSIGset")),
                                 column(2,selectInput("edgeSIG","Choose edge",as.character(1:5),selected='5',multiple=FALSE),offset=.33)
                             ),
                             fluidRow(
                               #column(4,uiOutput("choosesquareSIG"),offset=.33),
                               
                               column(4,selectInput("studySIG2","Choose study",studies,selected=studies[1],multiple=FALSE)),
                               column(4,uiOutput("selectSIGset2")),
                               column(2,selectInput("edgeSIG2","Choose edge",as.character(1:5),selected='5',multiple=FALSE),offset=.33)
                             ),
                             fluidRow(
                               column(6,sliderInput("defSig","define signature",50,8000,100),offset=.33),
                               column(2,selectInput("whichsig","choose subset",c("all","up","down"))),
                               column(2,selectInput("scale","scale",c("none","row","column")))
                             ),
                             fluidRow(
                                column(5,sliderInput("volcanoSIGlow","significance: lower bound",0,30,0,step=.05)),
                               column(5,sliderInput("volcanoSIG","significance: upper bound",0,30,0,step=.05))
                             ),
                             fluidRow(
                               column(5,sliderInput("volcanoFClow","Log2 fold change: lower bound",-5,5,0,step=0.05),offset=.33),
                               column(5,sliderInput("volcanoFC","Log2 fold change: upper bound",-5,5,0,step=0.05),offset=.33)
                               
                             ),fluidRow(
                                column(2,selectInput("vpcol","choose colour scheme",c("logfc","module"),selected="module",multiple=FALSE)),
                                
                                column(2,selectInput("plotcellmat","plot cellprop matrix",c("yes","no"),selected="yes",multiple=FALSE)),
                                column(3,sliderInput("scalefactor","scale second plot",0.7,1.3,1,step=0.005),offset=.33),
                                column(3,sliderInput("plotheight2","plot height",1000,2000,1000,step=100),offset=.33)
                             ),fluidRow(
                                column(3,uiOutput("selectModuleSIG"))
                                
                             )
                             ),
                    #INTERFACE ELEMENTS: Tabs
                  tabPanel("Modules",
                      tags$head(tags$style("#igraphSP{height:85vh !important;}")),
                      tags$head(tags$style("#d3graphSP{height:85vh !important;}")),
                      tags$head(tags$style("#subjectplot{height:85vh !important; overflow-y: scroll}")),
                      tags$head(tags$style("#module2d{height:85vh !important; overflow-y: scroll}")),
                      tags$head(tags$style("#modPheno{height:85vh !important; overflow-y: scroll}")),
                      tags$head(tags$style("#radar{height:85vh !important; overflow-y: scroll}")),
                      tags$head(tags$style("#intracor{height:85vh !important; overflow-y: scroll}")),
                      tags$head(tags$style("#MEvar{height:85vh !important;}")),
                      tags$head(tags$style("#projections{height:85vh !important;}")),
                      tags$head(tags$style("#sankey{height:85vh !important;}")),
                      tags$head(tags$style("#modcor{height:85vh !important;overflow-y: scroll}")),
                  h3("Datasets",style="color:#FF0000"),
                  fluidRow(
                    column(4,uiOutput("choosesquare"),offset=.33),
                    column(2,selectInput("edge","Choose edge",as.character(1:5),selected='5',multiple=FALSE),offset=.33),
                    column(2,sliderInput("factor","Font Magnification",0.2,5,1,0.1),offset=.33)
                  ),
                  h3("All modules",style="color:#FF0000"),
                  #fluidRow(
                  #column(3,uiOutput("menamegen2")))
                  #,
                  
                  h4("Bipartite graphs",style="color:#32C3EE"),
                  fluidRow(
                     column(3,selectInput("nodetype_proj","Choose nodetype 1",c("wgcna","baylor","reactomePW","ImmunePW","PalWangPW","cellEx","cellprop","pheno"),selected="wgnca",multiple=FALSE)),
                     column(3,selectInput("nodetype_proj2","Choose nodetype 2",c("wgcna","baylor","reactomePW","ImmunePW","PalWangPW","cellEx","cellprop","pheno","SYMBOL"),selected="baylor",multiple=FALSE)),
                     column(3,selectInput("plotwhich","Choose which plot",c("graph"=1,"proj1"=2,"proj2"=3),multiple=FALSE))
                  ),
                  fluidRow(column(3,uiOutput("chooseVST2"),offset=.33),column(3,uiOutput("chooseVT2"),offset=.33),column(3,selectInput("bygroupselect2","choose var. group",c("all","vst","vt"),multiple = FALSE),offset=.33)),
                  
                  h4("igraph controls",style="color:#0000FF"),
                  fluidRow(
                     column(3,selectInput("layout","Choose layout",c("layout_with_fr","layout_with_dh","layout_nicely","layout_with_gem","layout_as_bipartite"),selected="layout_with_fr",multiple=FALSE)),
                     column(2,sliderInput("vlc","Label size ",0.5,2,1,0.001)),
                     column(2,sliderInput("legendcex","Legend size ",0.5,2,1,0.001)),
                     column(2,sliderInput("vertexsize","Vertex size",2,30,15,0.1)),
                     column(2,sliderInput("edgeWidth","Edge width",0.5,5,0.5,0.1))
                  ),
                  fluidRow(
                     #column(2.5,sliderInput("phenoWeightLow","phenoWeightLow",0,1,0,0.05)),
                     #column(2.5,sliderInput("phenoWeightHigh","phenoWeightHigh",0,1,1,0.05)),
                     column(3,sliderInput("phenoWeightLow","RSQlow",0,1,0,0.05)),
                     column(3,sliderInput("phenoWeightHigh","RSQlow",0,1,1,0.05))
                  ),
                  
                  h4("Sankey plots",style="color:#32C3EE"),
                  fluidRow(
                     column(3,selectInput("nodetype_sankey1","Choose nodetype 1",c("wgcna","ImmunePW","cellEx"),selected="wgnca",multiple=FALSE)),
                     column(3,selectInput("nodetype_sankey2","Choose nodetype 2",c("ImmunePW","cellEx"),selected="cellEx",multiple=FALSE)),
                     column(3,selectInput("linkgroup","Choose link group",c("node1","node2")))),
                  
                  h3("Intracor",style="color:#FF0000"),
                  fluidRow(
                  #column(4,uiOutput("allsubsetnames")),
                  #column(4,selectInput("doIntracor","intracor?",c(FALSE,TRUE),selected=FALSE,multiple=FALSE)),
                  #column(4,selectInput("matchedset","Choose matchedset",c("Q_9_blood.PCF.defTB","Q_10_blood.PCF.probTB","Q_11_blood.PCF.defPC","Q_12_blood.PCF.probPC","Q_13_blood.PCF.HIVneg","Q_14_blood.PCF.HIVpos","Q_15_blood.PCF.HIVposHD"),selected="Q_9_blood.PCF.defTB",multiple=FALSE),offset=.33)
                  column(4,selectInput("matchedset","Choose matchedset",c(matchedsets),selected=matchedsets[1],multiple=FALSE),offset=.33),
                  column(3, uiOutput("phenonames2"))
                  #GBP
                  
                  ),
                  h3("Single module",style="color:#FF0000"),
                  h4("Select WGCNA module",style="color:#0000FF"),
                  fluidRow(
                    column(3,uiOutput("menamegen1")),
                    column(3,uiOutput("menamegen2")),
                    column(2, uiOutput("subsetdefs1")),
                    column(2,uiOutput("subsetdefs2"))
                  ),
                  
                  
                  h4("Module-phenotype correlation",style="color:#0000FF"),
                  fluidRow(
                    column(3, uiOutput("phenonames")),
                    column(2,selectInput("method","Method",c("pearson","spearman"),selected="pearson")),
                    column(2,selectInput("legpos","Legend",c("bottomleft","bottomright","topleft","topright")))#,
                    #column(2, uiOutput("subsetdefs1")),
                    #column(2,uiOutput("subsetdefs2"))
                  ),
                  
                  h4(paste('Modcor'),style="color:#0000FF"),
                  
                  fluidRow(
                    column(2,selectInput("studyModcor","Choose study",studies,selected=studies[1],multiple=FALSE)),
                    column(4,uiOutput("selectMCset")),
                    column(3,uiOutput("selectModuleModcor")),
                    column(2,selectInput("edgeModcor","select edge",as.character(1:5),selected='1',multiple=FALSE))
                  ),
                  
                  
                  
                  
                  ),
                  
                  tabPanel("Virtual Cells",
                           tags$head(tags$style("#cellcor{height:90vh !important;overflow-y: scroll}")),
                           tags$head(tags$style("#cellmatrix{height:90vh !important;overflow-y: scroll}")),
                           tags$head(tags$style("#gigabar{height:80vh !important;}")),
                 

                  h3(paste('Cell Matrix'),style="color:#FF0000"),
                  fluidRow(
                     column(4,selectInput("study","Choose study",studies,selected=studies[1],multiple=FALSE)),
                     column(4,uiOutput("selectCCset")),
                     column(2,selectInput("edges","select edge",as.character(1:5),selected='1',multiple=FALSE))
                  ),
                  fluidRow(
                    column(3,sliderInput("sliderpw","Pathway filter",0,3,0,0.001)),
                    column(3,sliderInput("slidercell","Cell filter",0,3,0,0.001),offset=0.33),
                    column(3,sliderInput("barfilter1","barfilter1",2,20,5,1),offset=0.33),
                    column(3,sliderInput("barfilter2","barfilter2",2,30,5,1),offset=0.33)
                  )
                  ,
                  fluidRow(
                    column(4,sliderInput("plotheightG1","plotheightG1",1000,6000,1000,500)),
                    column(4,selectInput("cellStudy","Select cell study",c("xp_1","xp_2","xp_3"),selected=1:3,multiple=TRUE))
                  ),
                  
                  h3(paste('Single Cell'),style="color:#FF0000"),
                  # fluidRow(
                  #   column(4,selectInput("study","Choose study",studies,selected=studies[1],multiple=FALSE)),
                  #   column(4,conditionalPanel(
                  #     condition = "input.study == 1",
                  #     selectInput("datasets1","Choose dataset",setlist1,selected=setlist1[1],multiple=FALSE)
                  #   ),
                  #   conditionalPanel(
                  #     condition = "input.study == 2",
                  #     selectInput("datasets2","Choose dataset",setlist2,selected=setlist2[1],multiple=FALSE)
                  #   ),
                  #   conditionalPanel(
                  #     condition = "input.study == 3",
                  #     selectInput("datasets3","Choose dataset",setlist3,selected=setlist3[1],multiple=FALSE)
                  #   ),
                  #   conditionalPanel(
                  #     condition = "input.study == 4",
                  #     selectInput("datasets4","Choose dataset",setlist4,selected=setlist4[1],multiple=FALSE)),offset=.5)
                  # ),
                  #new
                  
                  
                  #end new
                  fluidRow(
                    column(4,selectInput("cellname","Choose cell name",cellnames2))
                    
                  )
                  ),#endG1 tab,
                  
                  
                  
                  #submitButton(text = "Apply Changes", icon = NULL, width = NULL),
                  tabPanel("Meta-Analysis",
                  tags$head(tags$style("#ModuleMeta{height:40vh !important;}")),
                  tags$head(tags$style("#d3graphMM{height:80vh !important;}")),
                  tags$head(tags$style("#igraph{height:90vh !important;}")),
                  tags$head(tags$style("#gigamat2{height:90vh !important;}")),
                  h3('ModuleMeta',style="color:#FF0000"),
                  fluidRow(
                    column(12,DT::dataTableOutput('x1'))),
                    column(4,fluidRow(selectInput("subset","Choose biologic processes",names(subset.list),selected="all",multiple=FALSE)),
                           fluidRow(checkboxInput("allsame","Show all same",FALSE)),
                           fluidRow(checkboxInput("different","Show differences",FALSE)),
                           fluidRow(checkboxInput('returncsv', 'output csv of module annotations?', FALSE))
                           
                    ),
                  fluidRow(
                    column(3,selectInput("layout1","Choose layout",c("layout_with_fr","layout_with_dh","layout_nicely","layout_with_gem"),selected="layout_with_fr",multiple=FALSE)),
                    column(3,sliderInput("vlc1","label size adjustment",0.5,2,1,0.001)),
                    column(3,sliderInput("legendcex1","legend size adjustment",0.5,2,1,0.001)),
                    column(3,sliderInput("vertexsize1","vertex size",2,30,15,0.1))
                  ),
                  fluidRow(
                    column(4,sliderInput("plotheightMeta","plotheightMeta",500,1500,700,150)),
                    column(4,sliderInput("moduleThresh","moduleThresh",-0.05,1,-0.05,0.01))
                  ),
                  
                  h3(paste('CellMatrix:',project),style="color:#FF0000"),
                  #DT::dataTableOutput('x1'),
                  # fluidRow(
                  #   column(12,DT::dataTableOutput('x2'))),
                  fluidRow(
                    column(8,selectInput("cells","Choose cells",allcells,selected = allcells[c(1,5,10)],multiple=TRUE)),
                    column(4,sliderInput("plotheightG2","plotheightG2",1000,6000,1000,500))
                  )#,
                  #sliderInput("sliderpw","Pathway filter",0,3,0,0.001),
                  #sliderInput("slidercell","Cell filter",0,3,0,0.001),
                  #checkboxInput('returnpdf', 'output pdf?', FALSE),
                  #downloadButton("downloadPlot","Save plot")
                  ),
                  
                  
                  
                  tabPanel("Save",br(),fluidRow(
                     column(3,checkboxInput('returnpdf', 'output pdf?', FALSE)),
                     column(3,downloadButton("downloadPlot","Save plot"),offset=1),
                     column(3,sliderInput("downloadPlotWidth","Plot width",7,14,12,0.5),offset=1),
                     column(3,sliderInput("downloadPlotHeight","Plot height",7,30,7,0.5)
                            ),offset=1))
                  , selected = "Modules")#end tabset
                
                
                ),#end sidebar
                #mainPanel####
                mainPanel(
                  
                  
                    
                  tabsetPanel(
                    #design
                    tabPanel("Study Design",imageOutput("design")),
                    #pheno
                    tabPanel("Phenotype",
                             tabsetPanel(
                               tabPanel("Subject numbers",tableOutput('subjectnumbers')),
                               tabPanel("Table 1",tableOutput('table1')),
                               tabPanel("allData",DT::dataTableOutput('allphenodata')),
                               tabPanel("table_NUM",DT::dataTableOutput('tableNUM')),
                               tabPanel("table_NUM2",tableOutput('tableNUM2')),
                               tabPanel("table_CAT",DT::dataTableOutput('tableCAT')),
                               tabPanel("table_CAT2",tableOutput('tableCAT2')),
                               tabPanel("phenoDT",DT::dataTableOutput('phenodatatable')),
                               tabPanel("phenoPlot",uiOutput('phenoPlot')),
                               tabPanel("Summary",tableOutput("summary")),
                               tabPanel("STATS",tableOutput("stats")),
                               tabPanel("phenoGridPlot",uiOutput('phenoGridPlot')),
                               tabPanel("phenoGridPlot2",uiOutput('phenoGridPlot2')),
                               tabPanel("cellpropGridPlot",uiOutput('cellpropGridPlot')),
                               tabPanel("wgcnaMEGridPlot",uiOutput('wgcnaMEGridPlot')),
                               tabPanel("code",verbatimTextOutput("code"))
                             )
                    ),
                    #BXP
                    tabPanel("Gene",
                             #tabsetPanel("BXP",
                             tabsetPanel(             
                                         tabPanel("selected genes",tableOutput("genelist")),
                                         tabPanel("boxplot",plotOutput("boxplot")),
                                         tabPanel("detableBXP",DT::dataTableOutput('detableBXP'))
                             )
                             
                    ),
                    
                    #Inflammasomes
                    tabPanel("Inflammasomes",plotOutput("inflammasomes")),
                    
                    tabPanel("Signatures",
                              tabsetPanel(
                                          tabPanel("Signature table",DT::dataTableOutput("sigTable1")),
                                          tabPanel("Volcano Plot",plotOutput("volcano")),
                                          tabPanel("Heatmap",uiOutput("heatmap")),
                                          tabPanel("ModuleMaps",sankeyNetworkOutput("moduleSankey")),
                                          tabPanel("Enrichment table",DT::dataTableOutput("enrichTable1")),
                                          tabPanel("Barplot",plotOutput("barplot1")),
                                          tabPanel("Dotplot",plotOutput("dotplot1")),
                                          tabPanel("emapplot",plotOutput("emapplot1")),
                                          tabPanel("cnetplot",plotOutput("cnetplot1")),
                                          tabPanel("GSEA table",DT::dataTableOutput("gseaTable2")),
                                          tabPanel("GSEA emapplot",plotOutput("emapplot2")),
                                          tabPanel("gseaplot",plotOutput("gseaplot2"))
                                         )


                              ),
                    # 
                    
                    tabPanel("All modules",
                        #tabsetPanel("All modules",
                        tabsetPanel(
                          tabPanel("module2d",plotOutput("module2d")),
                            tabPanel("module stats",DT::dataTableOutput("dtwgcna")),
                            tabPanel("module table",htmlOutput("moduleTable")),
                            tabPanel("moduleHM",plotOutput("moduleHM")),
                            tabPanel("sigenrichHM",plotOutput("sigenrichHM")),
                             tabPanel("bipartite",plotOutput("projections")),
                             tabPanel("sankey",sankeyNetworkOutput("sankey"))
                             )),
                    
                    tabPanel("Intracor",
                             #tabsetPanel("All modules",
                             tabsetPanel(
                                tabPanel("intracor",plotOutput("intracor")),
                                tabPanel("intracorModules",DT::dataTableOutput("dtintracor")),
                                tabPanel("intracorPheno",DT::dataTableOutput("dtintracorpheno")),
                                tabPanel("intracorPhenoPlot",plotOutput("phenointracor"))
                             )),
                    
                    tabPanel("Single module",
                             #tabsetPanel("Single module",
                             tabsetPanel(
                                tabPanel("igraphSP",plotOutput("igraphSP")),
                                tabPanel("d3graphSP",forceNetworkOutput("d3graphSP",width="100%",height="400px")),
                                tabPanel("moduleNodesSP",DT::dataTableOutput('moduleNodesSP')),
                                tabPanel("moduleEdgesSP",DT::dataTableOutput('moduleEdgesSP')),
                                tabPanel("radarplots",plotOutput("radar")),
                                tabPanel("modcor",plotOutput("modcor")),
                                tabPanel("MEvar",plotOutput("MEvar")),
                                tabPanel("subjectPlot",plotOutput("subjectplot")),
                                tabPanel("modPheno",plotOutput("modPheno"))
                                #, selected = "d3graphSP")),
                             )),
                    
                    tabPanel("Virtual Cells",
                             #tabsetPanel("cellcor",
                             tabsetPanel(
                                         
                                         tabPanel("Cell Matrix",uiOutput("cellmatrix")),
                                         tabPanel("Cell and pathway summary",plotOutput("gigabar")),
                                         tabPanel("Single Cell",plotOutput("cellcor"))
                             )),
                    
                    
                    tabPanel("Meta-Analysis",
                             #tabsetPanel("Chaussabel",
                             tabsetPanel(            
                                         #tabPanel("ModuleMeta",plotOutput("pheat")),
                                         tabPanel("ModuleMeta",uiOutput("plot.pheat")),
                                         tabPanel("igraph",plotOutput("igraph")),
                                         tabPanel("d3graphMM",forceNetworkOutput("d3graphMM",width="100%",height="100%")),
                                         tabPanel("moduleNodesMM",DT::dataTableOutput('moduleNodesMM')),
                                         tabPanel("moduleEdgesMM",DT::dataTableOutput('moduleEdgesMM')),
                                         tabPanel("wgcna",plotOutput("wgcnacol")),
                                         tabPanel("Virtual Cells",uiOutput("cellmatrixMeta")),
                                         tabPanel("Cell type ratios",plotOutput("ratioFarm")))
                             
                     )#,
                    # tabPanel("Shiny Session",
                    #          verbatimTextOutput("clientdataText")),
                    # tabPanel("Output info",
                    #          textOutput("output$info "))
                    
                    
                    
                    #,selected="Single module")
                    )#delete this if uncommenting above
                )
                  
              )
)


# Define server logic ####
server <- shinyServer(function(input, output,session) {
   # Store in a convenience variable
   cdata <- session$clientData
   
   
   #libraries
  #library(RNeo4j)
   library(neo2R)
  #graph<-startGraph(graphstring)#!!!!! note the address!!
  graph<-startGraph(graphstring,database="animadb",username="neo4j",password="anima")
  #REACTIVE FUNCTIONS####
   #SP
   setfun<-reactive({input$set})
   setfun2<-reactive({input$matchedset})
   me1fun<-reactive({input$me1})
   me2fun<-reactive({input$me2})
   s1fun<-reactive({input$subset1})
   s2fun<-reactive({input$subset2})
   phenofun<-reactive({input$pheno})
   phenofunIC<-reactive({input$pheno2})
   methodfun<-reactive({input$method})
   edgefun<-reactive({input$edge})
   #buildfun<-reactive({input$build})
   nodefun<-reactive({input$nodetype_proj})
   nodefun2<-reactive({input$nodetype_proj2})
   plotwhichfun<-reactive({input$plotwhich})
   legposfun<-reactive({input$legpos})
   
  # new factor function
   factorfun<-reactive({input$factor})
   
   sankeynode1fun<-reactive({input$nodetype_sankey1})
   sankeynode2fun<-reactive({input$nodetype_sankey2})
   sankeycell1fun<-reactive({input$cellSet1})
   sankeycell2fun<-reactive({input$cellSet2})
   linkgroupfun<-reactive({input$linkgroup})
   #modcor
   studyModcorfun<-reactive({input$studyModcor})
   squareModcorfun<-reactive({c(input$datasetsMC)})
   edgeModcorfun<-reactive({input$edgeModcor})
   moduleModcorfun<-reactive({input$moduleSelect})
   #cellcor
   studyCellcorfun<-reactive({input$study})
   squareCellcorfun<-reactive({c(input$datasetsCC)})
   edgeCellcorfun<-reactive({input$edges})
   cellCfun<-reactive({input$cellname})
   
   #Signatures
   studySIGfun<-reactive({input$studySIG})
   squareSIGfun<-reactive({c(input$datasetsSIG)})
   edgeSIGfun<-reactive({input$edgeSIG})
   
   studySIGfun2<-reactive({input$studySIG2})
   squareSIGfun2<-reactive({c(input$datasetsSIG2)})
   edgeSIGfun2<-reactive({input$edgeSIG2})
   
   defSigfun<-reactive({input$defSig})
   
   whichSIGfun<-reactive({input$whichsig})
   scaleSIGfun<-reactive({input$scale})
   vpplotfun<-reactive({input$vpcol})
   cellmatfun<-reactive({input$plotcellmat})
   heatheightfun<-reactive({input$plotheight2})
   scalefactorfun<-reactive(input$scalefactor)
   
   whichModuleFun<-reactive({input$moduleSelectSIG})
   
   #Pheno
   squarephenofun<-metaReactive({..(input$PhenoSquare)})
   phenofun2<-reactive({input$phenochoice})
   phenoplotheightfun<-reactive({input$phenoplotheight})
   classesfun<-reactive({input$phenoclasses})
   phenovstfun<-reactive({input$phenovst})
   phenovtfun<-reactive({input$phenovt})
   fillvarfun<-reactive({input$fillvarselect})
   bygroupfun<-reactive({input$bygroupselect})
   tab1varfun<-reactive({input$table1vars})
   
   
   #Phenobipartite
   phenovstfun2<-reactive({input$phenovst2})
   phenovtfun2<-reactive({input$phenovt2})
   bygroupfun2<-reactive({input$bygroupselect2})
   
   #G1
   usestudy<-reactive({input$study})
   useset<-reactive({c(input$datasetsCC)})
   useedge<-reactive({input$edges})
   slider1<-reactive({input$sliderpw})
   slider2<-reactive({input$slidercell})
   barfilterfun1<-reactive({input$barfilter1})
   barfilterfun2<-reactive({input$barfilter2})
   plotheightCM<-reactive({input$plotheightG1})
   #BXP
   setfunBXP<-reactive({input$BXPset})
   orderPfun<-reactive({input$orderP})
   useregex<-reactive({input$regex})
   yminfun<-reactive({input$ymin})
   ymaxfun<-reactive({input$ymax})
   returncsvfun2<-reactive({input$returncsv2})
   #MM
   useallsame<-reactive({input$allsame})
   usediff<-reactive({input$different})
   plotheightM<-reactive({input$plotheightMeta})
   threshfun<-reactive({input$moduleThresh})
   
   #igraph_plotter controls single module and bipartite
   layoutfun<-reactive({input$layout})
   vlc<-reactive({input$vlc})
   legendcexfun<-reactive({input$legendcex})
   vertexsizefun<-reactive({input$vertexsize})
   ewf<-reactive({input$edgeWidth})
   phenoWeightLowfun<-reactive({input$phenoWeightLow})
   phenoWeightHighfun<-reactive({input$phenoWeightHigh})
   #igraph_plotter controls meta-analysis
   layoutfun1<-reactive({input$layout1})
   vlc1<-reactive({input$vlc1})
   legendcexfun1<-reactive({input$legendcex1})
   vertexsizefun1<-reactive({input$vertexsize1})
   
   usematrixfun<-reactive({input$x1_cells_selected})
   returncsvfun<-reactive({input$returncsv})
   #process dataset selection for UI
   output$x1 = DT::renderDataTable(
     datasetmatrix, server = FALSE,
     selection = list(target = 'cell', mode = 'multiple', selected = matrix(c(1,2,3,5,5,5),ncol=2)),
     options=list(lengthMenu = c(5, 10, 15,20,50), pageLength = 5),filter="top"
   ) 
   #G2
   #slider1<-reactive({input$sliderpw})
   #slider2<-reactive({input$slidercell})
   usematrixfun<-reactive({input$x1_cells_selected})
   # output$x2 = DT::renderDataTable(
   #   datasetmatrix, server = FALSE,
   #   selection = list(target = 'cell', mode = 'multiple', selected = matrix(c(1,2,3,5,5,5),ncol=2)),
   #   options=list(lengthMenu = c(5, 10, 15,20,50), pageLength = 5)
   # )
   cellfun<-reactive({input$cells})
   
   #Inflammasomes
   grxfun<-reactive({input$inflammasome})
   setfunMACHINE<-reactive({input$MACHINEset})
   
   downloadPlotWidthFun<-reactive({input$downloadPlotWidth})
   downloadPlotHeightFun<-reactive({input$downloadPlotHeight})
   
   #REACTIVE UI INPUT FUNCTIONS####
   
   #SP
   
   output$allsubsetnames<-renderUI({
     #prefix<-paste("/home/rstudio/output/build/",buildfun(),"/",sep="")
     prefix<-"/home/rstudio/output/build/"
     dirlisting<-system(paste("ls",prefix),intern=TRUE)
     allsubsetnames<-dirlisting[grep("Q_",dirlisting)]
     print("prefix:")
     print(prefix)
     print("dirlisting")
     print(dirlisting)
     print("allsubsetnames")
     print(allsubsetnames)
     #selectInput("build","Choose build",buildlisting,selected=buildlisting[length(buildlisting)],multiple=FALSE)
   })
   
   output$choosesquare<-renderUI({
     print("invoking input$set in renderUI function:output$choosesquare")
     #prefix<-paste("/home/rstudio/output/build/",buildfun(),"/",sep="")
     prefix<-"/home/rstudio/output/build/"
     dirlisting<-system(paste("ls",prefix),intern=TRUE)
     allsubsetnames<-dirlisting[grep("Q_",dirlisting)]
     selectInput("set","Choose square",allsubsetnames,selected="1",multiple=FALSE)
   })
   
   output$choosesquareBXP<-renderUI({
     print("invoking input$set in renderUI function:output$choosesquareBXP")
     #prefix<-paste("/home/rstudio/output/build/",buildfun(),"/",sep="")
     prefix<-"/home/rstudio/output/build/"
     dirlisting<-system(paste("ls",prefix),intern=TRUE)
     allsubsetnames<-dirlisting[grep("Q_",dirlisting)]
     selectInput("BXPset","Choose square",allsubsetnames,selected="1",multiple=FALSE)
   })
   
   output$choosesquareMACHINE<-renderUI({
     print("invoking input$set in renderUI function:output$choosesquareMACHINE")
     #prefix<-paste("/home/rstudio/output/build/",buildfun(),"/",sep="")
     prefix<-"/home/rstudio/output/build/"
     dirlisting<-system(paste("ls",prefix),intern=TRUE)
     allsubsetnames<-dirlisting[grep("Q_",dirlisting)]
     selectInput("MACHINEset","Choose square",allsubsetnames,selected="1",multiple=FALSE)
   })
   

   output$choosesquarePheno<-renderUI({
     print("invoking input$PhenoSquare in renderUI function:output$choosesquarePheno")
     query<-"MATCH (n:personPheno) RETURN DISTINCT n.square as square"
     res<-cypher(graph,query)
     squarevec<-res$square
     selectInput("PhenoSquare","Choose square",squarevec,selected="1",multiple=FALSE)
   })
   
   output$chooseVST<-renderUI({
      
      
      query<-paste("MATCH (p:personPheno {square:'",squarephenofun(),"'})-[r]-(vst:varsubtype)-[r2]-(vt:vartype) RETURN DISTINCT vst.name AS vst, vt. name as vt",sep="")
      groups<-cypher(graph,query)
      vst.list<-unique(groups$vst)
      vt.list<-unique(groups$vt)
      selectInput("phenovst","Choose VST",vst.list,selected="1",multiple=FALSE)
   })
   
   output$chooseTable1Data<-renderUI({
     
     
     query<-paste("MATCH (p:personPheno {square:'",squarephenofun(),"'})-[r]-(vst:varsubtype)-[r2]-(vt:vartype) RETURN DISTINCT vst.name AS vst, vt. name as vt",sep="")
     groups<-cypher(graph,query)
     vst.list<-unique(groups$vst)
     vt.list<-unique(groups$vt)
     selectInput("table1vars","Choose T1 vars",vst.list,selected="1",multiple = TRUE)
   })
   
   output$chooseVT<-renderUI({
     
      query<-paste("MATCH (p:personPheno {square:'",squarephenofun(),"'})-[r]-(vst:varsubtype)-[r2]-(vt:vartype) RETURN DISTINCT vst.name AS vst, vt. name as vt",sep="")
      groups<-cypher(graph,query)
      vst.list<-unique(groups$vst)
      vt.list<-unique(groups$vt)
      selectInput("phenovt","Choose VT",vt.list,selected="1",multiple=FALSE)
   })
   
   output$chooseVST2<-renderUI({
     squarelist<-unlist(strsplit(setfun(),"_"))
     if(length(squarelist)==3){
       square<-squarelist[3]
     }else{square<-paste(squarelist[c(3,4)],collapse="_")}
     
     query<-paste("MATCH (p:personPheno {square:'",square,"'})-[r]-(vst:varsubtype)-[r2]-(vt:vartype) RETURN DISTINCT vst.name AS vst, vt. name as vt",sep="")
     groups<-cypher(graph,query)
     vst.list2<-unique(groups$vst)
     vt.list2<-unique(groups$vt)
     selectInput("phenovst2","Choose VST",vst.list2,selected="1",multiple=FALSE)
   })
   
   output$chooseVT2<-renderUI({
     squarelist<-unlist(strsplit(setfun(),"_"))
     if(length(squarelist)==3){
       square<-squarelist[3]
     }else{square<-paste(squarelist[c(3,4)],collapse="_")}
     
     
     query<-paste("MATCH (p:personPheno {square:'",square,"'})-[r]-(vst:varsubtype)-[r2]-(vt:vartype) RETURN DISTINCT vst.name AS vst, vt. name as vt",sep="")
     groups<-cypher(graph,query)
     vst.list2<-unique(groups$vst)
     vt.list2<-unique(groups$vt)
     selectInput("phenovt2","Choose VT",vt.list2,selected="1",multiple=FALSE)
   })
   
   output$choosePheno<-renderUI({
     square<-squarephenofun()
     query<-paste("MATCH (n:personPheno {square:'",square,"'}) RETURN DISTINCT n.name as pheno",sep="")
     res<-cypher(graph,query)
     phenovec<-res$pheno
     selectInput("phenochoice","Choose pheno",phenovec,selected="1",multiple=FALSE)
   })
   
   output$menamegen1<-renderUI({
     print("computing menamegen1")
     #prefix<-paste("/home/rstudio/output/build/",buildfun(),"/",sep="")
     prefix<-"/home/rstudio/output/build/"
     print("prefix2:")
     print(prefix)
     print("setfun()")
     print(setfun())
     data0<-read.csv(paste(prefix,setfun(),"/5_WGCNA_4k_tv/results/ModuleEigengenes_all_samples_withPheno.csv",sep=""))
     #print("uiline1")
     data<-read.csv(paste(prefix,setfun(),"/5_WGCNA_4k_tv/results/ModuleEigengenes_all_samples.csv",sep=""))
     #print("uiline2")
     menames<-sort(colnames(data)[2:length(colnames(data))])
     selectInput("me1","Choose ME1",menames,selected=1,multiple=FALSE)
     
     })
   output$menamegen2<-renderUI({
     print("computing menamegen2")
     #prefix<-paste("/home/rstudio/output/build/",buildfun(),"/",sep="")
     prefix<-"/home/rstudio/output/build/"
     data0<-read.csv(paste(prefix,setfun(),"/5_WGCNA_4k_tv/results/ModuleEigengenes_all_samples_withPheno.csv",sep=""))
     #print("uiline1")
     data<-read.csv(paste(prefix,setfun(),"/5_WGCNA_4k_tv/results/ModuleEigengenes_all_samples.csv",sep=""))
     #print("uiline2")
     menames<-sort(colnames(data)[2:length(colnames(data))])
     selectInput("me2","Choose ME2",menames,selected=2,multiple=FALSE)
     
   })
   output$phenonames<-renderUI({
     print("computing phenonames")
     #print("line1 inside reactive")
     #prefix<-paste("/home/rstudio/output/build/",buildfun(),"/",sep="")
     prefix<-"/home/rstudio/output/build/"
     data0<-read.csv(paste(prefix,setfun(),"/5_WGCNA_4k_tv/results/ModuleEigengenes_all_samples_withPheno.csv",sep=""))
     #print("line2")
     data<-read.csv(paste(prefix,setfun(),"/5_WGCNA_4k_tv/results/ModuleEigengenes_all_samples.csv",sep=""))
     #print("line3")
     
     squarelist<-unlist(strsplit(setfun(),"_"))
     if(length(squarelist)==3){
       square<-squarelist[3]
       }else{square<-paste(squarelist[c(3,4)],collapse="_")}
     
     menames<-sort(colnames(data)[2:length(colnames(data))])
     phenovec<-colnames(data0[,which(!colnames(data0)%in%menames&!colnames(data0)=="X"&!colnames(data0)=="Row.names")])
     #print("phenovec ui function")
     #print(phenovec)
     selectInput("pheno","Choose pheno",phenovec,selected=2,multiple=FALSE)
   })
   
   output$phenonames2<-renderUI({
      print("computing phenonames")
      #print("line1 inside reactive")
      #prefix<-paste("/home/rstudio/output/build/",buildfun(),"/",sep="")
      prefix<-"/home/rstudio/output/build/"
      data0<-read.csv(paste(prefix,setfun2(),"/5_WGCNA_4k_tv/results/ModuleEigengenes_all_samples_withPheno.csv",sep=""))
      #print("line2")
      data<-read.csv(paste(prefix,setfun2(),"/5_WGCNA_4k_tv/results/ModuleEigengenes_all_samples.csv",sep=""))
      #print("line3")
      
      squarelist<-unlist(strsplit(setfun2(),"_"))
      if(length(squarelist)==3){
         square<-squarelist[3]
      }else{square<-paste(squarelist[c(3,4)],collapse="_")}
      
      menames<-sort(colnames(data)[2:length(colnames(data))])
      phenovec<-colnames(data0[,which(!colnames(data0)%in%menames&!colnames(data0)=="X"&!colnames(data0)=="Row.names")])
      #print("phenovec ui function")
      #print(phenovec)
      selectInput("pheno2","Choose intracor pheno",phenovec,selected=2,multiple=FALSE)
   })
   
   output$subsetdefs1<-renderUI({
     print("computing subsetdefs1")
     #print("line1 inside reactive")
     #prefix<-paste("/home/rstudio/output/build/",buildfun(),"/",sep="")
     prefix<-"/home/rstudio/output/build/"
     data0<-read.csv(paste(prefix,setfun(),"/5_WGCNA_4k_tv/results/ModuleEigengenes_all_samples_withPheno.csv",sep=""))
     #print("line2")
     data<-read.csv(paste(prefix,setfun(),"/5_WGCNA_4k_tv/results/ModuleEigengenes_all_samples.csv",sep=""))
     #print("line3")
     squarelist<-unlist(strsplit(setfun(),"_"))
     if(length(squarelist)==3){
       square<-squarelist[3]
     }else{square<-paste(squarelist[c(3,4)],collapse="_")}
     print("square");print(square)
     contrasts<-cypher(graph,paste("MATCH (n:wgcna {square:'",square,"'}) RETURN DISTINCT n.contrastvar AS contrastvars",sep=""))
     c1<-contrasts$contrastvars[1]
     c2<-contrasts$contrastvars[2]
     print(paste("c2:",c2,sep=""))
     
     classes<-cypher(graph,paste("MATCH (n:wgcna {square:'",square,"'}) RETURN DISTINCT n.contrast AS contrasts",sep=""))
     #cL1<-unlist(strsplit(classes$contrasts[1],"-"))[c(2,1)]
     #cL2<-unlist(strsplit(classes$contrasts[2],"-"))[c(2,1)]
     
     cL1<-sort(unlist(strsplit(classes$contrasts[1]," - ")))
     cL2<-sort(unlist(strsplit(classes$contrasts[2]," - ")))
     print(paste("cL1:",cL1,sep=""))
     print(paste("cL2:",cL2,sep=""))
     
     subset1vec<-eval(parse(text=paste("data0$",c1,sep="")))
     subset1classes<-sort(unique(subset1vec))
     subset2vec<-eval(parse(text=paste("data0$",c2,sep="")))
     subset2classes<-sort(unique(subset2vec))
     print("subset1classes")
     print(subset1classes)
     selectInput("subset1","Subset1",subset1classes,selected=c(1,2),multiple=TRUE)
   })
   
   output$subsetdefs2<-renderUI({
     print("computing subsetdefs2")
     print("line1 inside reactive")
     #print(setfun2())
     #prefix<-paste("/home/rstudio/output/build/",buildfun(),"/",sep="")
     prefix<-"/home/rstudio/output/build/"
     #edit debug replaced setfun2 with setfun
     data0<-read.csv(paste(prefix,setfun(),"/5_WGCNA_4k_tv/results/ModuleEigengenes_all_samples_withPheno.csv",sep=""))
     print("line2")
     data<-read.csv(paste(prefix,setfun(),"/5_WGCNA_4k_tv/results/ModuleEigengenes_all_samples.csv",sep=""))
     print("line3")
     squarelist<-unlist(strsplit(setfun(),"_"))
     print("squarelist")
     print(squarelist)
     if(length(squarelist)==3){
       square<-squarelist[3]
     }else{square<-paste(squarelist[c(3,4)],collapse="_")}
     
     contrasts<-cypher(graph,paste("MATCH (n:wgcna {square:'",square,"'}) RETURN DISTINCT n.contrastvar AS contrastvars",sep=""))
     c1<-contrasts$contrastvars[1]
     c2<-contrasts$contrastvars[2]
     
     classes<-cypher(graph,paste("MATCH (n:wgcna {square:'",square,"'}) RETURN DISTINCT n.contrast AS contrasts",sep=""))
     #cL1<-unlist(strsplit(classes$contrasts[1],"-"))[c(2,1)]
     #cL2<-unlist(strsplit(classes$contrasts[2],"-"))[c(2,1)]
     
     cL1<-sort(unlist(strsplit(classes$contrasts[1]," - ")))
     cL2<-sort(unlist(strsplit(classes$contrasts[2]," - ")))
     
     subset1vec<-eval(parse(text=paste("data0$",c1,sep="")))
     subset1classes<-sort(unique(subset1vec))
     subset2vec<-eval(parse(text=paste("data0$",c2,sep="")))
     subset2classes<-sort(unique(subset2vec))
     selectInput("subset2","Subset2",subset2classes,selected=c(1,2),multiple=TRUE)
   })
   
   #MCset 
   output$selectMCset<-renderUI({
     print("doing MCset")
     print(studyModcorfun())
     selectInput("datasetsMC","Choose dataset",setlistlist[[as.numeric(studyModcorfun())]],selected=setlistlist[[1]],multiple=FALSE)
   })
   #CCset
   output$selectCCset<-renderUI({
     selectInput("datasetsCC","Choose dataset",setlistlist[[as.numeric(studyCellcorfun())]],selected=setlistlist[[1]],multiple=FALSE)
   })
   output$selectSIGset<-renderUI({
     print("doing SIGset")
     print("studySIGfun()")
     print(studySIGfun())
     print(setlistlist[[as.numeric(studySIGfun())]])
     selectInput("datasetsSIG","Choose dataset",setlistlist[[as.numeric(studySIGfun())]],selected=setlistlist[[1]],multiple=FALSE)
   })
   output$selectSIGset2<-renderUI({
     print("doing SIGset2")
     print("studySIGfun2()")
     print(studySIGfun2())
     print(setlistlist[[as.numeric(studySIGfun2())]])
     selectInput("datasetsSIG2","Choose dataset",setlistlist[[as.numeric(studySIGfun2())]],selected=setlistlist[[1]],multiple=FALSE)
   })
   
   
   #modcor
   output$selectModuleSIG<-renderUI({
     print("doing moduleSIG")
     reacstudySIG<-as.numeric(studySIGfun())
     print("reacstudySIG")
     print(reacstudySIG)
     #reacsetMC<-as.numeric(squareModcorfun()[reacstudyMC])#debugging
     reacsetSIG<-as.numeric(squareSIGfun())
     print("squareSIGfun()")
     print(squareSIGfun())
     print("reacsetSIG")
     print(reacsetSIG)
     print("reacsetSIG")
     print(reacsetSIG)
     reacedgeSIG<-as.numeric(edgeSIGfun())
     print("edgeSIGfun")
     print(edgeSIGfun)

     studySIG<-names(studies)[reacstudySIG]
     print("studySIG")
     print(studySIG)

     setlistSIG<-setlistnameslist[[reacstudySIG]]
     print("setlistSIG")
     print(setlistSIG)
    
     squareSIG<-setlistSIG[[reacsetSIG]]
     print("squareSIG")
     print(squareSIG)
     
     moduleq<-paste("MATCH (n:wgcna {square:'",squareSIG,"'}) RETURN DISTINCT n.name AS modulename",sep="")
     modulenames<-cypher(graph,moduleq)
     modulenames<-modulenames$modulename[which(modulenames$modulename!="grey")]
     print(modulenames)
     modulenames<-c("All",modulenames)
     selectInput("moduleSelectSIG","Choose module SIG",modulenames,selected=1,multiple=FALSE)
   })
   
   output$selectModuleModcor<-renderUI({
      print("doing moduleModcor")
      reacstudyMC<-as.numeric(studyModcorfun())
      print("reacstudyMC")
      print(reacstudyMC)
      #reacsetMC<-as.numeric(squareModcorfun()[reacstudyMC])#debugging
      reacsetMC<-as.numeric(squareModcorfun())
      print("squareModcorfun()")
      print(squareModcorfun())
      print("reacsetMC")
      print(reacsetMC)
      print("reacsetMC")
      print(reacsetMC)
      reacedgeMC<-as.numeric(edgeModcorfun())
      print("reacedgeMC")
      print(reacedgeMC)
      
      studyMC<-names(studies)[reacstudyMC]
      print("studyMC")
      print(studyMC)
      
      setlistMC<-setlistnameslist[[reacstudyMC]]
      print("setlistMC")
      print(setlistMC)
      
      squareMC<-setlistMC[[reacsetMC]]
      print("squareMC")
      print(squareMC)
      moduleq<-paste("MATCH (n:wgcna {square:'",squareMC,"'}) RETURN DISTINCT n.name AS modulename",sep="")
      modulenames<-cypher(graph,moduleq)
      modulenames<-modulenames$modulename[which(modulenames$modulename!="grey")]
      print(modulenames)
      selectInput("moduleSelect","Choose module",modulenames,selected=1,multiple=FALSE)
      
   })
   
   
   #print("renderUI done")

#SESSION DATA OUTPUT####
   
   # # Values from cdata returned as text
   # output$clientdataText <- renderText({
   #    print("print(str(cdata))")
   #    #print(str(cdata))
   #    cnames <- names(cdata)
   #    
   #    allvalues <- lapply(cnames, function(name) {
   #       paste(name, cdata[[name]], sep = " = ")
   #    })
   #    
   #    res<-paste(allvalues, sep = ";;",collapse="\n")
   #    print(res)
   #    res
   # })
   # 
   # 
   # output$info <- renderText({
   #    info <- getCurrentOutputInfo()
   #    jsonlite::toJSON(
   #       list(
   #          bg = info$bg(),
   #          fg = info$fg(),
   #          accent = info$accent(),
   #          font = info$font()
   #       ),
   #       auto_unbox = TRUE
   #    )
   # })
   # 
#PLOT AND DATA OUTPUT####
   
   #Study design
   output$design <- renderImage({
     # When input$n is 3, filename is ./images/image3.jpeg
     filename <- file.path("/home/rstudio/source_data/design.pdf")
     
     # Return a list containing the filename and alt text
     list(src = filename)
     
   }, deleteFile = FALSE)
   
   #Phenotype####
   
   phenoqueryfun<-reactive({
     square<-squarephenofun()
     pheno<-phenofun2()
     query<-paste("MATCH (p:person)-[r]-(x) WHERE p.square =~ '",square,"' AND x.name ='",pheno,"' RETURN x.personName AS studyID, p.class1 AS class1, p.class2 AS class2, x.value as value",sep="")
     res<-cypher(graph,query)
     res$numeric<-as.numeric(res$value)
     res$factor<-as.factor(res$value)
     res$class1<-as.factor(res$class1)
     res$class2<-as.factor(res$class2)
     print(str(res))
     res
     res1<-res[which(res$class1==levels(res$class1)[1]),]
     class1avec<-res1$value
     res2<-res[which(res$class1==levels(res$class1)[2]),]
     class1bvec<-res2$value
     
     res3<-res[which(res$class2==levels(res$class2)[1]),]
     class2avec<-res3$value
     res4<-res[which(res$class2==levels(res$class2)[2]),]
     class2bvec<-res4$value
     
     print("debug splitting")
     print(class1avec)
     print(class1bvec)
     print(class2avec)
     print(class2bvec)
     print("class1avec==class1bvec")
     print(class1avec==class1bvec)
     print("sum(class1avec==class1bvec)")
     print(sum(class1avec==class1bvec))
     print("length(class1avec)")
     print(length(class1avec))
     if(sum(class1avec==class1bvec)==length(class1avec)){
        res<-res1
     }else if(sum(class2avec==class2bvec)==length(class2avec)){
        res<-res3
     }
     res
   })
   
   allphenofun<-reactive({
     #print(questions)
     square<-squarephenofun()
     #print(str(pData(eval(parse(text=square)))))
     print(paste("questions$",square,"$matrixPDname[[1]]",sep=""))
     filename<-eval(parse(text=(paste("questions$",square,"$matrixPDname[[1]]",sep=""))))
     mpd<-read.csv(file.path("/home/rstudio/source_data",filename),row.names=1,stringsAsFactors = FALSE)
     varclass<-mpd["varclass",]
     print(varclass)
     numerics<-which(varclass=="n")
     df<-pData(eval(parse(text=square)))
     df[,numerics]<-sapply(df[,numerics], as.numeric)
     df<-Filter(function(x) !(all(x=="")), df)
     
     DT <- as.data.table(df)
     DT[,which(unlist(lapply(DT, function(x)!all(is.na(x))))),with=F]
   })
   
   output$subjectnumbers<-metaRender(renderTable,{
     #square<-..(squarephenofun())
     ques<-eval(parse(text=paste("questions$",..(input$PhenoSquare),sep="")))
     out<-eval(parse(text=paste("table(",..(input$PhenoSquare),"$",ques$contrast_variable[[1]],",",..(input$PhenoSquare),"$",ques$contrast_variable[[2]],")",sep="")))
     out2<-as.data.frame.matrix(out)
     out2
     
     
   },striped=TRUE,rownames=TRUE,colnames=TRUE)
   
   output$allphenodata<-renderDataTable(allphenofun(),filter="top")
   
   
   table1fun<-reactive({
     #print(questions)
     square<-squarephenofun()
     #print(str(pData(eval(parse(text=square)))))
     print(paste("questions$",square,"$matrixPDname[[1]]",sep=""))
     filename<-eval(parse(text=(paste("questions$",square,"$matrixPDname[[1]]",sep=""))))
     mpd<-read.csv(file.path("/home/rstudio/source_data",filename),row.names=1,stringsAsFactors = FALSE)
     varclass<-mpd["varclass",]
     print(varclass)
     numerics<-which(varclass=="n")
     factors<-which(varclass=="c")
     df<-pData(eval(parse(text=square)))
     df[,numerics]<-sapply(df[,numerics], as.numeric)
     df[,factors]<-sapply(df[,factors], as.factor)
     
     #df<-Filter(function(x) !(all(x=="")), df)
     
     DT <- as.data.table(df)
     DT[,which(unlist(lapply(DT, function(x)!all(is.na(x))))),with=F]
     
     square<-squarephenofun()
     selectVSTs<-tab1varfun()
     print(square)
     print(selectVSTs)
     vstvec<-paste(selectVSTs,sep="','",collapse="','")
     print(vstvec)
     
     queryA<-paste("MATCH (n:wgcna {square:'",square,"',edge:'",1,"'}) RETURN DISTINCT n.contrastvar as cv1",sep="")
     print(queryA)
     res1<-cypher(graph,queryA)
     con1<-as.character(res1$cv1)
     
     queryB<-paste("MATCH (n:wgcna {square:'",square,"',edge:'",3,"'}) RETURN DISTINCT n.contrastvar as cv2",sep="")
     res2<-cypher(graph,queryB)
     con2<-as.character(res2$cv2)
     
     print(con1)
     print(con2)
     print(df[,1:12])
     queryC<-paste("MATCH (p:personPheno {square:'",square,"'})-[r]-(v1:varsubtype) WHERE (v1.name IN ['",vstvec,"']) RETURN DISTINCT p.name AS var",sep="")
     res3<-cypher(graph,queryC)
     vars<-as.character(res3$var)
     
     select<-c(con1,con2,vars)
     print("selecting table 1 variables and the two contrasts relevant to the square")
     print(select)
     print(colnames(df))
     datasubset<-df[,select]
     
    #  #query<-paste("MATCH (pr:person)-[r0]-(p:personPheno {square:'",square,"'})-[r]-(v1:varsubtype)-[r2]-(v2:vartype) WHERE (v1.name IN ['",vstvec,"'] AND p.personName =~ '.*?blood') RETURN p.personName AS name,p.name AS var,p.value AS value, v1.name AS vst, v2.name AS VT",sep="")
    #  query<-paste("MATCH (pr:person)-[r0]-(p:personPheno {square:'",square,"'})-[r]-(v1:varsubtype)-[r2]-(v2:vartype) WHERE (v1.name IN ['",vstvec,"']) RETURN p.personName AS name,p.name AS var,p.value AS value, v1.name AS vst, v2.name AS VT",sep="")
    #  
    #  print(query)
    #  res<-cypher(graph,query)
    #  print(res)
    #  df.wide <- pivot_wider(res[,1:4], names_from = var, values_from = value)
    #  
    #  print(paste("questions$",square,"$matrixPDname[[1]]",sep=""))
    #  filename<-eval(parse(text=(paste("questions$",square,"$matrixPDname[[1]]",sep=""))))
    #  mpd<-read.csv(file.path("/home/rstudio/source_data",filename),row.names=1,stringsAsFactors = FALSE)
    #  varclass<-mpd["varclass",]
    #  print(varclass)
    #  
    #  numerics<-which(varclass=="n")
    #  numericlist<-colnames(mpd)[numerics]
    #  print(numericlist)
    #  
    #  factors<-which(varclass=="c")
    #  faclist<-colnames(mpd)[factors]
    #  
    #  myNumerics<-colnames(df.wide)[which(colnames(df.wide)%in%numericlist)]
    #  myFactors<-colnames(df.wide)[which(colnames(df.wide)%in%faclist)]
    # 
    #  print(myNumerics)
    #  print(myFactors)
    #  
    #  df.wide[,myNumerics]<-lapply(df.wide[,myNumerics],as.numeric)
    #  
    #  df.wide[,myFactors]<-lapply(df.wide[,myFactors],as.factor)
    #  print(str(df.wide))
    #  
    df.wide<-datasubset
    print(df.wide)
    
    for (iterator in colnames(df.wide)) {
      print(iterator)
      expr<-paste("table1::label(df.wide$",as.character(iterator),")<-iterator",sep="")
      print(expr)
      eval(parse(text=expr))

      }
    #table1(df.wide,labels(df.wide))
    
    
    #classes
    byClass<-classesfun()
    if (byClass=="class1"){
      strat<-con1
    }else if (byClass=="class2"){
      strat<-con2
    }else{
      strat<-paste(con1,con2,sep="*")
    } 
    
    #t1expr<-paste("table1(~ ",paste(colnames(df.wide)[-c(1,2)],collapse=" + "), "| ", paste(con1,con2,sep=" + "), ", data=df.wide, overall=F, extra.col=list(`P-value`=pvalue),topclass='Rtable1-zebra')",sep="")
    print(df.wide)
    df.wide <- df.wide[,colSums(is.na(df.wide))<nrow(df.wide)]
    t1expr<-paste("table1(~ ",paste(colnames(df.wide)[-c(1,2)],collapse=" + "), "| ", strat, ", data=df.wide, overall=F, extra.col=list(`P-value`=pvalue),topclass='Rtable1-zebra')",sep="")
    
    
    print(t1expr)
    
    t1<-eval(parse(text=t1expr))
    t2<-t1kable(t1)
    print(str(t2))
    t2

    })
   
   output$table1<-renderText(table1fun())
   
   
   
   output$code<-renderPrint({
      expandChain(
         #Preparing the plot
         #output$phenoGridPlotPrepare(),
         #plot the plot
         #output$phenoGridPlot()
         output$subjectnumbers()
      )
   })
   
   
   tablefunNUM<-reactive({
      square<-squarephenofun()
      table1num<-data.frame("vartype"=character(),"varsubtype"=character(),"var"=character(),"All n"=integer(),"All median"=numeric(),"All IQR"=numeric(),
                            "grp1 n"=integer(),"grp1 median"=numeric(),"grp1 IQR"=numeric(),"grp2 n"=numeric(),"grp2 median"=numeric(),"grp2 IQR"=numeric(),"testname"=character(),"P-value"=numeric())
      
      table1cat<-data.frame("vartype"=character(),"varsubtype"=character(),"var"=character(),"All n"=integer(),"grp1"=character(),"grp2"=character(),"testname"=character(),"P-value"=numeric())
      
      count<-0
      count2<-0
      # 1. Get vartypes for square
      query<-paste("MATCH (pp:personPheno {square:'",square,"'})-[r]-(vst:varsubtype)-[r2]-(vt:vartype) RETURN DISTINCT(vt.name)",sep="")
      res<-cypher(graph,query)
      vartypes<-res[,1]
      
      # 2. For loop over vartypes
      for (vt in vartypes){
         # 3. Get varsubtypes for vartype
         print("vt:")
         print(vt)
         query<-paste("MATCH (pp:personPheno {square:'",square,"'})-[r]-(vst:varsubtype)-[r2]-(vt:vartype) WHERE vt.name = '",vt,"' RETURN DISTINCT(vst.name)",sep="")
         res2<-cypher(graph,query)
         varsubtypes<-res2[,1]
         
         # 4. For loop over varsubtypes
         for (vst in varsubtypes){
            print("vst:")
            print(vst)
            # 5. Get variables for varsubtype
            query<-paste("MATCH (pp:personPheno {square:'",square,"'})-[r]-(vst:varsubtype)-[r2]-(vt:vartype) WHERE vt.name = '",vt,"' AND vst.name = '",vst,"' RETURN DISTINCT(pp.name)",sep="")
            res3<-cypher(graph,query)
            vars<-res3[,1]
            print("vars")
            print(vars)
            # 6. For loop over variables
            for (var in vars){
               print("var:")
               print(var)
               #count<-count+1
               # 7. Get data for the variable
               query<-paste("MATCH (ps:person)-[r0]-(pp:personPheno {square:'",square,"'})-[r]-(vst:varsubtype)-[r2]-(vt:vartype) WHERE vt.name = '",vt,"' AND vst.name = '",vst,"' AND pp.name = '",var,"' RETURN 
                        pp.personName AS sampleID, ps.class1 AS class1, ps.class2 AS class2, pp.value AS value",sep="")
               data<-cypher(graph,query)
               
               print("DEBUG ON")
               print(square)
               print(vt)
               print(vst)
               print(var)
               print(head(data))
               
               print("DEBUG OFF")
               print("-------------------------------------------------------------------------------------------------")
               
               # 8. Check if the class1 subsets are identical. This is the case for matched samples (blood/fluid where the variable isn't different between blood and fluid, but not otherwise) select only data for one class1 type
               c1<-sort(data$value[which(data$class1==unique(data$class1)[1])])
               c2<-sort(data$value[which(data$class1==unique(data$class1)[2])])
               print("debug 0")
               print("c1")
               print(c1)
               print("c2")
               print(c2)
               
               needtosplit<-"no"
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
               print("needtosplit")
               print(needtosplit)
               #changeindex<-which(c1!=c2)
               # 9. Define variable class (numeric or character)
               datavec<-as.numeric(data$value)
               print("debug3")
               if (sum(is.na(datavec))==length(datavec)) {
                  data$value<-as.character(data$value)
               }else{
                  data$value<-as.numeric(data$value)
               }
               # 10. if numeric
               print("debug4")
               print("class(data$value)")
               print(class(data$value))
               print("needtosplit")
               print(needtosplit)
               
               print(class(data$value))
               if (class(data$value)=="numeric" & needtosplit == "no"){
                  print("debug5")
                  count<-count+1#remove NA
                  print(data)
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
                  stats<-"fail"
                  
                  try(stats<-kruskal.test(value~class2,data=data))
                  print("debug5.1 - stats output")
                  print(stats)
                  print(stats[1]!="fail")
                  print("class2names")
                  print(class2names)
                  print(length(class2names))
                  print("end debug 5.1")
                  
                  if(stats[1]!="fail"&length(class2names)==2){
                     table1num[count,]<-c(vt,vst,var,alln,allmed,alliqr,sum[1,2],sum[1,3],sum[1,4],sum[2,2],sum[2,3],sum[2,4],stats$method,stats$p.value)
                  }else if(length(class2names)==1) {
                     table1num[count,]<-c(vt,vst,var,alln,allmed,alliqr,sum[1,2],sum[1,3],sum[1,4],NA,NA,NA,"NA",NA)
                  }
                  
               }else if (class(data$value)=="numeric" & needtosplit == "yes"){
                  print("debug6")
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
                  if(stats[1]!="fail"){
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
                  if(stats[1]!="fail"){
                     table1num[count,]<-c(vt,vst,paste(var,splitclasses[2],sep="."),alln,allmed,alliqr,sum[1,2],sum[1,3],sum[1,4],sum[2,2],sum[2,3],sum[2,4],stats$method,stats$p.value)
                  }else{
                     table1num[count,]<-c(vt,vst,paste(var,splitclasses[2],sep="."),alln,allmed,alliqr,sum[1,2],sum[1,3],sum[1,4],sum[2,2],sum[2,3],sum[2,4],"NA",NA)
                  }
                  
                  print("hey!")
               }else if (class(data$value)=="character" & needtosplit == "no"){
                  print("new debug start")
                  #remove NA
                  data<-data[which(data$value!=""),]
                  
                  if (nrow(data)>0){
                     print("in nrow(data)>0 loop")
                     data$class2<-as.factor(data$class2)
                     alln<-nrow(data)
                     data$value<-as.factor(data$value)
                     ft<-ftable(class2~value,data=data)
                     #print(data)
                     #print(ft)
                     if(nrow(ft)>=2){
                        print("in nrow(ft) loop ")
                        ft[,1]<-ft[,1]/colSums(ft)[1]
                        ft[,2]<-ft[,2]/colSums(ft)[2]
                        statname="Chi-squared test"
                        
                        stats<-"fail"
                        
                        try(stats<-chisq.test(ft))
                        print("loop stats done")
                        ftdf<-as.data.frame((ft))
                        class2names<-as.character(levels(ftdf$class2))
                        ftdf1<-ftdf[which(ftdf$class2==class2names[1]),]
                        ftdf2<-ftdf[which(ftdf$class2==class2names[2]),]
                        
                        #grp1<-paste(paste(ftdf1$value,ftdf1$Freq,sep=":"),collapse="__")
                        #grp2<-paste(paste(ftdf2$value,ftdf2$Freq,sep=":"),collapse="__")
                        
                        grp1<-paste(paste(ftdf1$value,sprintf("%.2f",ftdf1$Freq),sep=":"),collapse="__")
                        grp2<-paste(paste(ftdf2$value,sprintf("%.2f",ftdf2$Freq),sep=":"),collapse="__")
                        
                        count2<-count2+1
                        if(stats[1]!="fail"&length(class2names)==2){
                           table1cat[count2,]<-c(vt,vst,var,alln,grp1,grp2,stats$method,stats$p.value)
                        }else if(length(class2names)==1){
                           table1cat[count2,]<-c(vt,vst,var,alln,grp1,grp2,"NA",NA)
                        }
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
      table1num
      
   })
   
   tablefunCAT<-reactive({
      square<-squarephenofun()
      table1num<-data.frame("vartype"=character(),"varsubtype"=character(),"var"=character(),"All n"=integer(),"All median"=numeric(),"All IQR"=numeric(),
                            "grp1 n"=integer(),"grp1 median"=numeric(),"grp1 IQR"=numeric(),"grp2 n"=numeric(),"grp2 median"=numeric(),"grp2 IQR"=numeric(),"testname"=character(),"P-value"=numeric())
      
      table1cat<-data.frame("vartype"=character(),"varsubtype"=character(),"var"=character(),"All n"=integer(),"grp1"=character(),"grp2"=character(),"testname"=character(),"P-value"=numeric())
      
      count<-0
      count2<-0
      # 1. Get vartypes for square
      query<-paste("MATCH (pp:personPheno {square:'",square,"'})-[r]-(vst:varsubtype)-[r2]-(vt:vartype) RETURN DISTINCT(vt.name)",sep="")
      res<-cypher(graph,query)
      vartypes<-res[,1]
      
      # 2. For loop over vartypes
      for (vt in vartypes){
         # 3. Get varsubtypes for vartype
         
         query<-paste("MATCH (pp:personPheno {square:'",square,"'})-[r]-(vst:varsubtype)-[r2]-(vt:vartype) WHERE vt.name = '",vt,"' RETURN DISTINCT(vst.name)",sep="")
         res2<-cypher(graph,query)
         varsubtypes<-res2[,1]
         
         # 4. For loop over varsubtypes
         for (vst in varsubtypes){
            # 5. Get variables for varsubtype
            query<-paste("MATCH (pp:personPheno {square:'",square,"'})-[r]-(vst:varsubtype)-[r2]-(vt:vartype) WHERE vt.name = '",vt,"' AND vst.name = '",vst,"' RETURN DISTINCT(pp.name)",sep="")
            res3<-cypher(graph,query)
            vars<-res3[,1]
            print("vars")
            print(vars)
            # 6. For loop over variables
            for (var in vars){
               #count<-count+1
               # 7. Get data for the variable
               query<-paste("MATCH (ps:person)-[r0]-(pp:personPheno {square:'",square,"'})-[r]-(vst:varsubtype)-[r2]-(vt:vartype) WHERE vt.name = '",vt,"' AND vst.name = '",vst,"' AND pp.name = '",var,"' RETURN 
                        pp.personName AS sampleID, ps.class1 AS class1, ps.class2 AS class2, pp.value AS value",sep="")
               data<-cypher(graph,query)
               
               print("DEBUG ON")
               print(square)
               print(vt)
               print(vst)
               print(var)
               print(head(data))
               
               print("DEBUG OFF")
               print("xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx")
               
               # 8. Check if the class1 subsets are identical. This is the case for matched samples (blood/fluid where the variable isn't different between blood and fluid, but not otherwise) select only data for one class1 type
               c1<-sort(data$value[which(data$class1==unique(data$class1)[1])])
               c2<-sort(data$value[which(data$class1==unique(data$class1)[2])])
               print("debug 0")
               print("c1")
               print(c1)
               print("c2")
               print(c2)
               
               needtosplit<-"no"
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
               print("needtosplit")
               print(needtosplit)
               #changeindex<-which(c1!=c2)
               # 9. Define variable class (numeric or character)
               datavec<-as.numeric(data$value)
               print("debug3")
               if (sum(is.na(datavec))==length(datavec)) {
                  data$value<-as.character(data$value)
               }else{
                  data$value<-as.numeric(data$value)
               }
               # 10. if numeric
               print("debug4")
               
               print("class(data$value)")
               print(class(data$value))
               print("needtosplit")
               print(needtosplit)
               
               if (class(data$value)=="numeric" & needtosplit == "no"){
                  print("debug5")
                  count<-count+1#remove NA
                  print(data)
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
                  stats<-"fail"
                  
                  try(stats<-kruskal.test(value~class2,data=data))#EDIT1
                  
                  if(stats[1]!="fail"&length(class2names)==2){
                     table1num[count,]<-c(vt,vst,var,alln,allmed,alliqr,sum[1,2],sum[1,3],sum[1,4],sum[2,2],sum[2,3],sum[2,4],stats$method,stats$p.value)
                  }else if(length(class2names)==1) {
                     table1num[count,]<-c(vt,vst,var,alln,allmed,alliqr,sum[1,2],sum[1,3],sum[1,4],NA,NA,NA,"NA",NA)
                  }
                  
               }else if (class(data$value)=="numeric" & needtosplit == "yes"){
                  print("debug6")
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
                  if(stats[1]!="fail"){
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
                  if(stats[1]!="fail"){
                     table1num[count,]<-c(vt,vst,paste(var,splitclasses[2],sep="."),alln,allmed,alliqr,sum[1,2],sum[1,3],sum[1,4],sum[2,2],sum[2,3],sum[2,4],stats$method,stats$p.value)
                  }else{
                     table1num[count,]<-c(vt,vst,paste(var,splitclasses[2],sep="."),alln,allmed,alliqr,sum[1,2],sum[1,3],sum[1,4],sum[2,2],sum[2,3],sum[2,4],"NA",NA)
                  }
                  
                  print("hey!")
               }else if (class(data$value)=="character" & needtosplit == "no"){
                  
                  #remove NA
                  data<-data[which(data$value!=""),]
                  
                  if (nrow(data)>0){
                     data$class2<-as.factor(data$class2)
                     alln<-nrow(data)
                     data$value<-as.factor(data$value)
                     ft<-ftable(class2~value,data=data)
                     if(nrow(ft)>=2){
                        ft[,1]<-ft[,1]/colSums(ft)[1]
                        ft[,2]<-ft[,2]/colSums(ft)[2]
                        statname="Chi-squared test"
                        
                        stats<-"fail"
                        
                        try(stats<-chisq.test(ft))
                        
                        ftdf<-as.data.frame((ft))
                        class2names<-as.character(levels(ftdf$class2))
                        ftdf1<-ftdf[which(ftdf$class2==class2names[1]),]
                        ftdf2<-ftdf[which(ftdf$class2==class2names[2]),]
                        
                        #grp1<-paste(paste(ftdf1$value,ftdf1$Freq,sep=":"),collapse="__")
                        #grp2<-paste(paste(ftdf2$value,ftdf2$Freq,sep=":"),collapse="__")
                        
                        grp1<-paste(paste(ftdf1$value,sprintf("%.2f",ftdf1$Freq),sep=":"),collapse="__")
                        grp2<-paste(paste(ftdf2$value,sprintf("%.2f",ftdf2$Freq),sep=":"),collapse="__")
                        
                        count2<-count2+1
                        if(stats[1]!="fail"&length(class2names)==2){
                           table1cat[count2,]<-c(vt,vst,var,alln,grp1,grp2,stats$method,stats$p.value)
                        }else if(length(class2names)==1){
                           table1cat[count2,]<-c(vt,vst,var,alln,grp1,grp2,"NA",NA)
                        }
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
      table1cat
      
   })
   
   output$tableNUM<-renderDataTable(tablefunNUM(),filter="top")
   
   output$tableNUM2<-renderTable({
     thetab<-tablefunNUM()
     xtable(thetab,digits=3)
     },striped=TRUE,bordered=TRUE)
   
   output$tableCAT<-renderDataTable(tablefunCAT(),filter="top")
   
   output$tableCAT2<-renderTable({
     thetab<-tablefunCAT()
     xtable(thetab,digits=3)
   },striped=TRUE,bordered=TRUE)
   
   # Downloadable csv of tableNUM
   output$tablenumcsv <- downloadHandler(
      filename = "tablenum.csv",
      content = function(file) {
         write.csv(tablefunNUM(), file, row.names = TRUE)
      }
   )
   # Downloadable csv of tableCAT
   output$tablecatcsv <- downloadHandler(
      filename = "tablecat.csv",
      content = function(file) {
         write.csv(tablefunCAT(), file, row.names = TRUE)
      }
   )
   
   output$phenodatatable<-renderDataTable(phenoqueryfun(),filter="top")
   
   output$phenoPlotPrepare<-renderPlot({
     data<-phenoqueryfun()
     print(data$numeric)
     print(is.na(data$numeric))
     
     
     ##here now
     if(length(unique(data$class1))==1){
        data$Categories<-data$class2
        
     }else if(length(unique(data$class2))==1){
        data$Categories<-data$class1
     }else{
        byClass<-classesfun()
        if (byClass=="class1"){
           data$Categories<-data$class1
        }else if (byClass=="class2"){
           data$Categories<-data$class2
        }else{
           data$Categories<-interaction(data$class1,data$class2)
        } 
        
     }
     
     
     
     
     
     
     if(sum(is.na(data$numeric))<length(data$numeric)){
       
       if(length(unique(data$class1))==2){
          print("numstats")
          numstat.kwt<-kruskal.test(numeric~class1,data=data)
          print("KWT:")
          print(numstat.kwt)
          #print(str(numstat.kwt))
          #print(data.frame(unlist(numstat.kwt,use.names = T)))
          KWT_class1<-data.frame(unlist(numstat.kwt,use.names = T))
          
          sol<-ddply(data,~class1,summarise,mean=mean(numeric,na.rm=TRUE),sd=sd(numeric,na.rm=TRUE),median=median(numeric,na.rm=TRUE),IQR=IQR(numeric,na.rm=TRUE))
          print(sol)
       }else{
          sol<-"not done"
          print("sol not done")
          
       }
       
        if(length(unique(data$class2))==2){
           print("numstats2")
           numstat.kwt2<-kruskal.test(numeric~class2,data=data)
           print("KWT2:")
           print(numstat.kwt2)
           KWT_class2<-data.frame(unlist(numstat.kwt2,use.names = T))
           
           sol2<-ddply(data,~class2,summarise,mean=mean(numeric,na.rm=TRUE),sd=sd(numeric,na.rm=TRUE),median=median(numeric,na.rm=TRUE),IQR=IQR(numeric,na.rm=TRUE))
           print(sol2)
        }else{
           sol2<-"not done"
           print("sol2 not done")
        }
        
       
       if(sol=="not done"){
          
          print(KWT_class2)
          KWT<-KWT_class2
          colnames(KWT)<-c("class2")
          colnames(sol2)[1]<-"Class"
          solcom<-sol2
          
       }else if(sol2=="not done"){
          print(KWT_class1)
          KWT<-KWT_class1
          colnames(KWT)<-c("class1")
          colnames(sol)[1]<-"Class"
          solcom<-sol
          
       }else{
          print(KWT_class1)
          print(KWT_class2)
          KWT<-cbind(KWT_class1,KWT_class2)
          colnames(KWT)<-c("class1","class2")
          
          colnames(sol)[1]<-"Class"
          colnames(sol2)[1]<-"Class"
          solcom<-rbind(sol,sol2)
       }
       
       
       
       output$summary<-renderTable(
         solcom,
         rownames=TRUE
         
       )
       
       output$stats<-renderTable(
         KWT,
         rownames=TRUE,
         digits=3
         #xtable(as.data.frame(numstat.kwt))
         
       )
       
       
       complevs<-eval(parse(text=paste("levels(as.factor(data$","Categories","))",sep="")))
       if (length(complevs)==4){
          my_comparisons <- list( c(complevs[1], complevs[2]),c(complevs[1], complevs[3]),c(complevs[2], complevs[4]),c(complevs[3], complevs[4]) )
       }else{
          my_comparisons <- list( c(complevs[1], complevs[2]) )
       }
       print("print(my_comparisons)")
       print(my_comparisons)
       
       # Basic dot plot
       # plotfig<-ggplot(data, aes(x=Categories, y=numeric, fill=Categories)) + 
       #   geom_boxplot(fill='seashell')+
       plotfig<-ggplot(data, aes(x=Categories, y=numeric, fill=Categories)) + 
         geom_boxplot(aes(alpha=0.3 )) +
         scale_fill_manual(values=c("lightsteelblue", "mistyrose", "moccasin","lightgreen")) +
         geom_dotplot(mapping=aes(x=Categories, y=numeric, color=Categories, fill = Categories),binaxis='y', stackdir='center',stackratio=1.2,method='dotdensity',inherit.aes = FALSE) +
         scale_color_manual(values=c("darkblue", "darkred", "orange","darkgreen")) +
         #geom_jitter() +
         #stat_summary(fun.y=median, geom="point", shape=18, color="red") +
         theme_minimal() +
         theme(axis.title.y = element_text(size = rel(2.5), angle = 90, color = "slategrey",margin = margin(t = 0, r = 50, b = 50, l = 50))) +
         theme(axis.title.x = element_text(size = rel(2.5), color = "slategrey",margin = margin(t = 50, r = 50, b = 0, l = 50))) + 
         theme(axis.text = element_text(size = rel(1.5))) +
         theme(legend.text = element_text(size = rel(1.2))) +
         theme(legend.title = element_text(size = rel(1.2))) +
         theme(legend.text=element_text(size=rel(1.4))) +
         theme(legend.key.size = unit(1,"cm")) +
         labs(y = phenofun2()) +
         ggtitle(paste(phenofun2(),"by class")) +
         theme(plot.title = element_text(size=rel(2.5),hjust = 0.5,color="maroon")) +
         stat_compare_means(comparisons = my_comparisons,method="wilcox.test",label.y.npc="center") + # Add pairwise comparisons p-value
         #stat_compare_means(comparisons = my_comparisons,method="kruskal.test",label.y.npc="center") + # Add pairwise comparisons p-value
         #stat_compare_means(comparisons = my_comparisons,method="anova",label.y.npc="center") + # Add pairwise comparisons p-value
         #stat_compare_means(label.y = 50) +     # Add global p-value
         scale_y_continuous(expand=expand_scale(mult=c(0,0.2)))
       
     }else{
       
       # ggplot(data=pd3cat,aes_string(x=names(pd3cat)[colu],fill=names(pd3cat)[colu]))+
       #   geom_bar(aes(y = (..count..)/sum(..count..)))+
       #   facet_grid(eval(parse(text=paste(contrast.variable[[1]],"~",contrast.variable[[2]],sep=""))))+
       #   theme_bw()+
       #   theme(axis.text.x = element_text(angle = 45, hjust = 1))+
       #   scale_x_discrete(name=paste("Proportion by",contrast.variable2,"\nP=",format(testres_c2$p.value,digits=3)))+
       #   scale_y_continuous(name=paste("Proportion by",contrast.variable1,"\nP=",format(testres_c1$p.value,digits=3)))+
       #   ggtitle(names(pd3cat)[[colu]])
       
       print("do newdat :)")
       newdat<-ftable(class1~value,data=data)
       print("newdat")
       print(class(newdat))
       print(str(newdat))
       print(newdat)
       testres<-chisq.test(newdat)
       print(testres)
       
       print("do newdat2 :)")
       newdat2<-ftable(class2~value,data=data)
       print("newdat2")
       print(newdat2)
       testres2<-chisq.test(newdat2)
       print(testres2)
       
       
       
       sumtab<-cbind(as.data.frame(newdat),as.data.frame(newdat2))
       
       print("dat frame conversion:")
       chisq1<-data.frame(unlist(testres,use.names = T))
       print(chisq1)
       chisq2<-data.frame(unlist(testres2,use.names = T))
       print(chisq2)
       chisqtab<-cbind(chisq1,chisq2)
       colnames(chisqtab)<-c("Class1","Class2")
       
       
       
       output$summary<-renderTable(
         sumtab,
         rownames=TRUE
         
       )
       
       output$stats<-renderTable(
         chisqtab[1:5,],
         rownames=TRUE,
         digits=3
         #xtable(as.data.frame(numstat.kwt))
         
       )
       
       #plotfig<-ggplot(data, aes(Categories, fill=factor)) + 
       plotfig<-ggplot(data, aes(value, fill=factor)) +   
         #geom_bar(aes(y = (..count..)/sum(..count..))) +
         geom_bar() +
         #scale_y_continuous(labels=percent) +
         facet_grid(class1~class2) +
         #geom_bar() +
         #stat_summary(fun.y=median, geom="point", shape=18, color="red") +
         #theme_minimal() +
         theme_light() +
         theme(axis.title.y = element_text(size = rel(2.5), angle = 90, color = "slategrey",margin = margin(t = 0, r = 50, b = 50, l = 50))) +
         theme(axis.title.x = element_text(size = rel(2.5), color = "slategrey",margin = margin(t = 50, r = 50, b = 0, l = 50))) + 
         theme(axis.text = element_text(size = rel(1.5))) +
         theme(strip.text = element_text(size = 14)) +
         theme(legend.text = element_text(size = rel(1.2))) +
         theme(legend.title = element_text(size = rel(1.2))) +
         theme(legend.text=element_text(size=rel(1.4))) +
         theme(legend.key.size = unit(1,"cm")) +
         labs(y = phenofun2()) +
         ggtitle(paste(phenofun2(),"by class")) +
         theme(plot.title = element_text(size=rel(2.5),hjust = 0.5,color="maroon"))
     }
     
     #ggplotly(plotfig,tooltip="studyID") 
     plotfig
   })
   
   output$phenoPlot <- renderUI({
     plotOutput("phenoPlotPrepare", height = phenoplotheightfun())
   })

   #############phenogridplot####
   output$phenoGridPlotPrepare<-renderPlot({
      # 1. select the square
      print("phenogrid start::###########")
      
      square<-squarephenofun()
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
      
      myData<-aggregate(data$numeric,by=list(cats=data$Categories),FUN=function(x){c(median=median(x,na.rm=T),iqr=IQR(x,na.rm=T),n=length(x))})
      myData<-do.call(data.frame,myData)
      colnames(myData)<-c("Categories","median","iqr","n")
      
      dodge <- position_dodge(width = 0.9)
       limits <- aes(ymax = myData$median + myData$iqr,
                     ymin = myData$median - myData$iqr)
      # p <- ggplot(data = myData, aes(x = names, y = mean, fill = names))
      # p + geom_bar(stat = "identity", position = dodge) +
      #    geom_errorbar(limits, position = dodge, width = 0.25) +
      
      plotfig<-ggplot(data, aes(x=eval(parse(text=fillvar)), y=numeric, fill=eval(parse(text=fillvar)))) + 
      #plotfig<-ggplot(myData, aes(x=eval(parse(text=fillvar)), y=median, fill=eval(parse(text=fillvar)))) + 
         geom_boxplot(aes(alpha=0.3 ),width=.25) +
         #geom_bar(stat = "identity") +
         #geom_col(stat = "median") +
         #geom_errorbar(limits,position=position_dodge(0.9),width = 0.25) +
         #geom_boxplot(width=.25) +
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
      
   })
   
   output$phenoGridPlot <- renderUI({
      plotOutput("phenoGridPlotPrepare", height = phenoplotheightfun())
   })
   #############phenogridplot2####
   output$phenoGridPlotPrepare2<-renderPlot({
      # 1. select the square
      print("phenogrid start::###########")
      
      square<-squarephenofun()
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
      plotfig<-ggplot(data, aes(x=eval(parse(text=fillvar)), fill=character)) + 
      #plotfig<-ggplot(data, aes(x=character, fill=eval(parse(text=fillvar)))) +   
         #geom_bar(aes(y = (..count..)/sum(..count..))) +
         geom_bar() +
         #scale_y_continuous(labels=percent) +
         #facet_grid(class1~class2) +
         #geom_bar() +
         #stat_summary(fun.y=median, geom="point", shape=18, color="red") +
         #theme_minimal() +
         theme_light() +
         theme(axis.title.y = element_text(size = rel(2.5), angle = 90, color = "slategrey",margin = margin(t = 0, r = 50, b = 50, l = 50))) +
         theme(axis.title.x = element_text(size = rel(2.5), color = "slategrey",margin = margin(t = 50, r = 50, b = 0, l = 50))) + 
         theme(axis.text = element_text(size = rel(1.5))) +
         theme(strip.text = element_text(size = 14)) +
         theme(legend.text = element_text(size = rel(1.2))) +
         theme(legend.title = element_text(size = rel(1.2))) +
         theme(legend.text=element_text(size=rel(1.4))) +
         theme(legend.key.size = unit(1,"cm")) +
         labs(y = phenofun2()) +
         ggtitle(paste(phenofun2(),"by class")) +
         theme(plot.title = element_text(size=rel(2.5),hjust = 0.5,color="maroon"))
      
      
      
      plotfig<-plotfig + facet_wrap(vars(pheno),scales="free",ncol=3)
      
      print(plotfig)
      
   })
   
   output$phenoGridPlot2 <- 
      renderUI({
         plotOutput("phenoGridPlotPrepare2", height = phenoplotheightfun())
      })
   
   #############cellpropgrid####
   output$cellpropGridPlotPrepare<-renderPlot({
      # 1. select the square
      print("phenogrid start::###########")
      
      square<-squarephenofun()
      print("square:")
      print(square)
      
      # 4. get the data for the group and square
      #select fill variable
      fillvar<-fillvarfun()
      # fillvar<-"Categories"
      # fillvar<-"class1"
      # fillvar<-"class2"
      
   
      query<-paste("MATCH (p:person)-[r0]-(pc:personCell {square:'",square,"'}) RETURN DISTINCT pc.personName as person, p.class1 as class1, p.class2 as class2, pc.name as cell, toFloat(pc.value) as value",sep="")
       
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
         geom_boxplot(aes(alpha=0.3 ),width=.25) +
         #geom_boxplot(width=.25) +
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
         labs(y = "value") +
         labs(x = "groupings") +
         ggtitle(paste("Celltype proportions by",fillvar)) +
         theme(plot.title = element_text(size=rel(2.5),hjust = 0.5,color="maroon")) +
         theme(strip.text.x = element_text(size = 12, colour = "blue", angle = 0)) +
         stat_compare_means(comparisons = my_comparisons,method="wilcox.test",label.y.npc="center") + # Add pairwise comparisons p-value
         #stat_compare_means(label.y = 50) +     # Add global p-value
         scale_y_continuous(expand=expand_scale(mult=c(0,0.2)))
      
      
      
      plotfig<-plotfig + facet_wrap(vars(cell),scales="free",ncol=2)
      
      print(plotfig)
      
   })
   
   output$cellpropGridPlot <- renderUI({
      plotOutput("cellpropGridPlotPrepare", height = phenoplotheightfun())
   })
   #####wgcnaMEgridplot####
   output$wgcnaMEGridPlotPrepare<-renderPlot({
      # 1. select the square
      print("phenogrid start::###########")
      
      square<-squarephenofun()
      print("square:")
      print(square)
      
      # 4. get the data for the group and square
      #select fill variable
      fillvar<-fillvarfun()
      # fillvar<-"Categories"
      # fillvar<-"class1"
      # fillvar<-"class2"
      
      
      query<-paste("MATCH (p:person)-[r0]-(pme:personME {square:'",square,"'}) RETURN DISTINCT pme.personName as person, p.class1 as class1, p.class2 as class2, pme.name as ME, toFloat(pme.value) as value",sep="")
      
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
         geom_boxplot(aes(alpha=0.3 ),width=.25) +
         #geom_boxplot(width=.25) +
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
         labs(y = "value") +
         labs(x = "groupings") +
         ggtitle(paste("Module eigengene by",fillvar)) +
         theme(plot.title = element_text(size=rel(2.5),hjust = 0.5,color="maroon")) +
         theme(strip.text.x = element_text(size = 12, colour = "blue", angle = 0)) +
         stat_compare_means(comparisons = my_comparisons,method="wilcox.test",label.y.npc="center") + # Add pairwise comparisons p-value
         #stat_compare_means(label.y = 50) +     # Add global p-value
         scale_y_continuous(expand=expand_scale(mult=c(0,0.2)))
      
      
      
      plotfig<-plotfig + facet_wrap(vars(ME),scales="free",ncol=3)
      
      print(plotfig)
      
   })
   
   output$wgcnaMEGridPlot <- renderUI({
      plotOutput("wgcnaMEGridPlotPrepare", height = phenoplotheightfun())
   })
   
   
   ############
   
   
   
   # output$code<-renderPrint({
   #    expandChain(
   #       #Preparing the plot
   #    output$phenoGridPlotPrepare(),
   #    #plot the plot
   #    output$phenoGridPlot()
   #    )
   # })
   
      
   observe({output$genelist<-renderTable({
      regexgl<-useregex()
      print(regexgl)
      genelistquery<-paste("MATCH (s:SYMBOL) WHERE s.name=~ '",regexgl,"' RETURN s.name as gene",sep="")
      genelist<-cypher(graph,genelistquery)
      genelist<-genelist[order(genelist),]
      print(genelist)
      genelist},priority=100)
   })
   
   #Gene####
   # observe({output$boxplot<-renderPlot({
   #   ##CODE
   #   squarelist<-unlist(strsplit(setfunBXP(),"_"))
   #   if(length(squarelist)==3){
   #     square<-squarelist[3]
   #   }else{square<-paste(squarelist[c(3,4)],collapse="_")}
   #   regex<-useregex()
   #   print(paste("regex:",regex))
   #   print("command")
   #   #print(paste("probe_boxplot4('",square,"','",regex,"',miny=",as.numeric(yminfun()),",maxy=",as.numeric(ymaxfun()),")",sep=""))#debug
   #   #eval(parse(text=paste("probe_boxplot4('",square,"','",regex,"',miny=",as.numeric(yminfun()),",maxy=",as.numeric(ymaxfun()),",orderP=",orderPfun(),")",sep="")))#debug
   #   eval(parse(text=paste("probe_boxplot5('",square,"','",regex,"',miny=",as.numeric(yminfun()),",maxy=",as.numeric(ymaxfun()),",orderP=",orderPfun(),")",sep="")))#debug
   #   if(input$returnpdf==TRUE){
   #     pdf("plot.pdf",width=16,height=as.numeric(input$plotheight)/100,onefile=FALSE)
   #     eval(parse(text=paste("probe_boxplot5('",square,"','",regex,"',miny=",as.numeric(yminfun()),",maxy=",as.numeric(ymaxfun()),",orderP=",orderPfun(),")",sep="")))
   #     dev.off()
   #   }
   #   
   #   dequery<-paste("MATCH (p:PROBE {square:'",square,"'})-[r]-(s:SYMBOL) WHERE s.name=~ '",regex,"' RETURN p.square as square,p.edge as edge, p.name as probe, s.name as gene, p.aveEXPR as expression, p.logfc as logfc, p.adjPVAL as qvalue",sep="")
   #   detable<-cypher(graph,dequery)
   #   print(detable)
   #   
   #   # genelistquery<-paste("MATCH (s:SYMBOL) WHERE s.name=~ '",regex,"' RETURN s.name as gene",sep="")
   #   # genelist<-cypher(graph,genelistquery)
   #   # print(genelist)
   #   # output$genelist<-renderTable(genelist)
   #   print(returncsvfun2())
   #   output$detableBXP<-DT::renderDataTable(detable,server = FALSE,options=list(lengthMenu = c(5, 10, 15,20,50), pageLength = 15),filter="top") 
   #   observe({if(returncsvfun2()){
   #     print("writing table")
   #     write.csv(detable,file.path(tabledir,paste("boxplotDEtable",Sys.time(),".csv",sep="")))}#
   #   })#end observe
   # }, height = as.numeric(input$plotheight))
   # 
   # },priority=80)
   
   observe({output$boxplot<-renderPlot({
      ##CODE
      fillvariable<-fillvarfun()
      squarelist<-unlist(strsplit(setfunBXP(),"_"))
      if(length(squarelist)==3){
         square<-squarelist[3]
      }else{square<-paste(squarelist[c(3,4)],collapse="_")}
      regex<-useregex()
      print(paste("regex:",regex))
      
      
      
      
      
      print("command")
      #print(paste("probe_boxplot4('",square,"','",regex,"',miny=",as.numeric(yminfun()),",maxy=",as.numeric(ymaxfun()),")",sep=""))#debug
      #eval(parse(text=paste("probe_boxplot4('",square,"','",regex,"',miny=",as.numeric(yminfun()),",maxy=",as.numeric(ymaxfun()),",orderP=",orderPfun(),")",sep="")))#debug
      eval(parse(text=paste("probe_boxplot5('",square,"','",regex,"',miny=",as.numeric(yminfun()),",maxy=",as.numeric(ymaxfun()),",orderP=",orderPfun(),")",sep="")))#debug
      #eval(parse(text=paste("probe_boxplot6('",square,"','",regex,"',miny=",as.numeric(yminfun()),",maxy=",as.numeric(ymaxfun()),",orderP=",orderPfun(),",fillvar='",fillvariable,"')",sep="")))#debug
      if(input$returnpdf==TRUE){
         pdf("plot.pdf",width=16,height=as.numeric(input$plotheight)/100,onefile=FALSE)
         eval(parse(text=paste("probe_boxplot5('",square,"','",regex,"',miny=",as.numeric(yminfun()),",maxy=",as.numeric(ymaxfun()),",orderP=",orderPfun(),")",sep="")))
         dev.off()
      }
      
      dequery<-paste("MATCH (p:PROBE {square:'",square,"'})-[r]-(s:SYMBOL) WHERE s.name=~ '",regex,"' RETURN p.square as square,p.edge as edge, p.name as probe, s.name as gene, p.aveEXPR as expression, p.logfc as logfc, p.adjPVAL as qvalue",sep="")
      detable<-cypher(graph,dequery)
      print(detable)
      
      # genelistquery<-paste("MATCH (s:SYMBOL) WHERE s.name=~ '",regex,"' RETURN s.name as gene",sep="")
      # genelist<-cypher(graph,genelistquery)
      # print(genelist)
      # output$genelist<-renderTable(genelist)
      print(returncsvfun2())
      output$detableBXP<-DT::renderDataTable(detable,server = FALSE,options=list(lengthMenu = c(5, 10, 15,20,50), pageLength = 15),filter="top") 
      observe({if(returncsvfun2()){
         print("writing table")
         write.csv(detable,file.path(tabledir,paste("boxplotDEtable",Sys.time(),".csv",sep="")))}#
      })#end observe
   }, height = as.numeric(input$plotheight))
   
   },priority=80)
   
   #Inflammasomes####
   plotInputMachines<-reactive({
     
     ##BEGIN INSERT
     print("begin machines")
     
     #query=paste("MATCH (s:SYMBOL)-[r]-(p:PROBE)-[r2]-(n:wgcna) WHERE p.square =~ '",labelsvec,"' AND p.edge IN [1,2,3,4] AND s.name =~ '",generegexlist[grxfun()],"' RETURN s.name AS gene, p.name as PROBE, p.square AS square, p.edge as edge, p.logfc as logfc,n.name AS wgcna",sep="")
     
      
      squarelist<-unlist(strsplit(setfunMACHINE(),"_"))
      print("debug-1")
      if(length(squarelist)==3){
         square<-squarelist[3]
      }else{square<-paste(squarelist[c(3,4)],collapse="_")}
      print(square)
      print("a square has been found")
     #query=paste("MATCH (s:SYMBOL)-[r]-(p:PROBE)-[r2]-(n:wgcna) WHERE n.square =~ '",square,"' AND p.edge IN [1,2,3,4] AND s.name =~ '",generegexlist[grxfun()],"' RETURN s.name AS gene, p.name as PROBE, p.square AS square, p.edge as edge, p.logfc as logfc,n.name AS wgcna",sep="")
     query=paste("MATCH (s:SYMBOL)-[r]-(p:PROBE)-[r2]-(n:wgcna) WHERE n.square =~ '",square,"' AND s.name =~ '",generegexlist[grxfun()],"' RETURN s.name AS gene, p.name as PROBE, p.square AS square, p.edge as edge, p.logfc as logfc,n.name AS wgcna",sep="")
     
     print("query made")
     print(Sys.time())
     
     res<-cypher(graph,query)
     print("cypher query executed made")
     print(Sys.time())
     
     res<-res[which(res$edge!='5'),]
     
     idmat<-unique(res[,c(1,2)])
     id<-as.character(apply(idmat,1,function(x){paste(x[1]," (",x[2],")",sep="")}))
     squaremat<-matrix(ncol=4,nrow=nrow(idmat),data=0)
     squaremat2<-cbind(idmat,squaremat)
     thisquare<-idmat
     print("results formatted")
     print(Sys.time())
     
     
     print(square)
     print("debug0")
     subset<-res[which(res$square==square),]
     subset$edge<-as.numeric(subset$edge)
     thismat<-squaremat2
     print("new debug 1")
     print("head(res)")
     print(head(res))
     print("subset")
     print(subset)
     for (row in 1:nrow(subset)){
       thismat[
          which(thismat$PROBE==subset[row,"PROBE"]),
          2+subset[row,"edge"]
          ]<-subset[row,"logfc"]
     }
     print("debug1")
     print("thismat")
     print(thismat)
     #individual squares
     #plotmat0<-data.frame(thismat[,3:ncol(thismat)])
     plotmat0<-data.frame(thismat[,3:ncol(thismat)])
     print("debug2")
     row.names(plotmat0)<-as.character(apply(thismat,1,function(x){paste(x[1]," (",x[2],")",sep="")}))
     print("debug3")
     genenames<-thismat[which(rowSums(plotmat0)!=0),][,1]
     print("debug4")
     probenames<-thismat[which(rowSums(plotmat0)!=0),][,2]
     print("debug5")
     plotmat0<-plotmat0[which(rowSums(plotmat0)!=0),]
     print("debug6")
     protFunction<-sapply(genenames,function(x){annot2[which(annot2$Symbol==x),2]})
     print("debug7")
     dict<-unique(subset[,c("gene","PROBE","wgcna")])
     print("debug8")
     wgcnaANNOT<-sapply(probenames,function(x){dict[which(dict$PROBE==x),3]})
     print("debug9")
     anndf<-data.frame(wgcnaANNOT,protFunction)
     print("debug10")
     rownames(anndf)<-row.names(plotmat0)
     print("debug11")
     edges<-paste("edge",1:4,sep="")
     colnames(plotmat0)<-paste("edge",1:4,sep="_")
     print("debug12")
     print("head(anndf)")
     print(head(anndf))
     anndf2<-data.frame(anndf[order(anndf$protFunction),])
     print("debug13")
     print(anndf2)
     plotmat0<-plotmat0[order(anndf$protFunction),]
     print("debug14")
     
     cur_dev <- dev.cur()
     
     heatplot<-pheatmap(plotmat0,cluster_rows = FALSE,cluster_cols = FALSE,scale = "none",col=blueWhiteRed(50),breaks=seq(-max(abs(plotmat0)),max(abs(plotmat0)),length.out=51),annotation_row = anndf2,main=paste(square,grxfun(),sep="_"))
     
     dev.set(cur_dev)
     print(heatplot)#pheatmap(plotmat0,cluster_rows = FALSE,cluster_cols = FALSE,scale = "none",col=blueWhiteRed(50),breaks=seq(-max(abs(plotmat)),max(abs(plotmat)),length.out=51),fontsize_row = 8,annotation_row = anndf2,main=paste(square,names(generegexlist)[grx],sep="_"),filename=file.path(figdir,paste(square,"_",names(generegexlist)[grx],".pdf",sep="")))
     
     #end individual squares
     
     #thisquare<-cbind(thisquare,thismat[,3:6])
     if(input$returnpdf==TRUE){
       pdf("plot.pdf",width=16,height=as.numeric(input$plotheightMACHINE)/100,onefile=FALSE)
       pheatmap(plotmat0,cluster_rows = FALSE,cluster_cols = FALSE,scale = "none",col=blueWhiteRed(50),breaks=seq(-max(abs(plotmat0)),max(abs(plotmat0)),length.out=51),annotation_row = anndf2,main=paste(square,grxfun(),sep="_"))
       dev.off()
     }
     #}
     
     
     # plotmat<-data.frame(thisquare[,3:ncol(thisquare)])
     # row.names(plotmat)<-id
     # genenames<-idmat[,1]
     # protFunction<-sapply(genenames,function(x){annot2[which(annot2$Symbol==x),2]})
     # anndf<-data.frame(protFunction)
     # rownames(anndf)<-row.names(plotmat)
     # edges<-paste("edge",1:4,sep="")
     # colnames(plotmat)<-as.vector(t(outer(labels, edges, paste, sep=".")))
     # anndf2<-data.frame(anndf[order(anndf$protFunction),])
     # plotmat<-plotmat[order(anndf$protFunction),]
     # pheatmap(plotmat,cluster_rows = FALSE,cluster_cols = FALSE,scale = "none",col=blueWhiteRed(50),breaks=seq(-max(abs(plotmat)),max(abs(plotmat)),length.out=51),fontsize_row = 8,gaps_col=seq(4,ncol(plotmat),4),annotation_row = anndf,main=names(generegexlist)[grx],filename=paste(file.path(figdir,names(generegexlist)[grx]),".pdf",sep=""))
     # pheatmap(plotmat,cluster_rows = FALSE,cluster_cols = FALSE,scale = "none",col=blueWhiteRed(50),breaks=seq(-max(abs(plotmat)),max(abs(plotmat)),length.out=51),fontsize_row = 8,gaps_col=seq(4,ncol(plotmat),4),annotation_row = anndf,main=names(generegexlist)[grx])
     # 
     # 
     # finalmat<-rbind(finalmat,plotmat)
     # nrowvector[grx+1]<-tail(nrowvector,n=1)+nrow(plotmat)
     
     
     
     ##END INSERT
   })
   
   output$inflammasomes<-renderPlot({
     plotInputMachines()
   },height=700)
   
   #Signatures####
   signatureQuery<-reactive({
     sliderval1<-input$volcanoFC
     sliderval2<-input$volcanoSIG
     sliderval3<-input$volcanoFClow
     sliderval4<-input$volcanoSIGlow
     reacstudySIG<-as.numeric(studySIGfun())
     print("reacstudySIG")
     print(reacstudySIG)
     #reacsetSIG<-as.numeric(squareModcorfun()[reacstudySIG])#debugging
     reacsetSIG<-as.numeric(squareSIGfun())
     print("squareSIGfun()")
     print(squareSIGfun())
     print("reacsetSIG")
     print(reacsetSIG)
     print("reacsetSIG")
     print(reacsetSIG)
     reacedgeSIG<-as.numeric(edgeSIGfun())
     print("reacedgeSIG")
     print(reacedgeSIG)
     
     
     
     studySIG<-names(studies)[reacstudySIG]
     print("studySIG")
     print(studySIG)
     
     setlistSIG<-setlistnameslist[[reacstudySIG]]
     print("setlistSIG")
     print(setlistSIG)
     
     squareSIG<-setlistSIG[[reacsetSIG]]
     print("squareSIG")
     print(squareSIG)
     
     sigsub<-whichSIGfun()
     if(sigsub=="up"){
       modifier<-" WHERE p.logfc > 0"
     }else if (sigsub=="down"){
       modifier<-" WHERE p.logfc < 0"
     }else{
       modifier<-""
     }
     
     print("debug new sig 1")
     siglength<-as.numeric(defSigfun())
     query<-paste("MATCH (s:SYMBOL)-[r]-(p:PROBE {square:'",squareSIG,"',edge:'",reacedgeSIG,"'})-[r2]-(n:wgcna) ",modifier," RETURN s.name as SYMBOL, p.name AS nuID, p.logfc as logfc, p.adjPVAL, n.name as module ORDER BY p.adjPVAL",sep="")
     
     sigresult<-cypher(graph,query)
     sigresult<-sigresult[!duplicated(sigresult$nuID), ]
     
     sigresult.trunc<-sigresult[1:siglength,]
     print(head(sigresult.trunc))
     sigresult.trunc
     
     
   })
   
   ##
   
   ##
   
   enrichfun<-reactive({
     sliderval1<-input$volcanoFC
     sliderval2<-input$volcanoSIG
     sliderval3<-input$volcanoFClow
     sliderval4<-input$volcanoSIGlow
     print("whichSIGfun()")
     print(whichSIGfun())
     
     sigresult.trunc<-signatureQuery()
     
     if(whichSIGfun()!="all"){
       sigresult.trunc<-subset(sigresult.trunc,logfc > sliderval3 & logfc < sliderval1 & -log10(p.adjPVAL) > sliderval4 & -log10(p.adjPVAL) < sliderval2 )
     }else{
       sigresult.trunc<-subset(sigresult.trunc, (logfc > abs(sliderval3) | logfc < -abs(sliderval3)) &
                                 (logfc < abs(sliderval1) & logfc > -abs(sliderval1)) & 
                                 -log10(p.adjPVAL) > sliderval4 &
                                 -log10(p.adjPVAL) < sliderval2 
                               
       )
     }
     if(whichModuleFun()!="All"){
        sigresult.trunc<-sigresult.trunc[which(sigresult.trunc$module==whichModuleFun()),]
     }
     
     #sigresult.trunc<-subset(sigresult.trunc, logfc < sliderval1 & -log10(p.adjPVAL) < sliderval2 & logfc > sliderval3 & -log10(p.adjPVAL) > sliderval4)
     head(sigresult.trunc)
     nuIDlist<-sigresult.trunc$nuID
     entrezlist<-nuID2EntrezID(nuIDlist,lib.mapping='lumiHumanIDMapping')
     print(entrezlist)
     enrich.result<-enrichPathway(entrezlist,readable=T)
     enrich.result
   })
   
   gseafun<-reactive({
     sliderval1<-input$volcanoFC
     sliderval2<-input$volcanoSIG
     sliderval3<-input$volcanoFClow
     sliderval4<-input$volcanoSIGlow
     sigresult.trunc<-signatureQuery()
     
     if(whichSIGfun()!="all"){
       sigresult.trunc<-subset(sigresult.trunc,logfc > sliderval3 & logfc < sliderval1 & -log10(p.adjPVAL) > sliderval4 & -log10(p.adjPVAL) < sliderval2 )
     }else{
       sigresult.trunc<-subset(sigresult.trunc, (logfc > abs(sliderval3) | logfc < -abs(sliderval3)) &
                      (logfc < abs(sliderval1) & logfc > -abs(sliderval1)) & 
                      -log10(p.adjPVAL) > sliderval4 &
                      -log10(p.adjPVAL) < sliderval2 
                    
       )
     }
     
     if(whichModuleFun()!="All"){
        sigresult.trunc<-sigresult.trunc[which(sigresult.trunc$module==whichModuleFun()),]
     }
     
     
     #sigresult.trunc<-subset(sigresult.trunc, logfc < sliderval1 & -log10(p.adjPVAL) < sliderval2 & logfc > sliderval3 & -log10(p.adjPVAL) > sliderval4)
     head(sigresult.trunc)
     nuIDlist<-sigresult.trunc$nuID
     entrezlist<-nuID2EntrezID(nuIDlist,lib.mapping='lumiHumanIDMapping')
     sigresult.trunc$entrez<-entrezlist
     print(head(sigresult.trunc))
     d<-sigresult.trunc[,c("entrez","logfc")]
     
     ## feature 1: numeric vector
     geneList = d[,2]
     ## feature 2: named vector
     names(geneList) = as.character(d[,1])
     ## feature 3: decreasing order
     geneList = sort(geneList, decreasing = TRUE)
     
     gsea.result<-gsePathway(geneList,nPerm=10000,pvalueCutoff=0.2,pAdjustMethod = "BH",verbose=F)
     gsea.result
   })
   
   
   output$sigTable1<-DT::renderDataTable({
     sliderval1<-input$volcanoFC
     sliderval2<-input$volcanoSIG
     sliderval3<-input$volcanoFClow
     sliderval4<-input$volcanoSIGlow
     res<-signatureQuery()
     
     if(whichSIGfun()!="all"){
       res2<-subset(res,logfc > sliderval3 & logfc < sliderval1 & -log10(p.adjPVAL) > sliderval4 & -log10(p.adjPVAL) < sliderval2 )
     }else{
       res2<-subset(res, (logfc > abs(sliderval3) | logfc < -abs(sliderval3)) &
                      (logfc < abs(sliderval1) & logfc > -abs(sliderval1)) & 
                      -log10(p.adjPVAL) > sliderval4 &
                      -log10(p.adjPVAL) < sliderval2 
                    
       )
     }
     if(whichModuleFun()!="All"){
        res2<-res2[which(res2$module==whichModuleFun()),]
     }
     #sigresult.trunc<-subset(sigresult.trunc, logfc < sliderval1 & -log10(p.adjPVAL) < sliderval2 & logfc > sliderval3 & -log10(p.adjPVAL) > sliderval4)
     datatable(res2,filter="top")},server = FALSE,options=list(lengthMenu = c(5, 10, 15,20,50), pageLength = 15) 
     )
   
   output$volcano<-renderPlot({
     reacstudySIG<-as.numeric(studySIGfun())
     reacsetSIG<-as.numeric(squareSIGfun())
     reacedgeSIG<-as.numeric(edgeSIGfun())
     studySIG<-names(studies)[reacstudySIG]
     setlistSIG<-setlistnameslist[[reacstudySIG]]
     squareSIG<-setlistSIG[[reacsetSIG]]
     
     sliderval1<-input$volcanoFC
     sliderval2<-input$volcanoSIG
     sliderval3<-input$volcanoFClow
     sliderval4<-input$volcanoSIGlow
     #res<-signatureQuery()
     
     res<-signatureQuery()
     res<-res[order(res$module),]
     print("###########################################################################################...")
     print("start debug 11.09.21")
     print("whichModuleFun")
     print(whichModuleFun())
     
     squareSIG<-squareSIGfun()
     print("squareSIG")
     print(squareSIG)
     moduleq<-paste("MATCH (n:wgcna {square:'",squareSIG,"'}) RETURN DISTINCT n.name AS modulename",sep="")
     modulenames<-cypher(graph,moduleq)
     modulenames<-modulenames$modulename[which(modulenames$modulename!="grey")]
     print(modulenames)
     modulenames<-c("All",modulenames)
     
     print("##########################################################################################444444////")
     
     if(whichModuleFun()!="All"){
        res<-res[which(res$module==whichModuleFun()),]
     }
     
     if(whichSIGfun()!="all"){
       res2<-subset(res,logfc > sliderval3 & logfc < sliderval1 & -log10(p.adjPVAL) > sliderval4 & -log10(p.adjPVAL) < sliderval2 )
     }else{
       res2<-subset(res, (logfc > abs(sliderval3) | logfc < -abs(sliderval3)) &
                         (logfc < abs(sliderval1) & logfc > -abs(sliderval1)) & 
                         -log10(p.adjPVAL) > sliderval4 &
                         -log10(p.adjPVAL) < sliderval2 
                         
                      )
     }
     
     #sigresult.trunc<-subset(sigresult.trunc, logfc < sliderval1 & -log10(p.adjPVAL) < sliderval2 & logfc > sliderval3 & -log10(p.adjPVAL) > sliderval4)
     
     if (vpplotfun()=="module"){
     vp<-ggplot(res,aes(logfc,-log10(p.adjPVAL)))+
       #geom_point(aes(col=logfc),size=9,alpha=0.4)+
       geom_point(aes(col=module),size=6,alpha=0.7)+ 
       #scale_color_gradient2(low="blue",mid="lightyellow",high="red") +
       scale_color_manual(values=unique(res$module))+
       #geom_text_repel(aes(label=SYMBOL),data=subset(res, logfc < sliderval1 & -log10(p.adjPVAL) < sliderval2 & logfc > sliderval3 & -log10(p.adjPVAL) > sliderval4 ),col="darkslategrey",size=6) +
       geom_text_repel(aes(label=SYMBOL),data=res2,col="darkslategrey",size=6) +
       geom_vline(xintercept=sliderval1, linetype="dashed", color = "red") +
       {if(whichSIGfun()=="all") geom_vline(xintercept=-sliderval1, linetype="dashed", color = "red")} +
       geom_vline(xintercept=sliderval3, linetype="dashed", color = "blue") +
       {if(whichSIGfun()=="all") geom_vline(xintercept=-sliderval3, linetype="dashed", color = "blue")} +
       geom_hline(yintercept=sliderval2, linetype="dashed", color = "red") +
       geom_hline(yintercept=sliderval4, linetype="dashed", color = "blue") +
       #theme_tufte() +
       theme_light() +
       #theme_minimal() +
       theme(axis.title.y = element_text(size = rel(2.5), angle = 90, color = "slategrey",margin = margin(t = 0, r = 50, b = 50, l = 50))) +
       theme(axis.title.x = element_text(size = rel(2.5), color = "slategrey",margin = margin(t = 50, r = 50, b = 0, l = 50))) + 
       theme(axis.text = element_text(size = rel(1.5))) +
       #theme(legend.text = element_text(size = rel(1.2))) +
       #theme(legend.title = element_text(size = rel(1.2))) +
       #theme(legend.text=element_text(size=rel(1.4))) +
       #theme(legend.key.size = unit(1,"cm")) +
       #labs(y = phenofun2()) +
       ggtitle(paste(setlistSIG[[reacsetSIG]],"edge:",reacedgeSIG)) +
       theme(plot.title = element_text(size=rel(2.5),hjust = 0.5,color="maroon"))
     }else{
        vp<-ggplot(res,aes(logfc,-log10(p.adjPVAL)))+
           geom_point(aes(col=logfc),size=9,alpha=0.4)+
           #geom_point(aes(col=module),size=9,alpha=0.5)+ 
           scale_color_gradient2(low="blue",mid="lightyellow",high="red") +
           #scale_color_manual(values=unique(res$module))+
           #geom_text_repel(aes(label=SYMBOL),data=subset(res, logfc < sliderval1 & -log10(p.adjPVAL) < sliderval2 & logfc > sliderval3 & -log10(p.adjPVAL) > sliderval4 ),col="darkslategrey",size=6) +
           geom_text_repel(aes(label=SYMBOL),data=res2,col="darkslategrey",size=6) +
           geom_vline(xintercept=sliderval1, linetype="dashed", color = "red") +
           {if(whichSIGfun()=="all") geom_vline(xintercept=-sliderval1, linetype="dashed", color = "red")} +
           geom_vline(xintercept=sliderval3, linetype="dashed", color = "blue") +
           {if(whichSIGfun()=="all") geom_vline(xintercept=-sliderval3, linetype="dashed", color = "blue")} +
           geom_hline(yintercept=sliderval2, linetype="dashed", color = "red") +
           geom_hline(yintercept=sliderval4, linetype="dashed", color = "blue") +
           #theme_tufte() +
           theme_light() +
           #theme_minimal() +
           theme(axis.title.y = element_text(size = rel(2.5), angle = 90, color = "slategrey",margin = margin(t = 0, r = 50, b = 50, l = 50))) +
           theme(axis.title.x = element_text(size = rel(2.5), color = "slategrey",margin = margin(t = 50, r = 50, b = 0, l = 50))) + 
           theme(axis.text = element_text(size = rel(1.5))) +
           #theme(legend.text = element_text(size = rel(1.2))) +
           #theme(legend.title = element_text(size = rel(1.2))) +
           #theme(legend.text=element_text(size=rel(1.4))) +
           #theme(legend.key.size = unit(1,"cm")) +
           #labs(y = phenofun2()) +
           ggtitle(paste(setlistSIG[[reacsetSIG]],"edge:",reacedgeSIG)) +
           theme(plot.title = element_text(size=rel(2.5),hjust = 0.5,color="maroon"))
     }
     
     vp
     
   },height=900)
   
   #output$heatmap<-renderPlot({
   output$heatmapMake<-renderPlot({
     
     scale<-scaleSIGfun()
     
     reacstudySIG<-as.numeric(studySIGfun())
     reacsetSIG<-as.numeric(squareSIGfun())
     reacedgeSIG<-as.numeric(edgeSIGfun())
     studySIG<-names(studies)[reacstudySIG]
     setlistSIG<-setlistnameslist[[reacstudySIG]]
     squareSIG<-setlistSIG[[reacsetSIG]]
     
     reacstudySIG2<-as.numeric(studySIGfun2())
     reacsetSIG2<-as.numeric(squareSIGfun2())
     reacedgeSIG2<-as.numeric(edgeSIGfun2())
     studySIG2<-names(studies)[reacstudySIG2]
     setlistSIG2<-setlistnameslist[[reacstudySIG2]]
     squareSIG2<-setlistSIG2[[reacsetSIG2]]
     
     sliderval1<-input$volcanoFC
     sliderval2<-input$volcanoSIG
     sliderval3<-input$volcanoFClow
     sliderval4<-input$volcanoSIGlow
     res<-signatureQuery()
     
     
     res<-signatureQuery()
     print("original res")
     print(head(res))
     
     if(whichSIGfun()!="all"){
       res<-subset(res,logfc > sliderval3 & logfc < sliderval1 & -log10(p.adjPVAL) > sliderval4 & -log10(p.adjPVAL) < sliderval2 )
     }else{
       res<-subset(res, (logfc > abs(sliderval3) | logfc < -abs(sliderval3)) &
                      (logfc < abs(sliderval1) & logfc > -abs(sliderval1)) & 
                      -log10(p.adjPVAL) > sliderval4 &
                      -log10(p.adjPVAL) < sliderval2 
                    
       )
     }
     if(whichModuleFun()!="All"){
        res<-res[which(res$module==whichModuleFun()),]
     }
     print("original res after if ")
     print(head(res))
     
     #res<-subset(res, logfc < sliderval1 & -log10(p.adjPVAL) < sliderval2 & logfc > sliderval3 & -log10(p.adjPVAL) > sliderval4 )
     
     #get wgcna modules for selected probes from visualised dataset
     probesvis<-res$nuID
     query2<-paste("MATCH (p:PROBE {square:'",squareSIG2,"',edge:'",reacedgeSIG2,"'})-[r]-(n:wgcna) RETURN p.name as nuID2, n.name as wgcna2",sep="")
     res2<-cypher(graph,query2)
     
     print("head(res) before merging")
     print(head(res))
     
     print("head(res2) before merging")
     print(head(res2))
     
     res<-merge(res,res2,by.x="nuID",by.y="nuID2",x.all=TRUE,y.all=FALSE)
     
     print("head(res) after merging")
     print(head(res))
     
     
     #get individuals
     query<-paste("MATCH (p:person {square:'",squareSIG2,"'}) RETURN p.name AS name, p.class1 AS class1, p.class2 AS class2",sep="")#debug
     people<-cypher(graph,query)
     
     print("people")
     print(people)
     
     #temp debug fix
     if(people$name[1]=="1"){
       people$name<-as.character(eval(parse(text=paste("pData(",squareSIG2,")$sampleID",sep=""))))
       
     }
     
     print(questions[[squareSIG2]]$pathway2[[1]])#debug
     edges=list(
       "edge1"=people[which(people$class2==questions[[squareSIG2]]$pathway2[[1]]),],
       "edge2"=people[which(people$class2==questions[[squareSIG2]]$pathway2[[2]]),],
       "edge3"=people[which(people$class1==questions[[squareSIG2]]$pathway1[[1]]),],
       "edge4"=people[which(people$class1==questions[[squareSIG2]]$pathway1[[2]]),],
       "edge5"=people[which(people$class1==questions[[squareSIG2]]$pathway1[[1]]|people$class1==questions[[squareSIG2]]$pathway1[[2]]),] 
     )
     
     print("edges")
     print(edges)
     
     #exprdata<-eval(parse(text=squareSIG))#debug
     exprdata<-eval(parse(text=squareSIG2))#debug
     sigvec<-res$nuID
     if(class(exprdata)=="LumiBatch"){
       data.v<-lumiT(exprdata,simpleOutput=FALSE)
       
     }else{
       data.v<-lumiT(exprdata,simpleOutput=FALSE,method="log2")
       
     }
     #Quantile normalise
     data.q<-lumiN(data.v,method="quantile")
     
     print("Normalisation done")
     
     heatdata.all<-data.q
     
     print("edge data selected")
     print("sigvec")
     print(sigvec)
     print(str(heatdata.all))
     
     print("head(exprs(heatdata.all))")
     print(head(exprs(heatdata.all)))
     
     try(colnames(heatdata.all)<-pData(heatdata.all)[,"sampleID"])
     
     try(print(pData(heatdata.all)[,"sample_name"]))
     
     try(colnames(heatdata.all)<-pData(heatdata.all)[,"sample_name"])
     
     print(edges[[reacedgeSIG2]])
     
     print(edges[[reacedgeSIG2]]$name)
     print("before")
     
     sigvecTV<-sigvec%in%rownames(exprs(heatdata.all))
     print("sigvectv")
     print(sigvecTV)
     print("newsigvec")
     print(length(sigvec))
     print(length(sigvecTV))
     newsigvec<-sigvec[sigvecTV]
     print(newsigvec)
     #print(rownames(exprs(heatdata.all)))
     print(colnames(exprs(heatdata.all)))
     print(edges[[as.numeric(reacedgeSIG2)]]$name)
     heatdata<-exprs(heatdata.all)[newsigvec,edges[[as.numeric(reacedgeSIG2)]]$name]
     
     
     print("signature data selected")
      annotation<-edges[[as.numeric(reacedgeSIG2)]]
      print("annotation")
      print(annotation)
      rownames(annotation)<-annotation$name
      annotation<-annotation[,2:3]
      print("annotation")
      print(annotation)
      #rowAnn<-as.data.frame(res$module[sigvecTV])
      #rowAnn<-as.data.frame(res$wgcna2[sigvecTV])
      rowAnn<-as.data.frame(res$wgcna2)
      
      colnames(rowAnn)<-"module"
      print("res")
      print(res)
      print("rowAnn")
      print(rowAnn)
      print(res$nuID)
      rownames(rowAnn)<-rownames(heatdata)
      print("rowAnn")
      print(rowAnn)
     print(str(heatdata))
     
     print(rownames(heatdata))
     print(res$nuID)
   
     cur_dev <- dev.cur()
     if(nrow(res)<100){boolrow<-TRUE}else{boolrow<-FALSE}
     
     print("debug")
     print(res$wgcna2)
     print(as.factor(res$wgcna2))
     print(levels(as.factor(res$wgcna2)))
     print(list("module"=levels(as.factor(res$wgcna2))))
     
     levs<-levels(as.factor(res$wgcna2))
     newlevs<-structure(levs, names=as.character(levs))
     print(levs)
     print(newlevs)
     
     c1l<-structure(levels(as.factor(c("#7b3294","#c2a5cf"))),names=as.character(levels(as.factor(annotation$class1))))
     c2l<-structure(levels(as.factor(c("#008837","#a6dba0"))),names=as.character(levels(as.factor(annotation$class2))))
     
     #ac<-list(module=newlevs)
     ac<-list(
        class1=c1l,
        class2=c2l,
        module=newlevs)
     
     ac2<-list(
        class1=c1l,
        class2=c2l)
     
     print("ac")
     print(ac)
     print("do heatmap")
     print("head(heatdata)")
     print(head(heatdata))
     #order heatdata
     heatdata<-heatdata[,order(colnames(heatdata))]
     
     if (cellmatfun()=="no"){
        heatm<-pheatmap(heatdata,annotation_col=annotation,annotation_row=rowAnn,show_rownames=boolrow,scale=scale,labels_row=res$SYMBOL,annotation_colors = ac,fontsize_row = 10)
        dev.set(cur_dev)
        print(heatm)
        
     }else if (cellmatfun()=="yes"){
        print("debug 0")
        heatm<-pheatmap(heatdata,annotation_col=annotation,annotation_row=rowAnn,show_rownames=boolrow,scale=scale,labels_row=res$SYMBOL,annotation_colors = ac,fontsize_row = 9)
        ##
        #Adding in plot of cellproptable
        
        #get the data
        
        #modquery<-paste("MATCH (c:cellEx)-[r]-(n:wgcna {square:'",squareSIG2,"'}) WHERE c.name =~ '.*_xp_3' RETURN DISTINCT c.name AS cell,n.name AS module,r.qvalue AS qval",sep="")
        
        #query<-"MATCH (p:personCell {square:'blood.PCF.defPC'}) RETURN p.personName as person,p.name AS cell, p.value AS value"
        query<-paste("MATCH (p:personCell {square:'",squareSIG2,"'}) RETURN p.personName as person,p.name AS cell, p.value AS value",sep="")
        res<-cypher(graph,query)
        
        persons<-unique(res$person)
        cells<-unique(res$cell)
        dfc<-data.frame()
        
        print("debug 1")
        
        for (person in persons) {
           res2<-res[which(res$person==person),]
           for (cell in cells){
              dfc[person,cell]<-as.numeric(res2[which(res2$cell==cell),"value"])
           }
           
        }
        
        print("debug 2")
        
        dfc<-dfc[edges[[as.numeric(reacedgeSIG2)]]$name,]
        dfc<-t(dfc[order(rownames(dfc)),])
        dfc<-dfc[,order.dendrogram(as.dendrogram(heatm$tree_col))]
        #dfc2<-melt(dfc)
        #colnames(dfc2)<-c("cell","person","value")
        ##HERE
        #bp<-barplot(as.matrix(dfc),las=2,col=brewer.pal(12,"Set3"),legend.text=rownames(dfc),beside=F,args.legend=c(x=ncol(dfc)*1.35,y=1,cex=.7))
        
        print("debug 3")
        
        bp2<-ggplot(res, aes(fill=cell, y=value, x=person)) + 
           geom_bar(position="stack", stat="identity")
        print("debug 3.1")
        
        print(head(dfc))
        print(head(annotation))
        print(head(rownames(dfc)))
        print(head(ac2))
        
        heatm2<-pheatmap(dfc,annotation_col=annotation,labels_row=rownames(dfc),cluster_cols=F,annotation_colors = ac2)
        
        #########-
        
        print("debug 4")
        
        query<-paste("MATCH (p:personME {square:'",squareSIG2,"'}) RETURN p.personName as person,p.name AS module, p.value AS value",sep="")
        res<-cypher(graph,query)
        
        persons<-unique(res$person)
        print(persons)
        MEs<-unique(res$module)
        print(MEs)
        dfm<-data.frame()
        
        print("debug 5")
        
        for (person in persons) {
           res2<-res[which(res$person==person),]
           print("debug 1")
           for (ME in MEs){
              dfm[person,ME]<-as.numeric(res2[which(res2$module==ME),"value"])
              print("debug 2")
           }
           
        }
        dfm<-dfm[edges[[as.numeric(reacedgeSIG2)]]$name,]
        
        print("debug 6")
        
        dfm<-t(dfm[order(rownames(dfm)),])
        print("debug 7")
        dfm<-dfm[,order.dendrogram(as.dendrogram(heatm$tree_col))]
        print("debug 8")
        #dfm2<-melt(dfm)
        print("debug 9")
        #colnames(dfm2)<-c("ME","person","value")
        print("debug 10")
        ##HERE
        #bp<-barplot(as.matrix(dfc),las=2,col=brewer.pal(12,"Set3"),legend.text=rownames(dfc),beside=F,args.legend=c(x=ncol(dfc)*1.35,y=1,cex=.7))
        bp3<-ggplot(res, aes(fill=module, y=value, x=person)) + 
           geom_bar(position="stack", stat="identity")
        heatm3<-pheatmap(dfm,annotation_col=annotation,labels_row=rownames(dfm),cluster_cols=F,annotation_colors = ac2)
        
        #########-
        print("debug 11")
        
        ##
        plot_list<-list(heatm[[4]],heatm2[[4]],heatm3[[4]])
        dev.set(cur_dev)
        #g<-do.call(grid.arrange,plot_list)
        #g<-grid.arrange(grobs=plot_list, ncol=1,nrow=2)
        #g<-plot_grid(heatm[[4]],heatm2[[4]],ncol=1,nrow=2,greedy=F,scale=c(1,scalefactorfun()))
        g<-plot_grid(heatm[[4]],heatm2[[4]],heatm3[[4]],ncol=1,align="v",scale=scalefactorfun(),rel_heights = c(2,1,1))
        #ggsave("g.pdf",g)
        if(input$returnpdf==TRUE){
           #pdf("plot.pdf",width=12,height=7,onefile=FALSE)
           ggsave("plot.pdf",g,device="pdf",width=downloadPlotWidthFun(),height=downloadPlotHeightFun(),bg="white")
           #dev.off()
        }
        print(g)
        #print(heatm)
        #barplot(as.matrix(dfc),las=2,col=brewer.pal(12,"Set3"),legend.text=rownames(dfc),beside=F,args.legend=c(x=ncol(dfc)*1.35,y=1,cex=.7))
        #par(mfrow=c(1,1))
     }
     
     
     
     
    
   })
   output$heatmap <- renderUI({
      plotOutput("heatmapMake", height = as.numeric(input$plotheight2))
      #plotOutput("gigamat", height = 7000)
   })
   
   
   #output:enrichPathway
   output$enrichTable1<-DT::renderDataTable(datatable(as.data.frame(enrichfun())),server = FALSE,options=list(lengthMenu = c(5, 10, 15,20,50), pageLength = 15),filter="top")
   
   #output:gsePathway
   output$gseaTable2<-DT::renderDataTable(datatable(as.data.frame(gseafun())),server = FALSE,options=list(lengthMenu = c(5, 10, 15,20,50), pageLength = 15),filter="top")
   
   output$barplot1<-renderPlot({
     plotdata<-enrichfun()
     head(as.data.frame(plotdata))
     barplot(plotdata)
     
   },height = 800)
   
   output$dotplot1<-renderPlot({
     plotdata<-enrichfun()
     head(as.data.frame(plotdata))
     dotplot(plotdata,font.size=20)
     
   },height = 800)
   
   output$emapplot1<-renderPlot({
     plotdata<-enrichfun()
     head(as.data.frame(plotdata))
     emapplot(enrichplot::pairwise_termsim(plotdata))
     
   },height = 800)
   
   output$cnetplot1<-renderPlot({
     plotdata<-enrichfun()
     head(as.data.frame(plotdata))
     cnetplot(enrichplot::pairwise_termsim(plotdata))
     
   },height = 800)
   
   output$emapplot2<-renderPlot({
     plotdata<-gseafun()
     head(as.data.frame(plotdata))
     emapplot(enrichplot::pairwise_termsim(plotdata))
     
   },height = 800)
   
   output$gseaplot2<-renderPlot({
     plotdata<-gseafun()
     head(as.data.frame(plotdata))
     gseaplot(enrichplot::pairwise_termsim(plotdata),geneSetID=1)
     
   },height = 800)
   
   output$moduleSankey<-renderSankeyNetwork({
      reacstudySIG<-as.numeric(studySIGfun())
      reacsetSIG<-as.numeric(squareSIGfun())
      reacedgeSIG<-as.numeric(edgeSIGfun())
      studySIG<-names(studies)[reacstudySIG]
      setlistSIG<-setlistnameslist[[reacstudySIG]]
      squareSIG<-setlistSIG[[reacsetSIG]]
      
      reacstudySIG2<-as.numeric(studySIGfun2())
      reacsetSIG2<-as.numeric(squareSIGfun2())
      reacedgeSIG2<-as.numeric(edgeSIGfun2())
      studySIG2<-names(studies)[reacstudySIG2]
      setlistSIG2<-setlistnameslist[[reacstudySIG2]]
      squareSIG2<-setlistSIG2[[reacsetSIG2]]
      
      set1name<-paste("mDat",squareSIG,sep="")
      set2name<-paste("mDat",squareSIG2,sep="")
      print(set1name)
      print(set2name)
      set1<-allmodules[[set1name]][[2]]
      set2<-allmodules[[set2name]][[2]]
      set1df<-data.frame(cbind(names(set1),set1))
      colnames(set1df)<-c("probe","module")
      set2df<-data.frame(cbind(names(set2),set2))
      colnames(set2df)<-c("probe","module")
      set1modules<-paste(set1df$module,"S1",sep="")
      set2modules<-paste(set2df$module,"S2",sep="")
      set1df$module<-set1modules
      set2df$module<-set2modules
      
      head(set1df)
      
      set1mlist<-unique(set1modules)
      set2mlist<-unique(set2modules)
      
      nodes<-c(set1mlist,set2mlist)
      nodes<-data.frame(nodes)
      colnames(nodes)<-"name"
      nodes<-cbind(seq(0,nrow(nodes)-1),nodes)
      colnames(nodes)<-c("index","name")
      
      links<-data.frame(matrix(ncol=6,nrow=0),stringsAsFactors = FALSE)
      colnames(links)<-c("source","value","target","node1","ngenes","node2")
      a=0
      b=length(set1mlist)+1
      for (md1 in set1mlist){
         b=length(set1mlist)+1
         for (md2 in set2mlist){
            probelist1<-set1df[which(set1df$module==md1),"probe"]
            #print(head(probelist1))
            probelist2<-set2df[which(set2df$module==md2),"probe"]
            #print(head(probelist2))
            count<-sum(probelist1%in%probelist2)
            links<-rbind(links,list(a,b,count,md1,count,md2))
            
            b=b+1
         }
         a=a+1
      }
      colnames(links)<-c("source","target","value","node1","ngenes","node2")
      head(links)
      head(nodes)
      links$source<-as.integer(links$source)
      links$target<-as.integer(links$target-1)
      links$value<-as.numeric(links$value)
      links$ngenes<-as.numeric(links$ngenes)
      
      sum(links$value)
      nrow(links)
      
      links2<-links[,1:3]
      nodes2<-data.frame(nodes[,2],stringsAsFactors = FALSE)
      colnames(nodes2)<-"name"
      
      #remove zero-count link table entries
      links<-links[which(links$ngenes>0),]
      
      data<-list("nodes"=nodes2,"links"=links)
      print("links and nodes done")
      
      print(sort(data$links$source))
      print(sort(data$links$target))
      
      sankeyNetwork(Links = data$links, Nodes = data$nodes, Source = 'source',Target = 'target',Value = 'value', NodeID = 'name',fontSize = 12,LinkGroup = 'node1')
      
      
      
      
   })
   
   #Single module####
   
   #igraphSP####
   output$igraphSP<-renderPlot({
     #get module
     #prefix<-paste("/home/rstudio/output/build/",buildfun(),"/",sep="")
     prefix<-"/home/rstudio/output/build/"
     #print("line1 inside reactive")
     data0<-read.csv(paste(prefix,setfun(),"/5_WGCNA_4k_tv/results/ModuleEigengenes_all_samples_withPheno.csv",sep=""))
     #print("line2")
     data<-read.csv(paste(prefix,setfun(),"/5_WGCNA_4k_tv/results/ModuleEigengenes_all_samples.csv",sep=""))
     #print("line3")
     squarelist<-unlist(strsplit(setfun(),"_"))
     if(length(squarelist)==3){
       square<-squarelist[3]
     }else{square<-paste(squarelist[c(3,4)],collapse="_")}
     #print("square")
     #print(square)
     menames<-sort(colnames(data)[2:length(colnames(data))])
     my.module<-me1fun()
     
  
     
     #query
     query.base<-paste("MATCH (c:CELL)-[r0]-(x1)-[r1]-(n:wgcna {square:'",square,"',edge:'",edgefun(),"',name:'",unlist(strsplit(me1fun(),"ME"))[2],"'})-[r2]-(x2) WHERE (x1:cellEx OR x1:cellprop) AND (x2:reactomePW OR x2:ImmunePW OR x2:pheno OR x2:PalWangPW)",sep="")
     print(query.base)
     nodelist<-c("c","x1","n","x2")
     edgetrips<-list(c("c","r0","x1"),c("x1","r1","n"),c("n","r2","x2"))
     
     
     igr<-igraph_plotter(query.base,nodelist,edgetrips,rimpar="diffME",plot=TRUE,csv=TRUE,prefix=cytodir,filename = "igraphWGCNAannot",vertexsize=vertexsizefun(),legendcex=legendcexfun(),plotd3=FALSE,vertex.label.cex=vlc(),lay_out=eval(parse(text=layoutfun())),edgewidth=ewf())
     
     # igr<-"FAIL"
     # igr<-try(igraph_plotter(query.base,nodelist,edgetrips,rimpar="diffME",plot=TRUE,csv=TRUE,prefix=cytodir,filename = "igraphWGCNAannot",vertexsize=vertexsizefun(),legendcex=legendcexfun(),plotd3=FALSE,vertex.label.cex=vlc(),lay_out=eval(parse(text=layoutfun())),edgewidth=ewf()))
     # #print(class(igr))
     # if(igr=="FAIL"){
     #    #query
     #    query.base<-paste("MATCH (n:wgcna {square:'",square,"',edge:'",edgefun(),"',name:'",unlist(strsplit(me1fun(),"ME"))[2],"'})-[r2]-(x2) WHERE (x2:cellEx OR x2:cellprop OR x2:reactomePW OR x2:ImmunePW OR x2:pheno OR x2:PalWangPW)",sep="")
     #    #print(query.base)
     #    nodelist<-c("n","x2")
     #    edgetrips<-list(c("n","r2","x2"))
     #    igr<-igraph_plotter(query.base,nodelist,edgetrips,rimpar="diffME",plot=FALSE,csv=TRUE,prefix=cytodir,filename = "igraphwgcna",plotd3=FALSE,return_graph=TRUE)
     # }
     
     
     if(input$returnpdf==TRUE){
       pdf("plot.pdf",width=14,height=12,onefile=FALSE)#,width=12,height=7,
       igr<-igraph_plotter(query.base,nodelist,edgetrips,rimpar="diffME",plot=TRUE,csv=TRUE,prefix=cytodir,filename = "igraphWGCNAannot",vertexsize=vertexsizefun(),legendcex=legendcexfun(),plotd3=FALSE,vertex.label.cex=vlc(),lay_out=eval(parse(text=layoutfun())),edgewidth=ewf())
       dev.off()
     }
   })
   
   #d3graphSP####
   output$d3graphSP<-renderForceNetwork({
     #get module
     #prefix<-paste("/home/rstudio/output/build/",buildfun(),"/",sep="")
     prefix<-"/home/rstudio/output/build/"
     #print("line1 inside reactive")
     data0<-read.csv(paste(prefix,setfun(),"/5_WGCNA_4k_tv/results/ModuleEigengenes_all_samples_withPheno.csv",sep=""))
     #print("line2")
     data<-read.csv(paste(prefix,setfun(),"/5_WGCNA_4k_tv/results/ModuleEigengenes_all_samples.csv",sep=""))
     #print("line3")
     squarelist<-unlist(strsplit(setfun(),"_"))
     if(length(squarelist)==3){
       square<-squarelist[3]
     }else{square<-paste(squarelist[c(3,4)],collapse="_")}
     print("square")
     print(square)
     menames<-sort(colnames(data)[2:length(colnames(data))])
     my.module<-me1fun()
     #query
     query.base<-paste("MATCH (c:CELL)-[r0]-(x1)-[r1]-(n:wgcna {square:'",square,"',edge:'",edgefun(),"',name:'",unlist(strsplit(me1fun(),"ME"))[2],"'})-[r2]-(x2) WHERE (x1:cellEx OR x1:cellprop) AND (x2:reactomePW OR x2:ImmunePW OR x2:pheno OR x2:PalWangPW)",sep="")
     print(query.base)
     nodelist<-c("c","x1","n","x2")
     edgetrips<-list(c("c","r0","x1"),c("x1","r1","n"),c("n","r2","x2"))
     #igraphplotter with igr output
     igr<-"FAIL"
     try(igr<-igraph_plotter(query.base,nodelist,edgetrips,rimpar="diffME",plot=FALSE,csv=TRUE,prefix=cytodir,filename = "igraphwgcna",plotd3=FALSE,return_graph=TRUE))
     print(class(igr))
     if(class(igr)!="character"){
       #query
       query.base<-paste("MATCH (n:wgcna {square:'",square,"',edge:'",edgefun(),"',name:'",unlist(strsplit(me1fun(),"ME"))[2],"'})-[r2]-(x2) WHERE (x2:cellEx OR x2:cellprop OR x2:reactomePW OR x2:ImmunePW OR x2:pheno OR x2:PalWangPW)",sep="")
       #print(query.base)
       nodelist<-c("n","x2")
       edgetrips<-list(c("n","r2","x2"))
       igr<-igraph_plotter(query.base,nodelist,edgetrips,rimpar="diffME",plot=FALSE,csv=TRUE,prefix=cytodir,filename = "igraphwgcna",plotd3=FALSE,return_graph=TRUE)
     }
     
     #get nodes and edges for DT
     currentnodes<-read.csv(file.path(cytodir,"igraphwgcna_nodes.csv"))
     print("currentnodes")
     print(currentnodes)
     output$moduleNodesSP<-DT::renderDataTable(currentnodes,server = FALSE,options=list(lengthMenu = c(5, 10, 15,20,50), pageLength = 50),filter="top") 
     currentedges<-read.csv(file.path(cytodir,"igraphwgcna_edges.csv"))
     dict<-currentnodes[,c(2,3)]
     for(row in 1:nrow(dict)){
       print(row)
       currentedges[currentedges==as.character(dict[row,"node_id"])]<-as.character(dict[row,"nodename"])
     }
     #print(currentedges)
     output$moduleEdgesSP<-DT::renderDataTable(currentedges,server = FALSE,options=list(lengthMenu = c(5, 10, 15,20,50), pageLength = 50),filter="top")
     
    
     V(igr)$name=V(igr)$nodename
     print("names before conversion")
     print(sort(V(igr)$name))
     V(igr)$name=make.names(V(igr)$name,unique=TRUE)# this becomes necessary where nodes of different type have the same name (e.g. pheno and CELL nodes named 'RBC')
     print('names after make.names before conversion')
     print(sort(V(igr)$name))
     #convert
     # Convert to object suitable for networkD3
     ig_d3 <- igraph_to_networkD3(igr, group = V(igr)$kind)
     print("conversion successful")
     print(ig_d3$links)
     print(ig_d3$nodes)
     # Create force directed network plot
     #R version 3.3.0 (D3 version 3)
     # forceNetwork(Links = ig_d3$links, Nodes = ig_d3$nodes,height=1200,width=1200,charge= -300,colourScale = JS("d3.scale.category20()"),
     #              linkDistance = 100,fontSize=14,opacity=.9,opacityNoHover = .8,radiusCalculation = JS(" Math.sqrt(300)+6"),
     #              Source = 'source', Target = 'target',
     #              NodeID = 'name', Group = 'group',legend=TRUE,zoom = TRUE,bounded=FALSE)
     #R version 3.3.3 (D3version 4)
     forceNetwork(Links = ig_d3$links, Nodes = ig_d3$nodes,height=1200,width=1200,charge= -300,colourScale = JS("d3.scaleOrdinal(d3.schemeCategory10)"),
                  linkDistance = 100,fontSize=14,opacity=.9,opacityNoHover = .8,radiusCalculation = JS(" Math.sqrt(300)+6"),
                  Source = 'source', Target = 'target',
                  NodeID = 'name', Group = 'group',legend=TRUE,zoom = TRUE,bounded=FALSE)
     
   })
   
   #moduleNodesSP + moduleEdgesSP####
   cytodirlist<-system(paste("ls",cytodir),intern=TRUE)
   print("igraphwgcna_nodes.csv%in%cytodirlist")
   print("igraphwgcna_nodes.csv"%in%cytodirlist)
   if("igraphwgcna_nodes.csv"%in%cytodirlist){
     currentnodes<-read.csv(file.path(cytodir,"igraphwgcna_nodes.csv"))
     print(currentnodes)
     output$moduleNodesSP<-DT::renderDataTable(currentnodes,server = FALSE,options=list(lengthMenu = c(5, 10, 15,20,50), pageLength = 15),filter="top") 
     currentedges<-read.csv(file.path(cytodir,"igraphwgcna_edges.csv"))
     print(currentedges)
     output$moduleEdgesSP<-DT::renderDataTable(currentedges,server = FALSE,options=list(lengthMenu = c(5, 10, 15,20,50), pageLength = 15),filter="top") 
   }
   #######Radarplots####
   output$radar<-renderPlot({
     #get module
     #prefix<-paste("/home/rstudio/output/build/",buildfun(),"/",sep="")
     prefix<-"/home/rstudio/output/build/"
     #print("line1 inside reactive")
     squarelist<-unlist(strsplit(setfun(),"_"))
     if(length(squarelist)==3){
       square<-squarelist[3]
     }else{square<-paste(squarelist[c(3,4)],collapse="_")}
     print("square")
     print(square)
     #menames<-sort(colnames(data)[2:length(colnames(data))])
     my.module<-me1fun()
     my.module<-substr(my.module,3,nchar(my.module))
     my.edge<-edgefun()
     #plot1###################
     
     
     #query<-paste("MATCH (n:wgcna {square:'nonLTBI.LTBI',name:'red'})-[r]-(x) WHERE (x:cellEx) RETURN DISTINCT n.name AS module,r.qvalue as qvalue, x.name as name, labels(x) as kind",sep="")
     query<-paste("MATCH (n:wgcna {square:'",square,"',name:'",my.module,"',edge:'",my.edge,"'})-[r]-(x) WHERE (x:cellEx) RETURN DISTINCT n.name AS module,toFloat(r.qvalue) as qvalue,x.name as name, labels(x) as kind",sep="")
     print(query)
     anq<-"MATCH (c:cellEx)-[r]-(c2:CELL) RETURN DISTINCT c.name as cells,c2.name AS group"
     cells<-cypher(graph,anq)
     
     
     
     cells$SIG<-rep(0,nrow(cells))
     
     res.cellex<-cypher(graph,query)
     #print("debug1")
     print("debug dups")
     print(res.cellex)
     
     res.cellex<-res.cellex[!duplicated(res.cellex[c(3)]),]
    
     print(res.cellex)
     print(is.null(res.cellex))
     if(!is.null(res.cellex)){
     res.cellex$SIG<--log10(res.cellex$qvalue)
     res.cellex$SIG[which(res.cellex$SIG>5)]<-5
     res.cellex2<-data.frame(as.numeric(res.cellex$SIG))
     rownames(res.cellex2)<-res.cellex$name
     colnames(res.cellex2)<-"sig"
     for (row in 1:nrow(res.cellex2)){
       cells$SIG[which(cells$cells==res.cellex$name[row])]<-res.cellex2$sig[row]
       
     }
     res.cellex3<-data.frame(as.numeric(cells$SIG))
     rownames(res.cellex3)<-cells$cells
     res.cellex3<-rbind(c(0,0),res.cellex3)
     rownames(res.cellex3)[1]<-"group"
     res.cellex4<-t(res.cellex3)
     
     #res.cellex4$group<-"cellex"
     
     
     plot1<-ggradar(res.cellex4,grid.max=5,values.radar = c("0","1","5"),plot.title = "cellEx sig",axis.label.size = 4)}else{plot1<-ggtext("NOTHING")}
     
     #plot2##################
     
     query<-paste("MATCH (n:wgcna {square:'",square,"',name:'",my.module,"',edge:'",my.edge,"'})-[r]-(x) WHERE (x:cellprop) RETURN DISTINCT n.name AS module,toFloat(r.qvalue) as qvalue,toFloat(r.weight) AS pr,toFloat(r.Rsq) AS Rsq,x.name as name, labels(x) as kind",sep="")
     anq<-"MATCH (c:cellprop)-[r]-(c2:CELL) RETURN DISTINCT c.name as cells,c2.name AS group"
     cells<-cypher(graph,anq)
     cells$SIG<-rep(0,nrow(cells))
     
     res.cellprop<-cypher(graph,query)
     #print("debug2")
     print(res.cellprop)
     print(is.null(res.cellprop))
     if(!is.null(res.cellprop)){
     res.cellprop$SIG<--log10(res.cellprop$qvalue)
     res.cellprop$SIG[which(res.cellprop$SIG>10)]<-10
     res.cellprop2<-data.frame(as.numeric(res.cellprop$SIG))
     rownames(res.cellprop2)<-res.cellprop$name
     colnames(res.cellprop2)<-"sig"
     for (row in 1:nrow(res.cellprop2)){
       cells$SIG[which(cells$cells==res.cellprop$name[row])]<-res.cellprop2$sig[row]
       
     }
     res.cellprop3<-data.frame(as.numeric(cells$SIG))
     rownames(res.cellprop3)<-cells$cells
     res.cellprop3<-rbind(c(0,0),res.cellprop3)
     rownames(res.cellprop3)[1]<-"group"
     
     res.cellprop4<-t(res.cellprop3)
     print(res.cellprop4)
     plot2a<-ggradar(res.cellprop4,grid.max=10,values.radar = c("0","1","10"),plot.title = "cellprop sig",axis.label.size = 4)
     print("debug3")
     res.cellprop3b<-res.cellprop[which(res.cellprop$pr!=0),]
     res.cellprop3c<-res.cellprop3b[order(res.cellprop3b$pr),]
     res.cellprop3c<-rbind(c(0,0),res.cellprop3c)
     rownames(res.cellprop3c)[1]<-"group"
     
     
     res.cellprop4<-t(res.cellprop3c$pr)
     res.cellprop4<-data.frame(res.cellprop4)
     colnames(res.cellprop4)<-res.cellprop3c$name
     print(res.cellprop4)
     if(ncol(res.cellprop4)>1){
     plot2<-ggradar(res.cellprop4,grid.max=1,grid.min=-1, grid.mid=0,values.radar = c("-1","0","1"),plot.title = "cellprop cor",axis.label.size = 4)}else{plot2<-ggtext("NOTHING")}
     
     }else{
       plot2<-ggtext("NOTHING")
       plot2a<-ggtext("NOTHING")
       }
     
     
     #plot3##################
     
     
     query<-paste("MATCH (n:wgcna {square:'",square,"',name:'",my.module,"',edge:'",my.edge,"'})-[r]-(x) WHERE (x:pheno) RETURN DISTINCT n.name AS module,toFloat(r.qvalue) as qvalue,toFloat(r.weight) as pr,toFloat(r.Rsq) AS Rsq,x.name as name, labels(x) as kind",sep="")
     
     anq<-"MATCH (c:pheno) RETURN DISTINCT c.name as phenos"
     phenos<-cypher(graph,anq)
     
     phenos$pr<-rep(0,nrow(phenos))
     phenos$qvalue<-rep(0,nrow(phenos))
     
     res.pheno<-cypher(graph,query)
     print("debug4")
     print(is.null(res.pheno))
     if(!is.null(res.pheno)){
       print(res.pheno)
     res.pheno$SIG<--log10(res.pheno$qvalue)
     res.pheno$SIG[which(res.pheno$SIG>10)]<-10
     res.pheno2<-data.frame(as.numeric(res.pheno$SIG))#
     rownames(res.pheno2)<-res.pheno$name
     colnames(res.pheno2)<-"qvalue"
     for (row in 1:nrow(res.pheno2)){
       phenos$qvalue[which(phenos$phenos==res.pheno$name[row])]<-res.pheno2$qvalue[row]
       
     }
     res.pheno3<-data.frame(as.numeric(phenos$qvalue))
     res.pheno3$phenos<-phenos$phenos
     rownames(res.pheno3)<-phenos$phenos
     colnames(res.pheno3)<-c("qvalue","phenos")
     
     res.pheno3b<-res.pheno3[which(res.pheno3$qvalue!=0),]
     res.pheno3c<-res.pheno3b[order(res.pheno3b$qvalue),]
     qorder<-order(res.pheno3b$qvalue)
     res.pheno3c<-rbind(c(0,0),res.pheno3c)
     rownames(res.pheno3c)[1]<-"group"
     
     
     res.pheno4<-t(res.pheno3c$qvalue)
     #res.pheno4<-data.frame(as.numeric(res.pheno4))
     colnames(res.pheno4)<-rownames(res.pheno3c)
     print(res.pheno4)
     print(class(res.pheno4))
     res.pheno4<-data.frame(res.pheno4)
     print(res.pheno4[,1])
     
     if(ncol(res.pheno4)>1){
     plot3a<-ggradar(res.pheno4,grid.max=10,grid.min=0, grid.mid=5,values.radar = c("0","5","10"),plot.title = "pheno sig",axis.label.size = 4)}else{plot3a<-ggtext("NOTHING")}
     
     }else{plot3a<-ggtext("NOTHING")}
     
     
     ###
     
     
     
     print("debug5")
     #res.pheno<-cypher(graph,query)
     
     print(is.null(res.pheno))
     if(!is.null(res.pheno)){
     res.pheno2<-data.frame(as.numeric(res.pheno$pr))#
     rownames(res.pheno2)<-res.pheno$name
     colnames(res.pheno2)<-"pr"
     for (row in 1:nrow(res.pheno2)){
       phenos$pr[which(phenos$phenos==res.pheno$name[row])]<-res.pheno2$pr[row]
       
     }
     res.pheno3<-data.frame(as.numeric(phenos$pr))
     res.pheno3$phenos<-phenos$phenos
     rownames(res.pheno3)<-phenos$phenos
     colnames(res.pheno3)<-c("pr","phenos")
     
     res.pheno3b<-res.pheno3[which(res.pheno3$pr!=0),]
     #res.pheno3c<-res.pheno3b[order(res.pheno3b$qvalue),]
     res.pheno3c<-res.pheno3b[qorder,]
     
     res.pheno3c<-rbind(c(0,0),res.pheno3c)
     rownames(res.pheno3c)[1]<-"group"
     
     res.pheno4<-t(res.pheno3c$pr)
     
     colnames(res.pheno4)<-rownames(res.pheno3c)
     print(res.pheno4)
     res.pheno4<-data.frame(res.pheno4)
     print(res.pheno4[,1])
     
     if(ncol(res.pheno4)>1){
     plot3<-ggradar(res.pheno4,grid.max=1,grid.min=-1, grid.mid=0,values.radar = c("-1","0","1"),plot.title = "pheno cor",axis.label.size = 4)}else{plot3<-ggtext("NOTHING")}
     
     
     }else{plot3<-ggtext("NOTHING")}
     
     
     
     #plot4##################
     
     query<-paste("MATCH (n:wgcna {square:'",square,"',name:'",my.module,"',edge:'",my.edge,"'})-[r]-(x) WHERE (x:ImmunePW) OR (x:PalWangPW) OR (x:reactomePW) RETURN DISTINCT n.name AS module,toFloat(r.qvalue) as qvalue ,x.name as name, labels(x) as kind",sep="")
     anq<-"MATCH (x) WHERE (x:ImmunePW) OR (x:PalWangPW) OR (x:reactomePW) RETURN DISTINCT x.name as pathways"
     pathways<-cypher(graph,anq)
     pathways$SIG<-rep(0,nrow(pathways))
     
     print("debug6")
     
     res.pw<-cypher(graph,query)
     print(is.null(res.pw))
     if(!is.null(res.pw)){
     
     res.pw$SIG<--log10(res.pw$qvalue)
     res.pw$SIG[which(res.pw$SIG>10)]<-10
     res.pw2<-data.frame(as.numeric(res.pw$SIG))
     rownames(res.pw2)<-res.pw$name
     colnames(res.pw2)<-"SIG"
     for (row in 1:nrow(res.pw2)){
       pathways$SIG[which(pathways$pathways==res.pw$name[row])]<-res.pw2$SIG[row]
       
     }
     res.pw3<-data.frame(as.numeric(pathways$SIG))
     res.pw3$pathways<-pathways$pathways
     rownames(res.pw3)<-pathways$pathways
     colnames(res.pw3)<-c("SIG","pathways")
     
     res.pw3b<-res.pw3[which(res.pw3$SIG!=0),]
     res.pw3c<-res.pw3b[order(res.pw3b$SIG),]
     res.pw3c<-rbind(c(0,0),res.pw3c)
     rownames(res.pw3c)[1]<-"group"
     
     res.pw4<-t(res.pw3c$SIG)
     
     colnames(res.pw4)<-rownames(res.pw3c)
     print("res.pw4")
     print(res.pw4)
     plot4<-ggradar(res.pw4,grid.max=10,grid.min=0, grid.mid=5,values.radar = c("0","5","10"),plot.title = "pathway sig",axis.label.size = 2)}else{plot4<-ggtext("NOTHING")}
     
     #####
     grid.arrange(plot1,plot2a,plot2,plot3a,plot3,plot4,nrow=2)
     
     
   })
   
   
   
   ######end radar
   
   #modcor####
   output$modcor<-renderPlot({
     print("START MODCOR PLOT")
     reacstudyMC<-as.numeric(studyModcorfun())
     #reacsetMC<-as.numeric(squareModcorfun()[reacstudyMC]) debugging
     reacsetMC<-as.numeric(squareModcorfun())
     reacedgeMC<-as.numeric(edgeModcorfun())
     
     print("reacstudyMC")
     print(reacstudyMC)
     print("reacsetMC")
     print(reacsetMC)
     print("reacedgeMC")
     print(reacedgeMC)
     
     studyMC<-names(studies)[reacstudyMC]
     setlistMC<-setlistnameslist[[reacstudyMC]]
     squareMC<-setlistMC[[reacsetMC]]
     
     print(studyMC)
     print(squareMC)
     print(moduleModcorfun())
     print(edgeModcorfun())
     mwat(studyMC,squareMC,moduleModcorfun(),edgeModcorfun())
     if(input$returnpdf==TRUE){
       pdf("plot.pdf",width=18,height=14,onefile=FALSE)
       mwat(studyMC,squareMC,moduleModcorfun(),edgeModcorfun())
       dev.off()
     }
     
   })
   
   #MEvar####
   output$MEvar<-renderPlot({
     print("begin MEvar")
     squarelist<-unlist(strsplit(setfun(),"_"))
     if(length(squarelist)==3){
       square<-squarelist[3]
     }else{square<-paste(squarelist[c(3,4)],collapse="_")}
     #for each ME
     
     #get ME vector
     #prefix<-paste("/home/rstudio/output/build/",buildfun(),"/",sep="")
     prefix<-"/home/rstudio/output/build/"
     data0<-read.csv(paste(prefix,setfun(),"/5_WGCNA_4k_tv/results/ModuleEigengenes_all_samples_withPheno.csv",sep=""))
     data<-read.csv(paste(prefix,setfun(),"/5_WGCNA_4k_tv/results/ModuleEigengenes_all_samples.csv",sep=""))
     #print("data is read")
     
     #get class subsets
     menames<-sort(colnames(data)[2:length(colnames(data))])
     
     contrasts<-cypher(graph,paste("MATCH (n:wgcna {square:'",square,"'}) RETURN DISTINCT n.contrastvar AS contrastvars",sep=""))
     c1<-contrasts$contrastvars[1]
     c2<-contrasts$contrastvars[2]
     
     data0<-data0[order(data0[,c1]),][order(data0[,c2]),]
     
     classes<-cypher(graph,paste("MATCH (n:wgcna {square:'",square,"'}) RETURN DISTINCT n.contrast AS contrasts",sep=""))
     
     #print(classes$contrasts[2])
     cL1<-unlist(strsplit(classes$contrasts[1]," - "))
     cL2<-unlist(strsplit(classes$contrasts[2]," - "))
     
     subset1vec<-eval(parse(text=paste("data0$",c1,sep="")))
     subset1classes<-sort(unique(subset1vec))
     axis1names<-sort(c(cL1[1],cL1[2]))
     print(axis1names)
     characterclasses1<-as.factor(subset1vec)
     levels(characterclasses1)<-axis1names
     print(as.character(characterclasses1))
     
     subset2vec<-eval(parse(text=paste("data0$",c2,sep="")))
     subset2classes<-sort(unique(subset2vec))
     axis2names<-sort(c(cL2[1],cL2[2]))
     characterclasses2<-as.factor(subset2vec)
     levels(characterclasses2)<-axis2names
     print(as.character(characterclasses2))
     
     #dev.new()
     sub<-which(subset1vec%in%s1fun()&subset2vec%in%s2fun())
     group1<-subset1vec[sub]#by TB
     group2<-subset2vec[sub]#by HIV
     #split in two or four
     
     print("group1")
     print(group1)
     print("group2")
     print(group2)
     MEvar<-me1fun()
     print("MEvar")
     print(MEvar)
     print(data[,MEvar])
     # #split
     # a1<-which(group1==1&group2==1)
     # a2<-which(group1==1&group2==2)
     # b1<-which(group1==2&group2==1)
     # b2<-which(group1==2&group2==2)
     #split
     a1<-which(characterclasses1==cL1[2]&characterclasses2==cL2[2])
     a2<-which(characterclasses1==cL1[1]&characterclasses2==cL2[2])
     b1<-which(characterclasses1==cL1[2]&characterclasses2==cL2[1])
     b2<-which(characterclasses1==cL1[1]&characterclasses2==cL2[1])
     print("subsetvectors")
     print(subset1vec)
     print(subset2vec)
     colvec<-character(length=nrow(data0))
     colvec[a1]<-"blue"
     colvec[a2]<-"green"
     colvec[b1]<-"yellow"
     colvec[b2]<-"red"
     
     #vectors
     print(paste(cL1[2],cL2[2]))
     print(var(data0[a1,MEvar]))
     print(paste(cL1[2],cL2[1]))
     print(var(data0[a2,MEvar]))
     print(paste(cL1[1],cL2[2]))
     print(var(data0[b1,MEvar]))
     print(paste(cL1[1],cL2[1]))
     print(var(data0[b2,MEvar]))
     
     par(mfrow=c(2,2))
     barplot(data0[b1,MEvar],col="yellow",main=paste(cL1[2],cL2[1],round(mad(data0[b1,MEvar]),3)),ylim=c(-max(abs(data0[,MEvar])),max(abs(data0[,MEvar]))))
     barplot(data0[b2,MEvar],col="red",main=paste(cL1[1],cL2[1],round(mad(data0[b2,MEvar]),3)),ylim=c(-max(abs(data0[,MEvar])),max(abs(data0[,MEvar]))))
     barplot(data0[a1,MEvar],col="blue",main=paste(cL1[2],cL2[2],round(mad(data0[a1,MEvar]),3)),ylim=c(-max(abs(data0[,MEvar])),max(abs(data0[,MEvar]))))
     barplot(data0[a2,MEvar],col="green",main=paste(cL1[1],cL2[2],round(mad(data0[a2,MEvar]),3)),ylim=c(-max(abs(data0[,MEvar])),max(abs(data0[,MEvar]))))
     
     #calculate variability
     
     #test variability to get at heterogeneity
     
     #look for drivers?
     
     
     
     
   })
   
   #modPheno####
   output$modPheno <- renderPlot({
      #new code
      #square
      #square<-set()
      squarelist<-unlist(strsplit(setfun(),"_"))
      if(length(squarelist)==3){
         square<-squarelist[3]
      }else{square<-paste(squarelist[c(3,4)],collapse="_")}
      #pheno variabe
      pheno<-phenofun()
      #module
      module<-me1fun()
      #module<-str_split(module,"ME")[[1]][2]
      module<-strsplit(module,"ME")[[1]][2]
      query<-paste("MATCH (p:personPheno {square:'",square,"',name:'",pheno,"'})-[r]-(p0:person)-[r1]-(p1:personME {name:'",module,"'}) RETURN p0.name AS id,p0.class1 as class1, p0.class2 AS class2, toFloat(p.value) as var, toFloat(p1.value) as ME",sep="")
      print(query)
      res<-cypher(graph,query)
      subset1vec<-factor(res$class1)
      subset1classes<-levels(subset1vec)
      subset2vec<-factor(res$class2)
      subset2classes<-levels(subset2vec)
      print("part 1 complete")
      
      
     # #prefix<-paste("/home/rstudio/output/build/",buildfun(),"/",sep="")
     # prefix<-"/home/rstudio/output/build/"
     # #print("line1 inside reactive")
     # data0<-read.csv(paste(prefix,setfun(),"/5_WGCNA_4k_tv/results/ModuleEigengenes_all_samples_withPheno.csv",sep=""))
     # 
     # #print("line2")
     # data<-read.csv(paste(prefix,setfun(),"/5_WGCNA_4k_tv/results/ModuleEigengenes_all_samples.csv",sep=""))
     # #print("line3")
     # squarelist<-unlist(strsplit(setfun(),"_"))
     # if(length(squarelist)==3){
     #   square<-squarelist[3]
     # }else{square<-paste(squarelist[c(3,4)],collapse="_")}
     # #print("square")
     # #print(square)
     # menames<-sort(colnames(data)[2:length(colnames(data))])
     # 
     # contrasts<-cypher(graph,paste("MATCH (n:wgcna {square:'",square,"'}) RETURN DISTINCT n.contrastvar AS contrastvars",sep=""))
     # c1<-contrasts$contrastvars[1]
     # c2<-contrasts$contrastvars[2]
     # 
     # data0<-data0[order(data0[,c1]),][order(data0[,c2]),]
     # 
     # classes<-cypher(graph,paste("MATCH (n:wgcna {square:'",square,"'}) RETURN DISTINCT n.contrast AS contrasts",sep=""))
     # #cL1<-unlist(strsplit(classes$contrasts[1],"-"))[c(2,1)]
     # #cL2<-unlist(strsplit(classes$contrasts[2],"-"))[c(2,1)]
     # #print("debug")
     # #print(classes$contrasts[1])
     # #print(classes$contrasts[2])
     # cL1<-unlist(strsplit(classes$contrasts[1]," - "))
     # cL2<-unlist(strsplit(classes$contrasts[2]," - "))
     # 
     # subset1vec<-eval(parse(text=paste("data0$",c1,sep="")))
     # subset1classes<-sort(unique(subset1vec))
     # 
     # subset2vec<-eval(parse(text=paste("data0$",c2,sep="")))
     # subset2classes<-sort(unique(subset2vec))
     # 
     # #dev.new()
     # sub<-which(subset1vec%in%s1fun()&subset2vec%in%s2fun())
     # group1<-subset1vec[sub]
     # group2<-subset2vec[sub]
     # #print("doing res")
     # #print(data0[sub,])
     # themethod<-methodfun()
     # res<-cor.test(eval(parse(text=paste("data0$",phenofun(),"[sub]",sep=""))),eval(parse(text=paste("data0$",me1fun(),"[sub]",sep=""))),method=themethod)
      
      themethod<-methodfun()
      
      #sub<-which(subset1vec%in%levels(subset1vec)[numeric(s1fun())]&subset2vec%in%levels(subset2vec)[numeric(s2fun())])
      s1sel<-as.numeric(s1fun())
      s2sel<-as.numeric(s2fun())
      toselect1<-levels(subset1vec)[s1sel]
      toselect2<-levels(subset2vec)[s2sel]
      print("s1sel")
      print(s1sel)
      print("s2sel")
      print(s2sel)
      print("toselect1")
      print(toselect1)
      print("toselect2")
      print(toselect2)
      
      sub<-which(subset1vec%in%toselect1&subset2vec%in%toselect2)
      
      print("subset1vec")
      print(subset1vec)
      print("subset2vec")
      print(subset2vec)
      print("s1fun")
      print(s1fun())
      print("s2fun")
      print(s2fun())
      
      print("sub")
      print(sub)
      print(res$var)
      print(res$ME)
   
     res.cor<-cor.test(res$var[sub],res$ME[sub],method=themethod)
     print("doing plot")
     
     #print(group1)
     #print(group2)
     #print(2*group1+group2)
     #print(2*data0[sub,c1]+data0[sub,c2])
     #print(paste(cL1[group1],cL2[group2]))
     
     #print(unique(2*group1+group2))
     #print(unique(paste(cL1[group1],cL2[group2])))
     #print(c(paste(cL1[1],cL2[1]),paste(cL1[1],cL2[2]),paste(cL1[2],cL2[1]),paste(cL1[2],cL2[2])))
     # plot(eval(parse(text=paste(phenofun(),"~",me1fun(),sep=""))),data=data0[sub,],pch=20,cex=2,col=unlist(strsplit(me1fun(),"ME"))[2],main=paste(me1fun(),"\n",themethod,"coef.=",sprintf("%.3f",res$estimate),"Pval=",sprintf("%.3f",res$p.value)),xlab=me1fun(),ylab=phenofun())
     # text(data0[sub,phenofun()]~data0[sub,me1fun()],labels=data0[sub,"Row.names"],col=(2*data0[sub,c1]+data0[sub,c2])-2,pos=2)
     # legend(legposfun(),sort(c(paste(cL1[1],cL2[1]),paste(cL1[1],cL2[2]),paste(cL1[2],cL2[1]),paste(cL1[2],cL2[2]))),fill=c(3,4,5,6)-2)
     # #legend("bottomleft",c(paste(cL1[1],cL2[1]),paste(cL1[1],cL2[2]),paste(cL1[2],cL2[1]),paste(cL1[2],cL2[2])),fill=c(3,4,5,6))
     # abline(lm(eval(parse(text=paste(phenofun(),"~",me1fun(),sep=""))),data=data0[sub,]), col="red")
     # 
     
     # if(input$returnpdf==TRUE){
     #   pdf("plot.pdf",width=9,height=9,onefile=FALSE)
     #   plot(eval(parse(text=paste(phenofun(),"~",me1fun(),sep=""))),data=data0[sub,],pch=20,cex=2,col=unlist(strsplit(me1fun(),"ME"))[2],main=paste(me1fun(),"\nPearsonR=",sprintf("%.3f",res$estimate),"Pval=",sprintf("%.3f",res$p.value)),xlab=me1fun(),ylab=phenofun())
     #   text(data0[sub,phenofun()]~data0[sub,me1fun()],labels=data0[sub,"Row.names"],col=(2*data0[sub,c1]+data0[sub,c2])-1,pos=2)
     #   legend(legposfun(),sort(c(paste(cL1[1],cL2[1]),paste(cL1[1],cL2[2]),paste(cL1[2],cL2[1]),paste(cL1[2],cL2[2]))),fill=c(3,4,5,6)-1)
     #   abline(lm(eval(parse(text=paste(phenofun(),"~",me1fun(),sep=""))),data=data0[sub,]), col="red")
     #   
     #   dev.off()
     # }
     #NEW
     # class1text<-character()
     # for (xx in 1:length(subset1vec)){
     #   class1text[xx]<-rev(cL1)[subset1vec[xx]]
     # }
     # class2text<-character()
     # for (xx in 1:length(subset2vec)){
     #   class2text[xx]<-rev(cL2)[subset2vec[xx]]
     # }
     # print(class1text)
     # print(class2text)
     # 
     # data0$classes<-interaction(class1text,class2text)
     # 
     
     # spdata<-data0[sub,]
     # print(data0$classes)
     res$classes<-interaction(res$class1,res$class2)
     #spdata<-res[sub,]
     spdata<-res[sub,]
     
     #sp<-ggplot(spdata,aes_string(ME,var))+
       sp<-ggplot(spdata,aes(ME,var))+
       geom_point(aes(col=classes),size=12,alpha=0.4)+
       geom_smooth(method=lm) + 
       geom_rug() +
       #geom_text_repel(aes(label=Row.names),data=data0[sub,],size=6) +
       geom_text_repel(aes(label=id),data=spdata,size=6) +
       #theme_tufte() +
       theme_light() +
       #theme_minimal() +
       labs(y = phenofun()) +
       #labs(x = me1fun()) +
       labs(x = module) +
       #labs(colour = paste("diffME edge:",edgefun())) +
       theme(axis.title.y = element_text(size = rel(2.5), angle = 90, color = "slategrey",margin = margin(t = 0, r = 50, b = 50, l = 50))) +
       theme(axis.title.x = element_text(size = rel(2.5), color = "slategrey",margin = margin(t = 50, r = 50, b = 0, l = 50))) + 
       theme(axis.text = element_text(size = rel(1.5))) +
       theme(aspect.ratio=1) +
       theme(legend.text = element_text(size = rel(1.2))) +
       #theme(legend.title = element_text(size = rel(1.2))) +
       #theme(legend.text=element_text(size=rel(1.4))) +
       #theme(legend.key.size = unit(1,"cm")) +
       #labs(y = phenofun2()) +
       ggtitle("Module eigengene/ phenotype correlation") +
       stat_cor(method = themethod,size=8,label.y.npc = "bottom",label.x.npc = "left",color="red") +
       theme(plot.title = element_text(size=rel(2.5),hjust = 0.5,color="maroon"))
     print(sp)
   #})
     
   })
   
   #subjectplot####
   output$subjectplot <- renderPlot({
     print("computing subjectplot")
     #prefix<-paste("/home/rstudio/output/build/",buildfun(),"/",sep="")
     prefix<-"/home/rstudio/output/build/"
     #print("line1 inside reactive")
     data0<-read.csv(paste(prefix,setfun(),"/5_WGCNA_4k_tv/results/ModuleEigengenes_all_samples_withPheno.csv",sep=""))
     #print("line2")
     data<-read.csv(paste(prefix,setfun(),"/5_WGCNA_4k_tv/results/ModuleEigengenes_all_samples.csv",sep=""),row.names = 1)
     #print("line3")
     squarelist<-unlist(strsplit(setfun(),"_"))
     print("setfun:")
     print(setfun())
     print("squarelist:")
     print(squarelist)
     if(length(squarelist)==3){
       square<-squarelist[3]
     }else{square<-paste(squarelist[c(3,4)],collapse="_")}
     print("square")
     print(square)
     menames<-sort(colnames(data)[2:length(colnames(data))])
     
     contrasts<-cypher(graph,paste("MATCH (n:wgcna {square:'",square,"'}) RETURN DISTINCT n.contrastvar AS contrastvars",sep=""))
     c1<-contrasts$contrastvars[1]
     c2<-contrasts$contrastvars[2]
     
     print("debug contrasts$contrastvars")
     print(c1)
     print(c2)
     
     classes<-cypher(graph,paste("MATCH (n:wgcna {square:'",square,"'}) RETURN DISTINCT n.contrast AS contrasts",sep=""))
     
     print("debug classes$contrasts")
     print(classes$contrasts[1])
     print(classes$contrasts[2])
     
     #cL1<-unlist(strsplit(classes$contrasts[1],"-"))[c(2,1)]
     #cL2<-unlist(strsplit(classes$contrasts[2],"-"))[c(2,1)]
     
     cL1<-unlist(strsplit(classes$contrasts[1]," - "))
     cL2<-unlist(strsplit(classes$contrasts[2]," - "))
     
     print(cL1)
     print(cL2)
     
     subset1vec<-eval(parse(text=paste("data0$",c1,sep="")))
     subset1classes<-sort(unique(subset1vec))
     subset2vec<-eval(parse(text=paste("data0$",c2,sep="")))
     subset2classes<-sort(unique(subset2vec))
     
     subset1fac<-factor(subset1vec,levels=c("1","2"),labels=sort(cL1))
     subset2fac<-factor(subset2vec,levels=c("1","2"),labels=sort(cL2))
     
     #dev.new()
     sub<-which(subset1vec%in%s1fun()&subset2vec%in%s2fun())
     group1<-subset1vec[sub]
     group2<-subset2vec[sub]
     #print(data0[sub,])
  #   plot(eval(parse(text=paste(me2fun(),"~",me1fun(),sep=""))),data=data0[sub,],pch=20,col="grey",main=paste("Subjects:",square),ylab=me2fun(),xlab=me1fun())
     #text(eval(parse(text=paste(me2fun(),"~",me1fun(),sep=""))),data=data0[sub,],labels=paste(cL1[group1],cL2[group2]),col=group1+2*group2)
  #   text(eval(parse(text=paste(me2fun(),"~",me1fun(),sep=""))),data=data0[sub,],labels=data0[sub,"Row.names"],col=(group1+2*group2)-2,pos=2)
  #   legend("bottomleft",sort(c(paste(cL1[1],cL2[1]),paste(cL1[1],cL2[2]),paste(cL1[2],cL2[1]),paste(cL1[2],cL2[2]))),fill=c(3,5,4,6)-2)
     #legend("bottomleft",sort(c(paste(cL1[1],cL2[1]),paste(cL1[2],cL2[1]),paste(cL1[1],cL2[2]),paste(cL1[2],cL2[2]))),fill=c(3,4,5,6))
     #legend("bottomleft",c(paste(cL1[1],cL2[1]),paste(cL1[1],cL2[2]),paste(cL1[2],cL2[1]),paste(cL1[2],cL2[2])),fill=c(3,4,5,6))
     
  #   abline(h=0.0,col="blue")
  #   abline(v=0.0,col="blue")
     # if(input$returnpdf==TRUE){
     #   pdf("plot.pdf",width=16,height=10,onefile=FALSE)
     #   plot(eval(parse(text=paste(me2fun(),"~",me1fun(),sep=""))),data=data0[sub,],pch=20,col="grey",main=paste("Subjects:",square),ylab=me2fun(),xlab=me1fun())
     #   text(eval(parse(text=paste(me2fun(),"~",me1fun(),sep=""))),data=data0[sub,],labels=data0[sub,"Row.names"],col=(group1+2*group2)-2,pos=2)
     #   legend("bottomleft",sort(c(paste(cL1[1],cL2[1]),paste(cL1[1],cL2[2]),paste(cL1[2],cL2[1]),paste(cL1[2],cL2[2]))),fill=c(3,5,4,6)-2)
     #   dev.off()
     # }
     
     print("debug new")
     print(data0[sub,])
     print(str(data0[sub,]))
     print("me2fun()")
     print(me2fun())
     print("me1fun()")
     print(me1fun())
     print("subset1vec")
     print(subset1vec)
     print(class(subset1vec))
     
     class1text<-character()
     for (xx in 1:length(subset1vec)){
       class1text[xx]<-rev(cL1)[subset1vec[xx]]
       }
     class2text<-character()
     for (xx in 1:length(subset2vec)){
       class2text[xx]<-rev(cL2)[subset2vec[xx]]
     }
     print("classtext")
     print(class1text)
     print(class2text)
     print("subsetfac")
     print(subset1fac)
     print(subset2fac)
     
     
     #data0$classes<-interaction(class1text,class2text)
     data0$classes<-interaction(subset1fac,subset2fac)
     spdata<-data0[sub,]
     print(data0$classes)
     print(data0[1:12,])
     sp<-ggplot(spdata,aes_string(me2fun(),me1fun()))+
       geom_point(aes(col=classes),size=12,alpha=0.4)+
       
       geom_text_repel(aes(label=Row.names),data=data0[sub,],size=6) +
       #theme_tufte() +
       theme_light() +
       #theme_minimal() +
       labs(x = me2fun()) +
       labs(y = me1fun()) +
       #labs(colour = paste("diffME edge:",edgefun())) +
       theme(axis.title.y = element_text(size = rel(2.5), angle = 90, color = "slategrey",margin = margin(t = 0, r = 50, b = 50, l = 50))) +
       theme(axis.title.x = element_text(size = rel(2.5), color = "slategrey",margin = margin(t = 50, r = 50, b = 0, l = 50))) + 
       theme(axis.text = element_text(size = rel(1.5))) +
       theme(aspect.ratio=1) +
       theme(legend.text = element_text(size = rel(1.2))) +
       #theme(legend.title = element_text(size = rel(1.2))) +
       #theme(legend.text=element_text(size=rel(1.4))) +
       #theme(legend.key.size = unit(1,"cm")) +
       #labs(y = phenofun2()) +
       ggtitle("Subjectplot") +
       theme(plot.title = element_text(size=rel(2.5),hjust = 0.5,color="maroon"))
     print(sp)
   })
   
   #All modules####
   
   #module2d####
   output$module2d<-renderPlot({
     #prefix<-paste("/home/rstudio/output/build/",buildfun(),"/",sep="")
     prefix<-"/home/rstudio/output/build/"
     datarec1<-read.csv(paste(prefix,setfun(),"/5_WGCNA_4k_tv/results/modules_as_classifiers_AUC_glm.csv",sep=""))
     squarelist<-unlist(strsplit(setfun(),"_"))
     if(length(squarelist)==3){
       square<-squarelist[3]
     }else{square<-paste(squarelist[c(3,4)],collapse="_")}
     #print("square")
     #print(square)
     contrast1<-cypher(graph,paste("MATCH (n:wgcna {square:'",square,"',edge:'1'}) RETURN DISTINCT n.contrastvar AS contrastvars",sep=""))
     contrast2<-cypher(graph,paste("MATCH (n:wgcna {square:'",square,"',edge:'3'}) RETURN DISTINCT n.contrastvar AS contrastvars",sep=""))
     
     print(contrast1)
     print(contrast2)
     c1<-contrast1$contrastvars[1]
     c2<-contrast2$contrastvars[1]
     print("c1")
     print(c1)
     print("c2")
     print(c2)
     query<-paste("MATCH (n:wgcna {square:'",square,"',edge:'",edgefun(),"'}) RETURN n.name as name, n.sigenrich as sigenrich, n.diffME as diffME, n.modAUC1 as modAUC1, n.modAUC2 as modAUC2, n.AIC1 as AIC1, n.AIC2 as AIC2, n.P_1 as P_1, n.P_2 as P_2, n.logit_est1 as logit_est1, n.logit_est2 as logit_est2",sep="")
     res<-cypher(graph,query)
     print(res)
     res2<-as.data.frame(lapply(res[,2:ncol(res)],as.numeric))
     res2$name<-res$name
     print(head(res2))
     print(class(res2$diffME))
     # plot(modAUC1~modAUC2,data=datarec1,pch=20,asp=1,col="grey",main=paste(square,": Modules differentiating ",c1," and ",c2,sep=""),xlab=paste(c2,"index"),ylab=paste(c1, "index"),height=900)#aspect ratio can be fixed with asp=1
     # shadowtext(datarec1$modAUC1~datarec1$modAUC2,labels=datarec1$X,col=datarec1$X,font=2,cex=1.0,xpd=TRUE,bg="lightgrey")
     # abline(h=0.75,col="blue")
     # abline(v=0.75,col="blue")
     
     #vp<-ggplot(res,aes(modAUC2,modAUC1))
     
     #######B
     vp<-ggplot(res2,aes(modAUC2,modAUC1))+
       geom_point(aes(col=as.numeric(diffME)),size=12,alpha=0.7)+
       #geom_point(size=12,alpha=0.4)+
       scale_color_gradient2(low="blue",mid="lightyellow",high="red") +
       #scale_color_gradient2() +
       geom_text_repel(aes(label=name),data=res2,size=7) +
       #theme_tufte() +
       #theme_light() +
       theme_minimal() +
       labs(x = paste(c2,"AUC")) +
       labs(y = paste(c1,"AUC")) +
       labs(colour = paste("diffME edge:",edgefun())) +
       theme(axis.title.y = element_text(size = rel(2.5), angle = 90, color = "skyblue4",margin = margin(t = 0, r = 50, b = 50, l = 50))) +
       theme(axis.title.x = element_text(size = rel(2.5), color = "skyblue4",margin = margin(t = 50, r = 50, b = 0, l = 50))) +
       theme(axis.text = element_text(size = rel(1.5))) +
       theme(aspect.ratio=1) +
       theme(legend.text = element_text(size = rel(1.2))) +
       theme(legend.title = element_text(size = rel(1.2))) +
       theme(legend.text=element_text(size=rel(1.4))) +
       theme(legend.key.size = unit(1,"cm")) +
       #labs(y = phenofun2()) +
       ggtitle(paste(square,": Modules differentiating ",c1," and ",c2,sep="")) +
       theme(plot.title = element_text(size=rel(2.5),hjust = 0.5,color="maroon"))
     #######E
     
     

     
     ##
     output$dtwgcna<-DT::renderDataTable({
         print("doing the data table for wgcna")
         datatable(res)%>%formatSignif(names(res)[2:ncol(res)],3)},server = FALSE,options=list(lengthMenu = c(5, 10, 15,20,50), pageLength = 50),filter="top")
     ##
     output$moduleTable<-renderText({
        query1<-paste("MATCH (n:wgcna {square:'",square,"',edge:'",edgefun(),"'}) RETURN n.name AS module",sep="")
        modules<-cypher(graph,query1)
        
        df <- data.frame(matrix(ncol = 6, nrow = 0))
        colnames(df)<-c("module","sigenrich","modAUC1","modAUC2","entity","qvalue")
        for (module in modules$module){
           for (entity in c("x:cellEx","x:cellprop","x:ImmunePW OR x:PalWangPW OR x:reactomePW")){
              query<-paste("MATCH (n:wgcna {square:'",square,"',edge:'",edgefun(),"',name:'",module,"'})
          WITH n 
          MATCH (n)-[r1]-(x) WHERE (",entity,")
          WITH * ORDER BY toFloat(r1.qvalue) LIMIT 3 
          RETURN DISTINCT n.name AS module, toFloat(n.diffME) as diffME, toFloat(n.sigenrich) as sigenrich, toFloat(n.modAUC1) as modAUC1, toFloat(n.modAUC2) as modAUC2, type(r1) AS reltype, labels(x) AS type, x.name AS entity, toFloat(r1.qvalue) as qv1, toFloat(r1.weight) AS PearsonR",sep="")
              
              res<-cypher(graph,query)
              df<-rbind(df,res)
           }
           
        }
        
        hx<-as_hux(df)
        #hx2<-merge_repeated_rows(hx,everywhere,c("module", "diffME", "sigenrich", "modAUC1", "modAUC2"))
        
        width(hx) <- 0.4
        #wrap(hx) <- TRUE
        hx<-set_wrap(hx, everywhere, "entity", TRUE )
        for (module in modules$module) {
           hx<-merge_repeated_rows(hx,hx$module==module,c("module", "diffME", "sigenrich", "modAUC1", "modAUC2"))
        }
        
        
        #hx2<-merge_repeated_rows(hx,everywhere,c("module"))
        hx3<-theme_article(hx)
        #hx3<-theme_plain(hx2)
        #hx3<-theme_mondrian(hx2)
        hx4<-set_background_color(hx3,evens, 6:8,"grey95")
        
        hx5<-hx4 %>% map_text_color(everywhere, c("diffME"),
                                    by_colorspace("blue","grey", "red"))
        
        hx6<-hx5 %>% map_text_color(everywhere, c("modAUC1","modAUC2"),
                                    by_colorspace("grey", "red"))
        
        
        #quick_html(hx6,file="/home/rstudio/output/project/tabular/moduleAnnotations.html")
        
        #htmlout<-read_html("/home/rstudio/output/project/tabular/moduleAnnotations.html")
        #htmlout<-readtext("/home/rstudio/output/project/tabular/moduleAnnotations.html")
        #htmlout$text
        print_notebook(hx6)
     })
     
     
     print(vp)
     
     # #Interesting, but not useful attempt to put wordclouds on module2d. extremely slow, result unreadable...
     # #plot(modAUC1~modAUC2,data=datarec1,pty = 's', type = 'n', xlim = c(0.5, 1), ylim = c(0.5, 1))
     # xx = grconvertX(x = datarec1$modAUC2, from = 'user', to = 'ndc')
     # print(xx)
     # #xx = datarec1$modAUC2
     # #print(xx)
     # yy = grconvertY(y = datarec1$modAUC1, from = 'user', to = 'ndc')
     # print(yy)
     # #yy = datarec1$modAUC1
     # #print(yy)
     # 
     # modules<-datarec1$X
     # plot(modAUC1~modAUC2,data=datarec1, pty = 's', type = 'p', xlim = c(0, 1), ylim = c(0, 1))
     # for (modu in 1:length(modules)){#debug length(modules)
     #   query<-paste("MATCH (n:wgcna {name:'",modules[modu],"'})-[r]-(x) WHERE (x:reactomePW OR x:PalWangPW OR x:ImmunePW OR x:cellEx OR x:pheno) return x.name AS annot",sep="")
     #   res<-cypher(graph,query)
     #   words<-c(res$annot,rep(modules[modu],5))
     #   print("done words")
     #   postscript("currentWC.ps",encoding="ISOLatin1")
     #   wordcloud(words,colors = modules[modu],min.freq = 3)
     #   dev.off()
     #   print("done ps")
     #   PostScriptTrace("currentWC.ps","currentWC.xml")
     #   print("done trace")
     #   currentWC<-readPicture("currentWC.xml")
     #   print("done read picture")
     #   #grid.symbols(currentWC, x = xx[modu], y = yy[modu],size=.3)
     #   #plot(modAUC1~modAUC2,data=datarec1, pty = 's', type = 'n', xlim = c(-1, 1), ylim = c(-1, 1))
     #   #xx = grconvertX(x = x, from = 'user', to = 'ndc')
     #   #yy = grconvertY(y = y, from = 'user', to = 'ndc')
     #   grid.picture(currentWC, x = xx[modu], y = yy[modu],width=.1,height=0.1)
     #  print("DONE insert") 
     # }
     
    
     if(input$returnpdf==TRUE){
       pdf("plot.pdf",width=12,height=7,onefile=FALSE)
       plot(modAUC1~modAUC2,data=datarec1,pch=20,col="grey",main=paste(square,": Modules differentiating ",c1," and ",c2,sep=""),xlab=paste(c2,"index"),ylab=paste(c1, "index"))
       shadowtext(datarec1$modAUC1~datarec1$modAUC2,labels=datarec1$X,col=datarec1$X,font=2,cex=1.0,xpd=TRUE,bg="lightgrey")
       abline(h=0.75,col="blue")
       abline(v=0.75,col="blue")
       dev.off()
     }
   })
   
   #moduleHM####
   
   output$moduleHM<-renderPlot({
     #get modules and diffME for square
     #prefix<-paste("/home/rstudio/output/build/",buildfun(),"/",sep="")
     prefix<-"/home/rstudio/output/build/"
     datarec1<-read.csv(paste(prefix,setfun(),"/5_WGCNA_4k_tv/results/modules_as_classifiers_AUC_glm.csv",sep=""))
     squarelist<-unlist(strsplit(setfun(),"_"))
     if(length(squarelist)==3){
       square<-squarelist[3]
     }else{square<-paste(squarelist[c(3,4)],collapse="_")}
     #print("square")
     #print(square)
     contrast1<-cypher(graph,paste("MATCH (n:wgcna {square:'",square,"',edge:1}) RETURN DISTINCT n.contrastvar AS contrastvars",sep=""))
     contrast2<-cypher(graph,paste("MATCH (n:wgcna {square:'",square,"',edge:3}) RETURN DISTINCT n.contrastvar AS contrastvars",sep=""))
     
     print(contrast1)
     print(contrast2)
     c1<-contrast1$contrastvars[1]
     c2<-contrast2$contrastvars[1]
     
     query<-paste("MATCH (n:wgcna {square:'",square,"'}) RETURN n.name as name, n.diffME as diffME, n.edge AS edge",sep="")
     res<-cypher(graph,query)
     print(res)
    
     hmmat<-matrix(nrow=5,ncol=(nrow(res)/5),data=NA)
     dimnames(hmmat)<-list(1:5,unique(res$name))
     print(hmmat)
     print(hmmat[,"black"])
     print("debug1")
     
     for (lines in 1:nrow(res)){
       ires<-res[lines,]
       #print(class(ires))
       #print(ires)
       hmmat[ires[[3]],as.character(ires[[1]])]<-as.numeric(ires[[2]])
       
     }
     datasetrow<-datasetmatrix[square,]
     rownames(hmmat)<-datasetrow
     cur_dev <- dev.cur()
     
     theheatmap<-pheatmap(hmmat,fontsize_row = 10 * factorfun(),fontsize_col = 10 * factorfun())
     dev.set(cur_dev)
     
     print(theheatmap)
     
     
   })
   
   #sigenrichHM####
   
   output$sigenrichHM<-renderPlot({
     #get modules and diffME for square
     #prefix<-paste("/home/rstudio/output/build/",buildfun(),"/",sep="")
     prefix<-"/home/rstudio/output/build/"
     datarec1<-read.csv(paste(prefix,setfun(),"/5_WGCNA_4k_tv/results/modules_as_classifiers_AUC_glm.csv",sep=""))
     squarelist<-unlist(strsplit(setfun(),"_"))
     if(length(squarelist)==3){
       square<-squarelist[3]
     }else{square<-paste(squarelist[c(3,4)],collapse="_")}
     #print("square")
     #print(square)
     contrast1<-cypher(graph,paste("MATCH (n:wgcna {square:'",square,"',edge:1}) RETURN DISTINCT n.contrastvar AS contrastvars",sep=""))
     contrast2<-cypher(graph,paste("MATCH (n:wgcna {square:'",square,"',edge:3}) RETURN DISTINCT n.contrastvar AS contrastvars",sep=""))
     
     print(contrast1)
     print(contrast2)
     c1<-contrast1$contrastvars[1]
     c2<-contrast2$contrastvars[1]
     
     query<-paste("MATCH (n:wgcna {square:'",square,"'}) RETURN n.name as name, n.sigenrich as sigenrich, n.edge AS edge",sep="")
     res<-cypher(graph,query)
     print(res)
     
     hmmat<-matrix(nrow=5,ncol=(nrow(res)/5),data=NA)
     dimnames(hmmat)<-list(1:5,unique(res$name))
     print(hmmat)
     print(hmmat[,"black"])
     print("debug1")
     
     for (lines in 1:nrow(res)){
       ires<-res[lines,]
       print(class(ires))
       print(ires)
       hmmat[ires[[3]],as.character(ires[[1]])]<-as.numeric(ires[[2]])
       
     }
     datasetrow<-datasetmatrix[square,]
     rownames(hmmat)<-datasetrow
     
     cur_dev <- dev.cur()
     
     #theheatmap<-pheatmap(hmmat,color = colorRampPalette(RColorBrewer::brewer.pal(n = 7, name = "Reds"))(100) ,scale = "none")
     theheatmap<-pheatmap(hmmat,color = colorRampPalette(c("white","yellow","orange","red","purple"))(100) ,scale = "none",fontsize_row = 10 * factorfun(),fontsize_col = 10 * factorfun())
     dev.set(cur_dev)
     
     print(theheatmap)
     
     
   })
   
   #intracor####
   if(doIntracor==TRUE){
   #
     vals1 <- reactiveValues() 
   intracorcalc<-reactive({  
     #ME intracorrelates
     #print("reading data")
     print("doing intracorcalc")
     squarelist<-unlist(strsplit(setfun2(),"_"))
     if(length(squarelist)==3){
       square<-squarelist[3]
     }else{square<-paste(squarelist[c(3,4)],collapse="_")}
     print("square")
     print(square)
     #prefix<-paste("/home/rstudio/output/build/",buildfun(),"/",sep="")
     prefix<-"/home/rstudio/output/build/"
     data0<-read.csv(paste(prefix,setfun2(),"/5_WGCNA_4k_tv/results/ModuleEigengenes_all_samples_withPheno.csv",sep=""))
     data<-read.csv(paste(prefix,setfun2(),"/5_WGCNA_4k_tv/results/ModuleEigengenes_all_samples.csv",sep=""))
     print("data is read")
     # #debug
     # data0<-read.csv("/home/rstudio/output/ModuleEigengenes_all_samples_withPheno.csv")
     # data<-read.csv("/home/rstudio/output/ModuleEigengenes_all_samples.csv")

     data<-data[order(data$X),]
     data.blood<-data[seq(1,nrow(data)-1,2),]
     data.fluid<-data[seq(2,nrow(data),2),]

     mes<-colnames(data.blood)[2:ncol(data.blood)]
     menames<-sort(colnames(data)[2:length(colnames(data))])
     print("menames")
     print(menames)
     intracor<-data.frame()
     for (name in mes){
       res<-cor.test(eval(parse(text=paste("data.blood$",name,sep=""))),eval(parse(text=paste("data.fluid$",name,sep=""))),method="p")
       # plot(eval(parse(text=paste("data.blood$",name,sep=""))),eval(parse(text=paste("data.fluid$",name,sep=""))),col=unlist(strsplit(name,"ME"))[2],main=paste(name,"\nPearsonR=",sprintf("%.3f",res$estimate),"Pval=",sprintf("%.3f",res$p.value)),pch=19,xlab="Blood",ylab="Fluid",xlim=c(-0.25,0.25),ylim=c(-0.25,0.25))
       # abline(lm(eval(parse(text=paste("data.blood$",name,sep="")))~eval(parse(text=paste("data.fluid$",name,sep="")))), col="red")
       # print(name)
       #print("here is a res")
       #print(res)
       intracor<-rbind(intracor,c(unlist(strsplit(name,"ME"))[2],res$estimate,res$p.value))
     }
     colnames(intracor)<-c("Pheno","PearsonR","P-value")
     write.csv(intracor,file.path(tabledir,paste(square,"intracor_ME.csv",sep="_")))
     vals1$square<-square
     vals1$intracor<-intracor
     vals1$menames<-menames
     
     
     })
     
   vals2<-reactiveValues()
   intracorcalc2<-reactive({
     intracorcalc()
     menames<-vals1$menames
     
     #ME intracorrelates
     #print("reading data")
     print("doing intracorcalc2")
     squarelist<-unlist(strsplit(setfun2(),"_"))
     if(length(squarelist)==3){
       square<-squarelist[3]
     }else{square<-paste(squarelist[c(3,4)],collapse="_")}
     print("square")
     print(square)
     #prefix<-paste("/home/rstudio/output/build/",buildfun(),"/",sep="")
     prefix<-"/home/rstudio/output/build/"
     data0<-read.csv(paste(prefix,setfun2(),"/5_WGCNA_4k_tv/results/ModuleEigengenes_all_samples_withPheno.csv",sep=""))
     data<-read.csv(paste(prefix,setfun2(),"/5_WGCNA_4k_tv/results/ModuleEigengenes_all_samples.csv",sep=""))
     print("data is read 2")
     phenovec<-colnames(data0[,which(!colnames(data0)%in%menames&!colnames(data0)=="X"&!colnames(data0)=="Row.names")])
     print(phenovec)
     pdata<-data0[order(data0$X),phenovec]
     pdatablood<-pdata[seq(1,nrow(pdata)-1,2),]
     pdatafluid<-pdata[seq(2,nrow(data),2),]
     
     
     
     intracor2<-data.frame()
     for (name in phenovec){
       print(name)
       res<-"Skip"
       bloodvec<-eval(parse(text=paste("pdatablood$",name,sep="")))
       fluidvec<-eval(parse(text=paste("pdatafluid$",name,sep="")))
       print("identical(bloodvec,fluidvec)")
       print(!identical(bloodvec,fluidvec))
       if(!identical(bloodvec,fluidvec)){
       try(res<-cor.test(bloodvec,fluidvec,method="p"))
       }
       # plot(eval(parse(text=paste("data.blood$",name,sep=""))),eval(parse(text=paste("data.fluid$",name,sep=""))),col=unlist(strsplit(name,"ME"))[2],main=paste(name,"\nPearsonR=",sprintf("%.3f",res$estimate),"Pval=",sprintf("%.3f",res$p.value)),pch=19,xlab="Blood",ylab="Fluid",xlim=c(-0.25,0.25),ylim=c(-0.25,0.25))
       # abline(lm(eval(parse(text=paste("data.blood$",name,sep="")))~eval(parse(text=paste("data.fluid$",name,sep="")))), col="red")
       # print(name)
       print("here is a res2")
       print(res)
       print(class(res))
       if(class(res)!="character"){
       intracor2<-rbind(intracor2,c(name,res$estimate,res$p.value))
       }
     }
     colnames(intracor2)<-c("Pheno","PearsonR","P-value")
     write.csv(intracor2,file.path(tabledir,paste(square,"intracor2_pheno.csv",sep="_")))
     print("intracorcalc2 is done")
     
      vals2$pdatablood<-pdatablood
       vals2$pdatafluid<-pdatafluid
       vals2$intracor2<-intracor2
     })
     
   output$intracor<-renderPlot({
     intracorcalc()
     intracor<-vals1$intracor
     square<-vals1$square
     
     #print("intracor is done")
     colnames(intracor)<-c("name","intracor","pval")
     query<-paste("MATCH (n:wgcna {square:'",square,"',edge:'",edgefun(),"'}) RETURN n.name as name, n.modAUC1 as modAUC1, n.modAUC2 as modAUC2, n.diffME as diffME, n.sigenrich as sigenrich",sep="")
     moduletable<-cypher(graph,query)
     #print("moduletable is done")
     moduletable2<-merge(moduletable,intracor,by.x="name",by.y="name")
     #print("moduletable2 is done")
     write.csv(moduletable2,file.path(tabledir,paste(square,"moduletable2.csv",sep="_")))
     output$dtintracor<-DT::renderDataTable(datatable(moduletable2)%>%formatSignif(names(moduletable2)[2:ncol(moduletable2)],3),server = FALSE,options=list(lengthMenu = c(5, 10, 15,20,50), pageLength = 50),filter="top") 
     
     # ####
     # plot(modAUC1~intracor,data=moduletable2,pch=20,col="grey",main=paste(square,": Compartmentalised vs. representative processes",sep=""),xlab="Representation Index",ylab="Compartmentalisation Index")
     # shadowtext(moduletable2$modAUC1~moduletable2$intracor,labels=moduletable2$name,col=moduletable2$name,font=2,cex=1.0,xpd=TRUE,bg="lightgrey")
     # ####
     
     # abline(h=0.75,col="blue")
     # #abline(v=0.5,col="blue")
     # abline(v=0.0,col="blue")
     print(str(moduletable2))
     moduletable2$intracor<-as.numeric(moduletable2$intracor)
     
     moduletable3<-as.data.frame(lapply(moduletable2[,2:ncol(moduletable2)],as.numeric))
     moduletable3$name<-moduletable2$name
     print(str(moduletable3))
     ####
     #######B
     vp<-ggplot(moduletable3,aes(intracor,modAUC1))+
       geom_point(aes(col=diffME),size=12,alpha=0.4)+
       scale_color_gradient2(low="blue",mid="lightyellow",high="red") +
       geom_text_repel(aes(label=name),data=moduletable3,size=6) +
       #scale_x_discrete() +
       #theme_tufte() +
       theme_light() +
       #theme_minimal() +
       labs(y = "Compartmentalisation index") +
       labs(x = "Representation index") +
       labs(colour = paste("diffME edge:",edgefun())) +
       theme(axis.title.y = element_text(size = rel(2.5), angle = 90, color = "slategrey",margin = margin(t = 0, r = 50, b = 50, l = 50))) +
       theme(axis.title.x = element_text(size = rel(2.5), color = "slategrey",margin = margin(t = 50, r = 50, b = 0, l = 50))) + 
       theme(axis.text = element_text(size = rel(1.5))) +
       theme(aspect.ratio=1) +
       theme(legend.text = element_text(size = rel(1.2))) +
       #theme(legend.title = element_text(size = rel(1.2))) +
       #theme(legend.text=element_text(size=rel(1.4))) +
       #theme(legend.key.size = unit(1,"cm")) +
       #labs(y = phenofun2()) +
       ggtitle(paste(square,": intracor",sep="")) +
       theme(plot.title = element_text(size=rel(2.5),hjust = 0.5,color="maroon"))
     #######
     print(vp)
     ####
     if(input$returnpdf==TRUE){
       pdf("plot.pdf",width=12,height=7,onefile=FALSE)
       plot(modAUC1~intracor,data=moduletable2,pch=20,col="grey",main=paste(square,": Compartmentalised vs. representative processes",sep=""),xlab="Representation Index",ylab="Compartmentalisation Index")
       shadowtext(moduletable2$modAUC1~moduletable2$intracor,labels=moduletable2$name,col=moduletable2$name,font=2,cex=1.5,xpd=TRUE,bg="lightgrey")
       #text(modAUC1~intracor,data=moduletable2,labels=moduletable2$name,col=moduletable2$name,font=2,cex=1.0,xpd=TRUE,bg=col2rgb(moduletable2$name)+5)
       abline(h=0.75,col="blue")
       abline(v=0.5,col="blue")
       dev.off()
     }
     
   })
   
     output$dtintracorpheno<-DT::renderDataTable({
         print("doing the data table for intracor2")
         intracorcalc2()
         intracor2<-vals2$intracor2
         print("the table should be ready")
         print(intracor2)
         
         #datatable(intracor2)%>%formatSignif(names(intracor2),3)},server = FALSE,options=list(lengthMenu = c(5, 10, 15,20,50), pageLength = 50))
         datatable(intracor2)%>%formatSignif(names(intracor2)[2:3],3)},server = FALSE,options=list(lengthMenu = c(5, 10, 15,20,50), pageLength = 50),filter="top")
     
     output$phenointracor<-renderPlot({
       intracorcalc2()
       intracor2<-vals2$intracor2
       pdatablood<-vals2$pdatablood
       
       pdatafluid<-vals2$pdatafluid
       
       #output$dtintracorpheno<-DT::renderDataTable(datatable(intracor2)%>%formatSignif(names(intracor2),3),server = FALSE,options=list(lengthMenu = c(5, 10, 15,20,50), pageLength = 50)) 
       
       #nameX<-phenofun2()
       nameX<-phenofunIC()
       
       bloodvecX<-eval(parse(text=paste("pdatablood$",nameX,sep="")))
       fluidvecX<-eval(parse(text=paste("pdatafluid$",nameX,sep="")))
       cordat<-data.frame(cbind(bloodvecX,fluidvecX))
       colnames(cordat)<-c("blood","fluid")
       phenointracorplot<-ggplot(cordat, aes(x=blood, y=fluid)) + 
         geom_point(shape=18, color="red",size=10) +
         geom_smooth(method=lm) + 
         stat_cor(method = "pearson",size=6,label.y.npc = "top",label.x.npc = "middle") +
         geom_rug() +
         theme_light() +
         theme(axis.title.y = element_text(size = rel(2.5), angle = 90, color = "slategrey",margin = margin(t = 0, r = 50, b = 50, l = 50))) +
         theme(axis.title.x = element_text(size = rel(2.5), color = "slategrey",margin = margin(t = 50, r = 50, b = 0, l = 50))) + 
         theme(axis.text = element_text(size = rel(1.5))) +
         theme(aspect.ratio=1) +
         ggtitle(nameX) +
         theme(plot.title = element_text(size=rel(2.5),hjust = 0.5,color="maroon"))
       phenointracorplot
     },height=900)
     
     
   }#end conditional intracor
   
   #sankey####
   output$sankey<-renderSankeyNetwork({
     squarelist<-unlist(strsplit(setfun(),"_"))
     if(length(squarelist)==3){
       square<-squarelist[3]
     }else{square<-paste(squarelist[c(3,4)],collapse="_")}
     
     #values
     node1<-sankeynode1fun()
     node2<-sankeynode2fun()
     cell1<-sankeycell1fun()
     cell2<-sankeycell2fun()
     
     #query builder
     if(node1=="wgcna"){
       query1=paste("MATCH (n:wgcna {square:'",square,"',edge:'5'})-[r0]-(p:PROBE)-[r1]-(s:SYMBOL)-[r2]-(x:",node2,") ",sep="")
     }else{
       query1=paste("MATCH (n:",node1,")-[r1]-(s:SYMBOL)-[r2]-(x:",node2,") ",sep="")  
       }
     
     # if(node1=="wgcna"){
     #   query2<-" RETURN n.name AS node1, COUNT(r0) AS nprobes, s.name AS SYMBOL, COUNT(r2) AS ngenes, x.name AS node2"
     # }else{
     #   query2<-" RETURN n.name AS node1, s.name AS SYMBOL, COUNT(r2) AS ngenes, x.name AS node2"
     # }
     
     query2<-" RETURN n.name AS node1, COUNT(r2) AS ngenes, x.name AS node2"
     
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
     
     #Scolors <- paste(links2$node1, collapse = '", "')
     #colorJS <- paste('d3.scaleOrdinal(["', Scolors, '"])')
     #colourScale = colorJS, 
     
     sankeyNetwork(Links = links2, Nodes = nodes, Source = 'source',Target = 'target',Value = 'value', NodeID = 'name',fontSize = 12, nodeWidth = 30,LinkGroup=linkgroupfun(), iterations = 128)
     
     # ##Begin inset
     # # #Cell
     # # query<-"MATCH (c:cellEx)-[r1]-(s:SYMBOL)-[r2]-(p:PROBE)-[r3]-(n:wgcna) WHERE p.square='blood.PCF.defTB' AND p.edge=5 RETURN c.name AS cell, COUNT(r1) AS ngenes, s.name AS SYMBOL, COUNT(r2) AS nprobes, p.name AS probename, p.logfc as logfc, n.name as wgcnamod, n.diffME AS diffME"
     # # res<-cypher(graph,query)
     # all<-unique(res[,c("node1","node2")])
     # colnames(all)<-c("source","target")
     # 
     # # maplist<-list()
     # # maplistnames<-c("wgcna","HaemAtlas","Abbas")
     # 
     # ##
     # list2<-all
     # nodes<-unique(c(list2$source,list2$target))
     # nodes<-data.frame(nodes)
     # colnames(nodes)<-"name"
     # indices<-cbind(seq(0,nrow(nodes)-1),nodes)
     # value<-numeric()
     # print("starting iterative query")
     # print(paste("need to do",nrow(all),"queries!"))
     # print(head(list2))
     # pb = txtProgressBar(min = 0, max = nrow(list2), initial = 0,style=3) 
     # for (row in 1:nrow(list2)){
     #   if (node1=="wgcna"){
     #     query<-paste("MATCH (c {name:'",list2[row,1],"',square:'",square,"'})-[r1]-(p:PROBE)-[r2]-(s:SYMBOL)-[r3]-(n {name:'",list2[row,2],"'}) RETURN COUNT(p) as value",sep="")
     #   }else{
     #     query<-paste("MATCH (c {name:'",list2[row,1],"'})-[r1]-(s:SYMBOL)-[r3]-(n {name:'",list2[row,2],"'}) RETURN COUNT(s) as value",sep="")
     #     }
     #   
     #   res<-cypher(graph,query)
     #   value<-c(value,res$value)
     #   setTxtProgressBar(pb,row)
     # }
     # close(pb)
     # print("done with queries")
     # list2.1<-cbind(list2,value)
     # for(index in 1:nrow(indices)){
     #   list2.1[list2.1==indices[index,2]]<-indices[index,1]
     # }
     # list2.1$source<-as.integer(list2.1$source)
     # list2.1$target<-as.integer(list2.1$target)
     # 
     # map<-list("nodes"=nodes,"links"=list2.1)
     # ##
     # # Plot
     # sankeyNetwork(Links = map$links, Nodes = map$nodes, Source = 'source',Target = 'target',Value = 'value', NodeID = 'name',fontSize = 12, nodeWidth = 30)

     # for(bloodset in 1:3){
     #   list1<-all[grep(paste("_",bloodset,sep=""),all$source),]
     #   nodes<-unique(c(list1$source,list1$target))
     #   nodes<-data.frame(nodes)
     #   colnames(nodes)<-"name"
     #   indices<-cbind(seq(0,nrow(nodes)-1),nodes)
     #   value<-numeric()
     #   for (row in 1:nrow(list1)){
     #     query<-paste("MATCH (c:cellEx {name:'",list1[row,1],"'})-[r1]-(s:SYMBOL)-[r2]-(p:PROBE)-[r3]-(n:wgcna {name:'",list1[row,2],"'}) WHERE p.square='blood.PCF.defTB' AND p.edge=5 RETURN COUNT(p) as value",sep="")
     #     res<-cypher(graph,query)
     #     value<-c(value,res$value)
     #     
     #   }
     #   list1.1<-cbind(list1,value)
     #   for(index in 1:nrow(indices)){
     #     list1.1[list1.1==indices[index,2]]<-indices[index,1]
     #   }
     #   list1.1$source<-as.integer(list1.1$source)
     #   list1.1$target<-as.integer(list1.1$target)
     #   
     #   map<-list("nodes"=nodes,"links"=list1.1)
     #   
     #   maplist[[maplistnames[bloodset]]]<-map
     #   
     # }
     # 
     #  # Plot
     # sankeyNetwork(Links = maplist[[1]]$links, Nodes = maplist[[1]]$nodes, Source = 'source',
     #               Target = 'target',Value = 'value', NodeID = 'name',
     #               fontSize = 12, nodeWidth = 30)
     # 
     # sankeyNetwork(Links = maplist[[2]]$links, Nodes = maplist[[2]]$nodes, Source = 'source',
     #               Target = 'target',Value = 'value', NodeID = 'name',
     #               fontSize = 12, nodeWidth = 30)
     # 
     # sankeyNetwork(Links = maplist[[3]]$links, Nodes = maplist[[3]]$nodes, Source = 'source',
     #               Target = 'target',Value = 'value', NodeID = 'name',
     #               fontSize = 12, nodeWidth = 30)
     ##End inset
     
     
     
     #map
     
     #plot
     
     
     
     
   })
   
   #projections####
   output$projections<-renderPlot({
     squarelist<-unlist(strsplit(setfun(),"_"))
     if(length(squarelist)==3){
       square<-squarelist[3]
     }else{square<-paste(squarelist[c(3,4)],collapse="_")}
     if(nodefun()%in%c("reactomePW","PalWangPW","ImmunePW","cellEx")){
       #query.base=paste("MATCH (a:",nodefun(),")-[r]-(b:",nodefun2(),")",sep="")
       query.base=paste("MATCH (a:",nodefun(),")-[r]-(b) WHERE (b:",nodefun2(),")",sep="")
     }else{
        #query.base=paste("MATCH (a:",nodefun()," {square:'",square,"',edge:'",edgefun(),"'})-[r]-(b:",nodefun2(),")",sep="")
        query.base=paste("MATCH (a:",nodefun()," {square:'",square,"',edge:'",edgefun(),"'})-[r]-(b) WHERE (b:",nodefun2(),")",sep="")
        }
     
     print(query.base)
     nodelist<-c("a","b")
     edgetrips<-list(c("a","r","b"))
     rimpar="diffME"
     nodeproj=nodefun()
     edgeproj=edgefun()
     if(nodefun2()!="pheno"&nodefun2()!="cellprop"){
       igraph_plotter(query.base,nodelist,edgetrips,rimpar,csv=F,prefix=cytodir,filename=paste("modnet",square,edgeproj,nodeproj,collapse="_"),vertexsize=vertexsizefun(),legendcex=legendcexfun(),optlabel="\n",optvalue="V(ig)$edge",optchar="V(ig)$square",lay_out=eval(parse(text=layoutfun())),plot_bipartite = TRUE,plotd3=FALSE,plotwhich=plotwhichfun(),vertex.label.cex=vlc(),edgewidth=ewf())
       #from igraphMM: igraph_plotter(query.base,nodelist,edgetrips,rimpar="diffEX",plot=TRUE,csv=TRUE,prefix=cytodir,filename = "igraphModuleMeta",vertexsize=vertexsizefun(),legendcex=legendcexfun(),plotd3=FALSE,vertex.label.cex=vlc(),lay_out=eval(parse(text=layoutfun())))
       
     }else if(nodefun2()=="cellprop"){
       query.base<-paste(query.base," AND toFloat(r.Rsq) > ",phenoWeightLowfun()," AND toFloat(r.Rsq) <= ",phenoWeightHighfun()," ",sep="")
       igraph_plotter(query.base,nodelist,edgetrips,rimpar,csv=F,prefix=cytodir,filename=paste("modnet",square,edgeproj,nodeproj,collapse="_"),vertexsize=vertexsizefun(),legendcex=legendcexfun(),optlabel="\n",optvalue="V(ig)$edge",optchar="V(ig)$square",lay_out=eval(parse(text=layoutfun())),plot_bipartite = TRUE,plotd3=FALSE,plotwhich=plotwhichfun(),vertex.label.cex=vlc(),edgewidth=ewf())
       #from igraphMM: igraph_plotter(query.base,nodelist,edgetrips,rimpar="diffEX",plot=TRUE,csv=TRUE,prefix=cytodir,filename = "igraphModuleMeta",vertexsize=vertexsizefun(),legendcex=legendcexfun(),plotd3=FALSE,vertex.label.cex=vlc(),lay_out=eval(parse(text=layoutfun())))
       
     }else{
       if(bygroupfun2()=="vt"){
         query.base=paste("MATCH (a:",nodefun()," {square:'",square,"',edge:'",edgefun(),"'})-[r]-(b:pheno)-[rv0]-(c:personPheno)-[rv1]-(v1:varsubtype)-[rv2]-(v2:vartype {name:'",phenovtfun2(),"'})",sep="")
         query.base<-paste(query.base," WHERE toFloat(r.Rsq) > ",phenoWeightLowfun()," AND toFloat(r.Rsq) <= ",phenoWeightHighfun()," ",sep="")
         igraph_plotter(query.base,nodelist,edgetrips,rimpar,csv=F,prefix=cytodir,filename=paste("modnet",square,edgeproj,nodeproj,collapse="_"),vertexsize=vertexsizefun(),legendcex=legendcexfun(),optlabel="\n",optvalue="V(ig)$edge",optchar="V(ig)$square",lay_out=eval(parse(text=layoutfun())),plot_bipartite = TRUE,plotd3=FALSE,plotwhich=plotwhichfun(),vertex.label.cex=vlc(),edgewidth=ewf())
         #from igraphMM: igraph_plotter(query.base,nodelist,edgetrips,rimpar="diffEX",plot=TRUE,csv=TRUE,prefix=cytodir,filename = "igraphModuleMeta",vertexsize=vertexsizefun(),legendcex=legendcexfun(),plotd3=FALSE,vertex.label.cex=vlc(),lay_out=eval(parse(text=layoutfun())))
         
       }else if(bygroupfun2()=="vst"){
         query.base=paste("MATCH (a:",nodefun()," {square:'",square,"',edge:'",edgefun(),"'})-[r]-(b:pheno)-[rv0]-(c:personPheno)-[rv1]-(v1:varsubtype {name:'",phenovstfun2(),"'})",sep="")
         query.base<-paste(query.base," WHERE toFloat(r.Rsq) > ",phenoWeightLowfun()," AND toFloat(r.Rsq) <= ",phenoWeightHighfun()," ",sep="")
         igraph_plotter(query.base,nodelist,edgetrips,rimpar,csv=F,prefix=cytodir,filename=paste("modnet",square,edgeproj,nodeproj,collapse="_"),vertexsize=vertexsizefun(),legendcex=legendcexfun(),optlabel="\n",optvalue="V(ig)$edge",optchar="V(ig)$square",lay_out=eval(parse(text=layoutfun())),plot_bipartite = TRUE,plotd3=FALSE,plotwhich=plotwhichfun(),vertex.label.cex=vlc(),edgewidth=ewf())
         #from igraphMM: igraph_plotter(query.base,nodelist,edgetrips,rimpar="diffEX",plot=TRUE,csv=TRUE,prefix=cytodir,filename = "igraphModuleMeta",vertexsize=vertexsizefun(),legendcex=legendcexfun(),plotd3=FALSE,vertex.label.cex=vlc(),lay_out=eval(parse(text=layoutfun())))
         
       }else{
         query.base=paste("MATCH (a:",nodefun()," {square:'",square,"',edge:'",edgefun(),"'})-[r]-(b:pheno)",sep="")
         query.base<-paste(query.base," WHERE toFloat(r.Rsq) > ",phenoWeightLowfun()," AND toFloat(r.Rsq) <= ",phenoWeightHighfun()," ",sep="")
         igraph_plotter(query.base,nodelist,edgetrips,rimpar,csv=F,prefix=cytodir,filename=paste("modnet",square,edgeproj,nodeproj,collapse="_"),vertexsize=vertexsizefun(),legendcex=legendcexfun(),optlabel="\n",optvalue="V(ig)$edge",optchar="V(ig)$square",lay_out=eval(parse(text=layoutfun())),plot_bipartite = TRUE,plotd3=FALSE,plotwhich=plotwhichfun(),vertex.label.cex=vlc(),edgewidth=ewf())
         #from igraphMM: igraph_plotter(query.base,nodelist,edgetrips,rimpar="diffEX",plot=TRUE,csv=TRUE,prefix=cytodir,filename = "igraphModuleMeta",vertexsize=vertexsizefun(),legendcex=legendcexfun(),plotd3=FALSE,vertex.label.cex=vlc(),lay_out=eval(parse(text=layoutfun())))
         
       }
       
       
     }
     
     
   })
   
  
  #Virtual Cells#### 
   
   #cellcor####
   output$cellcor<-renderPlot({
     print("reacstudyCC")
     reacstudyCC<-as.numeric(studyCellcorfun())
     print(reacstudyCC)
     
     print("reacsetCC")
     reacsetCC<-as.numeric(squareCellcorfun())#[reacstudyCC])
     print(reacsetCC)
     
     print("reacedgeCC")
     reacedgeCC<-as.numeric(edgeCellcorfun())
     print(reacedgeCC)
     
     studyCC<-names(studies)[reacstudyCC]
     setlistCC<-setlistnameslist[[reacstudyCC]]
     squareCC<-setlistCC[[reacsetCC]]
     
     print(studyCC)
     print(squareCC)
     print(edgeCellcorfun())
     
     cellCore4(study=studyCC,squareC=squareCC,edgeC=as.character(reacedgeCC),cellC=cellCfun(),pwm="other",PalWang=FALSE)
     #cellCore2(study=studyname,squareC=set[setcount],edgeC=as.character(edge),cellGroup=cell,pwm="other",PalWang=FALSE,plotCell=FALSE)
     # print(dev.cur())
     # cur_dev <- dev.cur()
     # #new
     # ccc<-cellCore3(study=studyCC,squareC=squareCC,edgeC=as.character(reacedgeCC),cellC=cellCfun(),pwm="other",PalWang=FALSE)
     # dev.set(cur_dev)
     # print(ccc)
     
   })
   
   #cellmatrix####
   plotInputGigamat1 <- reactive({
     #cellcordir<-file.path("/home/rstudio/output/build",buildfun(),"12_CELLCOR/figures")
     cellcordir<-file.path("/home/rstudio/output/build/CELLCOR/figures")
     print(cellcordir)
     load(file.path(cellcordir,"csl.RData"))
     #"/home/rstudio/output/build//version_2017-03-25_20_17_30//12_CELLCOR/figures/csl.RData"
     ##CODE
     
     reacstudy<-as.numeric(studyCellcorfun())
     reacset<-as.numeric(squareCellcorfun())#[reacstudy])
     reacedge<-as.numeric(edgeCellcorfun())
     print("print(c(reacstudy,reacset,reacedge))")
     print(c(reacstudy,reacset,reacedge))
     xc<-studycellscorelist[[reacstudy]][[reacset]][[reacedge]]
     
     #xc<-edgeset
     xcvec<-vector()
     for (xcm in xc){
       xcvec<-c(xcvec,names(xcm))
     }
     print("debug1")
     xcvec<-unique(xcvec)
     datamatrix<-matrix(nrow=length(xcvec),ncol=length(names(xc)),data=0)
     dimnames(datamatrix)<-list(xcvec,names(xc))
     for (xcm in 1:length(xc)){
       datamatrix[names(xc[xcm][[1]]),names(xc)[xcm]]<-as.numeric(unlist(xc[xcm]))
     }
     
     #more subsetting based on study
     print("input$cellStudy")
     print(input$cellStudy)
     print(str(input$cellStudy))
     print("grep(paste(input$cellStudy,collapse='|'),colnames(datamatrix))")
     print(grep(paste(input$cellStudy,collapse="|"),colnames(datamatrix)))
     datamatrix<-datamatrix[,grep(paste(input$cellStudy,collapse="|"),colnames(datamatrix))]
     
     print("debug2")
     
     datamatrix[is.na(datamatrix)]<-0
     breaksvals<-max(abs(datamatrix))
     tv2<-apply(datamatrix,1,function(ff){max(abs(ff),na.rm=TRUE)})>slider1()
     tv2[is.na(tv2)]<-0
     tv3<-apply(datamatrix,2,function(ff){max(abs(ff),na.rm=TRUE)})>slider2()
     tv3[is.na(tv3)]<-0
     plotMatrix<-datamatrix[as.logical(tv2),as.logical(tv3)]
     print("debug3")
     
     
     
     #cex calculation
     cexval1<-3*sqrt(28/ncol(plotMatrix))
     cexval2<-2*sqrt(120/nrow(plotMatrix))
     
    
     print("debug4")
     print(plotMatrix[1:5,1:5])
     
     print("debug5")
     print(dev.cur())
     cur_dev <- dev.cur()
     #new
     tp<-pheatmap(plotMatrix,col=blueWhiteRed(60),breaks=seq(floor(-breaksvals),ceiling(breaksvals),length.out=61),fontsize_col=9+cexval1,fontsize_row=9+cexval2,main=paste("Study:",names(studies)[reacstudy],"Set:",setlistnameslist[[reacstudy]][reacset],"Edge:",reacedge,"PW filter:",sprintf("%.3f",slider1()),"Cell filter:",sprintf("%.3f",slider2())))
     dev.set(cur_dev)
     print(tp)
     
     #pheatmap(datamatrix,col=blueWhiteRed(60),breaks=seq(floor(-breaksvals),ceiling(breaksvals),length.out=61))
     if(input$returnpdf==TRUE){
       pdf("plot.pdf",width=16,height=as.numeric(input$plotheightG1)/100,onefile=FALSE)
       pheatmap(plotMatrix,col=blueWhiteRed(60),breaks=seq(floor(-breaksvals),ceiling(breaksvals),length.out=61),fontsize_col=9+cexval1,fontsize_row=9+cexval2,main=paste("Study:",names(studies)[reacstudy],"Set:",setlistnameslist[[reacstudy]][reacset],"Edge:",reacedge,"PW filter:",sprintf("%.3f",slider1()),"Cell filter:",sprintf("%.3f",slider2())))
       #pheatmap(plotMatrix,col=blueWhiteRed(60),breaks=seq(floor(-breaksvals),ceiling(breaksvals),length.out=61),fontsize_col=9+cexval1,fontsize_row=9+cexval2,main=paste("study"))
       dev.off()
     }
     #end new
     
     

     
     
   })
   
   plotInputGigbar <- reactive({
     ##CODE
     #cellcordir<-file.path("/home/rstudio/output/build/",buildfun(),"/12_CELLCOR/figures")
     cellcordir<-file.path("/home/rstudio/output/build/CELLCOR/figures")
     load(file.path(cellcordir,"csl.RData"))
     reacstudy<-as.numeric(usestudy())
     reacset<-as.numeric(useset())#[reacstudy]
     reacedge<-as.numeric(useedge())
     xc<-studycellscorelist[[reacstudy]][[reacset]][[reacedge]]
     #xc<-edgeset
     xcvec<-vector()
     for (xcm in xc){
       xcvec<-c(xcvec,names(xcm))
     }
     xcvec<-unique(xcvec)
     datamatrix<-matrix(nrow=length(xcvec),ncol=length(names(xc)),data=0)
     dimnames(datamatrix)<-list(xcvec,names(xc))
     for (xcm in 1:length(xc)){
       datamatrix[names(xc[xcm][[1]]),names(xc)[xcm]]<-as.numeric(unlist(xc[xcm]))
     }
     datamatrix[is.na(datamatrix)]<-0
     rs<-apply(datamatrix,1,function(row){sum(abs(row))})
     cs<-apply(datamatrix,2,function(row){sum(abs(row))})
     
     rs.up<-apply(datamatrix,1,function(row){sum(row[which(row>0)])})
     cs.up<-apply(datamatrix,2,function(row){sum(row[which(row>0)])})
     
     rs.down<-apply(datamatrix,1,function(row){abs(sum(row[which(row<0)]))})
     cs.down<-apply(datamatrix,2,function(row){abs(sum(row[which(row<0)]))})
     
     allrows<-as.matrix(cbind(rs.up,rs.down))
     sumvec<-apply(allrows,1,sum)
     
     allcols<-as.matrix(cbind(cs.up,cs.down))
     sumvec2<-apply(allcols,1,sum)
     barfilter1<-barfilterfun1()
     barfilter2<-barfilterfun2()
     par(mar=c(5,30,2,2),xpd=T,mfrow=c(1,2),oma = c(0, 0, 2, 0))
     barplot(t(allrows[order(sumvec,decreasing=F),][(nrow(allrows)-20):nrow(allrows),][-1:-barfilter1,]),beside=F,horiz=T,las=2,col=c("red","blue"),xlab="Cumulative pathway score")
     barplot(t(allcols[order(sumvec2,decreasing=F),][-1:-barfilter2,]),beside=F,horiz=T,las=2,col=c("red","blue"),xlab="Cumulative pathway score")
     mtext(paste("Study:",names(studies)[reacstudy],"Set:",setlistnameslist[[reacstudy]][reacset],"Edge:",reacedge),outer=TRUE,cex=2)
     if(input$returnpdf==TRUE){
       pdf("plot.pdf",width=18,height=7,onefile=FALSE)
       par(mar=c(5,30,2,2),xpd=T,mfrow=c(1,2),oma = c(0, 0, 2, 0))
       barplot(t(allrows[order(sumvec,decreasing=F),][(nrow(allrows)-20):nrow(allrows),][-1:-barfilter1,]),beside=F,horiz=T,las=2,col=c("red","blue"),xlab="Cumulative pathway score")
       barplot(t(allcols[order(sumvec2,decreasing=F),][-1:-barfilter2,]),beside=F,horiz=T,las=2,col=c("red","blue"),xlab="Cumulative pathway score")
       mtext(paste("Study:",names(studies)[reacstudy],"Set:",setlistnameslist[[reacstudy]][reacset],"Edge:",reacedge),outer=TRUE,cex=2)
       dev.off()
     }
     
   })
   
   output$cellmatrix <- renderUI({
     plotOutput("gigamat", height = as.numeric(input$plotheightG1))
     #plotOutput("gigamat", height = 7000)
   })
   
   output$gigamat <- renderPlot({
     plotInputGigamat1()
   })

   output$gigabar<-renderPlot({
     plotInputGigbar()
     
   })
   
   #Meta-Analysis####
   
   #ModuleMeta####
   
   plotInput <- reactive({
     usematrix<-usematrixfun()
     usematrix<-usematrix[order(usematrix[,2]),]
     usematrix<-usematrix[order(usematrix[,1]),]
     print("begin loop")
     print(usematrix)
     
     
     #data$sets<-sort(data$sets)
     data$sets<-sort(unlist2(setlistnameslist))
     
     setlist=data$sets[unique(usematrix[,1])]
     edgelist=as.character(unique(usematrix[,2]))
     print(setlist)
     #print(usesets())
     print(edgelist)
     #print(useedges())
     #lst.all<-moduleMeta(type="sets",setlist=usesets(),edgelist=useedges(),extent=111)
     lst.all<-moduleMeta(type="sets",setlist=setlist,edgelist=edgelist,extent=111,plotsimple=FALSE)
     
     lst<-lst.all[[1]]
     
     print("lst dims")
     print(str(lst))
     
     allmat2<-c()
     allmatnames<-c()
     allmatnames_label<-c()
     #setlist<-usesets()
     #iteredges<-c(1,2,5)
     chosenSubset=subset.list[[input$subset]]
     for (iter in 1:length(lst)){
       print("iter")
       usematrixiter<-NA
       print("unique")
       print(unique(usematrix[,1][iter]))
       usematrixiter<-usematrix[which(usematrix[,1]==unique(usematrix[,1])[iter]),]
       print("usematrixiter")
       print(class(usematrixiter))
       print(usematrixiter)
       setstoUse<-NA
       print("new debug")
       print("class usematrixiter")
       print(class(usematrixiter))
       #ew code
       #if(class(usematrixiter)=="matrix"){
       if(class(usematrixiter)[1]=="matrix"){
         setstoUse<-as.character(usematrixiter[,2])
       }else if (class(usematrixiter)=="integer"){setstoUse<-as.character(usematrixiter[2])}
       print("setstoUse")
       print(class(setstoUse))
       print(setstoUse)
       smallmat=lst[[iter]]
       #rownames(smallmat)<-paste(setlist[iter],"_edge_",edgelist,sep="")
       smallmat=smallmat[which(rownames(smallmat)%in%setstoUse),]
       print(class(smallmat))
       #row.names(smallmat)<-paste(setlist[iter],"_edge_",setstoUse,sep="")
       #print(smallmat)
       allmat2<-rbind(allmat2,smallmat)
       allmatnames_label<-c(allmatnames_label,paste(setlist[iter],"_edge_",setstoUse,sep=""))
       datasetrow<-datasetmatrix[setlist[iter],]
       allmatnames<-c(allmatnames,paste(setlist[iter],"_",datasetrow[as.numeric(setstoUse)],sep=""))
       
       #print(allmat2[,1:4])
     }
     rownames(allmat2)<-allmatnames
     #optional allsame code
     print("debug A")
     print(str(allmat2))
     #print(max(abs(allmat2)))
     
     allmat3<-allmat2[,chosenSubset]
     print(str(allmat3))
     
     
     print("which(apply(allmat3,2,function(x){min(abs(x))>threshfun()}))")
     print(which(apply(allmat3,2,function(x){min(abs(x))>threshfun()})))
     print("which(apply(allmat3,2,function(x){min(abs(x))<threshfun()}))")
     print(which(apply(allmat3,2,function(x){min(abs(x))<threshfun()})))
     print("which(apply(allmat3,2,function(x){min(abs(x))==threshfun()}))")
     print(which(apply(allmat3,2,function(x){min(abs(x))==threshfun()})))
     print("which(apply(allmat3,2,function(x){min(abs(x))==threshfun()}))")
     print("which(apply(allmat3,2,function(x){is.na(x)}))")
     print(which(apply(allmat3,2,function(x){is.na(x)})))
      
     
     allmat3[!is.finite(allmat3)] <- 0 # the matrix contains non-finite elements like NA, Inf, NaN. Not sure why. replace these with zero. The threshhold function would otherwise remove them anyway.
     
     allmat3<-allmat3[,which(apply(allmat3,2,function(x){min(abs(x))>threshfun()}))]
     
     
     
     print(str(allmat3))
     print("debug B")
     allmat3.allup<-apply(allmat3,2,function(x){sum(x>0)==length(x)})
     
   
     allmat3.down<-apply(allmat3,2,function(x){sum(x<0)==length(x)})
     
 
     allmat3.allsame<-allmat3.allup|allmat3.down
    
    
     if(useallsame()==TRUE){allmat3<-allmat3[,which(allmat3.allsame)]}else if(usediff()==TRUE){allmat3<-allmat3[,which(!allmat3.allsame)]}
     
     print("debugB1.1")
     print(class(allmat3))
     print(max(abs(allmat3)))
           
           
     print(str(allmat3))
     #breaks<-seq(-max(abs(allmat2)),max(abs(allmat2)),length.out=51)
     breaks<-seq(-max(abs(allmat3)),max(abs(allmat3)),length.out=51)
     modulelist<-colnames(allmat3)
     #cex calculation
     cexval<-2*sqrt(258/ncol(allmat3))
     cexval<-cexval*as.numeric(vlc1())
     print("debug B2")
     
     print("chosen subset")
     print(chosenSubset)
     
     finalmods<-colnames(allmat3)
     print("finalmods")
     print(finalmods)
     
     print(head(sdf))
     print(head(chosenSubset))
     anncol<-sdf[chosenSubset,]
     print("anncol before")
     print(anncol)
     anncol<-anncol[which(rownames(anncol)%in%finalmods),]
     
     
     print("debug C")
     print(anncol)
     print(dim(anncol))
     print(str(anncol))
     print(rownames(anncol))
     print(colnames(anncol))
     #anncol<-sdf
     print(head(allmat3))
     print(dim(allmat3))
     print(str(allmat3))
     
     factor_colors = hsv((seq(0, 1, length.out = length(colnames(anncol)) + 1)[-1] + 0.2)%%1, 0.7, 0.95)
     
     ann_col = list()
     for (bb in 1:length(colnames(anncol))){
       #ann_col[[bb]]=c("FALSE"="white","TRUE"=colors()[sample(2:length(colors()),1)])
       ann_col[[colnames(anncol)[bb]]]=c("FALSE"="white","TRUE"=factor_colors[bb])
     }
     print("debug D")
     print(ann_col)
     
     
     cur_dev <- dev.cur()
     pheatout <- pheatmap(allmat3,col=blueWhiteRed(50),scale="none",breaks=breaks,annotation_col = anncol,fontsize_col=6+cexval,fontsize_row=9+cexval,main=paste(input$subset,c("all_same","alldiff")[c(useallsame(),usediff())],"n:",ncol(allmat3)),annotation_colors=ann_col,annotation_legend = F)
     pheatres <- allmat3[c(pheatout$tree_row[["order"]]),pheatout$tree_col[["order"]]]
     dev.set(cur_dev)
     print("printing heatmap")
     print(pheatout)
     print("have finished printing heatmap")
     print(pheatres)
     
     #the plotting
     #
     #pheatmap(allmat3,col=blueWhiteRed(50),scale="none",breaks=breaks,fontsize_col=6+cexval,fontsize_row=9+cexval,main=paste(input$subset,c("all_same","alldiff")[c(useallsame(),usediff())],"n:",ncol(allmat3)),filename = NA)
     #new
     #
     #print(tp)
     
     
     if(input$returnpdf==TRUE){
       pdf("plot.pdf",width=16,height=7,onefile=FALSE)
       pheatmap(allmat3,col=blueWhiteRed(50),scale="none",breaks=breaks,fontsize_col=6+cexval,fontsize_row=9+cexval,main=paste(input$subset,c("all_same","alldiff")[c(useallsame(),usediff())],"n:",ncol(allmat3)),width=16,height=7)
       dev.off()
     }
     if(returncsvfun()){
     mlv<-paste(modulelist,collapse="|")
     tablequery<-paste("MATCH (n:baylor)-[r]-(x) WHERE (x:cellEx OR x:reactomePW OR x:ImmunePW OR x:PalWangPW) AND n.name =~ '",mlv,"' RETURN n.name AS modname, r.qvalue as qval, labels(x) AS kind, x.name as name",sep="")
     tableres<-unique(cypher(graph,tablequery))
     tableres<-tableres[order(tableres$modname),]
     print(returncsvfun())
     
       print("writing table")
       write.csv(tableres,file.path(tabledir,paste("CMA",paste(allmatnames_label,collapse="_"),".csv",sep="")))}#Chaussabel Module Annotation
    
     
     
     return(modulelist) 
   })
  
   
  output$plot.pheat <- renderUI({
      #withSpinner(plotOutput("pheat", height = input$plotheightMeta))
      plotOutput("pheat", height = input$plotheightMeta)
    })
    
  output$pheat <- renderPlot({
      
      plotInput()
      
    })
 
  #igraph####
  output$igraph<-renderPlot({
      data$sets<-sort(data$sets)
      usematrix<-input$x1_cells_selected
      print(usematrix)
      setlist=data$sets[unique(usematrix[,1])]
      edgelist=as.character(unique(usematrix[,2]))
      print(setlist)
      #print(usesets())
      print(edgelist)
      #print(useedges())
      #lst.all<-moduleMeta(type="sets",setlist=usesets(),edgelist=useedges(),extent=111)
      lst.all<-moduleMeta(type="sets",setlist=setlist,edgelist=edgelist,extent=111,plotsimple=FALSE)
      
      lst<-lst.all[[1]]
      print("lst dims")
      print(str(lst))
      
      allmat2<-c()
      allmatnames<-c()
      #setlist<-usesets()
      #iteredges<-c(1,2,5)
      chosenSubset=subset.list[[input$subset]]
      for (iter in 1:length(lst)){
        print("iter")
        usematrixiter<-NA
        print("unique")
        print(unique(usematrix[,1][iter]))
        usematrixiter<-usematrix[which(usematrix[,1]==unique(usematrix[,1])[iter]),]
        print("usematrixiter")
        print(class(usematrixiter))
        print(usematrixiter)
        setstoUse<-NA
        if(class(usematrixiter)=="matrix"){
          setstoUse<-as.character(usematrixiter[,2])
        }else if (class(usematrixiter)=="integer"){setstoUse<-as.character(usematrixiter[2])}
        print("setstoUse")
        print(class(setstoUse))
        print(setstoUse)
        smallmat=lst[[iter]]
        #rownames(smallmat)<-paste(setlist[iter],"_edge_",edgelist,sep="")
        smallmat=smallmat[which(rownames(smallmat)%in%setstoUse),]
        print(class(smallmat))
        #row.names(smallmat)<-paste(setlist[iter],"_edge_",setstoUse,sep="")
        #print(smallmat)
        allmat2<-rbind(allmat2,smallmat)
        allmatnames<-c(allmatnames,paste(setlist[iter],"_edge_",setstoUse,sep=""))
        #print(allmat2[,1:4])
      }
      rownames(allmat2)<-allmatnames  
      #optional allsame code
      allmat3<-allmat2[,chosenSubset]
      
      allmat3<-allmat3[,which(apply(allmat3,2,function(x){min(abs(x))>threshfun()}))]
      
      allmat3.allup<-apply(allmat3,2,function(x){sum(x>0)==length(x)})
      #print("allmat3.allup")
      #print(allmat3.allup)
      allmat3.down<-apply(allmat3,2,function(x){sum(x<0)==length(x)})
      #print("allmat3.down")
      #print(allmat3.down)
      allmat3.allsame<-allmat3.allup|allmat3.down
      #print("allmat3.allsame")
      #print(allmat3.allsame)
      if(useallsame()==TRUE){allmat3<-allmat3[,which(allmat3.allsame)]}else if(usediff()==TRUE){allmat3<-allmat3[,which(!allmat3.allsame)]}
      breaks<-seq(-abs(max(allmat3)),abs(max(allmat3)),length.out=51)
      modulelist<-colnames(allmat3)
      moduleinput<-as.character(paste(modulelist,collapse="|"))
      print(moduleinput)
      query.base<-paste("MATCH (b:baylor {square:'",setlist[1],"' ,edge:'5'})-[r]-(p) WHERE (p:reactomePW OR p:ImmunePW OR p:PalWangPW OR p:cellEx) AND b.name =~ '(?i)",moduleinput,"'",sep="")
      print(query.base)
      nodelist<-c("b","p")
      edgetrips<-list(c("b","r","p"))
      igr<-igraph_plotter(query.base,nodelist,edgetrips,rimpar="diffEX",plot=TRUE,csv=TRUE,prefix=cytodir,filename = "igraphModuleMeta",vertexsize=vertexsizefun1(),legendcex=legendcexfun1(),plotd3=FALSE,vertex.label.cex=vlc1(),lay_out=eval(parse(text=layoutfun1())))
    })
    
  #d3graphMM####
  output$d3graphMM<-renderForceNetwork({
      data$sets<-sort(data$sets)
      #get igraph
      usematrix<-input$x1_cells_selected
      print(usematrix)
      setlist=data$sets[unique(usematrix[,1])]
      edgelist=as.character(unique(usematrix[,2]))
      print(setlist)
      #print(usesets())
      print(edgelist)
      #print(useedges())
      #lst.all<-moduleMeta(type="sets",setlist=usesets(),edgelist=useedges(),extent=111)
      lst.all<-moduleMeta(type="sets",setlist=setlist,edgelist=edgelist,extent=111,plotsimple=FALSE)
      
      lst<-lst.all[[1]]
      print("lst dims")
      print(str(lst))
      
      allmat2<-c()
      allmatnames<-c()
      #setlist<-usesets()
      #iteredges<-c(1,2,5)
      chosenSubset=subset.list[[input$subset]]
      for (iter in 1:length(lst)){
        print("iter")
        usematrixiter<-NA
        print("unique")
        print(unique(usematrix[,1][iter]))
        usematrixiter<-usematrix[which(usematrix[,1]==unique(usematrix[,1])[iter]),]
        print("usematrixiter")
        print(class(usematrixiter))
        print(usematrixiter)
        setstoUse<-NA
        if(class(usematrixiter)=="matrix"){
          setstoUse<-as.character(usematrixiter[,2])
        }else if (class(usematrixiter)=="integer"){setstoUse<-as.character(usematrixiter[2])}
        print("setstoUse")
        print(class(setstoUse))
        print(setstoUse)
        smallmat=lst[[iter]]
        #rownames(smallmat)<-paste(setlist[iter],"_edge_",edgelist,sep="")
        smallmat=smallmat[which(rownames(smallmat)%in%setstoUse),]
        print(class(smallmat))
        #row.names(smallmat)<-paste(setlist[iter],"_edge_",setstoUse,sep="")
        #print(smallmat)
        allmat2<-rbind(allmat2,smallmat)
        allmatnames<-c(allmatnames,paste(setlist[iter],"_edge_",setstoUse,sep=""))
        #print(allmat2[,1:4])
      }
      rownames(allmat2)<-allmatnames  
      #optional allsame code
      allmat3<-allmat2[,chosenSubset]
      
      allmat3<-allmat3[,which(apply(allmat3,2,function(x){min(abs(x))>threshfun()}))]
      
      allmat3.allup<-apply(allmat3,2,function(x){sum(x>0)==length(x)})
      #print("allmat3.allup")
      #print(allmat3.allup)
      allmat3.down<-apply(allmat3,2,function(x){sum(x<0)==length(x)})
      #print("allmat3.down")
      #print(allmat3.down)
      allmat3.allsame<-allmat3.allup|allmat3.down
      #print("allmat3.allsame")
      #print(allmat3.allsame)
      if(useallsame()==TRUE){allmat3<-allmat3[,which(allmat3.allsame)]}else if(usediff()==TRUE){allmat3<-allmat3[,which(!allmat3.allsame)]}
      breaks<-seq(-abs(max(allmat3)),abs(max(allmat3)),length.out=51)
      modulelist<-colnames(allmat3)
      moduleinput<-as.character(paste(modulelist,collapse="|"))
      print(moduleinput)
      query.base<-paste("MATCH (b:baylor {square:'",setlist[1],"' ,edge:'5'})-[r]-(p) WHERE (p:reactomePW OR p:ImmunePW OR p:PalWangPW OR p:cellEx) AND b.name =~ '(?i)",moduleinput,"'",sep="")
      print(query.base)
      nodelist<-c("b","p")
      edgetrips<-list(c("b","r","p"))
      #print("debug before igraph")
      ##NEW PASTE
      #igraphplotter with igr output
      
      igr<-igraph_plotter(query.base,nodelist,edgetrips,rimpar="diffEX",plot=FALSE,csv=TRUE,prefix=cytodir,filename = "igraphModuleMeta",plotd3=FALSE,return_graph=TRUE)
      
      V(igr)$name=V(igr)$nodename
      
      #new
      #cex calculation
      
      pheatout <- pheatmap(allmat3)
      pheatres <- allmat3[c(pheatout$tree_row[["order"]]),pheatout$tree_col[["order"]]]
      ordervec<-order(colnames(pheatres))
      print("ordervec")
      print(ordervec)
      print(pheatres)
      ordermat<-data.frame(cbind(colnames(pheatres),1:ncol(pheatres)))
      colnames(ordermat)<-c("module","order")
      
      #end new
      
      #get nodes and edges for DT
      currentnodes<-read.csv(file.path(cytodir,"igraphModuleMeta_nodes.csv"))
      #print(currentnodes)
      #currentnodes<-currentnodes[ordervec]
      output$moduleNodesMM<-DT::renderDataTable(currentnodes,server = FALSE,options=list(lengthMenu = c(5, 10, 15,20,50), pageLength = 50),filter="top")
      currentedges<-read.csv(file.path(cytodir,"igraphModuleMeta_edges.csv"))
      dict<-currentnodes[,c(2,3)]
      for(row in 1:nrow(dict)){
        #print(row)
        currentedges[currentedges==as.character(dict[row,"node_id"])]<-as.character(dict[row,"nodename"])
      }
      
      for(row in 1:nrow(ordermat)){
        #print(row)
        currentedges$ordering[currentedges$source==as.character(ordermat[row,"module"])]<-as.numeric(ordermat[row,"order"])
      }
      
      
      output$moduleEdgesMM<-DT::renderDataTable(currentedges,server = FALSE,options=list(lengthMenu = c(5, 10, 15,20,50), pageLength = 50),filter="top")
      
      #convert
      # Convert to object suitable for networkD3
      ig_d3 <- igraph_to_networkD3(igr, group = V(igr)$kind)
      
      # Create force directed network plot
      # forceNetwork(Links = ig_d3$links, Nodes = ig_d3$nodes, 
      #              Source = 'source', Target = 'target', 
      #              NodeID = 'name', Group = 'group',legend=TRUE,zoom = TRUE,bounded=TRUE)
      ##END NEW PASTE
      # #R version 3.3.0
      # forceNetwork(Links = ig_d3$links, Nodes = ig_d3$nodes,height=1200,width=1200,charge= -300,colourScale = JS("d3.scale.category20()"),
      #              linkDistance = 65,fontSize=14,opacity=.9,opacityNoHover = .7,radiusCalculation = JS(" Math.sqrt(300)+6"),
      #              Source = 'source', Target = 'target',
      #              NodeID = 'name', Group = 'group',legend=TRUE,zoom = TRUE,bounded=FALSE)
      
      #R version 3.3.3
      forceNetwork(Links = ig_d3$links, Nodes = ig_d3$nodes,height=1200,width=1200,charge= -300,colourScale = JS("d3.scaleOrdinal(d3.schemeCategory10)"),
                   linkDistance = 65,fontSize=14,opacity=.9,opacityNoHover = .7,radiusCalculation = JS(" Math.sqrt(300)+6"),
                   Source = 'source', Target = 'target',
                   NodeID = 'name', Group = 'group',legend=TRUE,zoom = TRUE,bounded=FALSE)
    })   
   
  #wgcnacol####

   output$wgcnacol<-renderPlot({
     usematrix<-input$x1_cells_selected
     setlist=data$sets[unique(usematrix[,1])]
     edgelist=as.character(unique(usematrix[,2]))
     lst.all<-moduleMeta(type="sets",setlist=setlist,edgelist=edgelist,extent=111)
     lst<-lst.all[[1]]
     lst.col<-lst.all[[2]]
     allmat2<-c()
     allmat2.col<-c()
     iteredges<-edgelist
     chosenSubset=subset.list[[input$subset]]
     for (iter in 1:length(lst)){
       smallmat=lst[[iter]]
       rownames(smallmat)<-paste(setlist[iter],"_edge_",iteredges,sep="")
       #rownames(smallmat)<-paste(setlist[iter],"_edge_",c(1,2,5),sep="")
       allmat2<-rbind(allmat2,smallmat)
       smallmat.col<-lst.col[[iter]]
       allmat2.col<-rbind(allmat2.col,smallmat.col)
     }
     print("wgcna debug")
     print(class(allmat2))
     print(str(allmat2))
     print(head(allmat2[,chosenSubset]))
     nmods<-length(allmat2[1,chosenSubset])
     print("wgcna debug hclust 1")
     #new
     allmat2[is.na(allmat2)]<-0
     
     den<-as.dendrogram(hclust(dist(t(allmat2[,chosenSubset]))))
     den.ord<-order.dendrogram(den)
     print("wgcna debug hclust 2")
     den2<-as.dendrogram(hclust(dist(allmat2[,chosenSubset])))
     den.ord2<-order.dendrogram(den2)
     #print(den.ord)
     allmat2.col<-allmat2.col[den.ord2,den.ord]
     print(dim(allmat2.col))
     
     
     yxt<-nrow(allmat2.col)
     xxt<-ncol(allmat2.col)
     #print(allmat2.col)
     # plot(c(0, 1.5*xxt), c(0, 1.5*yxt), type = "n", xlab = "", ylab = "",frame=FALSE,axes=FALSE)
     # axis(2, at=1:xxt-0.5, labels=colnames(allmat2.col),las=2,cex.axis=.6)
     # axis(1, at=1:yxt, labels=rownames(allmat2.col),cex.axis=.6,las=2)
     # i <-1:xxt
     # for(j in 1:yxt){
     #   rect(i, j, i+1, j+1, col=allmat2.col[j,],border="white",lwd=.25 )
     # }
     print(allmat2.col)
     image(t(allmat2[den.ord2,den.ord]),col=t(allmat2.col),axes=FALSE,xlab="",ylab="")
     axis(1, at = 1:ncol(allmat2.col), labels=colnames(allmat2.col),tick=FALSE)
     axis(4, at = 1:nrow(allmat2.col), labels=rownames(allmat2.col),tick=FALSE)
     
   })
   
  #Virtual Cells####
   plotInputGigamat2 <- reactive({
     print("debugG2_1")
     #cellcordir<-file.path("/home/rstudio/output/build/12_CELLCOR/figures")
     cellcordir<-file.path("/home/rstudio/output/build/CELLCOR/figures")
     print(cellcordir)
     load(file.path(cellcordir,"csl.RData"))
     ##CODE
     print("begin input gigamat")
     # slider1<-reactive({input$sliderpw})
     # slider2<-reactive({input$slidercell})
     # #new
     # usematrixfun<-reactive({input$x1_cells_selected})
     # cellfun<-reactive({input$cells})
     usematrix<-usematrixfun()
     #usematrix<-usematrix[order(usematrix[,2]),]
    #usematrix<-usematrix[order(usematrix[,1]),]
     print("begin loop")
     print(usematrix)
     
     #usematrix<-usematrix[order(usematrix[,2]),]
     #usematrix<-usematrix[order(usematrix[,1]),]
     print("begin loop")
     print("this is usematrix:")
     print(usematrix)
     studyvec<-numeric()
     for (sets in 1:length(setlist)){
       currentvec<-rep(sets,length(setlist[[sets]]))
       studyvec<-c(studyvec,currentvec)
     }
     studyvec<-studyvec[usematrix[,1]]
     
     #this needs work
     # if(project=="IMPI"){
     #   studyvec[which(usematrix[,1]<4)]<-1
     #   studyvec[which(usematrix[,1]>3)]<-2 
     #   
     # }else if (project=="SLE"){
     #   studyvec[which(usematrix[,1]==1)]<-1
     #   studyvec[which(usematrix[,1]==4)]<-2
     #   studyvec[which(usematrix[,1]%in%c(2,3,5,6))]<-3
     #   
     # }else if (project=="GBP"){
     #   studyvec[which(usematrix[,1]%in%c(1))]<-1#berry
     #   studyvec[which(usematrix[,1]%in%c(2,3,4))]<-3#eu
     #   studyvec[which(usematrix[,1]%in%c(5))]<-4#he
     #   studyvec[which(usematrix[,1]%in%c(6,7,8,9,10,11,12))]<-2#impi
     #   
     #   
     # }else if (project=="method"){
     #   studyvec[which(usematrix[,1]%in%c(1,2))]<-1
     #   studyvec[which(usematrix[,1]%in%c(3,4,5))]<-2
     #   studyvec[which(usematrix[,1]%in%c(6,7,8))]<-3
     #   
     # }
     
     print("studyvec")
     print(studyvec)
     print("usematrix")
     print(usematrix)
     usematrix2<-cbind(studyvec,usematrix)
     print("usematrix2")
     print(usematrix2)
     #squarenames<-rownames(datasetmatrix)[usematrix2[,2]]
     squarenames<-rownames(datasetmatrix)
     #print("squarenames")
     #print(squarenames)
     #\new
     cells<-cellfun()
     #print(cells)
     cellindexall<-1:length(allcells)
     #print("cellindexall")
     #print(cellindexall)
     cellindex<-cellindexall[which(allcells%in%cells)]
     #print("cellindex")
     #print(cellindex)
     #loop to make datamatrix
     relevantPathways<-character()
     colnamevector<-character()
     
     print("begin for loop1")
     for (thecell in 1:length(cellindex)){
       print("thecell")
       pastecell<-allcells[cellindex[thecell]]
       print(pastecell)
       for (row in 1:nrow(usematrix2)){
         #print("row")
         #print(row)
         #count<-count+1
         iterstudy<-usematrix2[row,1]
         #print("iterstudy")
         #print(iterstudy)
         #if(iterstudy==1){itersquare<-usematrix2[row,2]}else{itersquare<-usematrix2[row,2]-3}
         itersquare<-squarenames[usematrix2[row,2]]
         print("itersquare")
         print(itersquare)
         iteredge<-usematrix2[row,3]
         print("iteredge")
         print(iteredge)
         pathwayvector<-studycellscorelist[[iterstudy]][[itersquare]][[iteredge]][[pastecell]]
         
         relevantPathways<-c(relevantPathways,names(pathwayvector))
         print(names(pathwayvector))
         colnamevector<-c(colnamevector,paste(
           names(studycellscorelist)[iterstudy],
           #names(studycellscorelist[[iterstudy]])[itersquare],
           itersquare,
           names(studycellscorelist[[iterstudy]][[itersquare]])[iteredge],
           pastecell,
           sep="_"
         ))
         #print("internal colnamevector")
         #print(colnamevector)
       }#end usematrix for loop
     }#end cell for loop
     print("finished for loop1")
     
     relevantPathways<-sort(unique(relevantPathways))
     print("colnamevector")
     print(colnamevector)
     print("relevantPathways")
     print(relevantPathways)
     
     
     matrixdim<-c(length(relevantPathways),length(colnamevector))
     print("matrixdim")
     print(matrixdim)
     datamatrix<-matrix(nrow=matrixdim[1],ncol=matrixdim[2],data=0)
     dimnames(datamatrix)<-list(relevantPathways,colnamevector)
     colcount<-1
     print("begin for loop2")
     for (thecell in 1:length(cellindex)){
       print("thecell_round2")
       pastecell<-allcells[cellindex[thecell]]
       print(pastecell)
       for (row in 1:nrow(usematrix2)){
         print("row")
         print(row)
         iterstudy<-usematrix2[row,1]
         print("iterstudy")
         print(iterstudy)
         #if(iterstudy==1){itersquare<-usematrix2[row,2]}else{itersquare<-usematrix2[row,2]-3}
         itersquare<-squarenames[usematrix2[row,2]]
         print("itersquare")
         print(itersquare)
         iteredge<-usematrix2[row,3]
         print("iteredge")
         print(iteredge)
         print("class iteredge")
         print(class(iteredge))
         print("test1")
         print("testx")
         print(paste("iterstudy",iterstudy,"itersquare",itersquare,"iteredge",iteredge,"pastecell",pastecell))
         pathwayvectordata<-studycellscorelist[[iterstudy]][[itersquare]][[iteredge]][[pastecell]]
         print(pathwayvectordata)
         if(is.null(pathwayvectordata)){pathwayvectordata<-rep(0,length(relevantPathways));names(pathwayvectordata)<-relevantPathways}
         # print("test2")
         # print(relevantPathways%in%names(pathwayvectordata))
         # print("test3")
         # print(relevantPathways[which(relevantPathways%in%names(pathwayvectordata))])
         # print("test4")
         # print(relevantPathways)
         # print(names(pathwayvectordata))
         # print("test5")
         # print(sum(relevantPathways%in%names(pathwayvectordata)))
         # print("length(relevantPathways)")
         # print(length(relevantPathways))
         # print("length(pathwayvectordata)")
         # print(length(pathwayvectordata))
         # print("colcount")
         # print(colcount)
         datamatrix[which(relevantPathways%in%sort(names(pathwayvectordata))),colcount]<-pathwayvectordata[order(names(pathwayvectordata))]
         print("datamatrix[,colcount]")
         #print(dimnames(datamatrix))
         #print(datamatrix[,colcount])
         colcount<-colcount+1
       }
       
     }
     print("finished for loop2")
     
     print("nrow(datamatrix)")
     print("ncol(datamatrix)")
     print(nrow(datamatrix))
     print(ncol(datamatrix))
     
     datamatrix[is.na(datamatrix)]<-0
     breaksvals<-max(abs(datamatrix))
     print("breaksvals")
     print(breaksvals)
     tv2<-apply(datamatrix,1,function(ff){max(abs(ff),na.rm=TRUE)})>slider1()
     tv2[is.na(tv2)]<-0
     #tv3<-apply(datamatrix,2,function(ff){max(abs(ff),na.rm=TRUE)})>slider2()
     #tv3[is.na(tv3)]<-0
     plotMatrix<-datamatrix
     #plotMatrix<-datamatrix[as.logical(tv2),as.logical(tv3)]
     plotMatrix<-datamatrix[as.logical(tv2),]
     
     #cex calculation
     cexval1<-3*sqrt(28/ncol(plotMatrix))
     cexval2<-2*sqrt(120/nrow(plotMatrix))
     print("ready to plot")
     print(head(plotMatrix))
     cur_dev <- dev.cur()
     pheatout <- pheatmap(plotMatrix,col=blueWhiteRed(60),breaks=seq(floor(-breaksvals),ceiling(breaksvals),length.out=61),fontsize_col=9+cexval1,fontsize_row=9+cexval2)
     dev.set(cur_dev)
     print(pheatout)
     #pheatmap(datamatrix,col=blueWhiteRed(60),breaks=seq(floor(-breaksvals),ceiling(breaksvals),length.out=61))
     if(input$returnpdf==TRUE){
       pdf("plot.pdf",width=16,height=as.numeric(input$plotheightG2)/100,onefile=FALSE)
       pheatmap(plotMatrix,col=blueWhiteRed(60),breaks=seq(floor(-breaksvals),ceiling(breaksvals),length.out=61),fontsize_col=9+cexval1,fontsize_row=9+cexval2)
       #pheatmap(plotMatrix,col=blueWhiteRed(60),breaks=seq(floor(-breaksvals),ceiling(breaksvals),length.out=61),fontsize_col=9+cexval1,fontsize_row=9+cexval2,main=paste("study"))
       
       
       dev.off()
     }
   })
   
  output$cellmatrixMeta <- renderUI({
    plotOutput("Gigamat2", height = as.numeric(input$plotheightG2))
  })
  
  output$Gigamat2 <- renderPlot({
    plotInputGigamat2()
  })
  
  #ratioFarm ####
  output$ratioFarm<-renderPlot({
    usematrix<-usematrixfun()
    usematrix<-usematrix[order(usematrix[,2]),]
    usematrix<-usematrix[order(usematrix[,1]),]
    print("begin loop")
    print(usematrix)
    
    edgenames<-1:5
    squarenames<-rownames(datasetmatrix)
    cellnames3
    
    dfr<-matrix(ncol=length(cellnames3),nrow=nrow(usematrix),data=NA)
    colnames(dfr)<-cellnames3
    rownamesdfr<-c()
    for (row in 1:nrow(usematrix)){
      sq<-squarenames[usematrix[row,1]]
      ed<-edgenames[usematrix[row,2]]
      query=paste("MATCH (c:cellprop) WHERE c.square='",sq,"' AND c.edge='",ed,"' RETURN c.name AS name, toFloat(c.ratio) AS ratio",sep="")
      res<-cypher(graph,query)
      print(res)
      dfr[row,res$name]<-res$ratio
      #datasetrow<-datasetmatrix[sq,]
      rownamesdfr<-c(rownamesdfr,paste(sq,datasetmatrix[sq,as.numeric(ed)],sep="_"))
      
    
      
    }
    rownames(dfr)<-rownamesdfr
    dfr[is.na(dfr)]<-1
    dfr[dfr==0]<-0.000001
    print("debug dfr")
    print(str(dfr))
    dfr<-log2(dfr)
    print(dfr)
    dfr<-dfr[,colSums(dfr)!=0]
    print(dfr)
    dfr[dfr>3]<-3
    dfr[dfr<(-3)]<-(-3)
    breaksvals<-max(abs(dfr))
    print(dfr)
    #breaksvals<-
    
    cur_dev <- dev.cur()
    theres<-pheatmap(dfr,col=blueWhiteRed(60),breaks=seq(floor(-breaksvals),ceiling(breaksvals),length.out=61),fontsize_row = 10 * factorfun(),fontsize_col = 10 * factorfun(),display_numbers = TRUE,number_format="%.1f",fontsize_number = 10 * factorfun())
    dev.set(cur_dev)
    print(theres)
  })
  
    
   output$downloadPlot <- downloadHandler(
     filename = function() { paste("ANIMAplot",Sys.time(),".pdf",sep="") },
     content = function(file) {
       file.copy("plot.pdf",file)
     }
   )  
   
})

# Run the application 
shinyApp(ui = ui, server = server)

