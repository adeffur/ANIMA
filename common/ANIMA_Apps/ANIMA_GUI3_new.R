#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

#File paths
tabledir="~/output/project/tabular"
cytodir="~/output/project/cytoscape"
venndir="~/output/project/vennDiagrams"
figdir="~/output/project/linePlots"
igraphdir="~/output/project/igraphs"
datadir="~/output/RData"
#cellcordir="~/output/project/cellcor"
modcordir="~/output/project/modcor"

for (dir in c(tabledir, cytodir,venndir,figdir,igraphdir,datadir)) {
  if (!file.exists(dir)) {
    dir.create(dir, recursive=T,showWarnings=T)
  }
}



library(shiny)

library(networkD3)
library(htmlwidgets)
library(igraph)
library(WGCNA)
library(raster)

library(DT)
library(pheatmap)

library(illuminaHumanv3.db)
library(illuminaHumanv4.db)

library("lumi", lib.loc="/usr/local/lib/R/site-library")
library(Heatplus)
# library(wordcloud)
# library(tm)
# library(grImport)
scripts<-system(paste("ls","~/scripts"),intern=TRUE)
datafiles<-system(paste("ls",datadir),intern=TRUE)
excl<-grep("86|393",datafiles)
datafiles<-datafiles[!1:length(datafiles)%in%excl]
source(file.path("~/source_data/questions.R"),local=TRUE)
source("~/source_data/setlist.R",local=TRUE)
library(RNeo4j)
graph<-startGraph(graphstring)#!!!!! note the address!!

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

#R-scripts root
dir.R.files <- file.path(dir.common, "03_Functions")
## R utilities folder
dir.util <- file.path(dir.R.files, 'util')

#Data root
dir.data_root<-file.path(dir.root,'source_data')

#annotation and module data
dir.annot<-file.path(dir.root,"source_data")

#Individual subject distributions based on 2 modules

buildlisting<-system(paste("ls","~/output/build/"),intern=TRUE)

source(file.path(dir.R.files,"igraph_plotter_newer.R"),local=TRUE)
source(file.path(dir.R.files,"moduleMeta.R"),local=TRUE)
source(file.path("~/source_data/questions.R"),local=TRUE)
source(file.path(dir.R.files,"probe_boxplot4.R"),local=TRUE)
source(file.path(dir.R.files,"mwat.R"),local=TRUE)

##Required for G1

#new
source("~/source_data/setlist.R",local=TRUE)
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
all<-1:258
print("debug start")
print(data$sets)
neutromods<-cypher(graph,"MATCH (b:baylor {edge:5})-[r]-(c:cellEx)-[r2]-(c2:CELL) WHERE c2.name = 'neutrophil' RETURN DISTINCT b.name")$b.name
monomods<-cypher(graph,"MATCH (b:baylor {edge:5})-[r]-(c:cellEx)-[r2]-(c2:CELL) WHERE c2.name = 'monocyte' RETURN DISTINCT b.name")$b.name
bcellmods<-cypher(graph,"MATCH (b:baylor {edge:5})-[r]-(c:cellEx)-[r2]-(c2:CELL) WHERE c2.name = 'B-cell' RETURN DISTINCT b.name")$b.name
cd4mods<-cypher(graph,"MATCH (b:baylor {edge:5})-[r]-(c:cellEx)-[r2]-(c2:CELL) WHERE c2.name = 'CD4+ T-cell' RETURN DISTINCT b.name")$b.name
pltmods<-cypher(graph,"MATCH (b:baylor {edge:5})-[r]-(c:cellEx)-[r2]-(c2:CELL) WHERE c2.name = 'platelet' RETURN DISTINCT b.name")$b.name
nkmods<-cypher(graph,"MATCH (b:baylor {edge:5})-[r]-(c:cellEx)-[r2]-(c2:CELL) WHERE c2.name = 'NK-cell' RETURN DISTINCT b.name")$b.name
lymmods<-cypher(graph,"MATCH (b:baylor {edge:5})-[r]-(c:cellEx)-[r2]-(c2:CELL) WHERE c2.name = 'lymphocyte' RETURN DISTINCT b.name")$b.name
rbcmods<-cypher(graph,"MATCH (b:baylor {edge:5})-[r]-(c:cellEx)-[r2]-(c2:CELL) WHERE c2.name = 'RBC' RETURN DISTINCT b.name")$b.name

#functions and pathways
ifnmods<-cypher(graph,"MATCH (b:baylor {edge:5})-[r]-(pw) WHERE pw.name =~ '(?i).*interferon.*|.*IFN.*|.*ISG15.*' RETURN DISTINCT b.name")$b.name
lysosomemods<-cypher(graph,"MATCH (b:baylor {edge:5})-[r]-(pw) WHERE (pw:reactomePW OR pw:ImmunePW OR pw:PalWangPW) AND pw.name =~ '(?i).*lysosom[a,e].*' RETURN DISTINCT b.name")$b.name
phagocytmods<-cypher(graph,"MATCH (b:baylor {edge:5})-[r]-(pw) WHERE (pw:reactomePW OR pw:ImmunePW OR pw:PalWangPW) AND pw.name =~ '(?i).*phagocyt[o,i].*' RETURN DISTINCT b.name")$b.name
antprocmods<-cypher(graph,"MATCH (b:baylor {edge:5})-[r]-(pw) WHERE (pw:reactomePW OR pw:ImmunePW OR pw:PalWangPW) AND pw.name =~ '(?i).*antigen.*proc.*pres.*' RETURN DISTINCT b.name")$b.name
glycolmods<-cypher(graph,"MATCH (b:baylor {edge:5})-[r]-(pw) WHERE (pw:reactomePW OR pw:ImmunePW OR pw:PalWangPW) AND pw.name =~ '(?i).*glycoly[s,t].*' RETURN DISTINCT b.name")$b.name
oxphosmods<-cypher(graph,"MATCH (b:baylor {edge:5})-[r]-(pw) WHERE (pw:reactomePW OR pw:ImmunePW OR pw:PalWangPW) AND pw.name =~ '(?i).*oxidat.*phos.*' RETURN DISTINCT b.name")$b.name
ribosomemods<-cypher(graph,"MATCH (b:baylor {edge:5})-[r]-(pw) WHERE (pw:reactomePW OR pw:ImmunePW OR pw:PalWangPW) AND pw.name =~ '(?i).*ribosom[a,e].*' RETURN DISTINCT b.name")$b.name
lipidmods<-cypher(graph,"MATCH (b:baylor {edge:5})-[r]-(pw) WHERE (pw:reactomePW OR pw:ImmunePW OR pw:PalWangPW) AND pw.name =~ '(?i).*lipid.*|.*cholesterol.*' RETURN DISTINCT b.name")$b.name
prolifmods<-cypher(graph,"MATCH (b:baylor {edge:5})-[r]-(pw) WHERE (pw:reactomePW OR pw:ImmunePW OR pw:PalWangPW) AND pw.name =~ '(?i).*cell.*cycle*|.*mitosis*|.*checkpoint.*|.*cyclin.*' RETURN DISTINCT b.name")$b.name


subset.list<-list("All"=all,"Neutrophils"=neutromods,"Monocytes"=monomods,"CD4 T cells"=cd4mods,"B cells"=bcellmods,"NK cells"=nkmods,"Lymphocytes"=lymmods,
                  "Red blood cells"=rbcmods,"Platelets"=pltmods,"Interferon signaling"=ifnmods,"Lysosome"=lysosomemods,
                  "Phagocytosis"=phagocytmods,"Antigen proc. pres."=antprocmods,"Glycolysis"=glycolmods,"Oxidative phosphorylation"=oxphosmods,"Ribosome"=ribosomemods,"Lipids"=lipidmods,"Cell proliferation"=prolifmods)

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
  edgesChar<-cypher(graph,paste("MATCH (n:wgcna {square:'",set,"',edge:",edge,"}) RETURN DISTINCT n.contrast AS comparison",sep=""))
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
annot<-read.csv(file.path("~/source_data","inflammasomes.csv"))
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

#Required for cellcore3
source(file.path(dir.R.files,"cellCore3.R"),local=TRUE)

query<-"MATCH (c:cellEx) RETURN DISTINCT c.name AS cellname2"
cellnamelist2<-cypher(graph,query)
cellnames2<-cellnamelist2$cellname2

query<-"MATCH (c:cellprop) RETURN DISTINCT c.name AS cellname3"
cellnamelist3<-cypher(graph,query)
cellnames3<-cellnamelist3$cellname3


  
# Define UI 
ui<-fluidPage(
  #shinythemes::themeSelector(),
              titlePanel(paste("ANIMA REGO:",project)),
              sidebarLayout(
                sidebarPanel(
                  tabsetPanel(
                    tabPanel("Gene",
                      tags$head(tags$style("#boxplot{height:100vh !important; overflow-y: scroll}")),
                      h3(paste('boxplot:',project),style="color:#FF0000"),
                      fluidRow(column(8,uiOutput("choosesquareBXP"),offset=.33),
                               column(4,selectInput("orderP","Boxplot order",list("Group contrast 1"="c(1,3,2,4)","Group contrast2"="c(1,2,3,4)"),selected=1,multiple=FALSE))),
                      fluidRow(
                        column(8,textInput("regex","Write regex",value="STAT.*")),
                        column(4,sliderInput("plotheight","plot height",1000,6000,1000,500),offset=0.33)
                      ),
                      fluidRow(
                        column(4,sliderInput("ymin","y min",0,20,6,1)),
                        column(4,sliderInput("ymax","y max",0,20,14,1),offset=1)),
                      fluidRow(checkboxInput('returncsv2', 'output csv of differential expression?', FALSE))
                               ),
                  
                    tabPanel("Machines",
                             tags$head(tags$style("#inflammasomes{height:90vh !important;}")),
                             h3(paste('Inflammasomes:',project),style="color:#FF0000"),
                             fluidRow(column(4,uiOutput("choosesquareMACHINE"),offset=.33)),
                             fluidRow(
                               column(8,selectInput("inflammasome","Select inflammasome",names(generegexlist),selected = names(generegexlist)[1],multiple=FALSE)),
                               column(4,sliderInput("plotheightMachine","plotheightMachine",500,1500,700,150))
                             )
                    ),
                    #INTERFACE ELEMENTS: Tabs
                  tabPanel("Modules",
                      tags$head(tags$style("#igraphSP{height:85vh !important;}")),
                      tags$head(tags$style("#d3graphSP{height:85vh !important;}")),
                      tags$head(tags$style("#subjectplot{height:85vh !important; overflow-y: scroll}")),
                      tags$head(tags$style("#module2d{height:85vh !important; overflow-y: scroll}")),
                      tags$head(tags$style("#modPheno{height:85vh !important; overflow-y: scroll}")),
                      tags$head(tags$style("#intracor{height:85vh !important; overflow-y: scroll}")),
                      tags$head(tags$style("#projections{height:85vh !important;}")),
                      tags$head(tags$style("#MEvar{height:85vh !important;}")),
                      tags$head(tags$style("#sankey{height:85vh !important;}")),
                      tags$head(tags$style("#modcor{height:85vh !important;overflow-y: scroll}")),
                  h3("Datasets",style="color:#FF0000"),
                  fluidRow(
                  column(4,uiOutput("allsubsetnames")),
                  
                  column(4,selectInput("matchedset","Choose matchedset",c("Q_9_blood.PCF.defTB","Q_10_blood.PCF.probTB","Q_11_blood.PCF.defPC","Q_12_blood.PCF.probPC","Q_13_blood.PCF.HIVneg","Q_14_blood.PCF.HIVpos","Q_15_blood.PCF.HIVposHD"),selected="Q_9_blood.PCF.defTB",multiple=FALSE),offset=.33)
                  #column(4,selectInput("matchedset","Choose matchedset",c("Q_6_IMPI.TBPC.BF"),selected="Q_6_IMPI.TBPC.BF",multiple=FALSE),offset=.33)#GBP
                  
                  ),
                  fluidRow(
                    column(4,uiOutput("choosesquare"),offset=.33),
                    column(2,selectInput("edge","Choose edge",1:5,selected=5,multiple=FALSE),offset=.33)
                  ),
                  h3("Modules",style="color:#FF0000"),
                  fluidRow(
                    column(3, uiOutput("menamegen1")),
                    column(3,uiOutput("menamegen2"))
                  ),
                  h3("2D",style="color:#FF0000"),
                  fluidRow(
                    column(3, uiOutput("phenonames")),
                    column(3,selectInput("method","Correlation method",c("pearson","spearman"),selected="pearson")),
                    column(3,selectInput("legpos","Place legend",c("bottomleft","bottomright","topleft","topright")))
                  ),
                  fluidRow(
                    column(3, uiOutput("subsetdefs1")),
                    column(3,uiOutput("subsetdefs2"))
                  ),
                  h3("Bipartite",style="color:#FF0000"),
                  fluidRow(
                    column(3,selectInput("nodetype_proj","Choose nodetype 1",c("wgcna","baylor","reactomePW","ImmunePW","PalWangPW","cellEx","cellprop","pheno"),selected="wgnca",multiple=FALSE)),
                    column(3,selectInput("nodetype_proj2","Choose nodetype 2",c("wgcna","baylor","reactomePW","ImmunePW","PalWangPW","cellEx","cellprop","pheno"),selected="baylor",multiple=FALSE)),
                    column(3,selectInput("plotwhich","Choose which plot",c("graph"=1,"proj1"=2,"proj2"=3),multiple=FALSE))
                  ),
                  fluidRow(
                    column(3,selectInput("layout","Choose layout",c("layout_with_fr","layout_with_dh","layout_nicely","layout_with_gem","layout_as_bipartite"),selected="layout_with_fr",multiple=FALSE)),
                    column(3,sliderInput("vlc","label size adjustment",0.5,2,1,0.001)),
                    column(3,sliderInput("legendcex","legend size adjustment",0.5,2,1,0.001)),
                    column(3,sliderInput("vertexsize","vertex size",2,30,15,0.1))
                  ),
                  h3("Sankey plots",style="color:#FF0000"),
                  fluidRow(
                    column(3,selectInput("nodetype_sankey1","Choose nodetype 1",c("wgcna","ImmunePW","cellEx"),selected="wgnca",multiple=FALSE)),
                    column(3,selectInput("nodetype_sankey2","Choose nodetype 2",c("ImmunePW","cellEx"),selected="cellEx",multiple=FALSE)),
                    column(3,selectInput("linkgroup","Choose link group",c("node1","node2")))),
                  # fluidRow(  
                  #   column(4,conditionalPanel(condition = "input.nodetype_sankey1 == 'cellEx'",selectInput("cellSet1","Choose cellSet",c("xp_1","xp_2","xp_3"),selected=1,multiple=FALSE))),
                  #   column(4,conditionalPanel(condition = "input.nodetype_sankey2 == 'cellEx'",selectInput("cellSet2","Choose cellSet",c("xp_1","xp_2","xp_3"),selected=1,multiple=FALSE)))
                  # ),
                  
                  h3(paste('Modcor'),style="color:#FF0000"),
                  # fluidRow(
                  #   column(4,selectInput("studyModcor","Choose study",studies,selected=studies[1],multiple=FALSE)),
                  #   column(4,conditionalPanel(
                  #     condition = "input.studyModcor == 1",
                  #     selectInput("datasets1MC","Choose dataset",setlist1,selected=setlist1[1],multiple=FALSE)
                  #   ),
                  #   conditionalPanel(
                  #     condition = "input.studyModcor == 2",
                  #     selectInput("datasets2MC","Choose dataset",setlist2,selected=setlist2[1],multiple=FALSE)
                  #   ),
                  #   conditionalPanel(
                  #     condition = "input.studyModcor == 3",
                  #     selectInput("datasets3MC","Choose dataset",setlist3,selected=setlist3[1],multiple=FALSE)
                  #   ),
                  #   conditionalPanel(
                  #     condition = "input.studyModcor == 4",
                  #     selectInput("datasets4MC","Choose dataset",setlist4,selected=setlist4[1],multiple=FALSE)),offset=.5)
                  # ),
                  
                  #new
                  fluidRow(
                    column(4,selectInput("studyModcor","Choose study",studies,selected=studies[1],multiple=FALSE)),
                    column(4,uiOutput("selectMCset"))
                    ),
                    
                  #end new
                  
                  fluidRow(
                    column(4,uiOutput("selectModuleModcor")),
                    column(4,selectInput("edgeModcor","select edge",1:5,selected=1,multiple=FALSE))
                  )
                  ),#endSP tab
                  
                  tabPanel("Virtual Cells",
                           tags$head(tags$style("#cellcor{height:90vh !important;overflow-y: scroll}")),
                           tags$head(tags$style("#gigamat{height:90vh !important;overflow-y: scroll}")),
                           tags$head(tags$style("#gigabar{height:80vh !important;}")),
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
                  fluidRow(
                    column(4,selectInput("study","Choose study",studies,selected=studies[1],multiple=FALSE)),
                    column(4,uiOutput("selectCCset"))
                  ),
                  
                  #end new
                           fluidRow(
                             column(4,selectInput("cellname","Choose cell name",cellnames2)),
                             column(2,selectInput("edges","select edge",1:5,selected=1,multiple=FALSE))
                           ),

                  h3(paste('Cell Matrix'),style="color:#FF0000"),
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
                  
                  h3(paste('Gigamatrix2:',project),style="color:#FF0000"),
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
                  
                  
                  
                  tabPanel("Save",br(),fluidRow(column(4,checkboxInput('returnpdf', 'output pdf?', FALSE)),column(4,downloadButton("downloadPlot","Save plot")),offset=1))
                  , selected = "Modules")#end tabset
                
                
                ),#end sidebar
                if(doIntracor==TRUE){mainPanel(
                  tabsetPanel(
                    #BXP
                    tabPanel("Gene",
                             tabsetPanel("BXP",
                                         tabPanel("boxplot",plotOutput("boxplot")),
                                         tabPanel("detableBXP",DT::dataTableOutput('detableBXP'))
                             )
                             
                    ),
                    #Machines
                    tabPanel("Machines",plotOutput("inflammasomes")),
                    tabPanel("Single module",
                      tabsetPanel("Single module",
                                  tabPanel("igraphSP",plotOutput("igraphSP")),
                                  tabPanel("d3graphSP",forceNetworkOutput("d3graphSP",width="100%",height="400px")),
                                  tabPanel("moduleNodesSP",DT::dataTableOutput('moduleNodesSP')),
                                  tabPanel("moduleEdgesSP",DT::dataTableOutput('moduleEdgesSP')),
                                  tabPanel("modcor",plotOutput("modcor")),
                                  tabPanel("MEvar",plotOutput("MEvar")),
                                  tabPanel("modPheno",plotOutput("modPheno")),
                                  tabPanel("subjectPlot",plotOutput("subjectplot"))
                                  , selected = "d3graphSP")),
                      
                    tabPanel("All modules",
                        tabsetPanel("All modules",
                             tabPanel("module2d",plotOutput("module2d")),
                             tabPanel("intracor",plotOutput("intracor")),
                             tabPanel("sankey",sankeyNetworkOutput("sankey")),
                             tabPanel("projections",plotOutput("projections"))
                             )),
                    
                    tabPanel("Virtual Cells",
                             tabsetPanel("cellcor",
                                         tabPanel("Single Cell",plotOutput("cellcor")),
                                         tabPanel("Cell Matrix",plotOutput("gigamat")),
                                         tabPanel("Cell and pathway summary",plotOutput("gigabar"))
                             )),
                    
                    
                    tabPanel("Meta-Analysis",
                             tabsetPanel("Chaussabel",
                                         tabPanel("ModuleMeta",plotOutput("pheat")),
                                         tabPanel("igraph",plotOutput("igraph")),
                                         tabPanel("d3graphMM",forceNetworkOutput("d3graphMM",width="100%",height="100%")),
                                         tabPanel("moduleNodesMM",DT::dataTableOutput('moduleNodesMM')),
                                         tabPanel("moduleEdgesMM",DT::dataTableOutput('moduleEdgesMM')),
                                         tabPanel("wgcna",plotOutput("wgcnacol")),
                                         tabPanel("Virtual Cells",plotOutput("gigamat2")),
                                         tabPanel("Ratio Farm",plotOutput("ratioFarm")))
                             
                    )
                    
                    
                    
                    ,selected="Single module")
                )}else{mainPanel(
                  tabsetPanel(
                    tabPanel("SP",
                             tabsetPanel("SP",
                                         tabPanel("module2d",plotOutput("module2d")),
                                         tabPanel("subjectPlot",plotOutput("subjectplot")),
                                         tabPanel("modPheno",plotOutput("modPheno")),
                                         tabPanel("d3graphSP",forceNetworkOutput("d3graphSP",width="100%",height="100%")),
                                         tabPanel("moduleNodesSP",DT::dataTableOutput('moduleNodesSP')),
                                         tabPanel("moduleEdgesSP",DT::dataTableOutput('moduleEdgesSP')),
                                         #tabPanel("intracor",plotOutput("intracor")),
                                         tabPanel("projections",plotOutput("projections")),
                                         tabPanel("sankey",sankeyNetworkOutput("sankey")),
                                         tabPanel("MEvar",plotOutput("MEvar")))),
                    tabPanel("modcor",
                             tabsetPanel("modcor",
                                         tabPanel("modcor",plotOutput("modcor"))
                             )),
                    tabPanel("cellcor",
                             tabsetPanel("cellcor",
                                         tabPanel("cellcor",plotOutput("cellcor"))
                             )),
                    tabPanel("G1",
                             tabsetPanel("G1",
                                         tabPanel("Gigamat1",plotOutput("gigamat")),
                                         tabPanel("barplot",plotOutput("gigabar")))),
                    #BXP
                    tabPanel("BXP",
                             tabsetPanel("BXP",
                                         tabPanel("boxplot",plotOutput("boxplot")),
                                         tabPanel("detableBXP",DT::dataTableOutput('detableBXP'))
                             )
                             
                    ),
                    tabPanel("MM",
                             tabsetPanel("MM",
                                         tabPanel("ModuleMeta",plotOutput("pheat")),
                                         tabPanel("igraph",plotOutput("igraph")),
                                         tabPanel("d3graphMM",forceNetworkOutput("d3graphMM",width="100%",height="100%")),
                                         tabPanel("moduleNodesMM",DT::dataTableOutput('moduleNodesMM')),
                                         tabPanel("moduleEdgesMM",DT::dataTableOutput('moduleEdgesMM')),
                                         tabPanel("wgcna",plotOutput("wgcnacol")))),
                    #G2
                    tabPanel("G2",plotOutput("gigamat2")),
                    #Machines
                    tabPanel("Machines",plotOutput("inflammasomes"))
                  )
                )
                }
                  
              )
)


# Define server logic 
server <- shinyServer(function(input, output) {
  #libraries
  library(RNeo4j)
  graph<-startGraph(graphstring)#!!!!! note the address!!
  #REACTIVE FUNCTIONS
   #SP
   setfun<-reactive({input$set})
   setfun2<-reactive({input$matchedset})
   me1fun<-reactive({input$me1})
   me2fun<-reactive({input$me2})
   s1fun<-reactive({input$subset1})
   s2fun<-reactive({input$subset2})
   phenofun<-reactive({input$pheno})
   methodfun<-reactive({input$method})
   edgefun<-reactive({input$edge})
   buildfun<-reactive({input$build})
   nodefun<-reactive({input$nodetype_proj})
   nodefun2<-reactive({input$nodetype_proj2})
   plotwhichfun<-reactive({input$plotwhich})
   legposfun<-reactive({input$legpos})
   
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
   
   
   #G1
   usestudy<-reactive({input$study})
   useset<-reactive({c(input$datasetsCC)})
   useedge<-reactive({input$edges})
   slider1<-reactive({input$sliderpw})
   slider2<-reactive({input$slidercell})
   barfilterfun1<-reactive({input$barfilter1})
   barfilterfun2<-reactive({input$barfilter2})
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
   
   #igraph_plotter controls single module and bipartite
   layoutfun<-reactive({input$layout})
   vlc<-reactive({input$vlc})
   legendcexfun<-reactive({input$legendcex})
   vertexsizefun<-reactive({input$vertexsize})
   
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
     options=list(lengthMenu = c(5, 10, 15,20,50), pageLength = 5)
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
   
   #Machines
   grxfun<-reactive({input$inflammasome})
   setfunMACHINE<-reactive({input$MACHINEset})
   #REACTIVE UI INPUT FUNCTIONS
   #SP
   output$allsubsetnames<-renderUI({
     prefix<-paste("~/output/build/",buildfun(),"/",sep="")
     dirlisting<-system(paste("ls",prefix),intern=TRUE)
     allsubsetnames<-dirlisting[grep("Q_",dirlisting)]
     print("prefix:")
     print(prefix)
     print("dirlisting")
     print(dirlisting)
     print("allsubsetnames")
     print(allsubsetnames)
     selectInput("build","Choose build",buildlisting,selected=buildlisting[length(buildlisting)],multiple=FALSE)
   })
   
   output$choosesquare<-renderUI({
     print("invoking input$set in renderUI function:output$choosesquare")
     prefix<-paste("~/output/build/",buildfun(),"/",sep="")
     dirlisting<-system(paste("ls",prefix),intern=TRUE)
     allsubsetnames<-dirlisting[grep("Q_",dirlisting)]
     selectInput("set","Choose square",allsubsetnames,selected="1",multiple=FALSE)
   })
   
   output$choosesquareBXP<-renderUI({
     print("invoking input$set in renderUI function:output$choosesquare")
     prefix<-paste("~/output/build/",buildfun(),"/",sep="")
     dirlisting<-system(paste("ls",prefix),intern=TRUE)
     allsubsetnames<-dirlisting[grep("Q_",dirlisting)]
     selectInput("BXPset","Choose square",allsubsetnames,selected="1",multiple=FALSE)
   })
   
   output$choosesquareMACHINE<-renderUI({
     print("invoking input$set in renderUI function:output$choosesquare")
     prefix<-paste("~/output/build/",buildfun(),"/",sep="")
     dirlisting<-system(paste("ls",prefix),intern=TRUE)
     allsubsetnames<-dirlisting[grep("Q_",dirlisting)]
     selectInput("MACHINEset","Choose square",allsubsetnames,selected="1",multiple=FALSE)
   })
   
   output$menamegen1<-renderUI({
     print("computing menamegen1")
     prefix<-paste("~/output/build/",buildfun(),"/",sep="")
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
     prefix<-paste("~/output/build/",buildfun(),"/",sep="")
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
     prefix<-paste("~/output/build/",buildfun(),"/",sep="")
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
   
   output$subsetdefs1<-renderUI({
     print("computing subsetdefs1")
     #print("line1 inside reactive")
     prefix<-paste("~/output/build/",buildfun(),"/",sep="")
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
     selectInput("subset1","Choose subset1",subset1classes,selected=c(1,2),multiple=TRUE)
   })
   
   output$subsetdefs2<-renderUI({
     print("computing subsetdefs2")
     print("line1 inside reactive")
     #print(setfun2())
     prefix<-paste("~/output/build/",buildfun(),"/",sep="")
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
     selectInput("subset2","Choose subset2",subset2classes,selected=c(1,2),multiple=TRUE)
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
   #modcor
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

#PLOT AND DATA OUTPUT
   #SP
   output$subjectplot <- renderPlot({
     print("computing subjectplot")
     prefix<-paste("~/output/build/",buildfun(),"/",sep="")
     #print("line1 inside reactive")
     data0<-read.csv(paste(prefix,setfun(),"/5_WGCNA_4k_tv/results/ModuleEigengenes_all_samples_withPheno.csv",sep=""))
     #print("line2")
     data<-read.csv(paste(prefix,setfun(),"/5_WGCNA_4k_tv/results/ModuleEigengenes_all_samples.csv",sep=""))
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
     
     subset1vec<-eval(parse(text=paste("data0$",c1,sep="")))
     subset1classes<-sort(unique(subset1vec))
     subset2vec<-eval(parse(text=paste("data0$",c2,sep="")))
     subset2classes<-sort(unique(subset2vec))
     
     dev.new()
     sub<-which(subset1vec%in%s1fun()&subset2vec%in%s2fun())
     group1<-subset1vec[sub]
     group2<-subset2vec[sub]
     #print(data0[sub,])
     plot(eval(parse(text=paste(me2fun(),"~",me1fun(),sep=""))),data=data0[sub,],pch=20,col="grey",main=paste("Subjects:",square),ylab=me2fun(),xlab=me1fun())
     #text(eval(parse(text=paste(me2fun(),"~",me1fun(),sep=""))),data=data0[sub,],labels=paste(cL1[group1],cL2[group2]),col=group1+2*group2)
     text(eval(parse(text=paste(me2fun(),"~",me1fun(),sep=""))),data=data0[sub,],labels=data0[sub,"Row.names"],col=(group1+2*group2)-2,pos=2)
     legend("bottomleft",sort(c(paste(cL1[1],cL2[1]),paste(cL1[1],cL2[2]),paste(cL1[2],cL2[1]),paste(cL1[2],cL2[2]))),fill=c(3,5,4,6)-2)
     #legend("bottomleft",sort(c(paste(cL1[1],cL2[1]),paste(cL1[2],cL2[1]),paste(cL1[1],cL2[2]),paste(cL1[2],cL2[2]))),fill=c(3,4,5,6))
     #legend("bottomleft",c(paste(cL1[1],cL2[1]),paste(cL1[1],cL2[2]),paste(cL1[2],cL2[1]),paste(cL1[2],cL2[2])),fill=c(3,4,5,6))
     
     abline(h=0.0,col="blue")
     abline(v=0.0,col="blue")
     if(input$returnpdf==TRUE){
       pdf("plot.pdf",width=16,height=10,onefile=FALSE)
       plot(eval(parse(text=paste(me2fun(),"~",me1fun(),sep=""))),data=data0[sub,],pch=20,col="grey",main=paste("Subjects:",square),ylab=me2fun(),xlab=me1fun())
       text(eval(parse(text=paste(me2fun(),"~",me1fun(),sep=""))),data=data0[sub,],labels=data0[sub,"Row.names"],col=(group1+2*group2)-2,pos=2)
       legend("bottomleft",sort(c(paste(cL1[1],cL2[1]),paste(cL1[1],cL2[2]),paste(cL1[2],cL2[1]),paste(cL1[2],cL2[2]))),fill=c(3,5,4,6)-2)
       dev.off()
     }
     
     
   })
   
   output$module2d<-renderPlot({
     prefix<-paste("~/output/build/",buildfun(),"/",sep="")
     datarec1<-read.csv(paste(prefix,setfun(),"/5_WGCNA_4k_tv/results/modules_as_classifiers_AUC_glm.csv",sep=""))
     squarelist<-unlist(strsplit(setfun(),"_"))
     if(length(squarelist)==3){
       square<-squarelist[3]
     }else{square<-paste(squarelist[c(3,4)],collapse="_")}
     #print("square")
     #print(square)
     contrasts<-cypher(graph,paste("MATCH (n:wgcna {square:'",square,"'}) RETURN DISTINCT n.contrastvar AS contrastvars",sep=""))
     c1<-contrasts$contrastvars[1]
     c2<-contrasts$contrastvars[2]
     
     plot(modAUC1~modAUC2,data=datarec1,pch=20,col="grey",main=paste(square,": Modules differentiating ",c1," and ",c2,sep=""),xlab=paste(c2,"index"),ylab=paste(c1, "index"))#aspect ratio can be fixed with asp=1
     shadowtext(datarec1$modAUC1~datarec1$modAUC2,labels=datarec1$X,col=datarec1$X,font=2,cex=1.0,xpd=TRUE,bg="lightgrey")

     abline(h=0.75,col="blue")
     abline(v=0.75,col="blue")
     
     
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
   
   output$modPheno <- renderPlot({
     prefix<-paste("~/output/build/",buildfun(),"/",sep="")
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
     
     contrasts<-cypher(graph,paste("MATCH (n:wgcna {square:'",square,"'}) RETURN DISTINCT n.contrastvar AS contrastvars",sep=""))
     c1<-contrasts$contrastvars[1]
     c2<-contrasts$contrastvars[2]
     
     data0<-data0[order(data0[,c1]),][order(data0[,c2]),]
     
     classes<-cypher(graph,paste("MATCH (n:wgcna {square:'",square,"'}) RETURN DISTINCT n.contrast AS contrasts",sep=""))
     #cL1<-unlist(strsplit(classes$contrasts[1],"-"))[c(2,1)]
     #cL2<-unlist(strsplit(classes$contrasts[2],"-"))[c(2,1)]
     #print("debug")
     #print(classes$contrasts[1])
     #print(classes$contrasts[2])
     cL1<-unlist(strsplit(classes$contrasts[1]," - "))
     cL2<-unlist(strsplit(classes$contrasts[2]," - "))
     
     subset1vec<-eval(parse(text=paste("data0$",c1,sep="")))
     subset1classes<-sort(unique(subset1vec))
     
     subset2vec<-eval(parse(text=paste("data0$",c2,sep="")))
     subset2classes<-sort(unique(subset2vec))
     
     dev.new()
     sub<-which(subset1vec%in%s1fun()&subset2vec%in%s2fun())
     group1<-subset1vec[sub]
     group2<-subset2vec[sub]
     #print("doing res")
     #print(data0[sub,])
     themethod<-methodfun()
     res<-cor.test(eval(parse(text=paste("data0$",phenofun(),"[sub]",sep=""))),eval(parse(text=paste("data0$",me1fun(),"[sub]",sep=""))),method=themethod)
     #print("doing plot")

     #print(group1)
     #print(group2)
     #print(2*group1+group2)
     #print(2*data0[sub,c1]+data0[sub,c2])
     #print(paste(cL1[group1],cL2[group2]))

     #print(unique(2*group1+group2))
     #print(unique(paste(cL1[group1],cL2[group2])))
     #print(c(paste(cL1[1],cL2[1]),paste(cL1[1],cL2[2]),paste(cL1[2],cL2[1]),paste(cL1[2],cL2[2])))
     plot(eval(parse(text=paste(phenofun(),"~",me1fun(),sep=""))),data=data0[sub,],pch=20,cex=2,col=unlist(strsplit(me1fun(),"ME"))[2],main=paste(me1fun(),"\n",themethod,"coef.=",sprintf("%.3f",res$estimate),"Pval=",sprintf("%.3f",res$p.value)),xlab=me1fun(),ylab=phenofun())
     text(data0[sub,phenofun()]~data0[sub,me1fun()],labels=data0[sub,"Row.names"],col=(2*data0[sub,c1]+data0[sub,c2])-2,pos=2)
     legend(legposfun(),sort(c(paste(cL1[1],cL2[1]),paste(cL1[1],cL2[2]),paste(cL1[2],cL2[1]),paste(cL1[2],cL2[2]))),fill=c(3,4,5,6)-2)
     #legend("bottomleft",c(paste(cL1[1],cL2[1]),paste(cL1[1],cL2[2]),paste(cL1[2],cL2[1]),paste(cL1[2],cL2[2])),fill=c(3,4,5,6))
     abline(lm(eval(parse(text=paste(phenofun(),"~",me1fun(),sep=""))),data=data0[sub,]), col="red")
     
     
     if(input$returnpdf==TRUE){
       pdf("plot.pdf",width=9,height=9,onefile=FALSE)
       plot(eval(parse(text=paste(phenofun(),"~",me1fun(),sep=""))),data=data0[sub,],pch=20,cex=2,col=unlist(strsplit(me1fun(),"ME"))[2],main=paste(me1fun(),"\nPearsonR=",sprintf("%.3f",res$estimate),"Pval=",sprintf("%.3f",res$p.value)),xlab=me1fun(),ylab=phenofun())
       text(data0[sub,phenofun()]~data0[sub,me1fun()],labels=data0[sub,"Row.names"],col=(2*data0[sub,c1]+data0[sub,c2])-1,pos=2)
       legend(legposfun(),sort(c(paste(cL1[1],cL2[1]),paste(cL1[1],cL2[2]),paste(cL1[2],cL2[1]),paste(cL1[2],cL2[2]))),fill=c(3,4,5,6)-1)
       abline(lm(eval(parse(text=paste(phenofun(),"~",me1fun(),sep=""))),data=data0[sub,]), col="red")
       
       dev.off()
     }
     
     
   })
   
   if(doIntracor==TRUE){
   output$intracor<-renderPlot({
     #ME intracorrelates######
     #print("reading data")
     #print(setfun2())
     squarelist<-unlist(strsplit(setfun2(),"_"))
     if(length(squarelist)==3){
       square<-squarelist[3]
     }else{square<-paste(squarelist[c(3,4)],collapse="_")}
     #print("square")
     #print(square)
     prefix<-paste("~/output/build/",buildfun(),"/",sep="")
     data0<-read.csv(paste(prefix,setfun2(),"/5_WGCNA_4k_tv/results/ModuleEigengenes_all_samples_withPheno.csv",sep=""))
     data<-read.csv(paste(prefix,setfun2(),"/5_WGCNA_4k_tv/results/ModuleEigengenes_all_samples.csv",sep=""))
     #print("data is read")
     # #debug
     # data0<-read.csv("~/output/ModuleEigengenes_all_samples_withPheno.csv")
     # data<-read.csv("~/output/ModuleEigengenes_all_samples.csv")

     data<-data[order(data$X),]
     data.blood<-data[seq(1,nrow(data)-1,2),]
     data.fluid<-data[seq(2,nrow(data),2),]


     mes<-colnames(data.blood)[2:ncol(data.blood)]
     menames<-sort(colnames(data)[2:length(colnames(data))])
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
     write.csv(intracor,file.path(tabledir,paste(square,"intracor_pheno.csv",sep="_")))
     phenovec<-colnames(data0[,which(!colnames(data0)%in%menames&!colnames(data0)=="X"&!colnames(data0)=="Row.names")])
     #print(phenovec)
     pdata<-data0[order(data0$X),phenovec]
     pdatablood<-pdata[seq(1,nrow(pdata)-1,2),]
     pdatafluid<-pdata[seq(2,nrow(data),2),]
     intracor2<-data.frame()
     for (name in phenovec){
       #print(name)
       res<-"Skip"
       try(res<-cor.test(eval(parse(text=paste("pdatablood$",name,sep=""))),eval(parse(text=paste("pdatafluid$",name,sep=""))),method="p"))
       # plot(eval(parse(text=paste("data.blood$",name,sep=""))),eval(parse(text=paste("data.fluid$",name,sep=""))),col=unlist(strsplit(name,"ME"))[2],main=paste(name,"\nPearsonR=",sprintf("%.3f",res$estimate),"Pval=",sprintf("%.3f",res$p.value)),pch=19,xlab="Blood",ylab="Fluid",xlim=c(-0.25,0.25),ylim=c(-0.25,0.25))
       # abline(lm(eval(parse(text=paste("data.blood$",name,sep="")))~eval(parse(text=paste("data.fluid$",name,sep="")))), col="red")
       # print(name)
       #print("here is a res2")
       #print(res)
       if(res!="Skip"){
       intracor2<-rbind(intracor2,c(name,res$estimate,res$p.value))
       }
     }
     colnames(intracor2)<-c("Pheno","PearsonR","P-value")
     write.csv(intracor2,file.path(tabledir,paste(square,"intracor2_pheno.csv",sep="_")))


     #print("intracor is done")
     colnames(intracor)<-c("name","intracor","pval")
     query<-paste("MATCH (n:wgcna {square:'",square,"',edge:5}) RETURN n.name as name, n.modAUC1 as modAUC1, n.modAUC2 as modAUC2, n.diffME as diffME, n.sigenrich as sigenrich",sep="")
     moduletable<-cypher(graph,query)
     #print("moduletable is done")
     moduletable2<-merge(moduletable,intracor,by.x="name",by.y="name")
     #print("moduletable2 is done")
     write.csv(moduletable2,file.path(tabledir,paste(square,"moduletable2.csv",sep="_")))
     plot(modAUC1~intracor,data=moduletable2,pch=20,col="grey",main=paste(square,": Compartmentalised vs. representative processes",sep=""),xlab="Representation Index",ylab="Compartmentalisation Index")
     shadowtext(moduletable2$modAUC1~moduletable2$intracor,labels=moduletable2$name,col=moduletable2$name,font=2,cex=1.0,xpd=TRUE,bg="lightgrey")
     abline(h=0.75,col="blue")
     #abline(v=0.5,col="blue")
     abline(v=0.0,col="blue")
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
   }#end conditional intracor
   
   output$d3graphSP<-renderForceNetwork({
     #get module
     prefix<-paste("~/output/build/",buildfun(),"/",sep="")
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
     query.base<-paste("MATCH (c:CELL)-[r0]-(x1)-[r1]-(n:wgcna {square:'",square,"',edge:",edgefun(),",name:'",unlist(strsplit(me1fun(),"ME"))[2],"'})-[r2]-(x2) WHERE (x1:cellEx OR x1:cellprop) AND (x2:reactomePW OR x2:ImmunePW OR x2:pheno OR x2:PalWangPW)",sep="")
     print(query.base)
     nodelist<-c("c","x1","n","x2")
     edgetrips<-list(c("c","r0","x1"),c("x1","r1","n"),c("n","r2","x2"))
     #igraphplotter with igr output
     igr<-"FAIL"
     try(igr<-igraph_plotter(query.base,nodelist,edgetrips,rimpar="diffME",plot=FALSE,csv=TRUE,prefix=cytodir,filename = "igraphwgcna",plotd3=FALSE,return_graph=TRUE))
     print(class(igr))
     if(igr=="FAIL"){
       #query
       query.base<-paste("MATCH (n:wgcna {square:'",square,"',edge:",edgefun(),",name:'",unlist(strsplit(me1fun(),"ME"))[2],"'})-[r2]-(x2) WHERE (x2:cellEx OR x2:cellprop OR x2:reactomePW OR x2:ImmunePW OR x2:pheno OR x2:PalWangPW)",sep="")
       #print(query.base)
       nodelist<-c("n","x2")
       edgetrips<-list(c("n","r2","x2"))
       igr<-igraph_plotter(query.base,nodelist,edgetrips,rimpar="diffME",plot=FALSE,csv=TRUE,prefix=cytodir,filename = "igraphwgcna",plotd3=FALSE,return_graph=TRUE)
     }
     
     #get nodes and edges for DT
     currentnodes<-read.csv(file.path(cytodir,"igraphwgcna_nodes.csv"))
     #print(currentnodes)
     output$moduleNodesSP<-DT::renderDataTable(currentnodes,server = FALSE,options=list(lengthMenu = c(5, 10, 15,20,50), pageLength = 50)) 
     currentedges<-read.csv(file.path(cytodir,"igraphwgcna_edges.csv"))
     dict<-currentnodes[,c(2,3)]
     for(row in 1:nrow(dict)){
       print(row)
       currentedges[currentedges==as.character(dict[row,"node_id"])]<-as.character(dict[row,"nodename"])
     }
     #print(currentedges)
     output$moduleEdgesSP<-DT::renderDataTable(currentedges,server = FALSE,options=list(lengthMenu = c(5, 10, 15,20,50), pageLength = 50))
     
     V(igr)$name=V(igr)$nodename
     #convert
     # Convert to object suitable for networkD3
     ig_d3 <- igraph_to_networkD3(igr, group = V(igr)$kind)
     print("conversion successful")
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
   
   output$igraphSP<-renderPlot({
     #get module
     prefix<-paste("~/output/build/",buildfun(),"/",sep="")
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
     query.base<-paste("MATCH (c:CELL)-[r0]-(x1)-[r1]-(n:wgcna {square:'",square,"',edge:",edgefun(),",name:'",unlist(strsplit(me1fun(),"ME"))[2],"'})-[r2]-(x2) WHERE (x1:cellEx OR x1:cellprop) AND (x2:reactomePW OR x2:ImmunePW OR x2:pheno OR x2:PalWangPW)",sep="")
     #print(query.base)
     nodelist<-c("c","x1","n","x2")
     edgetrips<-list(c("c","r0","x1"),c("x1","r1","n"),c("n","r2","x2"))
     #igraphplotter with igr output
     igr<-igraph_plotter(query.base,nodelist,edgetrips,rimpar="diffME",plot=TRUE,csv=TRUE,prefix=cytodir,filename = "igraphWGCNAannot",vertexsize=vertexsizefun(),legendcex=legendcexfun(),plotd3=FALSE,vertex.label.cex=vlc(),lay_out=eval(parse(text=layoutfun())))
     if(input$returnpdf==TRUE){
       pdf("plot.pdf",width=14,height=12,onefile=FALSE)#,width=12,height=7,
       igr<-igraph_plotter(query.base,nodelist,edgetrips,rimpar="diffME",plot=TRUE,csv=TRUE,prefix=cytodir,filename = "igraphWGCNAannot",vertexsize=vertexsizefun(),legendcex=legendcexfun(),plotd3=FALSE,vertex.label.cex=vlc(),lay_out=eval(parse(text=layoutfun())))
       dev.off()
     }
   })
   
   
   output$projections<-renderPlot({
     squarelist<-unlist(strsplit(setfun(),"_"))
     if(length(squarelist)==3){
       square<-squarelist[3]
     }else{square<-paste(squarelist[c(3,4)],collapse="_")}
     if(nodefun()%in%c("reactomePW","PalWangPW","ImmunePW","cellEx")){
       query.base=paste("MATCH (a:",nodefun(),")-[r]-(b:",nodefun2(),")",sep="")
     }else{query.base=paste("MATCH (a:",nodefun()," {square:'",square,"',edge:",edgefun(),"})-[r]-(b:",nodefun2(),")",sep="")}
     
     print(query.base)
     nodelist<-c("a","b")
     edgetrips<-list(c("a","r","b"))
     rimpar="diffME"
     nodeproj=nodefun()
     edgeproj=edgefun()
     
     
     
     igraph_plotter(query.base,nodelist,edgetrips,rimpar,csv=T,prefix=cytodir,filename=paste("modnet",square,edgeproj,nodeproj,collapse="_"),vertexsize=vertexsizefun(),legendcex=legendcexfun(),optlabel="\n",optvalue="V(ig)$edge",optchar="V(ig)$square",lay_out=eval(parse(text=layoutfun())),plot_bipartite = TRUE,plotd3=FALSE,plotwhich=plotwhichfun(),vertex.label.cex=vlc())
     #from igraphMM: igraph_plotter(query.base,nodelist,edgetrips,rimpar="diffEX",plot=TRUE,csv=TRUE,prefix=cytodir,filename = "igraphModuleMeta",vertexsize=vertexsizefun(),legendcex=legendcexfun(),plotd3=FALSE,vertex.label.cex=vlc(),lay_out=eval(parse(text=layoutfun())))
     
     
   })
   
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
       query1=paste("MATCH (n:wgcna {square:'",square,"',edge:5})-[r0]-(p:PROBE)-[r1]-(s:SYMBOL)-[r2]-(x:",node2,") ",sep="")
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
   
   
   output$MEvar<-renderPlot({
     print("begin MEvar")
     squarelist<-unlist(strsplit(setfun(),"_"))
     if(length(squarelist)==3){
       square<-squarelist[3]
     }else{square<-paste(squarelist[c(3,4)],collapse="_")}
     #for each ME
     
     #get ME vector
     prefix<-paste("~/output/build/",buildfun(),"/",sep="")
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
   
   cytodirlist<-system(paste("ls",cytodir),intern=TRUE)
   print("igraphwgcna_nodes.csv%in%cytodirlist")
   print("igraphwgcna_nodes.csv"%in%cytodirlist)
   if("igraphwgcna_nodes.csv"%in%cytodirlist){
     currentnodes<-read.csv(file.path(cytodir,"igraphwgcna_nodes.csv"))
     print(currentnodes)
     output$moduleNodesSP<-DT::renderDataTable(currentnodes,server = FALSE,options=list(lengthMenu = c(5, 10, 15,20,50), pageLength = 15)) 
     currentedges<-read.csv(file.path(cytodir,"igraphwgcna_edges.csv"))
     print(currentedges)
     output$moduleEdgesSP<-DT::renderDataTable(currentedges,server = FALSE,options=list(lengthMenu = c(5, 10, 15,20,50), pageLength = 15)) 
   }
   
   #modcor
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
   
   #cellcor
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
     
     cellCore3(study=studyCC,squareC=squareCC,edgeC=as.character(reacedgeCC),cellC=cellCfun(),pwm="other",PalWang=FALSE)
   })
   
   #G1
   plotInputGigamat1 <- reactive({
     cellcordir<-file.path("~/output/build",buildfun(),"12_CELLCOR/figures")
     print(cellcordir)
     load(file.path(cellcordir,"csl.RData"))
     #"~/output/build//version_2017-03-25_20_17_30//12_CELLCOR/figures/csl.RData"
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
     
    
     
     
     #new
     pheatmap(plotMatrix,col=blueWhiteRed(60),breaks=seq(floor(-breaksvals),ceiling(breaksvals),length.out=61),fontsize_col=9+cexval1,fontsize_row=9+cexval2,main=paste("Study:",names(studies)[reacstudy],"Set:",setlistnameslist[[reacstudy]][reacset],"Edge:",reacedge,"PW filter:",sprintf("%.3f",slider1()),"Cell filter:",sprintf("%.3f",slider2())))
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
     cellcordir<-file.path("~/output/build/",buildfun(),"/12_CELLCOR/figures")
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
   
   observe({output$gigamat<-renderPlot({
     plotInputGigamat1()
     
   },height=as.numeric(input$plotheightG1))})
   
   
   output$gigabar<-renderPlot({
     plotInputGigbar()
     
   })
   
   #BXP
   observe({output$boxplot<-renderPlot({
     ##CODE
     squarelist<-unlist(strsplit(setfunBXP(),"_"))
     if(length(squarelist)==3){
       square<-squarelist[3]
     }else{square<-paste(squarelist[c(3,4)],collapse="_")}
     regex<-useregex()
     print(paste("regex:",regex))
     print("command")
     print(paste("probe_boxplot4('",square,"','",regex,"',miny=",as.numeric(yminfun()),",maxy=",as.numeric(ymaxfun()),")",sep=""))
     eval(parse(text=paste("probe_boxplot4('",square,"','",regex,"',miny=",as.numeric(yminfun()),",maxy=",as.numeric(ymaxfun()),",orderP=",orderPfun(),")",sep="")))
     if(input$returnpdf==TRUE){
       pdf("plot.pdf",width=16,height=as.numeric(input$plotheight)/100,onefile=FALSE)
       eval(parse(text=paste("probe_boxplot4('",square,"','",regex,"',miny=",as.numeric(yminfun()),",maxy=",as.numeric(ymaxfun()),",orderP=",orderPfun(),")",sep="")))
       dev.off()
     }
     dequery<-paste("MATCH (p:PROBE {square:'",square,"'})-[r]-(s:SYMBOL) WHERE s.name=~ '",regex,"' RETURN p.square as square,p.edge as edge, p.name as probe, s.name as gene, p.aveEXPR as expression, p.logfc as logfc, p.adjPVAL as qvalue",sep="")
     detable<-cypher(graph,dequery)
     print(detable)
     print(returncsvfun2())
     output$detableBXP<-DT::renderDataTable(detable,server = FALSE,options=list(lengthMenu = c(5, 10, 15,20,50), pageLength = 15)) 
     observe({if(returncsvfun2()){
       print("writing table")
       write.csv(detable,file.path(tabledir,paste("boxplotDEtable",Sys.time(),".csv",sep="")))}#
     })#end observe
     }, height = as.numeric(input$plotheight))
   
   })
   
   #MM
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
     breaks<-seq(-abs(max(allmat2)),abs(max(allmat2)),length.out=51)
     modulelist<-colnames(allmat3)
     #cex calculation
     cexval<-2*sqrt(258/ncol(allmat3))
     
     pheatout <- pheatmap(allmat3,col=blueWhiteRed(50),scale="none",breaks=breaks,fontsize_col=6+cexval,fontsize_row=9+cexval,main=paste(input$subset,c("all_same","alldiff")[c(useallsame(),usediff())],"n:",ncol(allmat3)))
     pheatres <- allmat3[c(pheatout$tree_row[["order"]]),pheatout$tree_col[["order"]]]
     print(pheatres)
     #the plotting
     pheatmap(allmat3,col=blueWhiteRed(50),scale="none",breaks=breaks,fontsize_col=6+cexval,fontsize_row=9+cexval,main=paste(input$subset,c("all_same","alldiff")[c(useallsame(),usediff())],"n:",ncol(allmat3)))
     if(input$returnpdf==TRUE){
       pdf("plot.pdf",width=16,height=7,onefile=FALSE)
       pheatmap(allmat3,col=blueWhiteRed(50),scale="none",breaks=breaks,fontsize_col=6+cexval,fontsize_row=9+cexval,main=paste(input$subset,c("all_same","alldiff")[c(useallsame(),usediff())],"n:",ncol(allmat3)),width=16,height=7)
       dev.off()
     }
     
     mlv<-paste(modulelist,collapse="|")
     tablequery<-paste("MATCH (n:baylor)-[r]-(x) WHERE (x:cellEx OR x:reactomePW OR x:ImmunePW OR x:PalWangPW) AND n.name =~ '",mlv,"' RETURN n.name AS modname, r.qvalue as qval, labels(x) AS kind, x.name as name",sep="")
     tableres<-unique(cypher(graph,tablequery))
     tableres<-tableres[order(tableres$modname),]
     print(returncsvfun())
     if(returncsvfun()){
       print("writing table")
       write.csv(tableres,file.path(tabledir,paste("CMA",paste(rownames(allmat3),collapse="_"),".csv",sep="")))}#Chaussabel Module Annotation
     return(modulelist) 
   })
   
   output$pheat <- renderPlot({ 
     plotInput()
   })
   
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
     nmods<-length(allmat2[1,chosenSubset])
     den<-as.dendrogram(hclust(dist(t(allmat2[,chosenSubset]))))
     den.ord<-order.dendrogram(den)
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
     breaks<-seq(-abs(max(allmat2)),abs(max(allmat2)),length.out=51)
     modulelist<-colnames(allmat3)
     moduleinput<-as.character(paste(modulelist,collapse="|"))
     print(moduleinput)
     query.base<-paste("MATCH (b:baylor {square:'",setlist[1],"' ,edge:5})-[r]-(p) WHERE (p:reactomePW OR p:ImmunePW OR p:PalWangPW OR p:cellEx) AND b.name =~ '(?i)",moduleinput,"'",sep="")
     print(query.base)
     nodelist<-c("b","p")
     edgetrips<-list(c("b","r","p"))
     igr<-igraph_plotter(query.base,nodelist,edgetrips,rimpar="diffEX",plot=TRUE,csv=TRUE,prefix=cytodir,filename = "igraphModuleMeta",vertexsize=vertexsizefun1(),legendcex=legendcexfun1(),plotd3=FALSE,vertex.label.cex=vlc1(),lay_out=eval(parse(text=layoutfun1())))
   })
   
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
     breaks<-seq(-abs(max(allmat2)),abs(max(allmat2)),length.out=51)
     modulelist<-colnames(allmat3)
     moduleinput<-as.character(paste(modulelist,collapse="|"))
     print(moduleinput)
     query.base<-paste("MATCH (b:baylor {square:'",setlist[1],"' ,edge:5})-[r]-(p) WHERE (p:reactomePW OR p:ImmunePW OR p:PalWangPW OR p:cellEx) AND b.name =~ '(?i)",moduleinput,"'",sep="")
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
     output$moduleNodesMM<-DT::renderDataTable(currentnodes,server = FALSE,options=list(lengthMenu = c(5, 10, 15,20,50), pageLength = 50))
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
     
     
     output$moduleEdgesMM<-DT::renderDataTable(currentedges,server = FALSE,options=list(lengthMenu = c(5, 10, 15,20,50), pageLength = 50))
     
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
   
   #G2
   plotInputGigamat2 <- reactive({
     print("debugG2_1")
     cellcordir<-file.path("~/output/build/",buildfun(),"/12_CELLCOR/figures")
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
     #this needs work
     if(project=="IMPI"){
       studyvec[which(usematrix[,1]<4)]<-1
       studyvec[which(usematrix[,1]>3)]<-2 
       
     }else if (project=="SLE"){
       studyvec[which(usematrix[,1]==1)]<-1
       studyvec[which(usematrix[,1]==4)]<-2
       studyvec[which(usematrix[,1]%in%c(2,3,5,6))]<-3
       
     }else if (project=="GBP"){
       studyvec[which(usematrix[,1]%in%c(1))]<-1#berry
       studyvec[which(usematrix[,1]%in%c(2,3,4))]<-3#eu
       studyvec[which(usematrix[,1]%in%c(5))]<-4#he
       studyvec[which(usematrix[,1]%in%c(6,7,8,9,10,11,12))]<-2#impi
       
       
     }else if (project=="method"){
       studyvec[which(usematrix[,1]%in%c(1,2))]<-1
       studyvec[which(usematrix[,1]%in%c(3,4,5))]<-2
       studyvec[which(usematrix[,1]%in%c(6,7,8))]<-3
       
     }
     
     print("studyvec")
     print(studyvec)
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
     pheatmap(plotMatrix,col=blueWhiteRed(60),breaks=seq(floor(-breaksvals),ceiling(breaksvals),length.out=61),fontsize_col=9+cexval1,fontsize_row=9+cexval2)
     #pheatmap(datamatrix,col=blueWhiteRed(60),breaks=seq(floor(-breaksvals),ceiling(breaksvals),length.out=61))
     if(input$returnpdf==TRUE){
       pdf("plot.pdf",width=16,height=as.numeric(input$plotheightG2)/100,onefile=FALSE)
       pheatmap(plotMatrix,col=blueWhiteRed(60),breaks=seq(floor(-breaksvals),ceiling(breaksvals),length.out=61),fontsize_col=9+cexval1,fontsize_row=9+cexval2)
       #pheatmap(plotMatrix,col=blueWhiteRed(60),breaks=seq(floor(-breaksvals),ceiling(breaksvals),length.out=61),fontsize_col=9+cexval1,fontsize_row=9+cexval2,main=paste("study"))
       
       
       dev.off()
     }
   })
   
  observe({output$gigamat2 <- renderPlot({
     plotInputGigamat2()
   },height=as.numeric(input$plotheightG2))
  })
  
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
      query=paste("MATCH (c:cellprop) WHERE c.square='",sq,"' AND c.edge=",ed," RETURN c.name AS name, c.ratio AS ratio",sep="")
      res<-cypher(graph,query)
      print(res)
      dfr[row,res$name]<-res$ratio
      rownamesdfr<-c(rownamesdfr,paste(sq,ed,sep="_"))
    }
    rownames(dfr)<-rownamesdfr
    dfr[is.na(dfr)]<-1
    dfr[dfr==0]<-0.000001
    dfr<-log2(dfr)
    print(dfr)
    dfr<-dfr[,colSums(dfr)!=0]
    print(dfr)
    dfr[dfr>3]<-3
    dfr[dfr<(-3)]<-(-3)
    breaksvals<-max(abs(dfr))
    print(dfr)
    #breaksvals<-
    pheatmap(dfr,col=blueWhiteRed(60),breaks=seq(floor(-breaksvals),ceiling(breaksvals),length.out=61))
    
  })
  
  plotInputMachines<-reactive({
    
    ##BEGIN INSERT
    print("begin machines")
    query=paste("MATCH (s:SYMBOL)-[r]-(p:PROBE)-[r2]-(n:wgcna) WHERE p.square =~ '",labelsvec,"' AND p.edge IN [1,2,3,4] AND s.name =~ '",generegexlist[grxfun()],"' RETURN s.name AS gene, p.name as PROBE, p.square AS square, p.edge as edge, p.logfc as logfc,n.name AS wgcna",sep="")
    res<-cypher(graph,query)
    idmat<-unique(res[,c(1,2)])
    id<-as.character(apply(idmat,1,function(x){paste(x[1]," (",x[2],")",sep="")}))
    squaremat<-matrix(ncol=4,nrow=nrow(idmat),data=0)
    squaremat2<-cbind(idmat,squaremat)
    thisquare<-idmat
    print("debug-2")
    #for (square in squares){
    
    squarelist<-unlist(strsplit(setfunMACHINE(),"_"))
    print("debug-1")
    if(length(squarelist)==3){
      square<-squarelist[3]
    }else{square<-paste(squarelist[c(3,4)],collapse="_")}
    
    
    print(square)
    print("debug0")
    subset<-res[which(res$square==square),]
    thismat<-squaremat2
    for (row in 1:nrow(subset)){
      thismat[which(thismat$PROBE==subset[row,"PROBE"]),2+subset[row,"edge"]]<-subset[row,"logfc"]
    }
    print("debug1")
    #individual squares
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
    anndf2<-data.frame(anndf[order(anndf$protFunction),])
    print("debug13")
    print(anndf2)
    plotmat0<-plotmat0[order(anndf$protFunction),]
    print("debug14")
    pheatmap(plotmat0,cluster_rows = FALSE,cluster_cols = FALSE,scale = "none",col=blueWhiteRed(50),breaks=seq(-max(abs(plotmat0)),max(abs(plotmat0)),length.out=51),annotation_row = anndf2,main=paste(square,grxfun(),sep="_"))
    #pheatmap(plotmat0,cluster_rows = FALSE,cluster_cols = FALSE,scale = "none",col=blueWhiteRed(50),breaks=seq(-max(abs(plotmat)),max(abs(plotmat)),length.out=51),fontsize_row = 8,annotation_row = anndf2,main=paste(square,names(generegexlist)[grx],sep="_"),filename=file.path(figdir,paste(square,"_",names(generegexlist)[grx],".pdf",sep="")))
    
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
   
  # observe({
  #   output$inflammasomes<-renderPlot({
  #   plotInputMachines()
  # },height=as.numeric(input$plotheightMACHINE))
  # 
  #   })
  
  
    output$inflammasomes<-renderPlot({
    plotInputMachines()
  },height=700)

    

   output$downloadPlot <- downloadHandler(
     filename = function() { paste("ANIMAplot",Sys.time(),".pdf",sep="") },
     content = function(file) {
       file.copy("plot.pdf",file)
     }
   )  
   
})

# Run the application 
shinyApp(ui = ui, server = server)

