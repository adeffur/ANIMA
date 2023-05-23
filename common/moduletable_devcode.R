install.packages('huxtable')
#install.packages('textreadr')

source("~/scripts/ANIMA_setup.R")
source(file.path(dir.data_root,"questions.R")) #this is specific to the project; each project will get its own set of questions tied to specific data.
source("~/source_data/setlist.R")
library(huxtable)

##### Moduletable

graph<-startGraph(graphstring,database="animadb",username="neo4j",password="anima")

query1<-"MATCH (n:wgcna {square:'blood.PCF.defPC',edge:'5'}) RETURN n.name AS module"
modules<-cypher(graph,query1)

#modules<-modules[which(modules$module%in%c("brown","blue","turquoise")),]
#modules<-modules[which(modules$module%in%c("green","darkgreen","turquoise")),]
modules<-modules[which(modules$module%in%c("brown","blue","pink")),]

df <- data.frame(matrix(ncol = 6, nrow = 0))
colnames(df)<-c("module","sigenrich","modAUC1","modAUC2","entity","qvalue")
#for (module in modules$module){
for (module in modules){
  for (entity in c("x:cellEx","x:cellprop","x:ImmunePW OR x:PalWangPW OR x:reactomePW")){
    query<-paste("MATCH (n:wgcna {square:'blood.PCF.defPC',edge:'5',name:'",module,"'})
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

width(hx) <- 0.6
#wrap(hx) <- TRUE
hx<-set_wrap(hx, everywhere, "entity", TRUE )
#for (module in modules$module) {
for (module in modules) {  
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


quick_html(hx6,file="/home/rstudio/output/moduleAnnotations.html")
#quick_latex(hx3,file="/home/rstudio/output/mthx.tex")
#write.csv(hx2,"/home/rstudio/output/mthx.csv")

##### Chaussabel annotation table for meta-analysis ####
#effect of organism burden
data<-read.csv("~/output/project/tabular/CMAprob.def.blood_edge_1_prob.def.blood_edge_2_prob.def.fluid_edge_1_prob.def.fluid_edge_2.csv")
hx<-as_hux(data[,2:5])
#hx2<-merge_repeated_rows(hx,everywhere,c("module", "diffME", "sigenrich", "modAUC1", "modAUC2"))


hx<-merge_repeated_rows(hx,everywhere,"modname")



#hx2<-merge_repeated_rows(hx,everywhere,c("module"))
hx3<-theme_article(hx)
#hx3<-theme_plain(hx2)
#hx3<-theme_mondrian(hx2)
#hx4<-set_background_color(hx3,evens, 2:4,"grey95")

hx5<-hx3 %>% map_text_color(everywhere, c("qval"),
                            by_colorspace("red","grey", "blue"))

quick_html(hx5,file="/home/rstudio/output/Burden_hx.html")

#effect of HIV
data<-read.csv("~/output/project/tabular/CMAblood.PCF.defPC_edge_3_blood.PCF.defPC_edge_4_blood.PCF.probPC_edge_3_blood.PCF.probPC_edge_4.csv")
hx<-as_hux(data[,2:5])
#hx2<-merge_repeated_rows(hx,everywhere,c("module", "diffME", "sigenrich", "modAUC1", "modAUC2"))


hx<-merge_repeated_rows(hx,everywhere,"modname")



#hx2<-merge_repeated_rows(hx,everywhere,c("module"))
hx3<-theme_article(hx)
#hx3<-theme_plain(hx2)
#hx3<-theme_mondrian(hx2)
#hx4<-set_background_color(hx3,evens, 2:4,"grey95")

hx5<-hx3 %>% map_text_color(everywhere, c("qval"),
                            by_colorspace("red","grey", "blue"))

quick_html(hx5,file="/home/rstudio/output/HIV_hx.html")

#effect of HD phenotype
data<-read.csv("~/output/project/tabular/CMAblood.PCF.HIVposHD_edge_3_blood.PCF.HIVposHD_edge_4.csv")
hx<-as_hux(data[,2:5])
#hx2<-merge_repeated_rows(hx,everywhere,c("module", "diffME", "sigenrich", "modAUC1", "modAUC2"))


hx<-merge_repeated_rows(hx,everywhere,"modname")



#hx2<-merge_repeated_rows(hx,everywhere,c("module"))
hx3<-theme_article(hx)
#hx3<-theme_plain(hx2)
#hx3<-theme_mondrian(hx2)
#hx4<-set_background_color(hx3,evens, 2:4,"grey95")

hx5<-hx3 %>% map_text_color(everywhere, c("qval"),
                            by_colorspace("red","grey", "blue"))

quick_html(hx5,file="/home/rstudio/output/HDpheno_hx.html")

#SAME_STORY
data<-read.csv("~/output/project/tabular/CMAberry.test_edge_5_blood.PCF.defPC_edge_5_blood.PCF.probPC_edge_5_nonLTBI.PTB_edge_5_nonLTBI.TBPC_edge_5.csv")
hx<-as_hux(data[,2:5])
#hx2<-merge_repeated_rows(hx,everywhere,c("module", "diffME", "sigenrich", "modAUC1", "modAUC2"))


hx<-merge_repeated_rows(hx,everywhere,"modname")



#hx2<-merge_repeated_rows(hx,everywhere,c("module"))
hx3<-theme_article(hx)
#hx3<-theme_plain(hx2)
#hx3<-theme_mondrian(hx2)
#hx4<-set_background_color(hx3,evens, 2:4,"grey95")

hx5<-hx3 %>% map_text_color(everywhere, c("qval"),
                            by_colorspace("red","grey", "blue"))

quick_html(hx5,file="/home/rstudio/output/SAMESTORY_hx.html")

#TB-PC in context
data<-read.csv("~/output/project/tabular/consensus_blood.csv")
hx<-as_hux(data[,2:5])
#hx2<-merge_repeated_rows(hx,everywhere,c("module", "diffME", "sigenrich", "modAUC1", "modAUC2"))


hx<-merge_repeated_rows(hx,everywhere,"modname")



#hx2<-merge_repeated_rows(hx,everywhere,c("module"))
hx3<-theme_article(hx)
#hx3<-theme_plain(hx2)
#hx3<-theme_mondrian(hx2)
#hx4<-set_background_color(hx3,evens, 2:4,"grey95")

hx5<-hx3 %>% map_text_color(everywhere, c("qval"),
                            by_colorspace("red","grey", "blue"))

quick_html(hx5,file="/home/rstudio/output/consensusBlood_hx.html")

#PTB_TB-PC
data<-read.csv("~/output/project/tabular/CMAPTB.TBPC_edge_1_PTB.TBPC_edge_2.csv")
hx<-as_hux(data[,2:5])
#hx2<-merge_repeated_rows(hx,everywhere,c("module", "diffME", "sigenrich", "modAUC1", "modAUC2"))


hx<-merge_repeated_rows(hx,everywhere,"modname")



#hx2<-merge_repeated_rows(hx,everywhere,c("module"))
hx3<-theme_article(hx)
#hx3<-theme_plain(hx2)
#hx3<-theme_mondrian(hx2)
#hx4<-set_background_color(hx3,evens, 2:4,"grey95")

hx5<-hx3 %>% map_text_color(everywhere, c("qval"),
                            by_colorspace("red","grey", "blue"))

quick_html(hx5,file="/home/rstudio/output/project/tabular/PTB-TBPC_hx.html")

