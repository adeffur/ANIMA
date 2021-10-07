#install.packages("table1")
library(table1)
#library(dplyr)
# 0. Initialize
options(scipen=999)
#install.packages("kableExtra")
library(kableExtra)

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

setwd("~/source_data/validation/master")
#get clinical data as df
pheno<-read.csv("SET2_validation.csv")
#make table 1 for select parameters
t1cols<-c(1,3,4,5,6,12,103,104,109,110,112,113,119,133,134,136,149,155)
phenot1<-pheno[,t1cols]
colnames(phenot1)
phenot1[, c(4,5,9,10,11,12,13,15)] <- sapply(phenot1[, c(4,5,9,10,11,12,13,15)], as.numeric)

phenot1[, -c(1,4,5,9,10,11,12,13,15)] <- sapply(phenot1[, -c(1,4,5,9,10,11,12,13,15)], as.factor)

table1(~ ETHNICITY + AGE_CALC + SEX + weight_kg + cxr_cardiomegaly + cxr_infiltrates + bl_CD4_cells_mcL + bl_hb_g_dl + bl_hb_g_dl + bl_wbc_10.9cells_L + bl_plt_10.9_L +pc_volume + pcf_lymPred  + pcf_ADA_U_L  + pcf_TBCUL  + sputum_TBCUL  + plf_TBCUL| rx_HIVstatus , data=phenot1, overall=F, extra.col=list(`P-value`=pvalue),topclass="Rtable1-zebra")

#make additional tables for HIV, TB, heart, pericardium, etc

#make stratificationdf
stratsdf<-phenot1[,c(1,6,16)]

#get neutro data
neutro<-read.csv("neutro_prot.csv")
neutro2<-merge(stratsdf,neutro,by="person_id",all.x=FALSE,all.y=TRUE)

#plot neutro data by HIV

Categories<-interaction(neutro2$rx_HIVstatus,neutro2$compartment)
neutro2$Categories<-Categories

data<-neutro2[,-c(1,2,3,4,5)]

df_long <- pivot_longer(data=data,cols=1:12,names_to="protein",values_to="value")

# 6. determine if categorical or numeric, if both split into 2 datasets

# 7. generate ggplot on grid; user may specify grid?; also a separate plot type for categorical and numeric, mixing two types on facet does not make sense

#need to identify from data:
# categorical vs numeric
# paired or unpaired
# number of classes that can be compared



complevs<-levels(data$Categories)
if (length(complevs)==4){
  my_comparisons <- list( c(complevs[1], complevs[2]),c(complevs[1], complevs[3]),c(complevs[2], complevs[4]),c(complevs[3], complevs[4]) )
}else{
  my_comparisons <- list( c(complevs[1], complevs[2]) )
}

plotfig<-ggplot(df_long, aes(x=Categories, y=value, fill=Categories)) + 
  geom_boxplot(aes(alpha=0.3 ),width=.25) +
  #geom_boxplot(width=.25) +
  scale_fill_manual(values=c("lightsteelblue", "mistyrose", "moccasin","lightgreen")) +
  #geom_dotplot(mapping=aes(x=Categories, y=value, color=Categories, fill = Categories),binaxis='y', stackdir='center',stackratio=0.7,method='dotdensity',inherit.aes = FALSE) +
  geom_jitter(mapping=aes(x=Categories, y=value),size=0.4,alpha=0.4,width=0.1) +
  
  #geom_jitter() +
  scale_color_manual(values=c("darkblue", "darkred", "orange","darkgreen")) +
  
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
  ggtitle("neutrophil proteins by category") +
  theme(plot.title = element_text(size=rel(2.5),hjust = 0.5,color="maroon")) +
  theme(strip.text.x = element_text(size = 12, colour = "blue", angle = 0)) +
  stat_compare_means(comparisons = my_comparisons,method="wilcox.test",label.y.npc="center") + # Add pairwise comparisons p-value
  #stat_compare_means(label.y = 50) +     # Add global p-value
  scale_y_continuous(expand=expand_scale(mult=c(0,0.2)))

plotfig<-plotfig + facet_wrap(vars(protein),scales="free",ncol=3)
print(plotfig)

#get NET data
#get neutro data
NETs<-read.csv("NETs.csv")
NETs2<-merge(stratsdf,NETs,by="person_id",all.x=FALSE,all.y=TRUE)

#plot NET data by HIV
Categories<-interaction(NETs2$rx_HIVstatus,NETs2$compartment)
NETs2$Categories<-Categories

data<-NETs2[,-c(1,2,3,4,5)]

df_long <- pivot_longer(data=data,cols=1:3,names_to="protein",values_to="value")

# 6. determine if categorical or numeric, if both split into 2 datasets

# 7. generate ggplot on grid; user may specify grid?; also a separate plot type for categorical and numeric, mixing two types on facet does not make sense

#need to identify from data:
# categorical vs numeric
# paired or unpaired
# number of classes that can be compared



complevs<-levels(data$Categories)
if (length(complevs)==4){
  my_comparisons <- list( c(complevs[1], complevs[2]),c(complevs[1], complevs[3]),c(complevs[2], complevs[4]),c(complevs[3], complevs[4]) )
}else{
  my_comparisons <- list( c(complevs[1], complevs[2]) )
}


plotfig<-ggplot(df_long, aes(x=Categories, y=value, fill=Categories)) + 
  geom_boxplot(aes(alpha=0.3 ),width=.25) +
  #geom_boxplot(width=.25) +
  scale_fill_manual(values=c("lightsteelblue", "mistyrose", "moccasin","lightgreen")) +
  #geom_dotplot(mapping=aes(x=Categories, y=value, color=Categories, fill = Categories),binaxis='y', stackdir='center',stackratio=0.7,method='dotdensity',inherit.aes = FALSE) +
  geom_jitter(mapping=aes(x=Categories, y=value),size=0.4,alpha=0.4,width=0.1) +
  
  #geom_jitter() +
  scale_color_manual(values=c("darkblue", "darkred", "orange","darkgreen")) +
  
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
  ggtitle("NETs by category") +
  theme(plot.title = element_text(size=rel(2.5),hjust = 0.5,color="maroon")) +
  theme(strip.text.x = element_text(size = 12, colour = "blue", angle = 0)) +
  stat_compare_means(comparisons = my_comparisons,method="wilcox.test",label.y.npc="center") + # Add pairwise comparisons p-value
  #stat_compare_means(label.y = 50) +     # Add global p-value
  scale_y_continuous(expand=expand_scale(mult=c(0,0.2)))

plotfig<-plotfig + facet_wrap(vars(protein),scales="free",ncol=3)
print(plotfig)

#get AB data
autoantibodies<-read.csv("autoantibodies.csv")
autoantibodies2<-merge(stratsdf,autoantibodies,by="person_id",all.x=FALSE,all.y=TRUE)


#plot AB data by HIV
Categories<-interaction(autoantibodies2$rx_HIVstatus,autoantibodies2$compartment)
autoantibodies2$Categories<-Categories

data<-autoantibodies2[,-c(1,2,3,4,5)]

df_long <- pivot_longer(data=data,cols=1:7,names_to="protein",values_to="value")

out<-df_long %>%
  group_by(protein,Categories) %>%
  summarise(n=n(),median=median(value),mad=mad(value))
out2<-as.data.frame(out)

out3<-knitr::kable(out2,"html")
save_kable(out3,file="autoantibody_sum.html")

hx<-as_hux(out2,scientific=FALSE)
hx<-merge_repeated_rows(hx)

quick_html(hx,file="hx_autoAB.html")
# 6. determine if categorical or numeric, if both split into 2 datasets

# 7. generate ggplot on grid; user may specify grid?; also a separate plot type for categorical and numeric, mixing two types on facet does not make sense

#need to identify from data:
# categorical vs numeric
# paired or unpaired
# number of classes that can be compared



complevs<-levels(data$Categories)
if (length(complevs)==4){
  my_comparisons <- list( c(complevs[1], complevs[2]),c(complevs[1], complevs[3]),c(complevs[2], complevs[4]),c(complevs[3], complevs[4]) )
}else{
  my_comparisons <- list( c(complevs[1], complevs[2]) )
}


plotfig<-ggplot(df_long, aes(x=Categories, y=value, fill=Categories)) + 
  geom_boxplot(aes(alpha=0.3 ),width=.25) +
  #geom_boxplot(width=.25) +
  scale_fill_manual(values=c("lightsteelblue", "mistyrose", "moccasin","lightgreen")) +
  #geom_dotplot(mapping=aes(x=Categories, y=value, color=Categories, fill = Categories),binaxis='y', stackdir='center',stackratio=0.7,method='dotdensity',inherit.aes = FALSE) +
  geom_jitter(mapping=aes(x=Categories, y=value),size=0.4,alpha=0.4,width=0.1) +
  
  #geom_jitter() +
  scale_color_manual(values=c("darkblue", "darkred", "orange","darkgreen")) +
  
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
  ggtitle("Autoantibodies by category") +
  theme(plot.title = element_text(size=rel(2.5),hjust = 0.5,color="maroon")) +
  theme(strip.text.x = element_text(size = 12, colour = "blue", angle = 0)) +
  stat_compare_means(comparisons = my_comparisons,method="wilcox.test",label.y.npc="center") + # Add pairwise comparisons p-value
  #stat_compare_means(label.y = 50) +     # Add global p-value
  scale_y_continuous(expand=expand_scale(mult=c(0,0.2)))

plotfig<-plotfig + facet_wrap(vars(protein),scales="free",ncol=3)
print(plotfig)
