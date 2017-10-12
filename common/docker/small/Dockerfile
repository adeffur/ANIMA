##Start with latest ROpenSci image
##This Dockerfile inspired by rocker/hadleyverse and rocker/ropensci Dockerfiles as we are extending ropensci

FROM rocker/ropensci:3.3.4

#This is me
MAINTAINER Armin Deffur <armin.deffur@absamail.co.za>

## Add RStudio binaries to PATH
ENV PATH /usr/lib/rstudio-server/bin/:$PATH 
ENV LANG en_US.UTF-8

#install CRAN packages
RUN rm -rf /tmp/*.rds \
&& install2.r --error \
    -r http://cran.rstudio.com \
    #-r http://www.omegahat.org/R \
    -r http://datacube.wu.ac.at \
    -r http://packages.ropensci.org \
    gplots \
    VennDiagram \
    pwr \
    mail \
    clusterRepro \
    doBy \
    gmodels \
    ROCR \
    aod \
    cowplot \
    heatmap.plus \
    dunn.test \
    DT \
    networkD3 \
&& rm -rf /tmp/downloaded_packages/ /tmp/*.rds

#bioconductor packages
RUN Rscript -e 'source("http://bioconductor.org/biocLite.R");biocLite();bio<-c("lumi","beadarray","annotate","lumiHumanAll.db","illuminaHumanv2.db","illuminaHumanv3.db","illuminaHumanv4.db","lumiHumanIDMapping","limma","WGCNA","gage","gageData","pathview","GOstats","goProfiles","ReactomePA","DOSE","clusterProfiler","RCytoscape","minet","csSAM","sva","pheatmap","flashClust","impute","hgu133a.db","hgu133b.db","Heatplus");biocLite(bio)' \
&& rm -rf /tmp/downloaded_packages/ /tmp/*.rds

#RNeo4j from GitHub
RUN Rscript -e 'devtools::install_github("nicolewhite/RNeo4j")' \
&& rm -rf /tmp/downloaded_packages/ /tmp/*.rds

#OOmpa packages
RUN Rscript -e 'source("http://bioinformatics.mdanderson.org/OOMPA/oompaLite.R");oompaLite()' \
&& rm -rf /tmp/downloaded_packages/ /tmp/*.rds

#CellMix
RUN Rscript -e 'source("http://bioconductor.org/biocLite.R");biocLite("CellMix", siteRepos = "http://web.cbio.uct.ac.za/~renaud/CRAN")' \
&& rm -rf /tmp/downloaded_packages/ /tmp/*.rds
