##Start with latest ROpenSci image
##This Dockerfile inspired by rocker/hadleyverse and rocker/ropensci Dockerfiles as we are extending ropensci

#FROM rocker/ropensci:3.4.1.
#FROM rocker/verse:3.4.3
#FROM rocker/verse:devel
FROM thomasweise/texlive:latest
#FROM rocker/verse:latest
FROM rocker/verse:4.2.0
#FROM rocker/verse:4.0.0
#FROM rocker/ml:latest

#This is me
MAINTAINER Armin Deffur <armin.deffur@absamail.co.za>

#Environment for timeout
ENV R_DEFAULT_INTERNET_TIMEOUT 1200


# ## Get Java (for h2o R package)
# RUN apt-get update -qq \
#   && apt-get -y --no-install-recommends install \
#     cmake \
#     default-jdk \
#     default-jre \
#   && R CMD javareconf

## h2o requires Java
#RUN install2.r -r http://r.adu.org.za/ h2o future coda
#RUN install2.r h2o future coda
#RUN install2.r greta
#RUN installGithub.r greta-dev/greta

#continue with anima libraries

## Add RStudio binaries to PATH
ENV PATH /usr/lib/rstudio-server/bin/:$PATH 
ENV LANG en_US.UTF-8


#install CRAN packages
RUN rm -rf /tmp/*.rds \
&& install2.r --error \
    -r https://cran.mirror.ac.za \
    #-r http://r.adu.org.za \
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
    #heatmap.plus \
    dunn.test \
    DT \
    networkD3 \
    tm \
    wordcloud \
    shinycssloaders \
    ggpubr \
    shinythemes \
    neo2R \
    psych \
    reactlog \
    UpSetR \
    vcfR \
&& rm -rf /tmp/downloaded_packages/ /tmp/*.rds

#udunits
RUN apt-get update

RUN apt-get install -y \
    libudunits2-0 \
    libudunits2-dev

#install more from CRAN
RUN rm -rf /tmp/*.rds \
&& install2.r --error \
    -r http://r.adu.org.za \ 
    -r http://cran.rstudio.com \
    udunits2 \
&& rm -rf /tmp/downloaded_packages/ /tmp/*.rds

# Bioconductor -----------------------------------------------------------------

#FROM bioconductor/bioconductor_docker:devel
ENV R_DEFAULT_INTERNET_TIMEOUT=1200

#bioconductor packages
RUN Rscript -e 'install.packages("BiocManager")'

RUN Rscript -e 'BiocManager::install();bio<-c("ReactomePA","Homo.sapiens","lumi","beadarray","annotate","lumiHumanAll.db","illuminaHumanv2.db","illuminaHumanv3.db","illuminaHumanv4.db","lumiHumanIDMapping","limma","WGCNA","gage","gageData","pathview","GOstats","goProfiles","DOSE","clusterProfiler","RCy3","minet","csSAM","sva","pheatmap","flashClust","impute","hgu133a.db","hgu133b.db","Heatplus","marray","RCy3");BiocManager::install(bio,timeout=1200)' 

#RUN Rscript -e 'source("http://bioconductor.org/biocLite.R");biocLite();bio<-c("Homo.sapiens","lumi","beadarray","annotate","lumiHumanAll.db","illuminaHumanv2.db","illuminaHumanv3.db","illuminaHumanv4.db","lumiHumanIDMapping","limma","WGCNA","gage","gageData","pathview","GOstats","goProfiles","ReactomePA","DOSE","clusterProfiler","RCy3","minet","csSAM","sva","pheatmap","flashClust","impute","hgu133a.db","hgu133b.db","Heatplus","marray","RCy3","maftools","VariantAnnotation","TVTB");biocLite(bio)' \
#&& rm -rf /tmp/downloaded_packages /tmp/*.rds

#install more from CRAN
#RUN rm -rf /tmp/*.rds \
#&& install2.r --error \
#    -r http://cran.rstudio.com \
#    CluMix \
#&& rm -rf /tmp/downloaded_packages/ /tmp/*.rds

#From GitHub
#RUN Rscript -e 'devtools::install_github("nicolewhite/RNeo4j")' \
#&& rm -rf /tmp/downloaded_packages/ /tmp/*.rds

RUN Rscript -e 'devtools::install_github("ricardo-bion/ggradar",dependencies=TRUE)' \
&& rm -rf /tmp/downloaded_packages/ /tmp/*.rds

RUN Rscript -e 'devtools::install_github("rstudio/shinymeta",dependencies=TRUE)' \
&& rm -rf /tmp/downloaded_packages/ /tmp/*.rds


#OOmpa packages
RUN Rscript -e 'source("http://silicovore.com/OOMPA/oompaLite.R");oompaLite()' \
&& rm -rf /tmp/downloaded_packages/ /tmp/*.rds

#CellMix
#RUN Rscript -e 'BiocManager::install("CellMix",siteRepos = "http://web.cbio.uct.ac.za/~renaud/CRAN")' \
#&& rm -rf /tmp/downloaded_packages/ /tmp/*.rds

RUN Rscript -e 'install.packages("BiocInstaller", repos="http://bioconductor.org/packages/2.13/bioc")' \
&& rm -rf /tmp/downloaded_packages/ /tmp/*.rds

RUN Rscript -e 'BiocManager::install("CellMix",site_repository = "http://web.cbio.uct.ac.za/~renaud/CRAN",update=FALSE)' \
&& rm -rf /tmp/downloaded_packages/ /tmp/*.rds







