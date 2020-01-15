##Start with latest ROpenSci image
##This Dockerfile inspired by rocker/hadleyverse and rocker/ropensci Dockerfiles as we are extending ropensci

#FROM rocker/ropensci:3.4.1.
#FROM rocker/verse:3.4.3
#FROM rocker/verse:devel
FROM rocker/verse:latest
#FROM rocker/ml:latest

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
    tm \
    wordcloud \
    shinycssloaders \
    ggpubr \
    psych \
&& rm -rf /tmp/downloaded_packages/ /tmp/*.rds

#install beamer
##### TEX related ########
# N.B: texlive-base may prompt a question about default paper size for applications
# setting DEBIAN_FRONTEND=noninteractive  force a default choice
RUN DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends \
  dvipng \
  ghostscript \
  texinfo \
  texlive-base \
  texlive-extra-utils \
  texlive-fonts-recommended \
  texlive-generic-recommended \
  texlive-latex-base \
  texlive-latex-recommended

## attempt at using custom texlive
#RUN apt-get install -y --no-install-recommends lmodern xzdec
# fix tlmgr which does not has correct exit status
RUN tlmgr init-usertree --usertree /usr/local/share/texmf/; true

RUN tlmgr --usertree /usr/local/share/texmf/ install \
      acronym background bigfoot blindtext cleveref enumitem etoolbox everypage \
      floatrow framed beamer translator \
      inconsolata pgf multirow preprint soul titling xstring && \
    texhash && updmap-sys


#udunits
RUN apt-get update

RUN apt-get install -y \
    libudunits2-0 \
    libudunits2-dev

#install more from CRAN
RUN rm -rf /tmp/*.rds \
&& install2.r --error \
    -r http://cran.rstudio.com \
    udunits2 \
&& rm -rf /tmp/downloaded_packages/ /tmp/*.rds

#bioconductor packages

RUN Rscript -e 'source("http://bioconductor.org/biocLite.R");biocLite();bio<-c("Homo.sapiens","lumi","beadarray","annotate","lumiHumanAll.db","illuminaHumanv2.db","illuminaHumanv3.db","illuminaHumanv4.db","lumiHumanIDMapping","limma","WGCNA","gage","gageData","pathview","GOstats","goProfiles","ReactomePA","DOSE","clusterProfiler","RCy3","minet","csSAM","sva","pheatmap","flashClust","impute","hgu133a.db","hgu133b.db","Heatplus","marray","RCy3");biocLite(bio)' \
&& rm -rf /tmp/downloaded_packages/ /tmp/*.rds

#install more from CRAN
#RUN rm -rf /tmp/*.rds \
#&& install2.r --error \
#    -r http://cran.rstudio.com \
#    CluMix \
#&& rm -rf /tmp/downloaded_packages/ /tmp/*.rds

#RNeo4j from GitHub
RUN Rscript -e 'devtools::install_github("nicolewhite/RNeo4j")' \
&& rm -rf /tmp/downloaded_packages/ /tmp/*.rds

RUN Rscript -e 'devtools::install_github("ricardo-bion/ggradar",dependencies=TRUE)' \
&& rm -rf /tmp/downloaded_packages/ /tmp/*.rds

#OOmpa packages
RUN Rscript -e 'source("http://silicovore.com/OOMPA/oompaLite.R");oompaLite()' \
&& rm -rf /tmp/downloaded_packages/ /tmp/*.rds

#CellMix
RUN Rscript -e 'source("http://bioconductor.org/biocLite.R");biocLite("CellMix",siteRepos = "http://web.cbio.uct.ac.za/~renaud/CRAN")' \
&& rm -rf /tmp/downloaded_packages/ /tmp/*.rds
