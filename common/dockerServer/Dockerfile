##Start with latest ROpenSci image
##This Dockerfile inspired by rocker/hadleyverse and rocker/ropensci Dockerfiles as we are extending ropensci

#FROM rocker/ropensci:3.4.1.
#FROM rocker/verse:3.4.3
#FROM rocker/verse:devel
FROM thomasweise/texlive:latest
#FROM rocker/verse:latest
FROM rocker/shiny-verse:3.6.3

#FROM rocker/ml:latest

#This is me
MAINTAINER Armin Deffur <armin.deffur@absamail.co.za>

#Add from rocker/tensorflow dockerfile

ENV WORKON_HOME /opt/virtualenvs
ENV PYTHON_VENV_PATH $WORKON_HOME/r-tensorflow

## Set up a user modifyable python3 environment
RUN apt-get update && apt-get install -y --no-install-recommends \
        libpython3-dev \
        python3-venv && \
    rm -rf /var/lib/apt/lists/*

RUN apt-get update && apt-get install -y --no-install-recommends apt-utils

RUN python3 -m venv ${PYTHON_VENV_PATH}


# Add rstudio user
RUN groupadd  rstudio \
&& useradd --gid rstudio --shell /bin/bash --create-home rstudio


RUN chown -R rstudio:rstudio ${WORKON_HOME}
ENV PATH ${PYTHON_VENV_PATH}/bin:${PATH}
## And set ENV for R! It doesn't read from the environment...
RUN echo "PATH=${PATH}" >> /usr/local/lib/R/etc/Renviron && \
    echo "WORKON_HOME=${WORKON_HOME}" >> /usr/local/lib/R/etc/Renviron && \
    echo "RETICULATE_PYTHON_ENV=${PYTHON_VENV_PATH}" >> /usr/local/lib/R/etc/Renviron

## Because reticulate hardwires these PATHs...
RUN ln -s ${PYTHON_VENV_PATH}/bin/pip /usr/local/bin/pip && \
    ln -s ${PYTHON_VENV_PATH}/bin/virtualenv /usr/local/bin/virtualenv

## install as user to avoid venv issues later
USER rstudio
RUN pip3 install \
    h5py==2.9.0 \
    pyyaml==3.13 \
    requests==2.21.0 \
    Pillow==5.4.1 \
    tensorflow==1.13.1 \
    tensorflow-probability==0.5.0 \
    keras==2.2.4 \
    --no-cache-dir
USER root
RUN install2.r -r http://r.adu.org.za/ reticulate tensorflow keras

#Add from rocker/ml Dockerfile

# Python Xgboost for CPU
USER rstudio
RUN pip3 --no-cache-dir install \
    xgboost==0.81 \
    wheel==0.33.0 \
    setuptools==40.8.0 \
    scipy==1.2.1
USER root

## Get Java (for h2o R package)
RUN apt-get update -qq \
  && apt-get -y --no-install-recommends install \
    cmake \
    default-jdk \
    default-jre \
  && R CMD javareconf

## h2o requires Java
RUN install2.r -r http://r.adu.org.za/ h2o future coda
#RUN install2.r greta
RUN installGithub.r greta-dev/greta

#continue with anima libraries

## Add RStudio binaries to PATH
ENV PATH /usr/lib/rstudio-server/bin/:$PATH 
ENV LANG en_US.UTF-8


#install psych package
#RUN rm -rf /tmp/*.rds \
#&& install2.r --error \
#    #-r https://cran.mirror.ac.za \
#    #-r http://r.adu.org.za \
#    #-r http://cran.us.r-project.org \
#    #-r http://cran.rstudio.com \
#    #-r http://www.omegahat.org/R \
#    #-r http://datacube.wu.ac.at \
#    #-r http://packages.ropensci.org \
#    psych \
#&& rm -rf /tmp/downloaded_packages/ /tmp/*.rds

RUN R -e "devtools::install_version('mnormt', '1.5-7')"

RUN R -e "devtools::install_version('psych', '1.9.12')"


#install CRAN packages
RUN rm -rf /tmp/*.rds \
&& install2.r --error \
    -r https://cran.mirror.ac.za \
    -r http://r.adu.org.za \
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
    shinythemes \
    neo2R \
&& rm -rf /tmp/downloaded_packages/ /tmp/*.rds

#install beamer
##### TEX related ########
# N.B: texlive-base may prompt a question about default paper size for applications
# setting DEBIAN_FRONTEND=noninteractive  force a default choice

#RUN DEBIAN_FRONTEND=noninteractive dpkg --purge --force-depends texlive-lang-other texlive-latex-extra tex-common texlive-fonts-recommended texlive-pictures texlive-metapost
#RUN DEBIAN_FRONTEND=noninteractive apt -y --fix-broken install
#RUN DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends -f texlive-lang-other texlive-latex-extra tex-common texlive-fonts-recommended texlive-pictures texlive-metapost

#RUN DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends -f texlive-base texlive-lang-other texlive-latex-extra tex-common texlive-fonts-recommended texlive-pictures texlive-metapost
#RUN DEBIAN_FRONTEND=noninteractive mktexlsr
#RUN DEBIAN_FRONTEND=noninteractive updmap-sys 

#RUN DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends \
#  texlive-base
  
#RUN DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends \  
  #dvipng \
  #ghostscript \
  #texinfo \
  #texlive-extra-utils \
  #texlive-fonts-recommended \
  #texlive-generic-recommended \
  #texlive-latex-base \
  #texlive-latex-recommended

## attempt at using custom texlive
#RUN apt-get install -y --no-install-recommends lmodern xzdec
# fix tlmgr which does not has correct exit status
#RUN tlmgr update --self

#RUN tlmgr init-usertree --usertree /usr/local/share/texmf/; true

#RUN tlmgr --usertree /usr/local/share/texmf/ install \
#      acronym background bigfoot blindtext cleveref enumitem etoolbox everypage \
#      floatrow framed beamer translator \
#      inconsolata pgf multirow preprint soul titling xstring && \
#    texhash && updmap-sys


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

#bioconductor packages

RUN Rscript -e 'install.packages("BiocManager")'

RUN Rscript -e 'BiocManager::install();bio<-c("Homo.sapiens","lumi","beadarray","annotate","lumiHumanAll.db","illuminaHumanv2.db","illuminaHumanv3.db","illuminaHumanv4.db","lumiHumanIDMapping","limma","WGCNA","gage","gageData","pathview","GOstats","goProfiles","ReactomePA","DOSE","clusterProfiler","RCy3","minet","csSAM","sva","pheatmap","flashClust","impute","hgu133a.db","hgu133b.db","Heatplus","marray","RCy3");BiocManager::install(bio)' 


#RUN Rscript -e 'source("http://bioconductor.org/biocLite.R");biocLite();bio<-c("Homo.sapiens","lumi","beadarray","annotate","lumiHumanAll.db","illuminaHumanv2.db","illuminaHumanv3.db","illuminaHumanv4.db","lumiHumanIDMapping","limma","WGCNA","gage","gageData","pathview","GOstats","goProfiles","ReactomePA","DOSE","clusterProfiler","RCy3","minet","csSAM","sva","pheatmap","flashClust","impute","hgu133a.db","hgu133b.db","Heatplus","marray","RCy3");biocLite(bio)' \
#&& rm -rf /tmp/downloaded_packages /tmp/*.rds

#install more from CRAN
#RUN rm -rf /tmp/*.rds \
#&& install2.r --error \
#    -r http://cran.rstudio.com \
#    CluMix \
#&& rm -rf /tmp/downloaded_packages/ /tmp/*.rds

#RNeo4j from GitHub
#RUN Rscript -e 'devtools::install_github("nicolewhite/RNeo4j")' \
#&& rm -rf /tmp/downloaded_packages/ /tmp/*.rds

RUN Rscript -e 'devtools::install_github("ricardo-bion/ggradar",dependencies=TRUE)' \
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



RUN apt-get update && apt-get install -y --no-install-recommends \
        libjpeg62-turbo-dev && \
    rm -rf /var/lib/apt/lists/*

RUN rm -rf /tmp/*.rds \
&& install2.r --error \
    -r https://cran.mirror.ac.za \
    -r http://r.adu.org.za \
    -r http://cran.rstudio.com \
    #-r http://www.omegahat.org/R \
    -r http://datacube.wu.ac.at \
    -r http://packages.ropensci.org \
    jpeg \
&& rm -rf /tmp/downloaded_packages/ /tmp/*.rds

RUN Rscript -e 'BiocManager::install();bio<-c("WGCNA");BiocManager::install(bio)' 

RUN apt-get update && apt-get install -y --no-install-recommends \
        libbz2-dev \
        liblzma-dev && \
    rm -rf /var/lib/apt/lists/*

RUN Rscript -e 'BiocManager::install();bio<-c("lumi");BiocManager::install(bio)' 




