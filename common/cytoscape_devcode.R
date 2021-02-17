source("/home/rstudio/scripts/ANIMA_setup.R")
source(file.path("/home/rstudio/source_data/questions.R"),local=TRUE)
source("/home/rstudio/source_data/setlist.R",local=TRUE)
library(RCy3)
cytoscapePing()
CytoscapeConnection(host="192.168.65.2",rpcPort=9000)
cytoscapePing(base.url = "loghost:1234")
cytoscapePing(base.url = "host.docker.internal:1234")

graphstring
