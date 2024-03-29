source('setUp.R')

#---OUT Dir
OUT    <- vector(length=3)
OUT[1] <- DIRS[grepl("GeneSets",DIRS)]
OUT[2] <- DIRS[grepl("Clustering",DIRS)]
OUT[3] <- DIRS[grepl("Graphs",DIRS)]

#---Check or create output dir
#eldir <- sprintf("%s/%s",dataDIR,subDIR[S])
#if( !file_test("-d",eldir) ){
#    dir.create(eldir)
#}

cldir <- sprintf("%s/%s",OUT[2],subDIR[S])
if( !file_test("-d",cldir) ){
    dir.create(cldir)
}

grdir <- sprintf("%s/%s",OUT[3],subDIR[S])
if( !file_test("-d",grdir) ){
    dir.create(grdir)
}
#---


#---READ IN GRAPH 
gg  <- igraph::read.graph(sprintf("%s/%s.gml",grdir,subDIR[S]),format="gml")
ids <- V(gg)$name;

#---
#  Run Spectral Clustering
#  Load the CDMSuite_0.1.0.tar.gz file on the command line using:
#> R CMD INSTALL CDMSuite_0.1.0.tar.gz
#---
library(CDMSuite)
Cnmin <- 1 #can set the minimum community size
el    <- as.data.frame(get.edgelist(gg,names=T))
spec  <- CDMSuite::spectral(DF=el, CnMIN=Cnmin)
mem   <- spec$K[match(V(gg)$name,spec$ID)]
gg    <- igraph::set.vertex.attribute(gg,"Spectral",V(gg), mem)
#---


#---
# Run lec Clustering
# No parameters to tune here
ptm <- proc.time()
lec     <- igraph::leading.eigenvector.community(gg) 
ll      <- igraph::leading.eigenvector.community(gg, start=membership(lec))
pet <- proc.time() - ptm
cat("lec \n")
cat(sprintf("time = %.3f \n", pet[[1]]))
cat(sprintf("mod=%.3f \n",max(lec$modularity)))
cat("---\n")
cc      <- matrix(NA, ncol=2, nrow=length(ids))
cc[,1]  <- as.character(ll$names)
cc[,2]  <- as.character(ll$membership)
cc      <- as.data.frame(cc)
outfile <- file(sprintf("%s/lec_communities.csv",cldir),"w")
cat("#communities",file=outfile,"\n")
write.table( cc, file=outfile, append=T, row.names=F, col.names=F, sep="\t", quote=F);
close(outfile);
#---
igraph::set.vertex.attribute(gg,"lec",V(gg),NA)
for( i in 1:length(ids) ){
  ind1 = which(cc[,1]==ids[i])
  Str <- "";
  if( length(ind1) != 0 ){if( Str == "" ){ Str <- as.character(cc[ind1[1],2]) }}
  V(gg)[i]$lec = as.integer(Str); 
}
#---

#---
# Run wt Clustering
#Parameters:
#Steps=4, optimal value as discussed in paper.
ptm <- proc.time()
wt  <- igraph::walktrap.community(gg)
pet <- proc.time() - ptm
cat("wt \n")
cat(sprintf("time = %.3f \n", pet[[1]]))
cat(sprintf("mod=%.3f \n",max(wt$modularity)))
cat("---\n")
cc      <- matrix(NA, ncol=2, nrow=length(ids))
cc[,1]  <- as.character(wt$names)
cc[,2]  <- as.character(wt$membership)
cc      <- as.data.frame(cc)
outfile <- file(sprintf("%s/wt_communities.csv",cldir),"w")
cat("#communities",file=outfile,"\n")
write.table( cc, file=outfile, append=T, row.names=F, col.names=F, sep="\t", quote=F);
close(outfile);
#---
igraph::set.vertex.attribute(gg,"wt",V(gg),NA)
for( i in 1:length(ids) ){
  ind1 = which(cc[,1]==ids[i])
  Str <- "";
  if( length(ind1) != 0 ){if( Str == "" ){ Str <- as.character(cc[ind1[1],2]) }}
  V(gg)[i]$wt = as.integer(Str); 
}
#---

#---
# Run fc Clustering
#No parameters needed to be tuned.
ptm <- proc.time()
fc  <- igraph::fastgreedy.community(gg)
pet <- proc.time() - ptm
cat("fc \n")
cat(sprintf("time = %.3f \n", pet[[1]]))
cat(sprintf("mod=%.3f \n",max(fc$modularity)))
cat("---\n")
cc      <- matrix(NA, ncol=2, nrow=length(ids))
cc[,1]  <- as.character(fc$names)
cc[,2]  <- as.character(fc$membership)
cc      <- as.data.frame(cc)
outfile <- file(sprintf("%s/fc_communities.csv",cldir),"w")
cat("#communities",file=outfile,"\n")
write.table( cc, file=outfile, append=T, row.names=F, col.names=F, sep="\t", quote=F);
close(outfile);
#---
igraph::set.vertex.attribute(gg,"fc",V(gg),NA)
for( i in 1:length(ids) ){
  ind1 = which(cc[,1]==ids[i])
  Str <- "";
  if( length(ind1) != 0 ){if( Str == "" ){ Str <- as.character(cc[ind1[1],2]) }}
  V(gg)[i]$fc = as.integer(Str); 
}
#---

#---
# Run louvain Clustering
ptm <- proc.time()
louvain <- igraph::cluster_louvain(gg)
pet <- proc.time() - ptm
cat("louvain \n")
cat(sprintf("time = %.3f \n", pet[[1]]))
cat(sprintf("mod=%.3f \n",max(louvain$modularity)))
cat("---\n")
cc       <- matrix(NA, ncol=2, nrow=length(ids))
cc[,1]   <- as.character(louvain$names)
cc[,2]   <- as.character(louvain$membership)
cc       <- as.data.frame(cc)
outfile  <- file(sprintf("%s/louvain_communities.csv",cldir),"w")
cat("#communities",file=outfile,"\n")
write.table( cc, file=outfile, append=T, row.names=F, col.names=F, sep="\t", quote=F);
close(outfile);
#---
#Save in graph 
igraph::set.vertex.attribute(gg,"louvain",V(gg),NA)
for( i in 1:length(ids) ){
  ind1 = which(cc[,1]==ids[i])
  Str <- "";
  if( length(ind1) != 0 ){if( Str == "" ){ Str <- as.character(cc[ind1[1],2]) }}
  V(gg)[i]$louvain = as.integer(Str); 
}
#---


#---
# Run infomap Clustering
ptm <- proc.time()
infomap <- igraph::cluster_infomap(gg)
pet <- proc.time() - ptm
cat("infomap \n")
cat(sprintf("time = %.3f \n", pet[[1]]))
cat("---\n")
cc       <- matrix(NA, ncol=2, nrow=length(ids))
cc[,1]   <- as.character(infomap$names)
cc[,2]   <- as.character(infomap$membership)
cc       <- as.data.frame(cc)
outfile  <- file(sprintf("%s/infomap_communities.csv",cldir),"w")
cat("#communities",file=outfile,"\n")
write.table( cc, file=outfile, append=T, row.names=F, col.names=F, sep="\t", quote=F);
close(outfile);
#---
#Save in graph 
igraph::set.vertex.attribute(gg,"infomap",V(gg),NA)
for( i in 1:length(ids) ){
  ind1 = which(cc[,1]==ids[i])
  Str <- "";
  if( length(ind1) != 0 ){if( Str == "" ){ Str <- as.character(cc[ind1[1],2]) }}
  V(gg)[i]$infomap = as.integer(Str); 
}
#---

#---
# Run sgG1 Clustering
#Parameters
#spins = 500, as used in paper, as big as the number of communities.
#gamma, optimal value is 1.0; this is when maximum modularity is reached. Bigger gamma value leads to more smaller communities. the stability of community structures varies considerably under the change of this parameter. Increasing gamma, tend to favour more edges between communities, i.e. lower the modularity value, rather than edges in communities.  

#gamma = 1.0 & spin=500, gives 25 communites, mod=0.38
#gamma = 2.0 & spin=500, gives 34 communites, mod=0.32
#gamma = 5.0 & spin=500, gives 85 communites, mod=0.25
    
ptm <- proc.time()
sg  <- igraph::spinglass.community(gg, spins=as.numeric(500),gamma=1)
pet <- proc.time() - ptm
cat("sgG1 \n")
cat(sprintf("time = %.3f \n", pet[[1]]))
cat(sprintf("mod=%.3f \n",sg$modularity))
cat("---\n")
cc      <- matrix(NA, ncol=2, nrow=length(sg$names))
cc[,1]  <- as.character(sg$names)
cc[,2]  <- as.character(sg$membership)
cc      <- as.data.frame(cc)
outfile <- file(sprintf("%s/sgG1_communities.csv",cldir),"w")
cat("#communities",file=outfile,"\n")
write.table( cc, file=outfile, append=T, row.names=F, col.names=F, sep="\t", quote=F);
close(outfile);
#---
igraph::set.vertex.attribute(gg,"sgG1",V(gg),0)
for( i in 1:length(ids) ){
  ind1 = which(cc[,1]==ids[i])
  Str <- "";
  if( length(ind1) != 0 ){if( Str == "" ){ Str <- as.character(cc[ind1[1],2]) }}
  V(gg)[i]$sgG1 = as.integer(Str); 
}
#---

##---Write .gml graph to file
igraph::write.graph(gg, sprintf("%s/%s.gml",grdir,subDIR[S]), "gml")
##---Write .graphml graph to file
igraph::write.graph(gg, sprintf("%s/%s.graphml",grdir,subDIR[S]), "graphml")

run=TRUE
if(run){

#---
# Run sgG2 Clustering    
ptm <- proc.time()
sg  <- igraph::spinglass.community(gg, spins=as.numeric(500),gamma=2)
pet <- proc.time() - ptm
cat("sgG2 \n")
cat(sprintf("time = %.3f \n", pet[[1]]))
cat(sprintf("mod=%.3f \n",sg$modularity))
cat("---\n")
cc      <- matrix(NA, ncol=2, nrow=length(sg$names))
cc[,1]  <- as.character(sg$names)
cc[,2]  <- as.character(sg$membership)
cc      <- as.data.frame(cc)
outfile <- file(sprintf("%s/sgG2_communities.csv",cldir),"w")
cat("#communities",file=outfile,"\n")
write.table( cc, file=outfile, append=T, row.names=F, col.names=F, sep="\t", quote=F);
close(outfile);
#---
igraph::set.vertex.attribute(gg,"sgG2",V(gg),0)
for( i in 1:length(ids) ){
  ind1 = which(cc[,1]==ids[i])
  Str <- "";
  if( length(ind1) != 0 ){if( Str == "" ){ Str <- as.character(cc[ind1[1],2]) }}
  V(gg)[i]$sgG2 = as.integer(Str); 
}
#---

#---
# Run sgG5 Clustering 
ptm <- proc.time()
sg  <- igraph::spinglass.community(gg, spins=as.numeric(500),gamma=5)
pet <- proc.time() - ptm
cat("sgG5 \n")
cat(sprintf("time = %.3f \n", pet[[1]]))
cat(sprintf("mod=%.3f \n",sg$modularity))
cat("---\n")
cc      <- matrix(NA, ncol=2, nrow=length(sg$names))
cc[,1]  <- as.character(sg$names)
cc[,2]  <- as.character(sg$membership)
cc      <- as.data.frame(cc)
outfile <- file(sprintf("%s/sgG5_communities.csv",cldir),"w")
cat("#communities",file=outfile,"\n")
write.table( cc, file=outfile, append=T, row.names=F, col.names=F, sep="\t", quote=F);
close(outfile);
#---
igraph::set.vertex.attribute(gg,"sgG5",V(gg),0)
for( i in 1:length(ids) ){
  ind1 = which(cc[,1]==ids[i])
  Str <- "";
  if( length(ind1) != 0 ){if( Str == "" ){ Str <- as.character(cc[ind1[1],2]) }}
  V(gg)[i]$sgG5 = as.integer(Str); 
}
#---

}#run

##---Write .gml graph to file
igraph::write.graph(gg, sprintf("%s/%s.gml",grdir,subDIR[S]), "gml")
##---Write .graphml graph to file
igraph::write.graph(gg, sprintf("%s/%s.graphml",grdir,subDIR[S]), "graphml")


