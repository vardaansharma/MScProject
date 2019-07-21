source('setUp.R')

library(CDMSuite)

setCnmin <- function( ALGN, N ){

    #set Cn min for Spectral method
    Cnmin = -1

    if( ALGN == "Spectral"      ){ Cnmin =  -1;  }
    if( ALGN == "Spectral1per"  ){ Cnmin =   1;  }
    if( ALGN == "Spectral25per" ){ Cnmin = 2.5;  }
    if( ALGN == "Spectral5per"  ){ Cnmin = 5.0;  }
    if( ALGN == "Spectral05per" ){ Cnmin = 0.5;  }
    if( ALGN == "Spectral025per"){ Cnmin = 0.25; }
    if( ALGN == "Spectral01per" ){ Cnmin = 0.1;  }
    
    if( Cnmin > 0 ){
        CnMIN = floor( (Cnmin*N)/100 )
        return(CnMIN)
    } 

    return(1)

}

removeVertexTerm <- function(GG,NAME){

    if( !is.null(igraph::get.vertex.attribute(GG,NAME)) ){
        GG <- igraph::remove.vertex.attribute(GG,name=NAME)
    }

    if( !is.null(igraph::get.vertex.attribute(GG,gsub("_","",NAME))) ){    
        GG <- igraph::remove.vertex.attribute(GG,name=gsub("_","",NAME))
    }

    return(GG)
    
}

#---Find Largest CC
findLCC <- function(GG){

    dec <- igraph::decompose.graph(GG)
    d=1
    CC=length(V(dec[[1]]))
    for( i in 1:length(dec) ){
        if(length(V(dec[[i]])) > CC){
            d=i
            CC=length(V(dec[[i]]))
        }
    }   

    GG  <- igraph::decompose.graph(GG)[[d]]
    return(GG)

}

#---Get all edges internal to a community
intraEdges <- function(GG, ALG, CC, INTRA=NULL, INTER=NULL){
  
  intra = NULL #edges in the community CC
  inter = NULL #edges going out from community CC   
  
  if( !is.null(igraph::get.vertex.attribute(GG,ALG)) ){
    # cat("\n in 1\n", ALG)
    coms <- get.vertex.attribute(GG,ALG)
    # cat(CC)
    if( length(which(coms == CC)) != 0 ){
      # cat("\n in 2 \n")
      ed_cc = E(GG)[inc(coms == CC)]
      
      all_edges_m <- get.edges(GG, ed_cc) #matrix representation
      
      inter = (ed_cc[!(all_edges_m[, 1] %in% V(GG)[coms == CC] & all_edges_m[, 2] %in% V(GG)[coms == CC])])
      
      intra = (ed_cc[(all_edges_m[, 1] %in% V(GG)[coms == CC] & all_edges_m[, 2] %in% V(GG)[coms == CC])])
      
    }
  }
  
  if( INTRA==TRUE && !is.null(intra) ){
    intra_m = get.edges(GG,intra)
    intra   = cbind(V(GG)$name[intra_m[,1]],V(GG)$name[intra_m[,2]])
    return(intra)
  }
  
  if( INTER==TRUE && !is.null(inter) ){
    inter_m = get.edges(GG,inter)
    inter   = cbind(V(GG)$name[inter_m[,1]],V(GG)$name[inter_m[,2]])
    return(inter)
  }
  
  return(NULL)
  
}

reCluster <- function( GG, ALGN, CnMAX, CnMIN ){
    
    # cat("\nin recluster function",ALGN)
    if( !is.null(igraph::get.vertex.attribute(GG,ALGN)) ){
        # cat("\n in the if\n")
        #--- algorithm clustering 1
        ALG1 <- get.vertex.attribute(GG,ALGN,V(GG))
        ALG1 <- cbind(V(GG)$name, ALG1)

        Cn <- table(as.numeric(ALG1[,2]))    
        cc <- names(Cn)[Cn > CnMAX]
        
        # if(length(cc) == 0){
        #   new_cnMax = CnMAX/2
        #   # cat(new_cnMax)
        #   cc <- names(Cn)[Cn > new_cnMax]
        # }
        # cat("cc:",length(cc))

        RES <- list()
        k=1
        for( i in 1:length(cc) ){
            # cat("in the for", cc)
            edCC = intraEdges(GG, ALGN, cc[i], INTRA=TRUE)

            if( !is.null(edCC) ){
                # cat("\n in the edcc not null \n")
                ggLCC    <- graph_from_data_frame(d=edCC, directed=F)
                
                #fc
                if( ALGN == "fc" ){                    
                    res      <- igraph::fastgreedy.community(ggLCC)
                    oo       <- cbind(res$names, res$membership)         
                }

                #lec
                if( ALGN == "lec" ){
                    lec     <- igraph::leading.eigenvector.community(ggLCC)
                    res     <- igraph::leading.eigenvector.community(ggLCC, start=membership(lec))
                    oo      <- cbind(res$names, res$membership)         
                }

                #wt
                if( ALGN == "wt" ){
                    res      <- igraph::walktrap.community(ggLCC)
                    oo       <- cbind(res$names, res$membership)         
                }

                
                #louvain
                if( ALGN == "louvain" ) {
                    res      <- igraph::cluster_louvain(ggLCC)
                    oo       <- cbind(res$names, res$membership)         
                }

                #infomap
                if( ALGN == "infomap" ){
                  # cat("\n in the if infomap \n")
                    res      <- igraph::cluster_infomap(ggLCC)
                    # cat(res)
                    oo       <- cbind(res$names, res$membership)         
                }

                #sgG1
                if( ALGN == "sgG1" ){
                    res      <- igraph::spinglass.community(findLCC(ggLCC), spins=as.numeric(500),gamma=1)
                    oo       <- cbind(res$names, res$membership)         
                }
                
                if( ALGN == 'sgG2'){
                  res      <- igraph::spinglass.community(findLCC(ggLCC), spins=as.numeric(500),gamma=2)
                  oo       <- cbind(res$names, res$membership)
                }
                
                if( ALGN == 'sgG2'){
                  res      <- igraph::spinglass.community(findLCC(ggLCC), spins=as.numeric(500),gamma=5)
                  oo       <- cbind(res$names, res$membership)
                }
                
                #Spectral
                if( grepl("Spectral", ALGN) ){        
                    res <- spectral(DF=edCC, CnMIN=CnMIN)
                    oo  <- cbind(res$ID, res$K)
                }
            
                RES[[k]]      <- oo
                names(RES)[k] <- cc[i]
                k=k+1
            }
            
        }#for


        if( length(RES) == 0 ){ 
          # cat("\n in len res is 0\n")
          return(NULL) }
        
        #--- algorithm clustering 2
        ALG2     <- cbind(ALG1, rep(-1, length(ALG1[,1])))
        indx     <- match(ALG2[,2],cc)
        indx     <- ifelse(is.na(indx),TRUE, FALSE)
        ALG2[,3] <- ifelse(indx, ALG2[,2], ALG2[,3])

        CCmax = max(as.numeric(ALG2[,3]))

        for( i in 1:length(RES) ){
 
            temp     <- RES[[i]]
            temp[,2] <- as.numeric(temp[,2]) + CCmax
        
            indx <- match(ALG2[,1],temp[,1])
            indx <- temp[indx,2]
        
            ALG2[,3] = ifelse(is.na(indx),ALG2[,3],indx)

            CCmax = max(as.numeric(ALG2[,3]))
        
        }

        #---reorder ALG2[,3]
        N = length(V(GG));
    
        temp    <- rep(-1, N)
        counter <- min(as.numeric(ALG2[,3]))
        Knew    <- 1;
        Kmax    <- max(as.numeric(ALG2[,3]))

        while( counter <= Kmax ){

            found=FALSE;
        
            for(v in 1:N ){
                if( as.numeric(ALG2[v,3]) == counter ){
                    temp[v] = Knew;
                    found=TRUE;
                }
            }

            if(found) Knew=Knew+1;
      
            counter=counter+1;
        }
        

        #---final 
        ALG3 <- cbind(ALG2, temp)
        return(ALG3)
    }

    return(NULL)
    
}



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

#---igraph's clustering algs.
#alg    <- vector(length=9)
#alg[1] <- "lourvain"
#alg[2] <- "infomap"
#alg[3] <- "fc"
#alg[4] <- "lec"
#alg[5] <- "sgG1"
#alg[6] <- "sgG2"
#alg[7] <- "sgG5"
#alg[8] <- "wt"
#alg[9] <- "spectral"

#---READ IN GRAPH 
gg  <- igraph::read.graph(sprintf("%s/%s.gml",grdir,subDIR[S]),format="gml")
ids <- V(gg)$name;

#set a maximum community size, i.e. 10% of the networks size
CnMAX <- floor((10*length(V(gg)))/100)

cat("> CnMAX = ", CnMAX,"\n")

#For example, lets runn over these clustering methods
ALGN <- ALGS[c(1,2,3,4,6,7,15,16,17,22)]
# ALGN <- ALGS[c(7,17)]


for( a in 1:length(ALGN) ){

    cat("reclustering ", ALGN[a], "...\n")

    #set Cn min for Spectral method
    CnMIN = setCnmin( ALGN[a], length(V(gg)) )

    cat("> CnMIN = ", CnMIN,"\n")
    
    oo = reCluster( gg, ALGN[a], CnMAX, CnMIN )

    if( !is.null(oo) ){
        cat("\n is reclutering \n")
        #--- reclustering name
        ALGN2 <- sprintf("%s_%s",ALGN[a],"2")

        removeVertexTerm(gg,ALGN2)
    
        if( is.null(igraph::get.vertex.attribute(gg,ALGN2)) ){
            gg <- igraph::set.vertex.attribute(gg,ALGN2,V(gg),as.numeric(oo[,4]))
        }

        cc      <- matrix(NA, ncol=2, nrow=length(oo[,1]))
        cc[,1]  <- as.character(oo[,1])
        cc[,2]  <- as.character(oo[,4])
        cc      <- as.data.frame(cc)
        # cat(ALGN2)
        outfile <- file(sprintf("%s/%s_communities.csv",cldir,ALGN2),"w")
        # outfile <- file(sprintf("%s/%s_communities_cnmax%s.csv",cldir,ALGN2,CnMAX),"w")
        cat("#communities",file=outfile,"\n")
        write.table( cc, file=outfile, append=T, row.names=F, col.names=F, sep="\t", quote=F);
        close(outfile);
        #---

    }
        
        cat("...done.\n")
    
}


##---Write .gml graph to file
igraph::write.graph(gg, sprintf("%s/%s.gml",grdir,subDIR[S]), "gml")
##---Write .graphml graph to file
igraph::write.graph(gg, sprintf("%s/%s.graphml",grdir,subDIR[S]), "graphml")


