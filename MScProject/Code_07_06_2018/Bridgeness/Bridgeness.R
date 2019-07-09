
source('../setUp.R');

#Calculate the Median absolute difference
MAD <- function( X, na.rm=F ){

    X <- as.numeric(X)

    if( na.rm ){
        X <- X[!is.na(X)]
    }
    
    Xmd <- median(X)

    MAD <- median(abs(X-Xmd))

    return(MAD)
}


#---Directories needed
OUT <- vector(length=3)
OUT[1] <- DIRS[grepl("Graphs",DIRS)]
OUT[2] <- DIRS[grepl("POWERlawFIT",DIRS)]
OUT[3] <- DIRS[grepl("SVI",DIRS)]

#---Check or create output dir
if( !file_test("-d",subDIR[S]) ){
    dir.create(subDIR[S])
}

#---Set Options
runBridge  <- vector(length=2)
runBridge[1] <- 0 #Calculate Bridgeness
runBridge[2] <- 1 #Plot Bridgeness


#---declare clustering algorithms in graph, and with a corresponding consensus matrix
alg    <- ALGS[c(1:10)]

#---Clustering results of the SVI algorithm
SVIFILES <- vector(length=2)
SVIFILES[1] <- "network.gml"
SVIFILES[2] <- "groups.txt"

#---VIP genes
VIP    <- vector(length=2)
VIP[1] <- "BrProteins_Pre.csv"
VIP[2] <- "BrProteins_PSD_reduced.csv"

vips=2;
if( grepl("Pre",subDIR[S]) ){
    vips=1;
}

#---read VIP gene list
VIPs   <- read.table(VIP[vips],sep="\t",header=F)[[1]]
VIPs   <- unique(VIPs)


#---load corresponding graph which was used to build the consensus matrices from 
gg <- igraph::read.graph(sprintf("%s/%s/%s.gml",OUT[1],subDIR[S],subDIR[S]),format="gml")

if(runBridge[1]){

    N    <- length(V(gg)$name)
    
    mm   <- readRDS(sprintf("%s/originalNtwrkMeasures.rds",OUT[2]))

    INDX <- which(names(mm)==subDIR[S])

    INDR <- match(V(gg)$name,mm[[INDX]][,1])
    
    CN  <- c('ENTREZ.ID','GENE.NAME','DEGREE','Closeness','Bet','CC','Cl','Clnorm','SP','PR',alg,sprintf("BRIDGE_%s",alg),sprintf("PROB_%s",alg),sprintf("MIX_%s",alg),'BR_CONSENSUS','BR_CONSENSUS_MAD',"BR_CONSENSUS_ADJ")
    
    meas     <- matrix(0, nrow=N, ncol=length(CN))
    colnames(meas) <- CN

    meas[,1] <- as.character(V(gg)$name)
    meas[,2] <- as.character(V(gg)$GeneName)

    meas[,3]  <- as.character(mm[[INDX]][INDR,2])
    meas[,4]  <- as.numeric(round(closeness(gg,mode="all",normalized=T),3))
    meas[,5]  <- as.numeric(mm[[INDX]][INDR,3])
    meas[,6]  <- as.numeric(mm[[INDX]][INDR,4])
    meas[,7]  <- as.numeric(mm[[INDX]][INDR,5])
    meas[,8]  <- as.character( (as.numeric(meas[,7]) - min(as.numeric(meas[,7])))/(max(as.numeric(meas[,7])) - min(as.numeric(meas[,7]))) )
    meas[,9]  <- as.numeric(mm[[INDX]][INDR,6])
    meas[,10] <- as.numeric(mm[[INDX]][INDR,7])


#START filling meas after PageRank column
FROM <- which(CN=="PR")

N <- length(V(gg))
M <- length(E(gg))

#run over each algorithm
for( a in 1:length(alg) ){

    cat("calculating Bridgeness for: ", alg[a], "\n")
    
    if( a != 10 ){
    
    #--- build reference matrix from the graph
    refin     <- matrix("",ncol=3,nrow=length(V(gg)$name))
    refin[,1] <- igraph::get.vertex.attribute(gg,"name",V(gg))
    refin[,2] <- igraph::get.vertex.attribute(gg,"name",V(gg))
    refin[,3] <- igraph::get.vertex.attribute(gg,alg[a],V(gg))

    #--- store alg. cluster results
    meas[,(FROM+a)] <- as.numeric(refin[,3])

    #Read in consensus matrix
    filein = read.csv(file=sprintf("%s/%s/%s/consensusmatrix.txt",rndDIR[1],subDIR[S],alg[a]), header=FALSE, sep=",");
    dimnames(filein)[2] <- dimnames(filein)[1]

    #format reference matrix
    #the reference matrix with the correct row.names(as a data.frame)
    refmat = as.matrix(refin);
    ref    = refmat[,2:3];
    rm           <- data.frame(ref);
    rownames(rm) <- rm$X1;
    rm$X1        <- NULL;
    names(rm)    <- 'cm';

    #format consensus matrix
    #the consensus matrix you may have made (as a numeric matrix) 
    conmat = as.matrix(filein);
    cm           <- data.frame(conmat);
    names(cm)    <- rownames(rm);
    rownames(cm) <- rownames(rm);
    cm           <- as.matrix(cm);

    rm(filein,refin)

    #get consensus matrix indices for each edge in edge list
    indA <- match(igraph::get.edgelist(gg)[,1],rownames(cm))
    indB <- match(igraph::get.edgelist(gg)[,2],rownames(cm))

    dat  <- data.frame(indA,indB)
    
    #get community numbers for each vertex in edge list for the algorithm
    elA    <- as.numeric(igraph::get.vertex.attribute(gg,alg[a],V(gg))[match(igraph::get.edgelist(gg)[,1],V(gg)$name)])
    elB    <- as.numeric(igraph::get.vertex.attribute(gg,alg[a],V(gg))[match(igraph::get.edgelist(gg)[,2],V(gg)$name)])
    
    ed      <- matrix(ncol=6,nrow=length(E(gg)))
    ed[,1]  <- igraph::get.edgelist(gg)[,1]
    ed[,2]  <- igraph::get.edgelist(gg)[,2]
    ed[,3]  <- elA
    ed[,4]  <- elB
    ed[,5]  <- apply(dat,1,function(x,mat) mat[x[1],x[2]], mat=cm)
    ed[,6]  <- abs(as.numeric(elA)-as.numeric(elB))
    
    Cmax  <- max(as.numeric(igraph::get.vertex.attribute(gg,alg[a],V(gg))))
    
    rm(cm)

    coms <- matrix(ncol=5,nrow=Cmax)

    for( m in 1:Cmax ){        

        Mtot = length(which(elA==m | elB==m))
        
        Min  = length(which(elA==m & elB==m))

        Mout = Mtot - Min

        Cn   = length(which(as.numeric(igraph::get.vertex.attribute(gg,alg[a],V(gg)))==m ))

        Kw=max(as.numeric(meas[,3]))
        
        coms[m,1] <- m
        coms[m,2] <- 2*Min
        coms[m,3] <- Mout
        coms[m,4] <- Cn
        coms[m,5] <- Kw
        
    }

    #container, to store probabilities of gene belonging to each community
    Vprobs           <- matrix(0,nrow=N,ncol=(2+Cmax))
    CNv              <- c("ENTREZ.ID","GENE.NAME",seq(1,Cmax,1))
    colnames(Vprobs) <- CNv   
        
    #loop over each vertex in the graph
    for( i in 1:length(V(gg)) ){

        #get edges belonging to the i'th veretx
        ind <- which(ed[,1] == V(gg)$name[i] | ed[,2] == V(gg)$name[i])

        #get community belonging to the i'th vertex       
        c <- as.numeric(igraph::get.vertex.attribute(gg,alg[a],V(gg)))[i]

        #reorder edge communities, so ed[,3] equals current community no: 'c'
        for( k in 1:length(ind) ){
            if( ed[ind[k],6] != 0 && ed[ind[k],4] == c ){
                ed[ind[k],4] <- ed[ind[k],3]
                ed[ind[k],3] <- c
            }
        }
        
        
        #prob of i'th vertex being found in community Cn relative to a random model
        #Lancichinetti et al, Statistical significance of communities in networks (2010). Physics Review E, 81, 046110.  
        K   <- length(ind)
        Kin <- 0
        indX <- which(ed[ind,3] == c & ed[ind,4] == c)
        if( length(indX) != 0 ){
            Kin <- length(indX)
        }
        Kout <- (K - Kin)

        #eqn 1 in (Lancichinetti et al, 2010)    
        t1 <- choose( (coms[c,3] - Kout), Kin )
        t2 <- choose( (2*M - sum(coms[c,2:3]) - K - coms[c,3] -Kout), (K-Kin) )
        t3 <- choose( (2*M - sum(coms[c,2:3]) - K), K)

        prob <- t1*t2/t3

        if( is.na(prob) ){ prob <- 0 }

        #PROB_
        meas[i,(FROM+2*length(alg)+a)] <- (1-prob)

        #MIX_
        meas[i,(FROM+3*length(alg)+a)] <- Kout/K

        #number of communities i'th vertex is connected too (via it's edges)
        ##cc <- unique(ed[ind,6])
        cc <- unique(ed[ind,4])
        
        #use sum of consensus values to calculate the likelihood of i'th
        #vertex beloning to to k'th community. 
        prob <- vector(length=length(cc))
        for( k in 1:length(cc) ){
            ##prob[k] = sum(as.numeric(ed[which(ed[ind,6]==cc[k]),5]))/length(ind)
            prob[k] = sum(as.numeric(ed[which(ed[ind,4]==cc[k]),5]))/length(ind)
        }

        #normalise
        prob <- prob/sum(prob)
        
        #calculate bridgeness of i'th vertex
        #Fuzzy communities and the concept of bridgeness in complex networks, T. Nepusz, arXiv, 2007
        b    <- sum( (prob - 1/Cmax) * (prob - 1/Cmax))

        Kzero <- Cmax - length(cc)
        b = b + sum(rep((1/(Cmax*Cmax)),times=Kzero))
        
        #store values
        #BRIDGE_
        meas[i,(FROM+length(alg)+a)]  <- 1-sqrt( Cmax/(Cmax-1) * b )

        #store vertex probs across communities
        INDC <- match(ed[ind,4],colnames(Vprobs))
        for( k in 1:length(INDC) ){
            Vprobs[i,INDC[k]] <- (as.numeric(ed[ind[k],5]) + as.numeric(Vprobs[i,INDC[k]]));
        }
        
            
    }#inner for

        #normalise vertex probs per cluster
        for(i in 1:length(V(gg)) ){
            Vprobs[i,3:length(colnames(Vprobs))] <- as.vector(Vprobs[i,3:length(colnames(Vprobs))])/sum(as.vector(Vprobs[i,3:length(colnames(Vprobs))]))
        }        
                
        #print out vertex probs for each community, using the non probalistic methods
        Vprobs[,1] <- as.character(V(gg)$name)
        Vprobs[,2] <- as.character(V(gg)$GeneName)
        write.table(Vprobs,sprintf("%s/%s_Vprobs.csv",subDIR[S],alg[a]),row.names=F,col.names=T,sep="\t",quote=F)
        

    } else {#run over SVI algorithm

        gg2 <- igraph::read.graph(sprintf("%s/%s/%s",OUT[3],subDIR[S],SVIFILES[1]),format="gml")

        gr  <- read.table(sprintf("%s/%s/%s",OUT[3],subDIR[S],SVIFILES[2]),sep="\t",header=F)
        
        #--- store alg. cluster results
        meas[,(FROM+a)] <- as.numeric(igraph::get.vertex.attribute(gg2,"group",V(gg2))[match(meas[,1],V(gg2)$extid)])+1

        #--- store bridgeness
        b <-  as.numeric(igraph::get.vertex.attribute(gg2,"bridgeness",V(gg2)))/as.numeric(igraph::get.vertex.attribute(gg2,"degree",V(gg2)))
        #BRIDGE_
        meas[,(FROM+length(alg)+a)] <- b[match(meas[,1],V(gg2)$extid)]
        

       #get community numbers for each vertex in edge list for the algorithm
       ed      <- matrix(ncol=6,nrow=length(E(gg)))
       ed[,1]  <- igraph::get.edgelist(gg)[,1]
       ed[,2]  <- igraph::get.edgelist(gg)[,2]
       ed[,3]  <- as.numeric(igraph::get.vertex.attribute(gg2,"group",V(gg2))[match(ed[,1],V(gg2)$extid)])+1
       ed[,4]  <- as.numeric(igraph::get.vertex.attribute(gg2,"group",V(gg2))[match(ed[,2],V(gg2)$extid)])+1
       ed[,6]  <- abs(as.numeric(ed[,3])-as.numeric(ed[,4]))
    
        Cmax  <- max(as.numeric(igraph::get.vertex.attribute(gg2,"group",V(gg2))))+1
    
        coms <- matrix(ncol=5,nrow=Cmax)

        for( m in 1:Cmax ){        

            Mtot = length(which(ed[,3]==m | ed[,4]==m))
        
            Min  = length(which(ed[,3]==m & ed[,4]==m))

            Mout = Mtot - Min

            Cn   = length(which(as.numeric(igraph::get.vertex.attribute(gg2,"group",V(gg2)))==(m-1) ))

            Kw=max(as.numeric(meas[,3]))

            coms[m,1] <- m
            coms[m,2] <- 2*Min
            coms[m,3] <- Mout
            coms[m,4] <- Cn
            coms[m,5] <- Kw
        
        }
    
        #loop over each vertex in the graph
        for( i in 1:length(V(gg2)) ){

            #get edges belonging to the i'th vertex
            ind <- which(ed[,1] == V(gg2)$extid[i] | ed[,2] == V(gg2)$extid[i])

            #get community belonging to the i'th vertex
            c <- as.numeric(igraph::get.vertex.attribute(gg2,"group",V(gg2))[i])+1

            #prob of i'th vertex being found in community Cn relative to a random model
            #Lancichinetti et al, Statistical significance of communities in networks (2010). Physics Review E, 81, 046110.  
            K   <- length(ind)
            Kin <- 0
            indX <- which(ed[ind,3] == c & ed[ind,4] == c)
            if( length(indX) != 0 ){
                Kin <- length(indX)
            }
            Kout <- (K - Kin)

            #use SVI max probability of i'th vertex being found in community Cn relative to a random model
            #PROB_     
            meas[i,(FROM+2*length(alg)+a)] <- as.numeric(max(gr[i,3:length(gr[1,])]))

            #MIX_
            meas[i,(FROM+3*length(alg)+a)] <- Kout/K
            
        
        }#for ith node

    }#ifelse
        
}#for


outfile <- file(sprintf("%s/%s_Measures.csv",subDIR[S],subDIR[S]),"w")
write.table(meas, file=outfile, append=T, row.names=F, col.names=T, sep="\t",quote=F);
close(outfile);

}


if( runBridge[2] ){

#Read-in all measures    
meas <- read.delim(sprintf("%s/%s_Measures.csv",subDIR[S],subDIR[S]),sep="\t",header=T)

#Build Consensus Bridgness 
indC   <- match(sprintf("BRIDGE_%s",alg),colnames(meas))
val    <- apply(meas[,indC],1,median, na.rm=T)
mad    <- apply(meas[,indC],1,MAD, na.rm=T)

XX     <- as.vector(unlist(meas[,indC]))
XXmed  <- median(XX, na.rm=T)
XXmad  <- MAD(XX, na.rm=T)
valAdj <- (val-XXmed)/XXmad
    
#Add Concensus Bridgeness to Measures file
meas[,grepl('BR_CONSENSUS',colnames(meas))]     <- val
meas[,grepl('BR_CONSENSUS_MAD',colnames(meas))] <- mad
meas[,grepl('BR_CONSENSUS_ADJ',colnames(meas))] <- valAdj
    
#Write Measures to file 
write.table(meas,sprintf("%s/%s_Measures.csv",subDIR[S],subDIR[S]), sep="\t", col.names=T, row.names=F,quote=F)

    
##---
## Plot Bridging (B) V. Semi-local Centrality (SL)
##---
    
##---Check or create output plot dir
plotDIR <- sprintf("%s/PLOTS/",subDIR[S])

if( !file_test("-d",plotDIR) ){
    dir.create(plotDIR)
}

##---Check or create output REGIONS dir
regDIR <- sprintf("%s/REGIONS/",subDIR[S])

if( !file_test("-d",regDIR) ){
    dir.create(regDIR)
}

    

#--- set rndm seed and for geom_text_repel
set.seed(42)

Force1    <- vector(length=3)
Force1[1] <- 10
Force1[2] <- 5
Force1[3] <- 10

Force3    <- vector(length=3)
Force3[1] <- 10
Force3[2] <- 1
Force3[3] <- 10
#----    
    

cons <- matrix(0,ncol=3,nrow=length(meas[,1]))
colnames(cons) <- c(colnames(meas)[c(1,2)],"reg")
cons[,1] <- meas[,1]
cons[,2] <- as.character(meas[,2])
cons[,3] <- rep(0,length(meas[,1]))

alg <- c(alg,"CONSENSUS")

#---B V. SL Plot Regions
for( a in 1:length(alg) ){

    dd <- data.frame(meas)
    
    str <- sprintf("BRIDGE_%s",alg[a])
    if( alg[a] == "CONSENSUS" ){
        str <- "BR_CONSENSUS"
    }
    indA    <- which(colnames(dd)==str)
    bridge  <- dd[,indA]

    indB <- which(colnames(dd)=="Cl")
    X    <- as.numeric(as.vector(dd[,indB]))
    X    <- (X-min(X))/(max(X)-min(X))
    Xlab <- "Local Centrality (Cl)"
    
    Y <- as.numeric(as.vector(bridge))    
    Ylab <- "Bridgeness (B)"            

    cat("Plotting B V. SL for: ", alg[a], "\n")
    
    xmin <- min(X)
    xmax <- max(X)
    
    ymin <- min(Y)
    ymax <- max(Y)


    cons <- fillRegions( cons, X, Y, 3 );   
    
    
    #---store Br. V. SL Region data for alg.
    oo <- data.frame(cons[,c(1,3)])
    outfile <- file(sprintf("%s/%s_REGION.csv",regDIR,alg[a]),"w")
    cat("#region",file=outfile,"\n")
    write.table(oo, file=outfile, append=T, row.names=F, col.names=F, sep="\t",quote=F);
    close(outfile);

    genes1 <- ifelse( cons[,3]==quad[1],cons[,2],"")
    genes2 <- ifelse( cons[,3]==quad[2],cons[,2],"")
    genes3 <- ifelse( cons[,3]==quad[3],cons[,2],"")
    genes4 <- ifelse( cons[,3]==quad[4],cons[,2],"")
    genes5 <- ifelse( cons[,3]==quad[5],cons[,2],"")
    genes6 <- ifelse( cons[,3]==quad[6],cons[,2],"")

    VIPsGN <- cons[match(VIPs,cons[,1]),2]
    
    genes1 <- ifelse(!is.na(match(genes1,VIPsGN)),genes1,"")
    genes2 <- ifelse(!is.na(match(genes2,VIPsGN)),genes2,"")
    genes3 <- ifelse(!is.na(match(genes3,VIPsGN)),genes3,"")
    genes4 <- ifelse(!is.na(match(genes4,VIPsGN)),genes4,"")
    genes5 <- ifelse(!is.na(match(genes5,VIPsGN)),genes5,"")
    genes6 <- ifelse(!is.na(match(genes6,VIPsGN)),genes6,"")
   

    GenesUL <- genes1
    GenesUR <- genes2
    GenesLL <- genes3
    GenesLR <- genes4

    if( Scheme == 2 ){

        GenesUL <- genes4
        GenesUR <- genes6
        GenesLL <- genes2
        GenesLR <- genes3
    }
    
    #plot
    gplot <- ggplot(dd,aes(x=as.numeric(as.vector(X)),y=as.numeric(as.vector(Y)) ))+geom_point(aes(alpha=X*Y),colour="magenta",show.legend=F)+

    #ul
    geom_label_repel(aes(label=as.vector(GenesUL)),fontface='bold',color='black',fill='white',box.padding=0.35,point.padding=0.5,segment.color='grey50',xlim = c(0.0,0.5), ylim = c(0.5, 1.0),force=Force1[S],size=rel(5),show.legend=F)+

    #ur
    geom_text_repel(aes(label=as.vector(GenesUR)),color='black',segment.color='grey50',xlim = c(0.5,1.0), ylim = c(0.5,1.0),force=Force3[S],size=rel(4.0),show.legend=F)+  

    #prim-loc-2
    geom_text_repel(aes(label=as.vector(GenesLL)),force=1, xlim=c(0.0,0.6),ylim=c(0.1,0.6),color="black",segment.color='grey50',size=rel(3.0),show.legend=F)+

    #lr
    geom_text_repel(aes(label=as.vector(GenesLR)),force=1, xlim=c(0.5,1.0),ylim=c(0.0,0.5),color="black",size=rel(4.0),show.legend=F)+ 
    
    labs(x=Xlab,y=Ylab,title=sprintf("%s",alg[a]))+
        xlim(c(0,1))+ylim(c(0,1))+
        theme(            
            axis.title.x=element_text(face="bold",size=rel(2.5)),
            axis.title.y=element_text(face="bold",size=rel(2.5)),
            legend.title=element_text(face="bold",size=rel(1.5)),
            legend.text=element_text(face="bold",size=rel(1.5)),
            legend.key=element_blank())+
        theme(panel.grid.major = element_line(colour = "grey40",size=0.2),
              panel.grid.minor = element_line(colour="grey40",size=0.1),
              panel.background = element_rect(fill = "white"),
              panel.border = element_rect(linetype="solid",fill=NA))+

        geom_vline(xintercept=0.5,colour="grey40",size=MainDivSize,linetype=2,show.legend=F)+
        geom_hline(yintercept=0.5,colour="grey40",size=MainDivSize,linetype=2,show.legend=F)+

        Gseg1 +
        Gseg2 +
        Gseg3 +

        GAnno1 +
        GAnno2 +
        GAnno3 +
        GAnno4 +
        GAnno5 +

       png(sprintf("%s/%s_%s_BrVCl.png",plotDIR,subDIR[S],alg[a]),width=WIDTH,height=HEIGHT,units="px")
       print(gplot)
       dev.off()

    
}

}


    
