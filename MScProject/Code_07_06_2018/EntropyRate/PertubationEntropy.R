# PPI network Permutation analysis
# Based on Figure 4C in: A. E. Teschendorff et al. INcreased signaling entropy in cancer requires the scale-free property of protein interaction networks, Scientific Reports, 5:9646.

source('../setUp.R')

#---Directories needed
OUT <- vector(length=3)
OUT[1] <- DIRS[grepl("Graphs",DIRS)]
OUT[2] <- DIRS[grepl("EntropyRate",DIRS)]

#---Check or create output dir
if( !file_test("-d",subDIR[S]) ){
    dir.create(subDIR[S])
}

maxLSi <- function( XX, BASE=0 ){

    XX  <- as.vector(as.numeric(XX))
    
    XXo <- XX[XX != 0]
    if( BASE == 2 ){
        return ( -sum( XXo * log2(XXo)) )
    }

    if( BASE == 10 ){
        return ( -sum( XXo * log10(XXo)) )
    }

    if( BASE == 0 ){
        return ( -sum( XXo * log(XXo)) )
    }
    
}

#---load corresponding graph which was used to build the consensus matrices from 
gg <- igraph::read.graph(sprintf("%s/%s/%s.gml",OUT[1],subDIR[S],subDIR[S]),format="gml")


#--- initial entropy rate
V    <- length(V(gg))
E    <- length(E(gg))
ki   <- as.vector(igraph::degree(gg))
Kbar <- mean(ki)

#--- get adjacency matrix for graph
A    <- get.adjacency(gg)

#--- get leading eigenvalue and vector
R     <- eigen(A)
Rindx <- which.max(R$values)
gamma <- R$values[Rindx]
nu    <- R$vectors[,Rindx]


#--- calculate max entropy rate, maxSr
Pij   <- (A * nu) / (gamma * nu)
Pi    <- ki/(2*E)
maxSr <- sum(as.vector(Pi * apply(Pij,1,maxLSi,BASE=0)))


#--- calculate initial configuration
Norm <- as.numeric(V*Kbar) #as.numeric(2*E)
SRo  <- as.numeric(1/Norm)*sum(ki*log(ki))


#--- perturbated each PPI node/gene

#--- expression values
xx    <- vector(length=2)
xx[1] <- 2  #active
xx[2] <- 16 #inactive

#--- perturbation expression values
lambda    <- vector(length=2)
lambda[1] <- 14   #active
lambda[2] <- -14  #inactive

#--- Norm for PI'
NORM      <- vector(length=2)
NORM[1]   <- 0
NORM[2]   <- 0

SRprime <- cbind(V(gg)$GeneName,ki,rep("",V),rep("",V))

for( v in 1:V ){

    #--- name of gene to perturb
    GN     <- as.character(SRprime[v,1])
    GNindx <- which(V(gg)$GeneName==GN)

    #--- PI'
    PIprime <- cbind( rep("",V), rep("",V) )

    #--- LS'
    LSprime <- cbind( rep("",V), rep("",V) )

    #--- reset NORM
    NORM[1] = 0; NORM[2] = 0;
    
    #--- calculate norm for PI'
    for( s in 1:length(lambda) ){
        X               <- rep(xx[s], V)
        X[GNindx[1]]    <- X[GNindx[1]] + lambda[s]
        NORM[s]         <- X %*% A %*% X
    }
    

    #--- find all neighors to v, i.e. N(v)
    Nv <- V(gg)$name[neighbors(gg,GNindx,mode="all")]

    oo <- cbind( ki, !(V(gg)$name %in% Nv) )
    

    #--- PI' when v is not N(v)
    for( s in 1:length(lambda) ){
        PIprime[,s] <- ifelse(oo[,2] == 1, (1/NORM[s] * xx[s] * xx[s] * as.numeric(oo[,1])), ".")
    }

    #--- PI' when v is v
    for( s in 1:length(lambda) ){

        X   <- as.numeric(xx[s])
        lam <- as.numeric(lambda[s])
        DEG <- as.numeric(oo[GNindx[1],1])

        PIprime[GNindx[1],s] <- ((X + lam) * DEG * X) / NORM[s] 
        
    }
    
    #--- PI' when v is N(v) 
    for( s in 1:length(lambda) ){
        PIprime[,s] <- ifelse(oo[,2] == 0, (1/NORM[s] * xx[s] * ( xx[s] + lambda[s] + (as.numeric(oo[,1]) - 1) * xx[s])),PIprime[,s])
    }


    #--- LS' when v is not N(v)
    for( s in 1:length(lambda) ){
        X <- as.numeric(xx[s])
        LSprime[,s] <- ifelse(oo[,2] == 1, (-log(X) + log(X*as.numeric(oo[,1]))),".")
    }

    #--- LS' when v is N(v)     
    Ni <- grep(0,oo[,2])

    for( i in 1:length(Ni) ){

        DEGi <- as.numeric(oo[Ni[i],1])
        SUM  <- DEGi-1     
        
        for( s in 1:length(lambda) ){

            X   <- as.numeric(xx[s])
            lam <- as.numeric(lambda[s])
                        
            dem <- X + lam + (DEGi -1) * X
            
            pij <- X / dem
            pi1 <- (X + lam) / dem

            LSi <- pij * log(pij)
            LSi <- - SUM * LSi - pi1 * log(pi1)
                    
            LSprime[Ni[i],s]  <- as.character(LSi)
    
        }

    }
    
    SRprime[v,3] <- sum( as.numeric(PIprime[,1]) * as.numeric(LSprime[,1]) )
    SRprime[v,4] <- sum( as.numeric(PIprime[,2]) * as.numeric(LSprime[,2]) )

}


#--- PLOT
subTIT    <- vector(length=3)
subTIT[1] <- ""#"Presynaptic"
subTIT[2] <- ""#"PSD"
subTIT[3] <- ""#"PSD"

colours <- c('lawngreen','firebrick2')

SRprime[,3] <- as.numeric(SRprime[,3])/maxSr
SRprime[,4] <- as.numeric(SRprime[,4])/maxSr

colnames(SRprime) <- c("GENE.NAME","DEGREE","SR_UP","SR_DOWN")

write.table(SRprime,sprintf("%s/SignalEntropyRate.csv",subDIR[S]),sep="\t",row.names=F,col.names=T,quote=F)

DF1 <- SRprime[,c(1,2,3)]
DF1 <- cbind(DF1,rep("UP",length(SRprime[,1])))

DF2 <- SRprime[,c(1,2,4)]
DF2 <- cbind(DF2,rep("DOWN",length(SRprime[,1])))

#--- Bottom 1% UP, i.e. OVER-EXPRESSED
XX  <- as.numeric(SRprime[,3])
MIN <- min(XX)
MAX <- max(XX)
XX2 <- (XX-MIN)/(MAX-MIN)
oo2 <- cbind(SRprime[,1],XX2)
oo2 <- oo2[order(as.numeric(oo2[,2])),]
ii  <- floor(1/100 * V)
GN  <- oo2[1:ii,1]
DF3 <- SRprime[match(GN,SRprime[,1]),c(1,2,3)]
DF3 <- cbind(DF3,rep("1%",length(GN)))

write.table(DF3,sprintf("%s/SignalEntropyRate_1percent_UP.csv",subDIR[S]),sep="\t",row.names=F,col.names=T,quote=F)
#---

DF  <- rbind(DF1,DF2)
colnames(DF) <- c("GENE.NAME","DEGREE","SR","GROUP")

DF <- as.data.frame(DF)

gplot <- ggplot(DF,aes(x=log(as.numeric(as.vector(DF$DEGREE))),y=as.numeric(as.vector(DF$SR)), colour=DF$GROUP) )+
    geom_point()+
    labs(x="log(k)",y="SR",title=subTIT[S])+
    theme(            
        axis.title.x=element_text(face="bold",size=rel(2)),
        axis.text.x =element_text(face="bold",size=rel(2)), 
        axis.title.y=element_text(face="bold",size=rel(2)),
        axis.text.y =element_text(face="bold",size=rel(2)), 
        legend.title=element_text(face="bold",size=rel(1.5)),
        legend.text=element_text(face="bold",size=rel(1.5)),
        legend.key=element_blank())+
    theme(panel.grid.major = element_line(colour = "grey40",size=0.2),
          panel.grid.minor = element_line(colour="grey40",size=0.1),
          panel.background = element_rect(fill = "white"),
          panel.border = element_rect(linetype="solid",fill=NA))+
     scale_color_manual("",breaks=levels(factor(DF$GROUP)),values=c(colours))+
   geom_hline(yintercept=SRo/maxSr,colour="grey40",size=2,linetype=2,show.legend=F)

png(sprintf("%s/SignalEntropyRate.png",subDIR[S]),width=WIDTH,height=HEIGHT,units="px");
print(gplot)
dev.off()
