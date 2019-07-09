
source('../setUp.R')

#Check or create out dir
if( !file_test("-d",subDIR[S]) ){
    dir.create(subDIR[S])
}
    
#---Clustering algorithms used 
alg    <- ALGS[c(1:10)]

#No consensus matrix for SVI
alg <- alg[!grepl("SVI",alg,fixed=T)]

#set CDF parameters
steps  = 100;
stepsi = 1.0/steps;
CDF <- vector(length=steps);

X   <- vector(length=steps)
for(x in 1:steps ){
    Xi=stepsi*x;
    X[x] <- Xi;
}      


#output results
oo <- matrix("",ncol=(1+length(alg)),nrow=steps)
colnames(oo) <- c("X",alg)
oo[,1]       <- X

#run over each algorithm
for( a in 1:length(alg) ){

    #print the algorithm we're running over
    cat("running over ", alg[a], "\n");

    #Read in consensus matrix
    CM = read.csv(file=sprintf("%s/%s/%s/consensusmatrix.txt",rndDIR[1],subDIR[S],alg[a]), header=FALSE, sep=",");

    #reset CDF vector
    CDF <- rep(0,length(CDF))

    #build consensus matrix CDF
    Fn  <- ecdf(CM[lower.tri(CM)])

    CDF <- Fn(X[1:steps])

    #store CDF result
    oo[,(1+a)] <- CDF

    #output file
    outfile <- file(sprintf("%s/%s_consensusCDFs.csv",subDIR[S],subDIR[S]),"w")
    write.table(oo, file=outfile, append=T, row.names=F, col.names=T, sep="\t",quote=F);
    close(outfile);

    rm(CM);
    
}

#---proportion of ambiguously clustered pairs (PAC)
# Senbabaoglu Y., Michailidis G. and Li J. Critical Limitations of consensus clustering in class discovery. Scientific Reports, 4, 6207, 2014. doi:10.1038/srep06207

#set parameters
X1  = 10 #0.1
X2  = 90 #0.9
PAC = -1

N   = length(oo[1,])

PAC <- as.numeric(oo[X2,2:N]) - as.numeric(oo[X1,2:N])
ALG <- colnames(oo)[2:N]

#print PAC results for each algorithm
cat("Alg. names  : ", ALG,"\n")
cat("PAC values  : ", PAC,"\n")

#find algorithm with minimum PAC value
PACmin <- which.min(PAC)

#print algorithm with min PAC value
cat("Algorithm ", ALG[PACmin], " has min PAC value of ", PAC[PACmin],"\n")


#---
# Plot consensus CDF for each algorithm
#---

CN <- colnames(oo)[2:N]
X  <- length(oo[,1])

#---set ploting colours
#colours    <- c('darkturquoise','lawngreen','dodgerblue','orange','firebrick2','slateblue2','forestgreen','midnightblue','deeppink2')
colours    <- c('royalblue','royalblue','grey70','grey70','firebrick2','royalblue','deeppink3','royalblue','royalblue')
#colours    <- c('firebrick2','blue2','deepskyblue3','cyan','magenta','darkorange3','orange','deeppink3','grey50');


#---reshape the oo (i.e. the output) dataframe to use in ggplot 
df <- data.frame(a=as.character(),b=as.character(),c=as.numeric(),d=as.numeric())

for( i in 2:N ){

    label <- as.character(rep(CN[i-1],X))
    col   <- as.character(rep(colours[i-1]),X)
    xx    <- as.numeric(oo[,1])
    cdf   <- as.numeric(oo[,i])

    df    <- rbind(df,data.frame(label,col,xx,cdf))
    
}

colnames(df) <- c("ALG","COL","X","CDF")

#---Check or create plot dir
plotDIR <- sprintf("%s/PLOTS/",subDIR[S])

if( !file_test("-d",plotDIR) ){
    dir.create(plotDIR)
}
  

#---Generate CDF plot
SIZE=2.0
gplot <- ggplot(df,aes(x=as.numeric(df$X),y=as.numeric(as.vector(df$CDF)),colour=df$ALG))+
    geom_line(size=SIZE,alpha=0.75)+
    labs(x="c",y="CDF(c)")+
    theme(axis.title.x=element_text(face="bold",size=rel(2.0)),
          axis.title.y=element_text(face="bold",size=rel(2.0)),
          legend.text=element_text(face="bold",size=rel(1.5)),
          legend.key=element_blank())+
    scale_color_manual("",breaks=c(levels(factor(df$ALG))),values=c(colours))+
    theme(panel.grid.major = element_line(colour = "grey40"),
    panel.grid.minor = element_line(colour="grey40",size=0.1),
    panel.background = element_rect(fill = "white"),
     panel.border = element_rect(linetype="solid",fill=NA))+
    geom_vline(xintercept=as.numeric(X1/steps),colour="grey10",size=SIZE,linetype=2,show.legend=F)+
    geom_vline(xintercept=as.numeric(X2/steps),colour="grey10",size=SIZE,linetype=2,show.legend=F)+
    guides(color = guide_legend(order=1),
           alpha = FALSE)


png(sprintf("%s/%sCDF.png",plotDIR,subDIR[S]),width=WIDTH,height=HEIGHT,units="px");
print(gplot)
dev.off()

gplot2 <- ggplot(df,aes(x=log10(as.numeric(df$X)),y=as.numeric(as.vector(df$CDF)),colour=df$ALG))+
    geom_line(size=SIZE,alpha=0.75)+
    labs(x="log(c)",y="CDF(c)")+
    theme(axis.title.x=element_text(face="bold",size=rel(2.0)),
          axis.title.y=element_text(face="bold",size=rel(2.0)),
          legend.text=element_text(face="bold",size=rel(1.5)),
          legend.key=element_blank())+
    scale_color_manual("",breaks=c(levels(factor(df$ALG))),values=c(colours))+
    theme(panel.grid.major = element_line(colour = "grey40"),
    panel.grid.minor = element_line(colour="grey40",size=0.1),
    panel.background = element_rect(fill = "white"),
     panel.border = element_rect(linetype="solid",fill=NA))+
    geom_vline(xintercept=log10(as.numeric(X1/steps)),colour="grey10",size=SIZE,linetype=2,show.legend=F)+
    geom_vline(xintercept=log10(as.numeric(X2/steps)),colour="grey10",size=SIZE,linetype=2,show.legend=F)+
    guides(color = guide_legend(order=1),
           alpha = FALSE)


png(sprintf("%s/%sCDFlog.png",plotDIR,subDIR[S]),width=WIDTH,height=HEIGHT,units="px");
print(gplot2)
dev.off()
