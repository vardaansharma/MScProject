#--------
# Once you've run "enrichmentAnalysis.R script, 
# to generate the output file "data_for_sigmoid_fit.csv"
# can run this script.
#--------

source('../setUp.R')

require("minpack.lm") 
require(cowplot)

### Sigmoid function ### create a function to generate sigmoid pattern
sigmoid <- function(pars, xx){

    a = as.numeric(pars[1])#lower asymptote,          ideal == 0
    b = as.numeric(pars[2])#upper asymptote,          ideal == 1
    c = as.numeric(pars[3])#gradiant, rate, or slope, ideal == -2
    d = as.numeric(pars[4])#inflextion point,         ideal == median(xx) = 3

    return( a + ((b-a)/(1+exp(-c*(xx-d)))) )
}

residFun <- function(pars, observed, xx) observed - sigmoid(pars,xx)

storeFitInfo <- function( ALG, LL ){

    tmp <- c(ALG) 
    
    for( i in 1:length(LL) ){
        tmp <- cbind(tmp, as.vector(LL)[i])
    }

    return(tmp)
    
}

storeParInfo <- function( ALG, LL ){

    tmp <- c(ALG) 
    
    for( i in 1:length(LL[,1]) ){
        tmp <- cbind(tmp, LL[i,1], LL[i,2])
    }

    return(tmp)
    
}


highlightRate <- function( rates, val=-2){

    indx <- which(rates==val)
    if( length(indx) == 0 ){
        indx = -1
    }

    return(indx)
    
}

plotSigmoid <- function( x, rates, model, alg="", pv=0 ){

    conf     <- NULL
    try(conf <- confint(model), FALSE)

    if( !is.null(conf) ){
        # adding the 95% confidence interval around the fitted coefficient
        lower = list(a=conf[1, 1], b=conf[2, 1], c=conf[3, 1], d=conf[4, 1])
        upper = list(a=conf[1, 2], b=conf[2, 2], c=conf[3, 2], d=conf[4, 2])
    }
        
    #fitted values
    y    <- model$m$lhs()
    yhat <- as.vector(fitted(model))

    if( !is.null(conf) ){
        ylower <- sigmoid(pars=lower, xx=x)
        yupper <- sigmoid(pars=upper, xx=x)
    } else {
        ylower <- rep(0, length(y))
        yupper <- rep(0, length(y))
    }
    #---

    #--- data.frame to plot
    df <- cbind( rep(alg,length(x)), x, y, yhat, ylower, yupper )

    #--- 'ideal' sigmoid curves, changing the rate values
    R = length(rates)
    Rsize = rep(1,R)
    Rcol  = rep("grey50",R)
    indx  = highlightRate( rates=rates, val=-2 )
    if( indx != -1 ){
        Rsize[indx[1]] = 2
        Rcol[indx[1]]  = "black"
    }
    for( r in 1:R ){
        pp = list(a=0, b=1, c=rates[r], d=round(median(x)) )
        yi <-  sigmoid(pars=pp, xx=x)
        df <- cbind(df,yi)
    }
    #---

    #--- data.frame to plot
    colnames(df) <- c("alg", "x", "y", "yhat", "ylower", "yupper", sprintf("yiR%f", seq(1,R,1)))
    df <- as.data.frame(df)

    #--- labels
    qq   = as.vector(quantile(as.numeric(df$x)))
    xval = qq
    xlab = as.character(qq)

    ylow <- ifelse( min(y) < 0, min(y), 0)
    yup  <- ifelse( max(y) > 1, max(y), 1)
    #---

    #--- should we plot the CI
    plotCI=TRUE
    if( is.null(conf) ){ plotCI = FALSE }

    #--- p.value from KS test for model against 'ideal' sigmoid with rate value = -2
    pv = as.numeric(pv)
    if( is.na(pv) ){ pv = 0 }
    #---

    #--- build the plot
    gplot <- ggplot(df, aes(as.numeric(df$x)))+
        geom_point(aes(y=as.numeric(as.vector(df$y))),   shape=1, size=2.5)+
        geom_line(aes(y=as.numeric(as.vector(df$yhat))), linetype="dashed", color="red", size=2)+
        geom_line(aes(y=as.numeric(as.vector(df$yiR1))), linetype="solid", color=Rcol[1], size=Rsize[1])+
        geom_line(aes(y=as.numeric(as.vector(df$yiR2))), linetype="solid", color=Rcol[2], size=Rsize[2])+
        geom_line(aes(y=as.numeric(as.vector(df$yiR3))), linetype="solid", color=Rcol[3], size=Rsize[3])+
        geom_line(aes(y=as.numeric(as.vector(df$yiR4))), linetype="solid", color=Rcol[4], size=Rsize[4])+
        geom_line(aes(y=as.numeric(as.vector(df$yiR4))), linetype="solid", color=Rcol[5], size=Rsize[5])+
        {if(plotCI)geom_line(aes(y=as.numeric(as.vector(df$ylower))), linetype="dashed", color="blue", size=2)}+
        {if(plotCI)geom_line(aes(y=as.numeric(as.vector(df$yupper))), linetype="dashed", color="blue", size=2)}+
        labs(x="log2(Fe)",y="Fraction of Enriched Communities",title=sprintf("%s, KS p.value = %3.e", alg, pv))+
        theme(axis.title.x=element_text(face="bold",size=rel(1.5)),
              axis.title.y=element_text(face="bold",size=rel(1.5)),
              legend.text=element_text(face="bold",size=rel(1.5)),
              plot.title=element_text(face="bold",size=rel(1.5)),
              legend.position="bottom")+
        coord_cartesian(xlim = c(min(x), max(x)), ylim=c(ylow, yup))+
        #scale_y_continuous(expand=c(0,0),limits=c(ylow,yup))+
        #scale_x_discrete(expand=c(0,0), limit=xval, labels=xlab)+
        theme(panel.grid.major = element_line(colour = "grey40"),
              panel.grid.minor = element_line(colour="grey40",size=0.1),
              panel.background = element_rect(fill = "white"),
              panel.border = element_rect(linetype="solid",fill=NA))+
        guides(color = FALSE,
               alpha = FALSE,
               size  = FALSE)
    #---
    
    return(list(gplot=gplot, df=df))
    
}

errorY <- function( x, N=100, SD=0.05){
    x = x+rnorm(N, 0, SD)
    return(sd(x)/sqrt(length(x)))
    
}

sigmaSqY <- function( x, N=100, SD=0.05){

    x = x+rnorm(N, 0, SD)

    N  = length(x)
    mu = mean(x)

    return( (1/N) * sum( (x-mu)^2 ) )
    
}

addNoise <- function( Y, MN=0, SD=0.05 ){
    return( Y+rnorm(length(Y),mean=MN, sd=SD) )
}


#goodness of fit test, KS
gofs <- function(x, rate, model, sigma2=NULL, countDATA=TRUE ){

    y    = model$m$lhs()
    yhat = fitted(model)

    R  = length(rate)
    KS = list()
    
    for( r in 1:R ){
    
        pp = list(a=0, b=1, c=rate[r], d=round(median(x)) )
        yi = sigmoid(pars=pp, xx=x)

        KS[[r]]      = ks.test(yhat, yi)
        names(KS)[r] = sprintf("rate_%f.1", rate[r]) 
    }

    return(KS)    
   
}


#---OUT Dir
OUT    <- vector(length=4)
OUT[1] <- DIRS[grepl("EnrichmentPackage",DIRS)]
OUT[2] <- DIRS[grepl("Graphs",DIRS)]
OUT[3] <- DIRS[grepl("parameterFiles",DIRS)]
OUT[4] <- DIRS[grepl("Annotations",DIRS)]

plotDIR <- sprintf("%s/PLOTS",OUT[1])
if( !file_test("-d",plotDIR) ){
    dir.create(plotDIR)
}
#---

tt <- read.delim("data_for_sigmoid_fit.csv",sep="\t", header=T, check.names=F)

N = length(colnames(tt))
x = as.numeric(colnames(tt)[3:N])

#Gaussian noise, the set of standard deviations test 
SDv   = c(0, 0.05, 0.1, 0.5)
SDlab = c("0","0.05","0.1","0.5")

#test fit against different 'idealised' sigmoid curves using KS gof statistic
rates = c(-10, -5, -2, -1, -0.5)

for( s in 1:length(SDv) ){

    models <- list()
    GPLOTS <- list()
    gof    <- list()

    #save the fit infomation for each run
    CNfit             <- c("alg", "isConv", "finTol", "stopCode", "stopMessage")
    #fitInfo           <- matrix("", ncol=length(CNfit), nrow=length(tt[,1]))
    #colnames(fitInfo) <- CNfit
    #fitInfo[,1]       <- tt[,1]

    fitInfo <- data.frame()
    parInfo <- data.frame()
    
    for( i in 1:length(tt[,1]) ){

        y = as.numeric(tt[i,3:N])/as.numeric(tt[i,2])
        y = addNoise(y, SD=SDv[s])#add gaussian noise to our data
    
        #starting parameter values for fit
        pp = list(a=0, b=round(max(y)), c=-2, d=round(median(x)) )
    
        #fit 
        m.s <- nlsLM(y ~ a + ((b - a)/(1 + exp(-c * (x - d)))), start = pp, trace = FALSE)
    
        models[[i]]      <- m.s
        names(models)[i] <- as.character(tt[i,1]) 

        rm(m.s)

        fitInfo <- rbind(fitInfo, storeFitInfo( names(models)[i], unlist(models[[i]]$convInfo) ) )

        parInfo <- rbind(parInfo, storeParInfo( names(models)[i], unlist(summary(models[[i]])$parameters[,1:2]) ))
        
    }

    #---save the fit infomation
    colnames(fitInfo) <- c("alg", "isConv", "finIter", "finTol", "stopCode", "stopMessage")
    write.table(fitInfo, sprintf("%s/FitInfo_sd_%s.csv",plotDIR, SDlab[s]), sep="\t", col.names=T, row.names=F, quote=F)

    colnames(parInfo) <- c("alg", c(rbind(names(models[[1]]$m$getPars()),sprintf("sd_%s",names(models[[1]]$m$getPars())))))
    write.table(parInfo, sprintf("%s/parInfo_sd_%s.csv",plotDIR, SDlab[s]), sep="\t", col.names=T, row.names=F, quote=F)
        
    #---print KS test results for each rate and noisy study here
    CN <- c("alg", sprintf("Rate_%.1f",rates))
    oo <- matrix("", ncol=length(CN), nrow=length(names(models)))
    colnames(oo) <- CN
    oo[,1] <- names(models)
    
    for( i in 1:length(names(models)) ){

        ks   <- gofs(x, rates, models[[i]])
        indx <- highlightRate( rates=rates, val=-2 )
        PV = as.numeric(ks[[1]]$p.value)
        if( indx != -1 ){
            PV = as.numeric(ks[[indx[1]]]$p.value)
        }
        tmp <- plotSigmoid( x=x, rates=rates, model=models[[i]], alg=names(models)[i], pv=PV)

        gof[[i]]    <- ks
        GPLOTS[[i]] <- tmp$gplot

        names(gof)[i]    <- names(models)[i]
        names(GPLOTS)[i] <- names(models)[i]

        for(j in 1:length(rates)){
            oo[i,(j+1)] = ks[[j]]$p.value
        }
    }

    p <- plot_grid(plotlist=GPLOTS, labels = "AUTO", label_size=50, label_fontface="bold")
    ggsave(sprintf("%s/Fitting_%s.png",plotDIR,SDlab[s]), p, width=20, height=20, device="png")

    write.table(oo, sprintf("%s/ks_pv_sd_%s.csv",plotDIR, SDlab[s]), sep="\t", col.names=T, row.names=F, quote=F)
    
    rm(models, GPLOTS, gof)

}
