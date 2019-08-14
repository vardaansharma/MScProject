rm(list=ls())

#Set the absolute path to this working directory
version="17_05_2019"
#mainDir <- sprintf("/afs/inf.ed.ac.uk/user/c/cmclean5/WORK/DATA/Human_synaptosome_%s",version)

#mainDir <- "C:/Users/homoe/OneDrive - University of Edinburgh/PROJECT/MScProject/PPI network/clustering"
mainDir <- "/afs/inf.ed.ac.uk/user/s18/s1882216/MSCPROJECT/MScProject/PPI network/clustering"


#Get Path to all top-level directories in folder
DIRS <- list.dirs(mainDir,recursive=F)

#---dataDIR, where to read PPI files
dataDIR <- DIRS[grepl("datasets",DIRS)]

#---subDIR for loading and storing PPI Graphs
subDIR <- list.files(path=dataDIR)
subDIR <- subDIR[grep(".csv",subDIR)]
subDIR <- unlist(strsplit(subDIR,".csv"))

#---Location for randomisation files
rndDIR    <- vector(length=1)
#rndDIR[1] <- sprintf("/disk/scratch/WORK/DATA/Human_synaptosome_%s",version)

#---C++ Dir
CPP    <- vector(length=4)
CPP[1] <- "Spectral"
CPP[2] <- "SPICi"
CPP[3] <- "SVI"
CPP[4] <- "Geodesic"

#---Declare all clustering algorithms
ALGS    <- vector(length=29)
ALGS[1]  <- "fc"#
ALGS[2]  <- "sgG1"#
ALGS[3]  <- "sgG2"
ALGS[4]  <- "sgG5"
ALGS[5]  <- "Spectral"#
ALGS[6]  <- "louvain"#
ALGS[7]  <- "infomap"#
ALGS[8]  <- "lec"
ALGS[9]  <- "wt"
ALGS[10] <- "SVI"
ALGS[11] <- "SPICi"
ALGS[12] <- "Geodesic"
ALGS[13] <- "CONSENSUS"
ALGS[14] <- "louvain2"
ALGS[15] <- "Spectral1per"
ALGS[16] <- "Spectral25per"#
ALGS[17] <- "Spectral5per"
ALGS[18] <- "fc2"
ALGS[19] <- "lec2"
ALGS[20] <- "wt2"
ALGS[21] <- "sgG12"
ALGS[22] <- "Spectral05per"
ALGS[23] <- "Spectral025per"
ALGS[24] <- "Spectral01per"
ALGS[25] <- "Spectral1per2"
ALGS[26] <- "Spectral25per2"
ALGS[27] <- "Spectral5per2"
ALGS[28] <- "Spectral05per2"
ALGS[29] <- "Spectral01per2"


#---HDO ID DISEASES of INTEREST
disn    <- vector(length=12);
disn[1]  <- "DOID:10652"#Alzheimer's_disease"
disn[2]  <- "DOID:3312"#bipolar_disorder"
disn[3]  <- "DOID:12849"#autistic_disorder"
disn[4]  <- "DOID:5419"#schizophrenia"
disn[5]  <- "DOID:0060041"#autism_spectrum_disorder
disn[6]  <- "DOID:1826"#epilepsy_syndrome
disn[7]  <- "DOID:1059"
disn[8]  <- "DOID:10763"
disn[9]  <- "DOID:12858"
disn[10] <- "DOID:14330"
disn[11] <- "DOID:9255"
disn[12] <- "DOID:2377"
#disn[12] <- "DOID:150" #Parent term for BI, SCH, ID, ASD, AUT
#disn[13] <- "DOID:331" #Parent term for AD, PD, HD, Epi


#---HDO Disease short names
dtype  <- vector(length=12);
dtype[1]   = "AD";
dtype[2]   = "BD";
dtype[3]   = "AUT";
dtype[4]   = "SCH";
dtype[5]   = "ASD";
dtype[6]   = "Epi";
dtype[7]   = "ID";
dtype[8]   = "HTN";
dtype[9]   = "HD";
dtype[10]  = "PD";
dtype[11]  = "FTD";
dtype[12]  = "MS";
#dtype[12]  = "DMH";
#dtype[13]  = "CNSD";

#---HDO Disease long names
disl  <- vector(length=12);
disl[1]   = "Alzheimer's_disease";
disl[2]   = "Bipolar_disorder";
disl[3]   = "Autistic_disorder";
disl[4]   = "Schizophrenia";
disl[5]   = "Autism_Spectrum_Disorder";
disl[6]   = "Epilepsy_syndrome";
disl[7]   = "Intellectual_Disability";
disl[8]   = "Hypertension";
disl[9]   = "Huntington's_disease";
disl[10]  = "Parkinson's_disease";
disl[11]  = "Frontotemporal_dementia";
disl[12]  = "Multiple_Sclerosis"
#disl[12]  = "disease_of_mental_health";
#disl[13]  = "central_nervous_system_disease";

#---parent HDO term
pHDO    <- vector(length=2)
pHDO[1] <- "DMH";
pHDO[2] <- "CNSD"; 

#---For Bonferroni correction
alpha <- vector(length=3)
alpha[1] <- 0.05
alpha[2] <- 0.01
alpha[3] <- 0.001

stars    <- vector(length=3)
stars[1] <- "*"
stars[2] <- "**"
stars[3] <- "***"

COLLAPSE <- vector(length=2)
COLLAPSE[1] <- ";"
COLLAPSE[2] <- "&"
c=1

#--- Set 
#--- Pre-load Graph of interest, stored in file 'graphs.csv'
pramFILES <- DIRS[grepl("parameterFiles",DIRS)]
Graph <- read.table(sprintf("%s/graphs.csv",pramFILES),header=F,sep=",",quote="")
S     <- as.vector(Graph[which(as.vector(Graph[,1]) == 1)[1],2])
S     <- match(S,subDIR)
#S     <- grep(S,subDIR)


#---All annotation types
Anno     <- vector(length=18)
Anno[1]  <- "topOnto_ovg"
Anno[2]  <- "SCH"
Anno[3]  <- "chua"
Anno[4]  <- "Domain"
Anno[5]  <- "Family"
Anno[6]  <- "GOMF"
Anno[7]  <- "GOBP"
Anno[8]  <- "GOCC"
Anno[9]  <- "CT_PMID27991900"
Anno[10] <- "CT_L2_PMID27991900"
Anno[11] <- "CT_PMID25700174"
Anno[12] <- "CT_PMID27471252"
Anno[13] <- "CT_PMID27716510"
Anno[14] <- "CT_PMID28602351"
Anno[15] <- "PW_PMID25599223"
Anno[16] <- "GOMF_PMID25599223"
Anno[17] <- "GOBP_PMID25599223"
Anno[18] <- "GOCC_PMID25599223"


#--- gene-disease annotations
# 1) ovg: ususal OMIN/Ensembl Var./geneRIF
# 2) ovPAPERS: OMIN/Ensembl Var. + p140Cap Gwas studies.
gdaDIR    <- vector(length=2)
gdaDIR[1] <- "ovg"
gdaDIR[2] <- "ovPAPERS" 

gdas <- 1

#--- TAX ID
TAXid    <- vector(length=3)
TAXid[1] <- "10090" #Mouse
TAXid[2] <- "9606"  #Human
TAXid[3] <- "7227"  #Fly

GOonto    <- vector(length=3)
GOonto[1] <- "MF"
GOonto[2] <- "BP"
GOonto[3] <- "CC"

#---WIDTH and HEIGHT for plots
WIDTH=480
HEIGHT=480

#set required R libraries
#source("http://www.bioconductor.org/biocLite.R")

library(igraph);
library(lattice);
library(ggplot2);
library(stringr);
library(gtable);
library(grid);
library(ggrepel);
library(scales);
# library(clusterCons);
library(plyr);
library(VennDiagram);
# library(Vennerable);
# library(biomaRt);
# library(latex2exp);
# library(knitr);
# library(poweRlaw);
# library(WriteXLS);
# library(gdata);
# library(methods);
# library(ggpubr);
# library(DBI);
# library(aricode);
# library(reshape2);

# library(org.Hs.eg.db)
# library(org.Mm.eg.db)
# library(org.Dm.eg.db)
# library(topGO)

#---Set Scheme in Regions.R for B V. SL Plot Regions in Bridgeness.R
#--- Select from 1: Basic scheme
#--- Select from 2: Deferential non-bridging, primary and secondray, bridging proteins

Scheme = 1
source(sprintf("%s/Regions.R",pramFILES))

#---Selected Regions for Venn, i.e. 1,2,3,4 
#   or combinations
Rsel <- c("1")


#set default R options 
options(stringsAsFactors=F)

printSessionInfo <- TRUE
#printSessionInfo <- FALSE

if( printSessionInfo ){
    print(sessionInfo())
}

#---Print Graph of interest
cat("\n")
cat("\n")
cat("***********\n")
cat("PPI Graph/file is: ", subDIR[S], " set in 'graphs.csv' \n")
cat("***********\n")
cat("\n")
cat("\n")


