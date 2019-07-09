#--- Merge output files in RESULTS dir ---#

library(igraph)

#---Script requires two inputs
# args[1] == The network file under study
#       1 == PPI_Presynaptic_Published
#       2 == PPI_PSP_clean_Published
#       3 == PPI_PSP_reduced
# args[2] == The number of iterations per job
args  <- commandArgs(TRUE);
ss    <- as.numeric(args[1])         #the network
PERMS_PER_JOB <- as.numeric(args[2]) #no: of iterations per job


#load corresponding graph
#---Directories
DIR    <- vector(length=2)
DIR[1] <- "../../datasets/"
DIR[2] <- "RESULTS"

#---Subdirectories, of a study, were consensus files are stored 
subDIR    <- vector(length=3)
subDIR[1] <- "PPI_Presynaptic_Published"
subDIR[2] <- "PPI_PSP_clean_Published"
subDIR[3] <- "PPI_PSP_reduced"

#---load corresponding graph which was used to build the consensus matrices from 
gg <- read.graph(sprintf("%s/%s.gml",DIR[1],subDIR[ss]),format="gml")

#---load batch parameters
NJOBS        = 0
NPERMS       = 0

#--- check output Dir
if( !file_test("-d",subDIR[ss]) ){
    dir.create(subDIR[ss])
}

#---HDO ID DISEASES of INTEREST
disn    <- vector(length=11);
disn[1]  <- "DOID:10652"#Alzheimer's disease
disn[2]  <- "DOID:3312" #bipolar disorder
disn[3]  <- "DOID:12849"#autistic disorder
disn[4]  <- "DOID:5419" #schizophrenia
disn[5]  <- "DOID:0060041"#autism_spectrum_disorder
disn[6]  <- "DOID:1826" #epilepsy_syndrome
disn[7]  <- "DOID:1059" #Intellectual Disability
disn[8]  <- "DOID:10763"#Hypertension
disn[9]  <- "DOID:12858"#Huntington's disease
disn[10] <- "DOID:14330"#parkinson's disease 
disn[11] <- "DOID:9255" #frontotemporal dementia
#disn[12] <- "DOID:150" #Parent term for BI, SCH, ID, ASD, AUT
#disn[13] <- "DOID:331" #Parent term for AD, PD, HD, Epi

#---HDO Disease short names
dtype  <- vector(length=11);
dtype[1]   = "AD";
dtype[2]   = "BI";
dtype[3]   = "AUT";
dtype[4]   = "SCH";
dtype[5]   = "ASD";
dtype[6]   = "Epi";
dtype[7]   = "ID";
dtype[8]   = "HTN";
dtype[9]   = "HD";
dtype[10]  = "PD";
dtype[11]  = "FTD";
#dtype[12]  = "DMH";
#dtype[13]  = "CNSD";


#---load gda's for the ppi network
GDA <- V(gg)$TopOntoOVG
NN  <- length(which(GDA!=""))

#---Remove Diseases with zero GDA's
remove <- c()
for( d in 1:length(dtype) ){
    IDS <- V(gg)$name[grepl(dtype[d],GDA,fixed=T)]
    if( length(IDS) == 0 ){
        cat(dtype[d], " => ", length(IDS),"\n")
	remove <- c(remove,d)
     }
}
if( length(remove) > 0 ){
   disn  <- disn [-remove]	
   dtype <- dtype[-remove]
}    
#--- 

oo <- matrix(".",nrow=NN,ncol=(length(dtype)+2))
colnames(oo) <- c("Gene.ID","Gene.Name",dtype)
oo[,1] <- V(gg)$name[GDA !=""]
oo[,2] <- V(gg)$GeneName[GDA !=""]

for( d in 1:length(dtype) ){

    IDS <- V(gg)$name[grepl(dtype[d],GDA,fixed=T)]

    N   <- length(IDS)
    ds  <- rep(0,N)
    for( i in 1:N ){

        ds[i] <- min(as.vector(shortest.paths(gg,IDS[i],IDS[-i])))

        indX <- which(oo[,1]==IDS[i])
        oo[indX,(2+d)] <- ds[i]
    }

}


#---Get all gene Entrez IDS
GNS <- V(gg)$name

tot <- matrix(0,nrow=length(GNS),ncol=(length(dtype)+2))
colnames(tot) <- c("Gene.ID","Gene.Name",dtype)
tot[,1] <- V(gg)$name
tot[,2] <- V(gg)$GeneName


#---loop through output files in RESULTS dir
files    <- vector(length=4)
files[1] <- "gene_disease_separation.csv"  
files[2] <- "mean_disease_separation.csv"
files[3] <- "sd_disease_separation.csv" 
files[4] <- "sAB_random_separation.csv"

subdirs  = list.files(path=sprintf("../%s/",DIR[2]));
nstudies = length(subdirs);

for( s in 1:nstudies ){
  
  st1 = sprintf("../%s/%s/%s",DIR[2],subdirs[s],files[1]);
  
  if( file.exists(st1) && file.info(st1)$size!=0 ){

     NJOBS=NJOBS+1

     tb = read.table(st1,header=T,sep="\t");
     
     for( d in 1:length(dtype) ){
     	  tot[,(2+d)] <- as.numeric(tot[,(2+d)]) + as.numeric(tb[,(2+d)])
     }    

  }
}

NPERMS       = NJOBS * PERMS_PER_JOB

stats <- matrix(0,ncol=2,nrow=3)
stats[1,1] <- as.character("Network")
stats[1,2] <- subDIR[ss]
stats[2,1] <- as.character("NJOBS")
stats[2,2] <- as.character(NJOBS)
stats[3,1] <- as.character("PERMS_PER_JOBS")
stats[3,2] <- as.character(PERMS_PER_JOB)

outfile <- file(sprintf("%s/stats.csv",subDIR[ss]),"w");
write.table( stats, file=outfile, append=T, row.names=F, col.names=T, sep="\t", quote=F);
close(outfile);


#Normalise Random ds
for( i in 1:length(dtype) ){
     tot[,(2+i)] <- as.numeric(as.vector(tot[,(2+i)]))/NPERMS
}

outfile <- file(sprintf("%s/random_%s",subDIR[ss],files[1]),"w");
write.table( tot, file=outfile, append=T, row.names=F, col.names=T, sep="\t", quote=F);
close(outfile);

Nn    <- length(dtype)
NELE  <- Nn*(Nn+1)/2

#---build output containers for script
#raw sAB of disease-disease overlap
RAWsAB <- matrix(0,nrow=NELE, ncol=(2+NPERMS))

#raw mean of disease-disease overlap
RAWm <- matrix(0,nrow=NELE, ncol=(2+nstudies))

#raw mean of disease-disease overlap
RAWsd <- matrix(0,nrow=NELE, ncol=(2+nstudies))

#mean of Disease-disease overlap
DABm <- matrix(0,nrow=NELE, ncol=3)

#SD of Disease-disease overlap
DABsd <- matrix(0,nrow=NELE, ncol=3)

for( k in 0:(NELE-1) ){

  #--- indexing for symmetric matrix
  i = floor( (2*Nn+1 - sqrt( (2*Nn+1)*(2*Nn+1) - 8*k ))/2 );
  j = k - Nn*i + i*(i-1)/2;

  i = i + 1;
  j = j + i;

  RAWsAB[(k+1),1] <- dtype[i]
  RAWsAB[(k+1),2] <- dtype[j]

  RAWm[(k+1),1] <- dtype[i]
  RAWm[(k+1),2] <- dtype[j]

  RAWsd[(k+1),1] <- dtype[i]
  RAWsd[(k+1),2] <- dtype[j]

  DABm[(k+1),1] <- dtype[i]
  DABm[(k+1),2] <- dtype[j]

  DABsd[(k+1),1] <- dtype[i]
  DABsd[(k+1),2] <- dtype[j]

 }

k=2
Mnorm=0
SDnorm=0

#---loop over all jobs in RESULTS dir.
for( s in 1:nstudies ){
  
  st1 = sprintf("../%s/%s/%s",DIR[2],subdirs[s],files[2]);
  st2 = sprintf("../%s/%s/%s",DIR[2],subdirs[s],files[3]);
  st3 = sprintf("../%s/%s/%s",DIR[2],subdirs[s],files[4]);
  
  if( (file.exists(st1) && file.info(st1)$size!=0) && 
      (file.exists(st2) && file.info(st2)$size!=0) && 
      (file.exists(st3) && file.info(st3)$size!=0) ){

     tb1 = read.table(st1,header=T,sep="\t");
     tb2 = read.table(st2,header=T,sep="\t");
     tb3 = read.table(st3,header=T,sep="\t");	

     ncols <- length(colnames(tb3))-2

     #raw random sAB vlaues 
     for( i in 1:ncols ){
        k=k+1
        RAWsAB[,k] <- as.numeric(tb3[,(2+i)])
     }

     #raw mean
     RAWm[,(2+s)] <- as.numeric(tb1[,3])

     #raw sd
     RAWsd[,(2+s)] <- as.numeric(tb2[,3])

     #pooled mean     
     DABm[,3] <- as.numeric(DABm[,3]) + (PERMS_PER_JOB) * as.numeric(as.vector(tb1[,3]))

     Mnorm = Mnorm + PERMS_PER_JOB

     #pooled sd     
     DABsd[,3] <- as.numeric(DABsd[,3]) + (PERMS_PER_JOB -1)*as.numeric(tb2[,3])

     SDnorm = SDnorm + (PERMS_PER_JOB -1)

     }    

  }

DABm[,3] <- as.numeric(DABm[,3])/Mnorm

DABsd[,3] <- as.numeric(DABsd[,3])/SDnorm

outfile <- file(sprintf("%s/RAW_random_%s",subDIR[ss],files[2]),"w");
write.table( RAWm, file=outfile, append=T, row.names=F, col.names=T, sep="\t", quote=F);
close(outfile);

outfile <- file(sprintf("%s/RAW_random_%s",subDIR[ss],files[3]),"w");
write.table( RAWsd, file=outfile, append=T, row.names=F, col.names=T, sep="\t", quote=F);
close(outfile);

outfile <- file(sprintf("%s/random_%s",subDIR[ss],files[2]),"w");
write.table( DABm, file=outfile, append=T, row.names=F, col.names=T, sep="\t", quote=F);
close(outfile);

outfile <- file(sprintf("%s/random_%s",subDIR[ss],files[3]),"w");
write.table( DABsd, file=outfile, append=T, row.names=F, col.names=T, sep="\t", quote=F);
close(outfile);

outfile <- file(sprintf("%s/RAW_random_%s",subDIR[ss],files[4]),"w");
write.table( RAWsAB, file=outfile, append=T, row.names=F, col.names=T, sep="\t", quote=F);
close(outfile);
