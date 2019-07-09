#--- first run source("runClustercons.R")

#--- convert resultmatrix back into a matrix
tt <- as.matrix(mr$resultmatrix@mrm);

#--- get max number of communities from refin
coms <- max(refin$V3);

dp     <- matrix(ncol=2,nrow=length(refin$V2)+coms);
dp[,2] <- -1;

k <- 1;
for( i in 1:coms ){

  dp[k,1] <- i;
  k = k + 1;
  
  #--- get ids of proteins in community coms
  ids <- refin$V2[which(refin$V3==as.integer(i))];

  for( j in 1:length(ids) ){ 
    rr  <- which(row.names(tt)==as.character(ids[j]));
    val <- tt[rr,i];
    dp[k,1] <- as.integer(ids[j]); 
    if(is.na(val) ){
      dp[k,2] <- as.double (0.0);
    }else{
      dp[k,2] <- as.double (val);
    }
    k = k + 1;
  }

}


#--- write outfile
dp  <- dp[dp[,2] != -1,];
dpt <- as.data.frame(dp);
#--- format output values
dpt[,2] <- sprintf("%.3f",dpt[,2]);
outfile <- file("memroblist.txt","w");
cat("#memroblist\n",file=outfile);
write.table(dpt,file=outfile,row.names=F,col.names=F,sep="\t",quote=F);
close(outfile);

#--- write outfile
dpt      <- as.data.frame(cr);
#--- format output values
dpt[,1] <- sprintf("%.3f",dpt[,1]);
outfile2 <- file("clusterrobustness.txt","w");
cat("#clusterrobustness\n",file=outfile2);
write.table(dpt,file=outfile2,row.names=T,col.names=F,sep="\t",quote=F);
close(outfile2);
