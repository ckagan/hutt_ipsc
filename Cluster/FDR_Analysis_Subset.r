setwd("/mnt/gluster/data/internal_supp/hutt_ipsc/Genotypes/ByChr")
library(plyr)
library(stringr)

#Get PC names
all.pvals <- list.files(path = "/mnt/gluster/data/internal_supp/hutt_ipsc/Genotypes/ByChr/",pattern=".imputed.1Mb.bonferroni.regressPCs.gemma.eqtls.txt",full.names=T)
pval.pcs.names = str_split(all.pvals,"[.]")
pval.pcs = c()
for(i in 1:length(pval.pcs.names)){
  pval.pcs[i] = str_split(pval.pcs.names[[i]][2],"PC")[[1]][2]
}
pcs.uni = unique(pval.pcs)
pcs = sort(as.numeric(pcs.uni), decreasing=FALSE)
print("Starting PC analysis")
##Analysis for loop
qcounts = c()
for(i in 1:length(pcs)){
  j = as.character(pcs[i])
  master <- read.table(paste("/mnt/gluster/data/internal_supp/hutt_ipsc/Genotypes/ByChr/PC",j,"_eQTLResults.LocalFDR.txt"),header=T)
  ##Create 1% FDR files
  #subset = master[master$corrected < 0.01,]
  #write.table(subset, paste("/mnt/gluster/data/internal_supp/hutt_ipsc/Genotypes/ByChr/PC",j,"_eQTLResults.LocalFDR.1%.txt", sep=""), quote=F, row.names=F, sep='\t')
  #Get the number of gene/SNP pairs at a 1% FDR
  qcount = sum(master$corrected < 0.01)
  qcounts[i] = qcount
}

pdf('PCAnalysis.pdf')
plot(1:length(pcs),qcounts,xlab="Number of PCs Removed",ylab="No. of eQTLs at Local FDR 1%",type="b",pch=20)
dev.off()
print("Finished PC analysis")