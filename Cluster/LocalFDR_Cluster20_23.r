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

##Analysis for loop
qcounts = c()
#for(i in 1:length(pcs)){
for(i in 20:length(pcs)){
  j = as.character(pcs[i])
  PC.pvals <- list.files(path = "/mnt/gluster/data/internal_supp/hutt_ipsc/Genotypes/ByChr/",pattern = paste("PC", j, ".imputed.1Mb.bonferroni.regressPCs.gemma.eqtls.txt", sep=""),full.names=T)
  #For each PC combine all chrs
  master=do.call("rbind", lapply(PC.pvals, read.table, header = FALSE))
  #Get a gene list
  genes = master$V1
  genelist= unique(genes)
  all=c()
  ##For each gene do an fdr correction  
  for(k in 1:length(genelist)){
    gene.temp=genelist[[k]]
    rows.temp = which(master$V1 == gene.temp)
    gene.only = master[rows.temp,]
    corrected = p.adjust(gene.only$V3,method="fdr")
    total = cbind(gene.only,corrected)
    ##Create an all file that containes the fdr correction (plus gene/SNP and origina pval)
    all = rbind(all,total)
  }
  #Get the number of gene/SNP pairs at a 1% FDR
  qcount = sum(all$corrected < 0.01)
  qcounts[i] = qcount
  #Write out results for future use if needed
  write.table(all, paste("PC",j,"_eQTLResults.LocalFDR.txt"), sep='\t', quote=F)
  #If you want to create a file in R 
  #assign(paste("PC",j, sep=""), all)
  #Create a pdf of all p-value distributions for each PC
  pdf("eQTLHistogram.pdf")
      hist(all$corrected,xlab="P-value",main=paste(j," PCs Removed"))
      rm(all)
      rm(master)
}
dev.off()