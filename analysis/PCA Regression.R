setwd("C:/Users/Courtney/Dropbox/LCL-iPSC/GEMMA eQTLs/Hutt iPSCs")
##Load in expression file (normalized and good probes only - by gene)
expr_gene = read.table('OriginGeneExpression_Normalized.txt', header=T, as.is=T, sep='\t', row.names=1)
samplenames = read.table('Covars.txt', header=T, sep ='\t')
#Re-order samplenames based on array location
samplenames = samplenames[order(samplenames$Order),]

hearts = c(56:60,67:72)
samplenames = samplenames[-hearts,]
colnames(expr_gene) = samplenames$Sample
goodprobes= read.table('ht12_probes_snps_ceu_hg19_af_0.05_map_37.txt', header=T, sep='\t')


batch.f = as.factor(samplenames$Batch)
expr_gene=as.matrix(expr_gene)
batch.residual = matrix(nrow= nrow(expr_gene), ncol = ncol(expr_gene))
rownames(batch.residual) = rownames(expr_gene)
colnames(batch.residual) = colnames(expr_gene)
for (i in 1:nrow(expr_gene)) {
  model= lm(expr_gene[i,]~ batch.f)
  batch.residual[i,] = resid(model)
}

x.pca = prcomp((t(expr_gene)), center = T,scale=T)
x.pca.sum = summary(x.pca)


#Converted categorical covariates to a factor so they are levels 
batch.f = as.factor(samplenames$Batch)
pas.f = as.factor(samplenames$Passage)
matbot.f = as.factor(samplenames$Matrigel.Bottle)
sex.f = as.factor(samplenames$Sex)

#Converted numerical covariates to a numeric so they are continuous
pluri.num =as.numeric(samplenames$Pluri)
novel.num = as.numeric(samplenames$Novelty)
age.num = as.numeric(samplenames$Age)
date.num = as.numeric(samplenames$Num.Date.Banked)

covars = list(batch.f,sex.f,age.num,pluri.num ,novel.num,matbot.f,pas.f, date.num)


#Loop to create data table
lmPCA_heatmap = function(pca, covars, npcs)
{
  results<-c()
  for (f in covars) {
    for (i in 1:npcs)
    {
      s = summary(lm(pca$x[,i]~f));
      results<-c(results,pf(s$fstatistic[[1]],
                            s$fstatistic[[2]],s$fstatistic[[3]], lower.tail = FALSE))
    }
  }
  resultsM<-matrix(nrow = length(covars), ncol = npcs, data =
                     results, byrow = TRUE)
  resultsM
  
}
pcaresults.heatmap = lmPCA_heatmap(x.pca,covars,73)
rownames(pcaresults.heatmap) = c("batch", "sex","age","pluri","novel","matrigel bottle", "passage", "date banked")
myCol <- c("red", "orange", "yellow", "grey")
# Defining breaks for the color scale
myBreaks <- c(0, 1e-5, 1e-2, 0.05, 1)
library(gplots)
heatmap.2(pcaresults.heatmap, col=myCol, breaks=myBreaks, margins=c(5,7),trace="none",main="PCA No Covariates Regressed", key=F, keysize=1.5,density.info="none",cexCol=0.7, cexRow = 0.8, Rowv = FALSE, Colv = FALSE, scale="none")
legend("topleft", fill = myCol, cex=0.8,
       legend = c("<1e-5", "1e-5 to 0.01", "0.01 to 0.05", ">0.05"))

## After regressing out batch
x.pca.batch = prcomp((t(batch.residual)), center = T,scale=T)
x.pca.sum.batch = summary(x.pca.batch)
pcaresults.heatmap.batch = lmPCA_heatmap(x.pca.batch,covars,73)
rownames(pcaresults.heatmap.batch) = c("batch", "sex","age","pluri","novel","matrigel bottle", "passage", "date banked")
heatmap.2(pcaresults.heatmap.batch, col=myCol, breaks=myBreaks, margins=c(5,7),trace="none",main="Batch Regressed", key=F, keysize=1.5,density.info="none",cexCol=0.7, cexRow = 0.8, Rowv = FALSE, Colv = FALSE, scale="none")
legend("topleft", fill = myCol, cex=0.8,
       legend = c("<1e-5", "1e-5 to 0.01", "0.01 to 0.05", ">0.05"))

##Pluri is next most significant so regress that out
batch.residual.pluri = matrix(nrow= nrow(expr_gene), ncol = ncol(expr_gene))
rownames(batch.residual.pluri) = rownames(expr_gene)
colnames(batch.residual.pluri) = colnames(expr_gene)
for (i in 1:nrow(batch.residual)) {
  model= lm(batch.residual[i,]~ pluri.num)
  batch.residual.pluri[i,] = resid(model)
}

x.pca.batch.pl = prcomp((t(batch.residual.pluri)), center = T,scale=T)
x.pca.sum.batch.pl = summary(x.pca.batch.pl)
pcaresults.heatmap.batch.pl = lmPCA_heatmap(x.pca.batch.pl,covars,73)
rownames(pcaresults.heatmap.batch.pl) = c("batch", "sex","age","pluri","novel","matrigel bottle", "passage", "date banked")
heatmap.2(pcaresults.heatmap.batch.pl, col=myCol, breaks=myBreaks, margins=c(5,7),trace="none",main="Pluri Regressed", key=F, keysize=1.5,density.info="none",cexCol=0.7, cexRow = 0.8, Rowv = FALSE, Colv = FALSE, scale="none")
legend("topleft", fill = myCol, cex=0.8,
       legend = c("<1e-5", "1e-5 to 0.01", "0.01 to 0.05", ">0.05"))

##Matrigel Bottle is next most significant so regress that out
batch.residual.pluri.mb = matrix(nrow= nrow(expr_gene), ncol = ncol(expr_gene))
rownames(batch.residual.pluri.mb) = rownames(expr_gene)
colnames(batch.residual.pluri.mb) = colnames(expr_gene)
for (i in 1:nrow(batch.residual.pluri)) {
  model= lm(batch.residual.pluri[i,]~ matbot.f)
  batch.residual.pluri.mb[i,] = resid(model)
}

x.pca.batch.pl.mb = prcomp((t(batch.residual.pluri.mb)), center = T,scale=T)
x.pca.sum.batch.pl.mb = summary(x.pca.batch.pl.mb)
pcaresults.heatmap.batch.pl.mb = lmPCA_heatmap(x.pca.batch.pl.mb,covars,73)
rownames(pcaresults.heatmap.batch.pl.mb) = c("batch", "sex","age","pluri","novel","matrigel bottle", "passage", "date banked")
heatmap.2(pcaresults.heatmap.batch.pl.mb, col=myCol, breaks=myBreaks, margins=c(5,7),trace="none",main="Matrigel Bottle Regressed", key=F, keysize=1.5,density.info="none",cexCol=0.7, cexRow = 0.8, Rowv = FALSE, Colv = FALSE, scale="none")
legend("topleft", fill = myCol, cex=0.8,
       legend = c("<1e-5", "1e-5 to 0.01", "0.01 to 0.05", ">0.05"))
