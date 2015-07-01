setwd("C:/Users/Courtney/Dropbox/LCL-iPSC/GEMMA eQTLs/Hutt iPSCs")

#Pull in ck iPSC gene expression
expr_gene = read.table('OriginGeneExpression_Normalized.txt', header=T, as.is=T, sep='\t', row.names=1)
exprnames = as.matrix(row.names(expr_gene))
gene.bed.ordered= read.table('ensemblCAGETSS_ipsc_sorted.bed', header=F, sep='\t')
nosex = gene.bed.ordered[1:10250,]

samplenames = read.table('Covars.txt', header=T, sep ='\t')
##Re-order samplenames based on array location
samplenames = samplenames[order(samplenames$Order),]
hearts = c(56:60,67:72)
samplenames = samplenames[-hearts,]

colnames(expr_gene)= samplenames$Findiv
ck = read.table("hutt.imputed.73subset.fam")[,2]
overlap.indiv = match(ck,colnames(expr_gene))
overlap.exprs = match(nosex$V4, rownames(expr_gene))
exprs.ordered = as.matrix(expr_gene[overlap.exprs,overlap.indiv])

#Make sure to use the PC regressed data
x.pca = prcomp((t(exprs.ordered)), center = T,scale=T)
PC1_13.residual = matrix(nrow= nrow(exprs.ordered), ncol = ncol(exprs.ordered))
rownames(PC1_13.residual) = rownames(exprs.ordered)
colnames(PC1_13.residual) = colnames(exprs.ordered)
for (i in 1:nrow(exprs.ordered)) {
  model= lm(exprs.ordered[i,]~ x.pca$x[,1:13])
  PC1_13.residual[i,] = resid(model)
}

exprs.o.t = t(PC1_13.residual)

geno.i = read.table('C:/Users/Courtney/Dropbox/LCL-iPSC/GEMMA eQTLs/Final Data/ipsc.eqtlgenotypes.raw', header=T, as.is=T, sep='\t', row.names=1)
geno.i.t = t(geno.i)

library(beeswarm)

pdf("Opposite_Direction_eQTLs.pdf")
eqtl.loc = as.matrix(exprs.o.t[,which(colnames(exprs.o.t) == "ENSG00000164048")])
eqtl.geno = as.matrix(geno.i.t[which(rownames(geno.i.t) == "rs12495221_C"),])
boxplot(eqtl.loc~eqtl.geno, main = "Expression in iPSCs")
beeswarm(eqtl.loc~eqtl.geno, add=T, col="black", vertical=T,pch=20)


dc.exprs = read.table("qqnorm.500ht.gccor.newcovcor.ordered.bimbam")
genes = read.table('qqnorm.500ht.gccor.newcovcor.genenames.txt', as.is=T)
colnames(dc.exprs)=genes$V1
geno.l = read.table('C:/Users/Courtney/Dropbox/LCL-iPSC/GEMMA eQTLs/Final Data/LCLall.eqtlgenotypes.raw', header=T, as.is=T, sep='\t')
geno.l.t = t(geno.l)

eqtl.loc = as.matrix(dc.exprs[,which(colnames(dc.exprs) == "ENSG00000164048")])
eqtl.geno = as.matrix(geno.l.t[which(rownames(geno.l.t) == "rs12495221_C"),])
boxplot(eqtl.loc~eqtl.geno, main = "Expression in LCLs")
beeswarm(eqtl.loc~eqtl.geno, add=T, col= "black", vertical=T,pch=20)


eqtl.loc = as.matrix(exprs.o.t[,which(colnames(exprs.o.t) == "ENSG00000146833")])
eqtl.geno = as.matrix(geno.i.t[which(rownames(geno.i.t) == "rs2571997_A"),])
boxplot(eqtl.loc~eqtl.geno, main = "Expression in iPSCs")
beeswarm(eqtl.loc~eqtl.geno, add=T, col="black", vertical=T,pch=20)

eqtl.loc = as.matrix(dc.exprs[,which(colnames(dc.exprs) == "ENSG00000146833")])
eqtl.geno = as.matrix(geno.l.t[which(rownames(geno.l.t) == "rs2571997_A"),])
boxplot(eqtl.loc~eqtl.geno, main = "Expression in LCLs")
beeswarm(eqtl.loc~eqtl.geno, add=T, col= "black", vertical=T,pch=20)
dev.off()
