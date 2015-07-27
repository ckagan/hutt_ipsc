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
  PC1_13.residual[i,] = resid(model) +model$coefficients[1]
}

exprs.o.t = t(PC1_13.residual)

#library(gplots)
#heatmap.2(cor(as.matrix(PC1_13.residual), use = "complete"), margins=c(7,8),trace="none",main="Gene Expression Correlation", key=T, keysize=1.5,density.info="none",cexCol=0.9)

##Load in iPSC geno
geno.i = read.table('C:/Users/Courtney/Dropbox/LCL-iPSC/GEMMA eQTLs/Final Data/ipsc.eqtlgenotypes.raw', header=T, as.is=T, sep='\t', row.names=1)
geno.i.t = t(geno.i)


#Load in dc data
dc = read.table("qqnorm.500ht.gccor.newcovcor.ordered.bimbam")
genes = read.table('qqnorm.500ht.gccor.newcovcor.genenames.txt', as.is=T)
colnames(dc)=genes$V1
dc.exprs.2 = t(dc)
d.pca = prcomp(dc, center = T,scale=T)
#Regress out 62 PCs
PC1_62.residual = matrix(nrow= nrow(dc.exprs.2), ncol = ncol(dc.exprs.2))
rownames(PC1_62.residual) = rownames(dc.exprs.2)
for (i in 1:nrow(dc.exprs.2)) {
  model= lm(dc.exprs.2[i,]~ d.pca$x[,1:62])
  PC1_62.residual[i,] = resid(model) +model$coefficients[1]
}

dc.exprs = t(PC1_62.residual)

mean(var(PC1_13.residual))
mean(var(PC1_62.residual))

#Load in dc geno
geno.l = read.table('C:/Users/Courtney/Dropbox/LCL-iPSC/GEMMA eQTLs/Final Data/LCLall.eqtlgenotypes.raw', header=T, as.is=T, sep='\t')
geno.l.t = t(geno.l)

library(beeswarm)
##Create pdfs
pdf("Opposite_Direction_eQTLs.pdf", family="Garamond")
eqtl.loc = as.matrix(exprs.o.t[,which(colnames(exprs.o.t) == "ENSG00000164048")])
eqtl.geno = as.matrix(geno.i.t[which(rownames(geno.i.t) == "rs12495221_C"),])
boxplot(eqtl.loc~eqtl.geno, main = "ZNF589 Expression in iPSCs",col=(c("orange", "lightgreen", "lightblue")))
beeswarm(eqtl.loc~eqtl.geno, add=T, col="black", vertical=T,pch=20)


eqtl.loc = as.matrix(dc.exprs[,which(colnames(dc.exprs) == "ENSG00000164048")])
eqtl.geno = as.matrix(geno.l.t[which(rownames(geno.l.t) == "rs12495221_C"),])
boxplot(eqtl.loc~eqtl.geno, main = "ZNF589 Expression in LCLs",col=(c("orange", "lightgreen", "lightblue")))
beeswarm(eqtl.loc~eqtl.geno, add=T, col= "black", vertical=T,pch=20)


eqtl.loc = as.matrix(exprs.o.t[,which(colnames(exprs.o.t) == "ENSG00000146833")])
eqtl.geno = as.matrix(geno.i.t[which(rownames(geno.i.t) == "rs2571997_A"),])
boxplot(eqtl.loc~eqtl.geno, main = expression(paste(italic("TRIM4")," Expression in iPSCs")),col=(c("orange", "lightgreen", "lightblue")), xlab = "Genotype (Number of Minor Alleles)", ylab = "Gene Expression", ylim=c(-2,3))#, xlim =c(0,7), at = 1:3*1, xaxt = "n", ylim = c(-1.5,1.5))
beeswarm(eqtl.loc~eqtl.geno, add=T, col="black", vertical=T,pch=20)

eqtl.loc = as.matrix(dc.exprs[,which(colnames(dc.exprs) == "ENSG00000146833")])
eqtl.geno = as.matrix(geno.l.t[which(rownames(geno.l.t) == "rs2571997_A"),])
boxplot(eqtl.loc~eqtl.geno, main = expression(paste(italic("TRIM4")," Expression in LCLs")),col=(c("orange", "lightgreen", "lightblue")),xlab = "Genotype (Number of Minor Alleles)", ylab = "Gene Expression", ylim = c(-2,3))#, add=T, at = 4:6*1, xaxt = "n")
beeswarm(eqtl.loc~eqtl.geno, add=T, col= "black", vertical=T,pch=20) #, at = 4:6*1)
dev.off()

#Most significant iPSC eQTL ENSG00000011405  rs17560341

eqtl.loc = as.matrix(exprs.o.t[,which(colnames(exprs.o.t) == "ENSG00000011405")])
eqtl.geno = as.matrix(geno.i.t[which(rownames(geno.i.t) == "rs17560341_A"),])
boxplot(eqtl.loc~eqtl.geno, main = "Expression in iPSCs")
beeswarm(eqtl.loc~eqtl.geno, add=T, col="black", vertical=T,pch=20)

#Largest effect size eQTL ENSG00000143632  rs867325
eqtl.loc = as.matrix(exprs.o.t[,which(colnames(exprs.o.t) == "ENSG00000143632")])
eqtl.geno = as.matrix(geno.i.t[which(rownames(geno.i.t) == "rs867325_A"),])
boxplot(eqtl.loc~eqtl.geno, main = "Expression in iPSCs")
beeswarm(eqtl.loc~eqtl.geno, add=T, col="black", vertical=T,pch=20)


boxplot(as.matrix(exprs.o.t[,which(colnames(exprs.o.t) == "ENSG00000181449")]), main = "SOX2 Expression")
boxplot(as.matrix(exprs.o.t[,which(colnames(exprs.o.t) == "ENSG00000204531")]), main = "OCT3/4 Expression")
boxplot(as.matrix(exprs.o.t[,which(colnames(exprs.o.t) == "ENSG00000215271")]), main = "HOMEZ Expression",xlim =c(0,1), at = .2:.2*1, xaxt = "n", ylim = c(-.25,.25))
boxplot(as.matrix(exprs.o.t[,which(colnames(exprs.o.t) == "ENSG00000187140")]), add=T, at = .8:.8*1, xaxt = "n") #main = "FOXD3 Expression",
var(as.matrix(exprs.o.t[,which(colnames(exprs.o.t) == "ENSG00000187140")]))
var(as.matrix(exprs.o.t[,which(colnames(exprs.o.t) == "ENSG00000215271")]))

boxplot(as.matrix(dc.exprs[,which(colnames(dc.exprs) == "ENSG00000215271")]), main = "HOMEZ LCL Expression") 
boxplot(as.matrix(dc.exprs[,which(colnames(dc.exprs) == "ENSG00000187140")]), main = "FOXD3 LCL Expression") 
var(as.matrix(dc.exprs[,which(colnames(dc.exprs) == "ENSG00000215271")]))
var(as.matrix(dc.exprs[,which(colnames(dc.exprs) == "ENSG00000187140")]))
