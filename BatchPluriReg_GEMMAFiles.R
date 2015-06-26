expr_gene = read.table('GeneExpression_Normalized_BatchPluriReg_GEMMA.txt', header=T, as.is=T, sep='\t', row.names=1)
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
exprs.ordered = expr_gene[overlap.exprs,overlap.indiv]
exprs.o.t = t(exprs.ordered)
write.table(exprs.o.t,"Expr.allgenes.ordered.batchplurreg.bimbam",row.names=F,col.names=F,quote=F,sep="\t")

## Create PCA file based on re-ordered expression with all expressed genes
exprs = read.table("Expr.allgenes.ordered.batchplurreg.bimbam",header=F,sep="\t")
htpca = prcomp(exprs,scale.=TRUE)
xhtpca = htpca$x
pca_results = summary(htpca)
pca_table = pca_results$importance
write.table(xhtpca,"hutt.allgenes.ordered.batchplurreg.pcs",col.names=F,row.names=F,quote=F,sep="\t")
