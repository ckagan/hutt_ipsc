setwd("C:/Users/Courtney/Desktop/GEMMATest")
master = read.table('hutt.imputed.1Mb.chrmspecific.mastercols.txt', header=F, sep='\t')
genelist = rownames(expr_gene)
genelist = as.matrix(genelist)
colnames(genelist) = c("GeneID")
colnames(master) = c("Gene", "rs", "Num", "Count", "Chr")
#Merge actual list with reference information from good probe document
chrlist.gene = merge(genelist, master, by.x = "GeneID", by.y = "Gene", all.x = T, all.y = F, sort=F)
clean.chrlist.gene = chrlist.gene[!duplicated(chrlist.gene$GeneID), ]
which(is.na(chrlist.gene[,2]))
chrlist.gene[2292650,]


ensg = cbind(goodprobes[,1:3], goodprobes[,7], goodprobes[,5:6])
gene.bed = merge(genelist, ensg, by.x = "GeneID", by.y = "goodprobes[, 7]", all.x = T, all.y = F, sort=F)
which(is.na(gene.bed[,2]))
test = unique(gene.bed)
# No unique values
write.table(gene.bed, 'ExpressedGeneList.txt', sep='\t', row.names=F, quote=F)

gene.bed.ordered= read.table('ExpressedGeneList.txt', header=T, sep='\t')
masterlist= read.table('ENSG_CAGETSS.txt', header=T, sep='\t')
CAGE.bed = merge(gene.bed.ordered, masterlist, by.x = "GeneID", by.y = "ENSG", all.x = T, all.y = F, sort=F)
which(is.na(CAGE.bed[,7]))
nas = c(13700:16436)
CAGE.bed.DConly = CAGE.bed[-nas,]
which(is.na(CAGE.bed.DConly[,7]))
write.table(CAGE.bed.DConly, 'Hutt.CAGE.Overlap.txt', sep='\t', quote=F)

##To directly filter mastercolum file
ENSG.ordered = as.matrix(CAGE.bed$GeneID)
mastercol= read.table('hutt.imputed.1Mb.mastercols.txt', header=F, sep='\t')
colnames(mastercol) = c("ENSG", "rs", "Gene", "Pos", "Chr")
mastercol.CK = merge(mastercol, ENSG.ordered, by.x = "ENSG", by.y = "V1", all.x = F, all.y = F, sort=F)

### Re-order Expression and filter out only DC expressed
exprnames = as.matrix(row.names(expr_gene))
idcoefs = read.table("/mnt/lustre/home/cusanovich/500HT/addSNP.coef.3671.square")
## In excel remove the duplicate genes
matcher = read.table("Hutt.CAGE.Overlap.txt", header =F, sep='\t')
ordered.ENSG = matcher$V4
matcherids = matcher$V4
matcherind = match(matcherids,exprnames[,1])
exprs.o = expr_gene[matcherind,]

## Create PCA file based on re-ordered expression with only Darren's expressed genes
htpca = prcomp(exprs.o,scale.=TRUE)
xhtpca = htpca$x
write.table(xhtpca,"hutt.DConly.Ordered.pcs.txt",col.names=F,row.names=F,quote=F,sep="\t")


pcs.o = pcs[matcherind,]
idcoefs.o = idcoefs[matcherind,matcherind]
write.table(ordered.ENSG,"ENSGList.DConly.Ordered.txt",row.names=F,col.names=F,quote=F,sep="\t")
write.table(exprs.o,"Expr.DConly.Ordered.txt",row.names=F,col.names=F,quote=F,sep="\t")

x.pca.sum = summary(htpca)
x.pca.sum$importance[2,1]
htpca.unsc = prcomp(exprs.o,scale.=F)
