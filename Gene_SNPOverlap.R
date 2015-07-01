setwd("C:/Users/Courtney/Dropbox/LCL-iPSC/GEMMA eQTLs/Hutt iPSCs")


###In excel remove duplicate ENSG and sort ExpressedGeneList Based on Chromosome and then start##
gene.bed.ordered= read.table('ExpressedGeneList.txt', header=T, sep='\t')
masterlist= read.table('ENSG_CAGETSS.txt', header=T, sep='\t')
CAGE.bed.DConly = merge(gene.bed.ordered, masterlist, by.x = "GeneID", by.y = "ENSG", all.x = F, all.y = F, sort=F)
which(is.na(CAGE.bed.DConly[,7]))
# nas = c(9203:11237)
# CAGE.bed.DConly = CAGE.bed[-nas,]
# which(is.na(CAGE.bed.DConly[,7]))
#write.table(CAGE.bed.DConly, 'Hutt.CAGE.Overlap.txt', sep='\t', quote=F,row.names=F)

setwd("C:/Users/Courtney/Dropbox/LCL-iPSC/GEMMA eQTLs/Hutt iPSCs")
##To directly filter mastercolum file - didn't work
# ENSG.ordered = as.matrix(CAGE.bed$GeneID)
# mastercol= read.table('hutt.imputed.1Mb.mastercols.txt', header=F, sep='\t')
# colnames(mastercol) = c("ENSG", "rs", "Gene", "Pos", "Chr")
# mastercol.CK = merge(mastercol, ENSG.ordered, by.x = "ENSG", by.y = "V1", all.x = F, all.y = F, sort=F)

### Re-order Expression and filter out only DC expressed
expr_gene = read.table('OriginGeneExpression_Normalized.txt', header=T, as.is=T, sep='\t', row.names=1)
exprnames = as.matrix(row.names(expr_gene))

## In excel remove the duplicate genes
matcher = read.table("Hutt.CAGE.Overlap.txt", header =T, sep='\t')
matcher = CAGE.bed
matcherids = matcher$GeneID
matcherind = match(matcherids,exprnames[,1])
exprs.o = expr_gene[matcherind,]
samplenames = read.table('Covars.txt', header=T, sep ='\t')
##Re-order samplenames based on array location
samplenames = samplenames[order(samplenames$Order),]
hearts = c(56:60,67:72)
samplenames = samplenames[-hearts,]
colnames(exprs.o)= samplenames$Findiv
ck = read.table("hutt.imputed.73subset.fam")[,2]
overlap.indiv = match(ck,colnames(exprs.o))
exprs.ordered = exprs.o[,overlap.indiv]
exprs.o.t = t(exprs.ordered)
write.table(exprs.o.t,"Expr.DConly.ordered.txt",row.names=F,col.names=F,quote=F,sep="\t")

## Create PCA file based on re-ordered expression with only Darren's expressed genes
exprs = read.table("Expr.DConly.ordered.txt",header=F,sep="\t")
htpca = prcomp(exprs,scale.=TRUE)
xhtpca = htpca$x
pca_results = summary(htpca)
pca_table = pca_results$importance
write.table(pca_table,"PC_importance_HuttiPSCs.txt",col.names=T,row.names=T,quote=F,sep="\t")
write.table(xhtpca,"hutt.DConly.ordered.pcs",col.names=F,row.names=F,quote=F,sep="\t")

reordercovar = match(ck,samplenames$Findiv)
covar.ordered = samplenames[reordercovar,]

pcs.o = pcs[matcherind,]
idcoefs = read.table("/mnt/lustre/home/cusanovich/500HT/addSNP.coef.3671.square")
idcoefs.o = idcoefs[matcherind,matcherind]
write.table(ordered.ENSG,"ENSGList.DConly.Ordered.txt",row.names=F,col.names=F,quote=F,sep="\t")


x.pca.sum = summary(htpca)
x.pca.sum$importance[2,1]
htpca.unsc = prcomp(exprs.o,scale.=F)

###Create Kinship Square
rawmat = read.table("addSNP1415.coef.3671",header=T)
orderlist = read.table("hutt.imputed.73subset.fam")[,2]

ordermat = matrix(NA,length(orderlist),length(orderlist))
for (i in 1:length(orderlist)) {
  for (j in 1:length(orderlist)) {
    row_index = rawmat[, 1] == orderlist[i] & rawmat[, 2] == orderlist[j] |
                rawmat[, 2] == orderlist[i] & rawmat[, 1] == orderlist[j]
    stopifnot(sum(row_index) == 1)
    ordermat[i, j] = rawmat[row_index, 3]
    ordermat[j, i] = rawmat[row_index, 3]  
  }
}
write.table(ordermat,"addSNP.1415.ordered.txt",row.name=F,col.names=F,quote=F)

# for(i in 1:dim(rawmat)[1]){
#   inds1 = match(rawmat[i,1],orderlist)
#   inds2 = grep(rawmat[i,2],orderlist)
#   ordermat[inds1,inds2] = rawmat[i,3]
#   ordermat[inds2,inds1] = rawmat[i,3]
# }

setwd("C:/Users/Courtney/Dropbox/LCL-iPSC/GEMMA eQTLs/Hutt iPSCs")
##Pulling out overlapping findivs
ck = read.table("hutt.imputed.73subset.fam")[,2]
dc = read.table("hutt.imputed.500ht.fam")[,2]
overlap = match(ck,dc)
# sum(is.na(overlap))
# #18 individual not in Darren's study
# fillers = sample(1:431, 18)
##Generated until no overlap was left
fillers = read.table('FindivList.newck.txt')[,2]
confirm = match(fillers, ck)
# match(overlap,fillers)
# ck.only = na.omit(overlap)
# ck_fillers = append(ck.only,fillers)
# newck = sort(ck_fillers, decreasing=FALSE)
newck = match(fillers,dc)
# dc.fam = read.table("hutt.imputed.500ht.fam")
# new.ck.fam = dc.fam[newck,]
# write.table(new.ck.fam, 'hutt.imputed.newck.fam',row.name=F,col.names=F,quote=F)

##Read in new fam from Plink after filtering
orderlist = read.table('hutt.imputed.newck.fam')[,2]
rawmat = read.table("addSNP1415.coef.3671",header=T)
ordermat = matrix(NA,length(orderlist),length(orderlist))
for (i in 1:length(orderlist)) {
  for (j in 1:length(orderlist)) {
    row_index = rawmat[, 1] == orderlist[i] & rawmat[, 2] == orderlist[j] |
      rawmat[, 2] == orderlist[i] & rawmat[, 1] == orderlist[j]
    stopifnot(sum(row_index) == 1)
    ordermat[i, j] = rawmat[row_index, 3]
    ordermat[j, i] = rawmat[row_index, 3]  
  }
}
write.table(ordermat,"addSNP.1415.newck.ordered.txt",row.name=F,col.names=F,quote=F)

##Create new gene expression and PC file
new.ck.fam = read.table('hutt.imputed.newck.fam')[,2]
dc = read.table("hutt.imputed.500ht.fam")[,2]
newck = match(new.ck.fam,dc)
ck.genes = read.table("ENSGList.DConly.Ordered.txt")
dc.genes= read.table('qqnorm.500ht.gccor.newcovcor.genenames.txt')
gene.overlap = match(ck.genes$V1, dc.genes$V1)
sum(is.na(gene.overlap))
newck.genes = na.omit(gene.overlap)
dc.exprs = read.table("qqnorm.500ht.gccor.newcovcor.ordered.bimbam")
exprs = dc.exprs[newck,newck.genes]

dim(exprs)

write.table(exprs,"qqnorm.newck.gccor.newcovcor.bimbam",col.names=F,row.names=F,quote=F,sep="\t")

htpca = prcomp(exprs,scale.=TRUE)
xhtpca = htpca$x
pca_results = summary(htpca)
pca_table = pca_results$importance
write.table(pca_table,"PC_importance_HuttLCLs.txt",col.names=T,row.names=T,quote=F,sep="\t")

write.table(xhtpca,"qqnorm.newck.gccor.newcovcor.pcs",col.names=F,row.names=F,quote=F,sep="\t")

#Create new gene expr and pc of DC data, but only on 73 indivs (all dc genes)
new.ck.fam = read.table('hutt.imputed.newck.fam')[,2]
dc = read.table("hutt.imputed.500ht.fam")[,2]
newck = match(new.ck.fam,dc)
dc.exprs = read.table("qqnorm.500ht.gccor.newcovcor.ordered.bimbam")
exprs = dc.exprs[newck,]

dim(exprs)

write.table(exprs,"qqnorm.newck.allgenes.gccor.newcovcor.bimbam",col.names=F,row.names=F,quote=F,sep="\t")

htpca = prcomp(exprs,scale.=TRUE)
xhtpca = htpca$x
pca_results = summary(htpca)
pca_table = pca_results$importance
write.table(pca_table,"PC_importance_HuttLCLs_allgenes.txt",col.names=T,row.names=T,quote=F,sep="\t")

write.table(xhtpca,"qqnorm.newck.allgenes.gccor.newcovcor.pcs",col.names=F,row.names=F,quote=F,sep="\t")


##Create files for iPSC all genes
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
exprs.ordered = expr_gene[overlap.exprs,overlap.indiv]
exprs.o.t = t(exprs.ordered)
write.table(exprs.o.t,"Expr.allgenes.ordered.bimbam",row.names=F,col.names=F,quote=F,sep="\t")

## Create PCA file based on re-ordered expression with all expressed genes
exprs = read.table("Expr.allgenes.ordered.bimbam",header=F,sep="\t")
htpca = prcomp(exprs,scale.=TRUE)
xhtpca = htpca$x
pca_results = summary(htpca)
pca_table = pca_results$importance
write.table(pca_table,"PC_importance_HuttiPSCs_allgenes.txt",col.names=T,row.names=T,quote=F,sep="\t")
write.table(xhtpca,"hutt.allgenes.ordered.pcs",col.names=F,row.names=F,quote=F,sep="\t")
write.table(nosex$V4,"ENSGList.allgenes.Ordered.txt",row.names=F,col.names=F,quote=F,sep="\t")

for (i in 1:22) {
  subset = nosex[nosex$V1 == paste("chr", i, sep=""),]
  write.table(subset$V4, paste("Genes.all.Ordered.chr", i,".genes.txt", sep=""), col.names=F, row.names=F,quote=F, sep='\t')
}
write.table(nosex,"Hutt.CAGE.allgenes.bed",row.names=F,col.names=F,quote=F,sep="\t")

###BROKEN
### Make Venn of overlapping genes between LCLs and Hutt iPSC #######
genes <- data.frame(subset$ENSG)
names(genes) <- "genes"
make.venn.quad <- function(geneset1, geneset2,geneset1.label, geneset2.label, univ){
  univ$g1 <- univ$probes %in% geneset1
  univ$g2 <- univ$probes %in% geneset2
  venn.placeholder <- draw.quad.venn(length(geneset1),length(geneset2), dim(univ[univ$g1 == T & univ$g2 == T , ])[1], c(geneset1.label, geneset2.label), fill=c("goldenrod", "plum4"), alpha=c(0.5, 0.5),col=NA, euler.d=T)
 }

dev.off()
# Make venn of full
pdf(file = "VennDiagram_GenesExpressed.pdf")
make.venn.quad(rownames(expr_gene), genes, "iPSC 11,237 genes", "LCLs 13,965 genes", genes)
dev.off()
test = rownames(expr_gene)
overlap1 <- genes$genes %in% rownames(expr_gene)
overlap2 <- genes$genes %in% genes
dim(genes[overlap1 == T & overlap2 == T , ])
draw.quad.venn(11237,13965, 9202, c("iPSC 11,237 genes", "LCLs 13,965 genes"),col=NA, euler.d=T)

matcherfindiv = match(orderlist,rawmat$Ind2)
rawmat$Ind2[14]
orderlist[3]
