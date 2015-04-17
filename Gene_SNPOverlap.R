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
# No unique values - however multiple entries for ~5k genes
#write.table(gene.bed, 'ExpressedGeneList.txt', sep='\t', row.names=F, quote=F)

###In excel remove duplicate ENSG and sort ExpressedGeneList Based on Chromosome and then start##
gene.bed.ordered= read.table('ExpressedGeneList.txt', header=T, sep='\t')
masterlist= read.table('ENSG_CAGETSS.txt', header=T, sep='\t')
CAGE.bed = merge(gene.bed.ordered, masterlist, by.x = "GeneID", by.y = "ENSG", all.x = T, all.y = F, sort=F)
which(is.na(CAGE.bed[,7]))
nas = c(9203:11237)
CAGE.bed.DConly = CAGE.bed[-nas,]
which(is.na(CAGE.bed.DConly[,7]))
write.table(CAGE.bed.DConly, 'Hutt.CAGE.Overlap.txt', sep='\t', quote=F,row.names=F)

setwd("C:/Users/Courtney/Dropbox/LCL-iPSC/GEMMA eQTLs/Hutt iPSCs")
##To directly filter mastercolum file
ENSG.ordered = as.matrix(CAGE.bed$GeneID)
mastercol= read.table('hutt.imputed.1Mb.mastercols.txt', header=F, sep='\t')
colnames(mastercol) = c("ENSG", "rs", "Gene", "Pos", "Chr")
mastercol.CK = merge(mastercol, ENSG.ordered, by.x = "ENSG", by.y = "V1", all.x = F, all.y = F, sort=F)

### Re-order Expression and filter out only DC expressed
exprnames = as.matrix(row.names(expr_gene))
idcoefs = read.table("/mnt/lustre/home/cusanovich/500HT/addSNP.coef.3671.square")
## In excel remove the duplicate genes
matcher = read.table("Hutt.CAGE.Overlap.txt", header =T, sep='\t')
ordered.ENSG = matcher$GeneID
matcherids = matcher$GeneID
matcherind = match(matcherids,exprnames[,1])
exprs.o = expr_gene[matcherind,]
exprs.o.t = t(exprs.o)

## Create PCA file based on re-ordered expression with only Darren's expressed genes
exprs = read.table("Expr.DConly.Ordered.txt",header=F,sep="\t")
htpca = prcomp(exprs,scale.=TRUE)
xhtpca = htpca$x
write.table(xhtpca,"hutt.DConly.Ordered.pcs.txt",col.names=F,row.names=F,quote=F,sep="\t")


pcs.o = pcs[matcherind,]
idcoefs.o = idcoefs[matcherind,matcherind]
write.table(ordered.ENSG,"ENSGList.DConly.Ordered.txt",row.names=F,col.names=F,quote=F,sep="\t")
write.table(exprs.o.t,"Expr.DConly.Ordered.txt",row.names=T,col.names=F,quote=F,sep="\t")
##In excel sort by indiv and remove row names(individuals)
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


##Pulling out overlapping findivs
ck = read.table("hutt.imputed.73subset.fam")[,2]
dc = read.table("hutt.imputed.500ht.fam")[,2]
overlap = match(ck,dc)
sum(is.na(overlap))
#18 individual not in Darren's study
fillers = sample(1:431, 18)
##Generated until no overlap was left
match(overlap,fillers)
ck.only = na.omit(overlap)
ck_fillers = append(ck.only,fillers)
newck = sort(ck_fillers, decreasing=FALSE)
dc.fam = read.table("hutt.imputed.500ht.fam")
new.ck.fam = dc.fam[newck,]
write.table(new.ck.fam, 'hutt.imputed.newck.fam',row.name=F,col.names=F,quote=F)

orderlist = new.ck.fam$V2
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
dc.exprs = read.table("qqnorm.500ht.gccor.newcovcor.bimbam.gz")
exprs = dc.exprs[newck,newck.genes]

dim(exprs)

write.table(exprs,"qqnorm.newck.gccor.newcovcor.bimbam.gz",col.names=F,row.names=F,quote=F,sep="\t")

htpca = prcomp(exprs,scale.=TRUE)
xhtpca = htpca$x
write.table(xhtpca,"qqnorm.newck.gccor.newcovcor.pcs",col.names=F,row.names=F,quote=F,sep="\t")

###BROKEN
### Make Venn of overlapping genes between LCLs and Hutt iPSC #######
genes <- data.frame(masterlist$ENSG)
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
