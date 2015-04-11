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
