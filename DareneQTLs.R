setwd("C:/Users/Courtney/Desktop/eQTLs")
eqtls= read.table('master.PC62.imputed.1Mb.bonferroni.gccor.newcovcor.regressPCs.gemma.eqtls.txt')
chosen = read.table('master.PC62.imputed.1Mb.bonferroni.gccor.newcovcor.regressPCs.gemma.chosen.txt')

overlap = match(ipsc.genes$Gen_ID,eqtls$V1)
over = unique(overlap)
sub = na.omit(overlap)
subset = eqtls[sub,]
fdr = p.adjust(subset$V3,method="BH")
dcmatrix.all = cbind(subset,fdr)

all =c()
genelist = as.data.frame(ipsc.genes$Gen_ID)
for(i in 1:dim(ipsc.genes)[1]){
  gene.temp= as.character(ipsc.genes[i,])
  rows.temp = which(as.character(eqtls$V1) == gene.temp)
  gene.only = eqtls[rows.temp,]
  ##Create an all file that containes the fdr correction (plus gene/SNP and origina pval)
  all = rbind(all,gene.only)
}
qs.all = p.adjust(all$V3,method="BH")
new.dc = cbind(all, qs.all)
colnames(new.dc) = c("ENSG", "Variant", "Unc", "FDR")
write.table(new.dc, "Darren_eQTLs_all.corrected.txt", sep='\t', row.names=F, quote=F)

#Just working with the bonf corrected
boverlap = match(chosen$V1, ipsc.genes$Gen_ID)
bover = unique(boverlap)
bsub = na.omit(bover)
bsubset = chosen[bsub,]

qs = p.adjust(bsubset$V4,method="BH")
dcmatrix = cbind(bsubset,qs)
qcount = sum(dcmatrix$qs < 0.05)
sig.dc = dcmatrix[dcmatrix$qs < 0.05,]
test = ipsc$V1 %in% sig.dc$V1
sum(test ==T)
test2 = sig.dc$V1 %in% ipsc$V1 
sum(test2==T)

shared = match(sig.dc$V1, ipsc$V1)
shared.g = na.omit(shared)
sharedeqtls = sig.dc[shared.g,]
write.table(sharedeqtls,"SharedeQTLs.txt", sep='\t', quote=F, row.names=F)
all.eqtls = append(as.character(sig.dc$V1), as.character(ipsc$V1))
write.table(all.eqtls, 'AlleQTLs.txt', quote=F, row.names =F, sep='\t')

make.venn.pair <- function(geneset1, geneset2, geneset1.label,
                           geneset2.label,universe){
 # pdf(file=paste(prefix, ".pdf", sep=""), width=7, height=7)
  venn.placeholder <-
    draw.pairwise.venn(length(geneset1),length(geneset2),
                       length(which(geneset1 %in% geneset2) == TRUE), c(paste(geneset1.label,
                                                                              length(geneset1), sep="\n"), paste(geneset2.label, length(geneset2),
                                                                                                                 sep="\n")), fill=c("goldenrod", "plum4"), alpha=c(0.5, 0.5),col=NA,
                       euler.d=T)
  complement.size <- dim(universe)[1] - length(geneset1) -
    length(geneset2) + length(which(geneset1 %in% geneset2) == TRUE)
  grid.text(paste(complement.size, " not eQTLs", sep=""),
            x=0.1, y=0.1)
  #dev.off()
  print(paste("Probes in a: ", length(geneset1), sep=""))
  print(paste("Probes in b: ", length(geneset2), sep=""))
  print(paste("Common probes: ", length(which(geneset1 %in% geneset2)
                                        == TRUE), sep=""))
}

make.venn.pair(dcmatrix[dcmatrix$qs < 0.05,]$V1,
               ipsc$V1,  "LCL eQTLs",
               "iPSC eQTLs", ipsc.genes)

# make.venn.pair(pluri$PluripotencyGenes,
#                ipsc.genes$Gen_ID,  "Pluripotency Genes",
#                "iPSC genes", ipsc.genes)

