setwd("C:/Users/Courtney/Desktop/eQTLs")
ck = read.table('master_PC13_eQTLResults.chosen.txt', header=T)
pvalmenp = read.table('master_PC13_eQTLResults.all.txt', header=T)
dc = read.table('master.PC62.imputed.1Mb.bonferroni.gccor.newcovcor.regressPCs.gemma.eqtls.txt', header=F)
ck.qc = p.adjust(ck$V4,method="BH")
ck = cbind(ck,ck.qc)
e = -log10(ppoints(length(ck$ck.qc)))
o = -log10(ck$ck.qc)
qqplot(e, o, main = "QQ Plot",lwd=1, xlab="Expected (-logP)", ylab="Observed (-logP)")
abline(0,1, col ='red')

qc = p.adjust(dc$V3,method="BH")
dc = cbind(dc,qc)
e = -log10(ppoints(length(dc$qc)))
o = -log10(dc$qc)
qqplot(e, o, main = "QQ Plot")
abline(0,1)
dev.off()



#QQ Plot
#pdf = ('/mnt/gluster/data/internal_supp/hutt_ipsc/AllGenes/eQTLs/QQplot.pdf')
pdf = ('QQplot_CK.all.pdf')
e = -log10(ppoints(length(pvalmenp$qs)))
o = -log10(pvalmenp$qs)
qqplot(e, o, main = "QQ Plot")
abline(0,1)
dev.off()



pvals = pvalmenp$qs
observed <- sort(pvals)
lobs <- -(log10(observed))

expected <- c(1:length(observed)) 
lexp <- -(log10(expected / (length(expected)+1)))

# pdf("QQplot_FDRLine_FDR5DE.pdf", width=6, height=6)
plot(c(0,7), c(0,7), col="red", lwd=1, type="l", main = "QQ Plot",xlab="Expected (-logP)", ylab="Observed (-logP)", xlim=c(0,7), ylim=c(0,7), las=1, xaxs="i", yaxs="i", bty="l")
points(lexp, lobs, pch=19, cex=.6, col="black") 