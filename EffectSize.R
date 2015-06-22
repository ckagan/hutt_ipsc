#Effect size

#Looking at Darren's subseted data from LCLs

#cnames = c("gene", "SNP", "pval", "beta", "SE")
#colnames(data) = cnames
#hist(data$beta)
data = read.table('/mnt/gluster/data/internal_supp/hutt_ipsc/DCAllGenes/eQTLs/master.PC8.imputed.1Mb.bonferroni.gccor.newcovcor.regressPCs.gemma.eqtls.txt', as.is=T)
data2 = read.table('/mnt/gluster/data/internal_supp/hutt_ipsc/AllGenes/eQTLs/master_PC13_eQTLResults.all.txt', as.is=T)
data3 = read.table('/mnt/gluster/data/internal_supp/hutt_ipsc/DCFull/eQTLs/master.PC62.imputed.1Mb.bonferroni.gccor.newcovcor.regressPCs.gemma.eqtls.txt', as.is=T)

pdf('EffectSizesLCL.pdf')
plot( p1, col=rgb(0,0,1,1/4), xlim=c(-2,2))  # first histogram
plot( p2, col=rgb(1,0,0,1/4), xlim=c(-2,2), add=T)  # second
plot( p3, col=rgb(0,1,0,1/4), xlim=c(-2,2), add=T)  # second
dev.off()


#and combine into your new data frame vegLengths
effectsizes <- rbind(LCL,iPSC)

#now make your lovely plot
ggplot(effectsizes, aes(length, fill = V4)) + geom_density(alpha = 0.2)