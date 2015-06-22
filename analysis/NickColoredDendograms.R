# Get methylation data from 1K most variable probes
setwd("C:/Users/Courtney/Google Drive/Origin Project/Plots from Courtney")
meth.final= read.table('1K_most_variable_iPSC_probes_2.txt', as.is=T, header=F)
inames = c("Ind1 F-iPSC", "Ind1 L-iPSC A", "Ind1 L-iPSC B", "Ind1 L-iPSC C", "Ind2 F-iPSC", "Ind2 L-iPSC A", "Ind2 L-iPSC B", "Ind2 L-iPSC C", "Ind3 F-iPSC", "Ind3 L-iPSC A", "Ind3 L-iPSC B", "Ind3 L-iPSC C", "Ind4 F-iPSC", "Ind4 L-iPSC A", "Ind4 L-iPSC B", "Ind4 L-iPSC C")
colnames(meth.final)=inames
# Make dendrogram

library("ClassDiscovery", lib.loc="~/R/win-library/3.0")

setwd("C:/Users/Courtney/Google Drive/Origin Project/Plots from Courtney")
cols = c("Navy", "Navy", "Navy", "Navy", "darkorange1", "darkorange1", "darkorange1", "darkorange1", "Black", "Black", "Black", "Black", "lightcoral", "lightcoral", "lightcoral", "lightcoral")

pdf(file = "Dendrogram_1K_most_variable_methprobes_withsex.pdf")
euc = hclust(dist(t(meth.final)))
plotColoredClusters(euc, colnames(meth.final), cols, cex = 1,lwd = 3, lty = 1,main = "Complete linkage clustering dendogram of 1,000 most variable DNA methylation loci", line = -1, xlab="", sub="")
dev.off()
