PC1_13.residual = matrix(nrow= nrow(expr_gene), ncol = ncol(expr_gene))
rownames(PC1_13.residual) = rownames(expr_gene)
colnames(PC1_13.residual) = colnames(expr_gene)
for (i in 1:nrow(expr_gene)) {
  model= lm(expr_gene[i,]~ x.pca$x[,1:13])
  PC1_13.residual[i,] = resid(model)
}

#Loop to create data table
lmPCA_heatmap = function(pca, covars, npcs)
{
  results<-c()
  for (f in covars) {
    for (i in 1:npcs)
    {
      s = summary(lm(pca$x[,i]~f));
      results<-c(results,pf(s$fstatistic[[1]],
                            s$fstatistic[[2]],s$fstatistic[[3]], lower.tail = FALSE))
    }
  }
  resultsM<-matrix(nrow = length(covars), ncol = npcs, data =
                     results, byrow = TRUE)
  resultsM
  
}
pcaresults.heatmap = lmPCA_heatmap(x.pca,covars,73)
rownames(pcaresults.heatmap) = c("batch", "sex","age","pluri","novel","matrigel bottle", "passage", "date banked")
myCol <- c("red", "orange", "yellow", "grey")
# Defining breaks for the color scale
myBreaks <- c(0, 1e-5, 1e-2, 0.05, 1)
library(gplots)
heatmap.2(pcaresults.heatmap, col=myCol, breaks=myBreaks, margins=c(5,7),trace="none",main="PCA No Covariates Regressed", key=F, keysize=1.5,density.info="none",cexCol=0.7, cexRow = 0.8, Rowv = FALSE, Colv = FALSE, scale="none")
legend("topleft", fill = myCol, cex=0.8,
       legend = c("<1e-5", "1e-5 to 0.01", "0.01 to 0.05", ">0.05"))

## After regressing out PCs 1-13
x.pca.PC1_13 = prcomp((t(PC1_13.residual)), center = T,scale=T)
x.pca.sum.PC1_13 = summary(x.pca.PC1_13)
pcaresults.heatmap.PC1_13 = lmPCA_heatmap(x.pca.PC1_13,covars,73)
rownames(pcaresults.heatmap.PC1_13) = c("batch", "sex","age","pluri","novel","matrigel bottle", "passage", "date banked")
heatmap.2(pcaresults.heatmap.PC1_13, col=myCol, breaks=myBreaks, margins=c(5,7),trace="none",main="PC 1-13 Regressed", key=F, keysize=1.5,density.info="none",cexCol=0.7, cexRow = 0.8, Rowv = FALSE, Colv = FALSE, scale="none")
legend("topleft", fill = myCol, cex=0.8,
       legend = c("<1e-5", "1e-5 to 0.01", "0.01 to 0.05", ">0.05"))
