## Pluritest Code for generating files ###

setwd("~/Arrays/Hutt iPSCs")
load("C:/Users/Courtney/Dropbox/LCL-iPSC/Pluritestsub_REnvironment.unk")
sample.names = 'Sample_names.txt'
arraydata = 'YGilad-CK-Mar6-15-ProbeLevelData-NotNormalized-NoBGSubtracted-FinalReport.txt'
# load packages
require(lumi)
require(xtable)
require(GO.db)
# set the stage
samplenames <- readLines(file(sample.names))
#IMPORT RAW DATA WITH LUMI
working.lumi<-lumiR(arraydata, convertNuID = FALSE, annotationColumn="PROBE_ID")
sampleNames(working.lumi) <- samplenames
#identify probe ID column
probeind <-match("PROBE_ID", colnames(fData(working.lumi)))
fData(working.lumi)[,probeind]<-gsub("\"","",fData(working.lumi)[,probeind])
#RSN NORMALIZATION OF THE DATA
A <- fData(working.lumi)[,probeind]  #matches on ILMN_Ids for lumi/RSN
B <- fData(H9targetArray)[,1] #for matches on ILMN_Ids for lumi/RSN
sel.match <- match(B,A)
working.lumi <- working.lumi[sel.match[!is.na(sel.match)],] #subsets the exprSet user.lumi down to subsetuser.lumi to match H9 ref array rows
working.lumi<-lumiN(working.lumi,method="rsn",target=H9targetArray[is.na(sel.match)==FALSE,])
#RSN NORMALIZATION OF THE DATA CONTINUED
# assume that illumina Probe Ids are in fData[,1]
A <- fData(working.lumi)[,probeind]
sel.match<-match(rownames(W15),A)
working.lumi.exprs <- working.lumi@assayData$exprs
rownames(working.lumi.exprs) <- A
#CALCULATION OF SCORES
try(
{
  sel<-match(rownames(W15),A)
  coef<-c(-1.267095e+02  ,4.567437e-03 , 4.377068e-03  ,1.043193e-03)
  working.lumi.int <- working.lumi.exprs[sel[!is.na(sel)],]
  H15.new<-predictH(working.lumi.int, W15[!is.na(sel),])
  H12.new<-predictH(working.lumi.int, W12[!is.na(sel),])
  rss.new<-apply((working.lumi.int - W12[!is.na(sel),]%*%H12.new)^2,2,sum)
  RMSE.new<-sqrt(rss.new/sum(!is.na(sel)))
  novel.new<-apply((working.lumi.int - W12[!is.na(sel),]%*%H12.new)^8,2,sum)
  novel.new<-(novel.new/sum(!is.na(sel)))^(1/8)
  s.new<-drop(coef[1] +coef[2:4]%*%H15.new[c(1,14,13),])
  print(s.new)
}
)
# plot MULTICLASS PLURITEST & overview
table.results<-matrix(,nrow=ncol(exprs(working.lumi)),ncol=5)
rownames(table.results)<-colnames(exprs(working.lumi))
colnames(table.results)<-c("pluri-raw","pluri logit-p","novelty","novelty logit-p","RMSD")

#Populate with data
table.results[,1]<-round(s.new,3)
table.results[,2]<-round(exp(s.new)/(1+exp(s.new)),3)
table.results[,3]<-round(novel.new,3)
table.results[,5]<-round(RMSE.new,3)

pdf("PluritestPlot.pdf")
par(mar=c(5,4,4,2))
par(xaxt='n')
plot(s.new,main="Pluritest",xlab="",ylab="Pluripotency Score",ylim=c(-130,70), cex = .8, pch = 16)
abline(h=20,lty="dashed",col="red")
abline(h=-28.92,lty="dashed",col="blue")
par(xaxt='s')
axis(1,at=c(1:length(s.new)),labels=names(s.new),las=2, cex.axis = .5)
dev.off()

pdf("PluritestHeatmap.pdf")
plot(s.new~novel.new,cex=.8,main="Overview", pch=16, ylab = "Pluripotency Score", xlab = "Novelty")
abline(v=1.67,col="black")
abline(h=20,col="black")
dev.off()

pdf("Novelty.pdf")
par(mar=c(5,4,4,2))
par(xaxt='n')
barplot(novel.new,main = "Novelty Score",names.arg=c(1:length(novel.new)),xlab="",xlim=c(0,length(novel.new)),width=.9,space=(1/9),ylim=c(0,4), ylab="Novelty")
par(xaxt='s')
axis(1,at=c(1:nrow(  table.results))-.4,labels=names(s.new),las=2, cex.axis=.5)
abline(h=1.67,col="black", lty="dashed")
dev.off()

# Save CSV FILE for TABLE
table.results[,5]<-round(RMSE.new,3)
write.csv(table.results,file="PluritestResults.csv")