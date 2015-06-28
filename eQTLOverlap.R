setwd("C:/Users/Courtney/Dropbox/LCL-iPSC/GEMMA eQTLs/GEMMATest/Enrichment")
liver = read.table('Adipose.eqtl.txt', header=T)
art = read.table('Artery_Aorta.eqtl.txt', header=T)
tib = read.table('Artery_Tibial.eqtl.txt', header=T)
ipsc = read.table('PC13_eQTLResults.corrected.chosen.txt', header=T)
CAGETSS = read.table('ENSG_CAGETSS.txt', header=T)
eso = read.table('Esophagus_Mucosa.eqtl.txt', header=T)
eso2 = read.table('Esophagus_Muscularis.eqtl.txt', header=T)
heart = read.table('Heart_Left_Ventricle.eqtl.txt', header=T)
lung = read.table('Lung.eqtl.txt', header=T)
musc = read.table('Muscle_Skeletal.eqtl.txt', header=T)
nerve = read.table('Nerve_Tibial.eqtl.txt', header=T)
skin = read.table('Skin_Sun_Exposed_Lower_leg.eqtl.txt', header=T)
stom = read.table('Stomach.eqtl.txt', header=T)
thy = read.table('Thyroid.eqtl.txt', header=T)
wb = read.table('Whole_Blood.eqtl.txt', header=T)
ipsc2 = as.data.frame(ipsc$V1)
colnames(ipsc2) = c("Gen_ID")
pluri = read.table('Plurigenes.txt', header=T)

ipsc.genes = read.table('ENSGList.allgenes.Ordered.txt')
ipsc.genes.tested = read.table('iPSCGenesTested.txt')
LCL.genes.tested = read.table('LCLGenesTested.txt')
LCL.genes = read.table('LCLExpressedGenes.txt')
colnames(ipsc.genes) = c("Gen_ID")

all=c()
all = c(as.character(CAGETSS$ENSG),as.character(liver$Gen_ID), as.character(ipsc2$Gen_ID),as.character(art$Gen_ID),as.character(tib$Gen_ID))
all = unique(all)

univ <- data.frame(all)
library(VennDiagram)
names(univ) <- "probes"

make.venn.quad <- function(geneset1, geneset2, geneset3, geneset4, geneset1.label, geneset2.label, geneset3.label, geneset4.label, univ){
  univ$g1 <- univ$probes %in% geneset1
  univ$g2 <- univ$probes %in% geneset2
  univ$g3 <- univ$probes %in% geneset3 
  univ$g4 <- univ$probes %in% geneset4 
  #pdf(file=paste(prefix, ".pdf", sep=""), width=7, height=7)
  venn.placeholder <- draw.quad.venn(length(geneset1),length(geneset2), length(geneset3), length(geneset4), dim(univ[univ$g1 == T & univ$g2 == T , ])[1], dim(univ[univ$g1 == T & univ$g3 == T , ])[1], dim(univ[univ$g1 == T & univ$g4 == T , ])[1], dim(univ[univ$g2 == T & univ$g3 == T , ])[1], dim(univ[univ$g2 == T & univ$g4 == T , ])[1], dim(univ[univ$g3 == T & univ$g4 == T , ])[1], dim(univ[univ$g1 == T & univ$g2 == T & univ$g3 == T , ])[1], dim(univ[univ$g1 == T & univ$g2 == T & univ$g4 == T , ])[1], dim(univ[univ$g1 == T & univ$g3 == T & univ$g4 == T , ])[1], dim(univ[univ$g2 == T & univ$g3 == T & univ$g4 == T , ])[1],  dim(univ[univ$g1 == T & univ$g2 == T & univ$g3 == T & univ$g4 == T , ])[1], c(geneset1.label, geneset2.label, geneset3.label, geneset4.label), fill=c("goldenrod", "plum4", "steelblue3", "darkolivegreen3"), alpha=c(0.5, 0.5, 0.5, 0.5),col=NA, euler.d=T)
  complement.size <- dim(univ[univ$g1 == F & univ$g2 == F & univ$g3 == F & univ$g4 == F , ])[1]
  grid.text(paste(complement.size, " not DE in any", sep=""), x=0.2, y=0.08)
  #dev.off()
}

#dev.off()
# Make venn of full
pdf(file = "VennDiagram_DE_FDR5.pdf")
make.venn.quad(as.character(liver$Gen_ID), as.character(ipsc2$Gen_ID),as.character(art$Gen_ID),as.character(tib$Gen_ID),"liver","ipsc","art","tib" ,univ)
dev.off()

res=matrix(data=NA, nrow=15418, ncol=13)
labs = c("liver", "art", "tib", "eso","eso2","heart","lung","musc","nerve","skin","stom","wb","ipsc2")
tissues = c(liver, art, tib, eso,eso2,heart,lung,musc,nerve,skin,stom,wb,ipsc2)
  for (i in 1:length(tissues)){
    m = tissues[i]
    test = univ$probes %in% m$Gen_ID
    res[,i] = test
  }
res = as.data.frame(res)
dim(univ[res$V1 == T & res$V2 == T & res$V3 == T & res$V4 == T & res$V5 == T & res$V6 == T & res$V7 == T & res$V8 == T & res$V9 == T & res$V10 == T & res$V11 == T & res$V12 == T & res$V13 == T, ])[1]



labs = c("liver", "art", "tib", "eso","eso2","heart","lung","musc","nerve","skin","stom","wb","ipsc2")
tissues = c(liver, art, tib, eso,eso2,heart,lung,musc,nerve,skin,stom,wb,ipsc2)
for (k in 1:length(tissues)){
  table= matrix(data=NA, nrow=13, ncol=3)
  props=c()
  for (i in 1:length(tissues)){
    j = tissues[i]
    m = tissues[k]
    test = m$Gen_ID %in% j$Gen_ID
    tr = sum(test == T)
    total = length(j$Gen_ID)
    prop = tr/total
    props[i]=prop
    table[i,1]=total
    table[i,2] =tr
  }
  plot(1:12, props[1:12], main = paste("Tissue Overlap with", labs[k]))
}


#Finding how many from their data is tested in mine and getting a proportion only for that
props.ipsconly=c()
for (i in 1:12){
  j = tissues[i]
  m = ipsc.genes
  subs = match(j$Gen_ID,m$Gen_ID)
  subs2 = na.omit(subs)
  subset = as.data.frame(ipsc.genes[subs2,])
  colnames(subset) = c("Gen_ID")
  test = ipsc2$Gen_ID %in% subset$Gen_ID
  tr = sum(test == T)
  total = length(subset$Gen_ID)
  table[i,2] = total
  table[i,3] = tr
  prop = tr/total
  props.ipsconly[i]=prop
}
counts = table(props.ipsconly)
plot(1:12, props.ipsconly, main = "Tissue Overlap with iPSCs - Overlapping genes only")


subset.size= matrix(data=NA, nrow=12, ncol=2)
for (i in 1:12){
  j = tissues[i]
  m = ipsc.genes
  subs = match(j$Gen_ID,m$Gen_ID)
  subs2 = na.omit(subs)
  subset = as.data.frame(ipsc.genes[subs2,])
  colnames(subset) = c("Gen_ID")
  assign(paste("subset.",labs[i],sep=""),subset)
  total = length(subset$Gen_ID)
  subset.size[i,1] = labs[i]
  subset.size[i,2] = total
}

tissues.subset = c(subset.liver, subset.art, subset.tib, subset.eso,subset.eso2,subset.heart,subset.lung,subset.musc,subset.nerve,subset.skin,subset.stom,subset.wb,ipsc2)
##1000 perms per tissue type
perm.output= matrix(data=NA, nrow=12, ncol=100000)
for(i in 1:100000) {
perm = sample(1:10250, 1028, replace=F)
perm.data = as.data.frame(ipsc.genes[perm,])
colnames(perm.data)=c("Gen_ID")
for (k in 1:12){
  j = tissues.subset[k] 
  m = perm.data
  test = m$Gen_ID %in% j$Gen_ID
  tr = sum(test == T)
  perm.output[k,i] = tr
}
}

for(i in 1:12){
print(sum(perm.output[i,] > table[i,3]))
}


ipsc = read.table('iPSC.PC13.gemma.chosen.txt', header=F)
LCL = read.table('LCL.PC8.gemma.chosen.txt', header=F)

##Take only 5% BH filtered eQTLs from the Bonf corrected
ipsc.qs = p.adjust(ipsc$V4,method="BH")
ipsc.cor = cbind(ipsc,ipsc.qs)
ipsc.eqtls = ipsc.cor[ipsc.cor$ipsc.qs <.05,]

LCL.qs = p.adjust(LCL$V4,method="BH")
LCL.cor = cbind(LCL,LCL.qs)
LCL.eqtls = LCL.cor[LCL.cor$LCL.qs <.05,]

##Find the shared eQTLs and those only for each cell type
shared = which(ipsc.eqtls$V1 %in% LCL.eqtls$V1)
shared.eqtls = ipsc.eqtls[shared,]
#23 genes with shared eQTLs

ipsc.only.eqtls = ipsc.eqtls[-shared,]
#1005 iPSC specific
shared.2 = which(LCL.eqtls$V1 %in% ipsc.eqtls$V1)
LCL.only.eqtls = LCL.eqtls[-shared.2,]
#128 LCL specific

tested.both = which(ipsc.only.eqtls$V1 %in% LCL.genes.tested$V1)
ipsc.only.tested.eqtls = ipsc.only.eqtls[tested.both,]
#920 iPSC speific and tested in LCLs
expr.both = which(ipsc.only.tested.eqtls$V1 %in% LCL.genes$V1)

ipsc.only.tested.expressed.eqtls = ipsc.only.tested.eqtls[expr.both,]
#920 eQTLs specific to iPSCs and expressed and tested in LCLs

tested.both = which(LCL.only.eqtls$V1 %in% ipsc.genes.tested$V1)
LCL.only.tested.eqtls = LCL.only.eqtls[tested.both,]
#94 LCL speific and tested in iPSCs
expr.both = which(LCL.only.tested.eqtls$V1 %in% ipsc.genes$V1)
#61 LCL specific and expressed in iPSCs

LCL.only.tested.expressed.eqtls = LCL.only.tested.eqtls[expr.both,]

##P-value of overlap between LCL eQTLs and iPSC eQTLs
perm.output= matrix(data=NA, nrow=100000, ncol=1)
for(i in 1:100000) {
  temp = sample(1:13483, 151, replace=F)
  perm.data = as.data.frame(LCL[temp,])
  test = ipsc.eqtls$V1 %in% perm.data$V1
    tr = sum(test == T)
    perm.output[i,1] = tr
  }
sum(perm.output[,1] > 22)
#29 so pval is 0.00029

##Do permutation for overall overlap
#First subset only genes that were tested and expressed in both tissues
tested.both = which(LCL$V1 %in% ipsc$V1)
LCL.both = LCL[tested.both,]
tested.both = which(ipsc$V1 %in% LCL$V1)
ipsc.both = ipsc[tested.both,]
#8,887 shared genes

perm.both.output= matrix(data=NA, nrow=100000, ncol=1)
for(i in 1:100000) {
  temp = sample(1:8887, 84, replace=F)
  perm.data = as.data.frame(LCL.both[temp,])
  temp.i = sample(1:8887,943, replace=F)
  perm.data.i= as.data.frame(ipsc.both[temp.i,])
  test = perm.data.i$V1 %in% perm.data$V1
  tr = sum(test == T)
  perm.both.output[i,1] = tr
}
sum(perm.both.output[,1] > 22)
#3 so pval 3e-05

##Look for GTEx overlap
gtex = read.table('GTExList.txt', header=F)
ipsc.eqtls.f = ipsc.only.tested.expressed.eqtls

colnames(ipsc.eqtls.f) = c("ENSG", "Variant", "p", "Bonf", "Qval")
colnames(gtex)= c("Gene", "SNP")
ipsc.gtex = merge(ipsc.eqtls.f, gtex, by.x = c("ENSG", "Variant"), by.y = c("Gene", "SNP"))
#7 iPSC eQTLs overlap with GTEx list

LCL.eqtls.f = LCL.only.tested.expressed.eqtls
colnames(LCL.eqtls.f) = c("ENSG", "Variant", "p", "Bonf", "Qval")
LCL.gtex = merge(LCL.eqtls.f, gtex, by.x = c("ENSG", "Variant"), by.y = c("Gene", "SNP"))
#0 LCL eQTLs overlap with GTEx list

colnames(shared.eqtls) = c("ENSG", "Variant", "p", "Bonf", "Qval")
shared.gtex = merge(shared.eqtls, gtex, by.x = c("ENSG", "Variant"), by.y = c("Gene", "SNP"))
#0 shared eQTLs overlap with GTEx list