Hutterite iPSC Panel
========================================================

This is the workflow used for QC, pluritest, and initial analysis on the 73 Hutterite iPSC samples and 11 heart samples (dox experiment). To generate final covariates including Pluri score and novelty we start with running pluritest. Start with script Clean_Pluritest_iPSCArray.R


```r
setwd("~/Arrays/Hutt iPSCs")
load("C:/Users/Courtney/Dropbox/LCL-iPSC/Pluritestsub_REnvironment.unk")
sample.names = 'Sample_names.txt'
arraydata = 'YGilad-CK-Mar6-15-ProbeLevelData-NotNormalized-NoBGSubtracted-FinalReport.txt'
# load packages
require(lumi)
```

```
## Warning: replacing previous import 'image' when loading 'graphics'
```

```
## Warning: replacing previous import 'nleqslv' when loading 'nleqslv'
```

```r
require(xtable)
```

```
## Warning: package 'xtable' was built under R version 3.0.3
```

```r
require(GO.db)
```

```
## Warning: package 'DBI' was built under R version 3.0.3
```

```r
# set the stage
samplenames <- readLines(file(sample.names))
#IMPORT RAW DATA WITH LUMI
working.lumi<-lumiR(arraydata, convertNuID = FALSE, annotationColumn="PROBE_ID")
```

```
## Warning: closing unused connection 5 (Sample_names.txt)
```

```
## Perform Quality Control assessment of the LumiBatch object ...
```

```r
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
```

```
## Perform rsn normalization ...
## 2015-03-11 16:30:15 , processing array  1 
## 2015-03-11 16:30:15 , processing array  2 
## 2015-03-11 16:30:15 , processing array  3 
## 2015-03-11 16:30:15 , processing array  4 
## 2015-03-11 16:30:16 , processing array  5 
## 2015-03-11 16:30:16 , processing array  6 
## 2015-03-11 16:30:16 , processing array  7 
## 2015-03-11 16:30:16 , processing array  8 
## 2015-03-11 16:30:16 , processing array  9 
## 2015-03-11 16:30:16 , processing array  10 
## 2015-03-11 16:30:17 , processing array  11 
## 2015-03-11 16:30:17 , processing array  12 
## 2015-03-11 16:30:17 , processing array  13 
## 2015-03-11 16:30:17 , processing array  14 
## 2015-03-11 16:30:17 , processing array  15 
## 2015-03-11 16:30:17 , processing array  16 
## 2015-03-11 16:30:17 , processing array  17 
## 2015-03-11 16:30:17 , processing array  18 
## 2015-03-11 16:30:18 , processing array  19 
## 2015-03-11 16:30:18 , processing array  20 
## 2015-03-11 16:30:18 , processing array  21 
## 2015-03-11 16:30:18 , processing array  22 
## 2015-03-11 16:30:18 , processing array  23 
## 2015-03-11 16:30:18 , processing array  24 
## 2015-03-11 16:30:19 , processing array  25 
## 2015-03-11 16:30:19 , processing array  26 
## 2015-03-11 16:30:19 , processing array  27 
## 2015-03-11 16:30:19 , processing array  28 
## 2015-03-11 16:30:19 , processing array  29 
## 2015-03-11 16:30:19 , processing array  30 
## 2015-03-11 16:30:19 , processing array  31 
## 2015-03-11 16:30:20 , processing array  32 
## 2015-03-11 16:30:20 , processing array  33 
## 2015-03-11 16:30:20 , processing array  34 
## 2015-03-11 16:30:20 , processing array  35 
## 2015-03-11 16:30:20 , processing array  36 
## 2015-03-11 16:30:20 , processing array  37 
## 2015-03-11 16:30:21 , processing array  38 
## 2015-03-11 16:30:21 , processing array  39 
## 2015-03-11 16:30:21 , processing array  40 
## 2015-03-11 16:30:21 , processing array  41 
## 2015-03-11 16:30:21 , processing array  42 
## 2015-03-11 16:30:21 , processing array  43 
## 2015-03-11 16:30:21 , processing array  44 
## 2015-03-11 16:30:22 , processing array  45 
## 2015-03-11 16:30:22 , processing array  46 
## 2015-03-11 16:30:22 , processing array  47 
## 2015-03-11 16:30:22 , processing array  48 
## 2015-03-11 16:30:22 , processing array  49 
## 2015-03-11 16:30:22 , processing array  50 
## 2015-03-11 16:30:22 , processing array  51 
## 2015-03-11 16:30:23 , processing array  52 
## 2015-03-11 16:30:23 , processing array  53 
## 2015-03-11 16:30:23 , processing array  54 
## 2015-03-11 16:30:23 , processing array  55 
## 2015-03-11 16:30:23 , processing array  56 
## 2015-03-11 16:30:23 , processing array  57 
## 2015-03-11 16:30:24 , processing array  58 
## 2015-03-11 16:30:24 , processing array  59 
## 2015-03-11 16:30:24 , processing array  60 
## 2015-03-11 16:30:24 , processing array  61 
## 2015-03-11 16:30:24 , processing array  62 
## 2015-03-11 16:30:24 , processing array  63 
## 2015-03-11 16:30:24 , processing array  64 
## 2015-03-11 16:30:25 , processing array  65 
## 2015-03-11 16:30:25 , processing array  66 
## 2015-03-11 16:30:25 , processing array  67 
## 2015-03-11 16:30:25 , processing array  68 
## 2015-03-11 16:30:25 , processing array  69 
## 2015-03-11 16:30:25 , processing array  70 
## 2015-03-11 16:30:25 , processing array  71 
## 2015-03-11 16:30:26 , processing array  72 
## 2015-03-11 16:30:26 , processing array  73 
## 2015-03-11 16:30:26 , processing array  74 
## 2015-03-11 16:30:26 , processing array  75 
## 2015-03-11 16:30:26 , processing array  76 
## 2015-03-11 16:30:26 , processing array  77 
## 2015-03-11 16:30:27 , processing array  78 
## 2015-03-11 16:30:27 , processing array  79 
## 2015-03-11 16:30:27 , processing array  80 
## 2015-03-11 16:30:27 , processing array  81 
## 2015-03-11 16:30:27 , processing array  82 
## 2015-03-11 16:30:27 , processing array  83 
## 2015-03-11 16:30:27 , processing array  84 
## 2015-03-11 16:30:28 , processing array  85
```

```r
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
```

```
##              16              21              25              43 
##        24.26793        25.75893        27.55310        29.51692 
##              46              61              62              67 
##        27.15599        25.40323        26.82581        26.75365 
##              73              69              22              23 
##        24.30677        23.78727        21.59470        22.17941 
##              20              18              33              36 
##        26.32808        26.18923        27.81019        23.98461 
##              44              45               5              54 
##        26.38043        24.15399        28.00584        24.35104 
##              64              10              12               8 
##        28.57498        26.64190        28.58306        27.98423 
##              19              27              41              38 
##        23.86254        26.82777        28.09090        31.16365 
##              55              52              57              63 
##        27.08759        22.45290        27.96159        23.00697 
##              59              74              70              26 
##        20.48502        27.24430        25.03894        20.90056 
##              30              40              37              50 
##        27.67938        22.54953        30.31327        29.96805 
##               1               4              58              65 
##        24.14134        24.80727        27.71237        25.75511 
##               7              68              11              15 
##        22.57494        22.82942        23.18641        27.66332 
##              35              32              39              48 
##        33.30140        34.77234        33.37310        35.08889 
##              47              72              24  HT14- 48h 1.25 
##        32.44501        30.90049        28.95333       -36.35994 
## HT14- 48H 0.625 HT14- 48H cont.     HT14- 48H 5   HT14- 48H 2.5 
##       -52.86083       -43.63621       -61.07164       -60.52280 
##              17              14               6              56 
##        29.41440        29.52387        27.90030        30.31478 
##              60              66    HT14- 48H 10 HT14- 96H cont. 
##        32.31204        33.27983       -59.97656       -44.67679 
## HT14- 96H 0.625  HT14- 96H 1.25   HT14- 96H 2.5     HT14- 96H 5 
##       -45.81746       -65.97621       -61.62268       -64.34667 
##              29              34              31              42 
##        33.25744        28.91509        30.79135        31.69517 
##              49               3               2              51 
##        28.44490        24.60073        30.06947        32.29281 
##              53              75               9              71 
##        28.51089        24.24579        29.90078        27.05281
```

```r
# plot MULTICLASS PLURITEST & overview
table.results<-matrix(,nrow=ncol(exprs(working.lumi)),ncol=5)
rownames(table.results)<-colnames(exprs(working.lumi))
colnames(table.results)<-c("pluri-raw","pluri logit-p","novelty","novelty logit-p","RMSD")
try(
{
print(s.new)
pdf("Origin pluritest Plot.pdf")
par(mar=c(12,4,4,2))
par(xaxt='n')
plot(s.new,main="pluripotency",xlab="",ylab="pluripotency",ylim=c(-130,70), cex = .6, pch = 20)
abline(h=25.25,lty="dashed",col="red")
abline(h=59.95,lty="dashed",col="red")
abline(h=-28.92,lty="dashed",col="lightblue")
abline(h=-130,lty="dashed",col="lightblue")
par(xaxt='s')
axis(1,at=c(1:length(s.new)),labels=names(s.new),las=2, cex.axis = .5)
dev.off()
}
)
```

```
##              16              21              25              43 
##        24.26793        25.75893        27.55310        29.51692 
##              46              61              62              67 
##        27.15599        25.40323        26.82581        26.75365 
##              73              69              22              23 
##        24.30677        23.78727        21.59470        22.17941 
##              20              18              33              36 
##        26.32808        26.18923        27.81019        23.98461 
##              44              45               5              54 
##        26.38043        24.15399        28.00584        24.35104 
##              64              10              12               8 
##        28.57498        26.64190        28.58306        27.98423 
##              19              27              41              38 
##        23.86254        26.82777        28.09090        31.16365 
##              55              52              57              63 
##        27.08759        22.45290        27.96159        23.00697 
##              59              74              70              26 
##        20.48502        27.24430        25.03894        20.90056 
##              30              40              37              50 
##        27.67938        22.54953        30.31327        29.96805 
##               1               4              58              65 
##        24.14134        24.80727        27.71237        25.75511 
##               7              68              11              15 
##        22.57494        22.82942        23.18641        27.66332 
##              35              32              39              48 
##        33.30140        34.77234        33.37310        35.08889 
##              47              72              24  HT14- 48h 1.25 
##        32.44501        30.90049        28.95333       -36.35994 
## HT14- 48H 0.625 HT14- 48H cont.     HT14- 48H 5   HT14- 48H 2.5 
##       -52.86083       -43.63621       -61.07164       -60.52280 
##              17              14               6              56 
##        29.41440        29.52387        27.90030        30.31478 
##              60              66    HT14- 48H 10 HT14- 96H cont. 
##        32.31204        33.27983       -59.97656       -44.67679 
## HT14- 96H 0.625  HT14- 96H 1.25   HT14- 96H 2.5     HT14- 96H 5 
##       -45.81746       -65.97621       -61.62268       -64.34667 
##              29              34              31              42 
##        33.25744        28.91509        30.79135        31.69517 
##              49               3               2              51 
##        28.44490        24.60073        30.06947        32.29281 
##              53              75               9              71 
##        28.51089        24.24579        29.90078        27.05281
```

```
## png 
##   2
```

```r
table.results[,1]<-round(s.new,3)
table.results[,2]<-round(exp(s.new)/(1+exp(s.new)),3)
table.results[,3]<-round(novel.new,3)
table.results[,5]<-round(RMSE.new,3)
```
Pluritest figure:


```r
par(mar=c(12,4,4,2))
par(xaxt='n')
plot(s.new,main="pluripotency",xlab="",ylab="pluripotency",ylim=c(-130,70), cex = .6, pch = 20)
abline(h=25.25,lty="dashed",col="red")
abline(h=59.95,lty="dashed",col="red")
abline(h=-28.92,lty="dashed",col="lightblue")
abline(h=-130,lty="dashed",col="lightblue")
par(xaxt='s')
axis(1,at=c(1:length(s.new)),labels=names(s.new),las=2, cex.axis = .5)
```

![plot of chunk unnamed-chunk-2](figure/unnamed-chunk-2-1.png) 

Why is table.results[,2]<-round(exp(s.new)/(1+exp(s.new)),3) the code to get the pluriscore 0/1? How is this being calculated? Is there a true cut-off value?
