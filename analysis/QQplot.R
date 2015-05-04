#QQ Plot
e = -log10(ppoints(length(pvalmenp)))
o = -log10(pvalmenp)
qqplot(e, o, main = "QQ Plot with final expression")
abline(0,1)