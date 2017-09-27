

# Building manuals and namespace using roxygen2
library(devtools)
document()


# Remove .Rhistory
unlink(".Rhistory")
unlink("R/.Rhistory")


# # # Force the formatR style
# library(formatR)
# getOption("width",80)
# for(rfile in dir("R",full.names=TRUE)){
#     tidy_source(source=rfile,file=rfile,width.cutoff=40,indent=4)
# }


cat("\nimport(zoo)\n",file="NAMESPACE",append=TRUE)
cat("import(ChIPpeakAnno)\n",file="NAMESPACE",append=TRUE)
cat("import(locfit)\n",file="NAMESPACE",append=TRUE)
cat("import(TxDb.Hsapiens.UCSC.hg19.knownGene)\n",file="NAMESPACE",append=TRUE)


cat("importFrom('GenomicRanges','GRanges')\n",file="NAMESPACE",append=TRUE)
cat("importFrom('S4Vectors','Rle')\n",file="NAMESPACE",append=TRUE)
cat("importFrom('gplots','colorpanel')\n",file="NAMESPACE",append=TRUE)
cat("importFrom('caTools','runmean')\n",file="NAMESPACE",append=TRUE)
cat("importFrom('csaw','correlateReads')\n",file="NAMESPACE",append=TRUE)
cat("importFrom('wordcloud','wordlayout')\n",file="NAMESPACE",append=TRUE)
cat("importFrom('DESeq','newCountDataSet','estimateSizeFactors',\n    'estimateDispersions','varianceStabilizingTransformation')\n",
    file="NAMESPACE",append=TRUE)
cat("importFrom('Biobase','exprs')\n",file="NAMESPACE",append=TRUE)
cat("importFrom('viper','rowTtest','msviper','msviperAnnot','ttestNull')\n",file="NAMESPACE",append=TRUE)
cat("importFrom('DiffBind','dba','dba.count')\n",file="NAMESPACE",append=TRUE)




cat("importFrom('graphics', 'abline', 'grid', 'layout', 'legend', 'lines',\n    'par', 'plot', 'points', 'rect', 'text')\n",
    file="NAMESPACE",append=TRUE)
cat("importFrom('stats', 'ks.test', 'pchisq', 'pnorm', 'pt', 'qnorm', 'qt',\n    'rnorm', 'setNames', 'quantile', 'var')\n",
    file="NAMESPACE",append=TRUE)
cat("importFrom('utils', 'read.csv', 'setTxtProgressBar', 'txtProgressBar', 'relist')\n",
    file="NAMESPACE",append=TRUE)


