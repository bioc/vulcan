

# Building manuals and namespace using roxygen2
library(devtools)
document()

cat("\nimport(zoo)\n",file="NAMESPACE",append=TRUE)
cat("import(IRanges)\n",file="NAMESPACE",append=TRUE)
cat("import(viper)\n",file="NAMESPACE",append=TRUE)
cat("import(ChIPpeakAnno)\n",file="NAMESPACE",append=TRUE)
cat("import(DESeq)\n",file="NAMESPACE",append=TRUE)
cat("import(DiffBind)\n",file="NAMESPACE",append=TRUE)
cat("import(TxDb.Hsapiens.UCSC.hg19.knownGene)\n",file="NAMESPACE",append=TRUE)
cat("import(csaw)\n",file="NAMESPACE",append=TRUE)
cat("import(gplots)\n",file="NAMESPACE",append=TRUE)
cat("import(wordcloud)\n",file="NAMESPACE",append=TRUE)
cat("import(caTools)\n",file="NAMESPACE",append=TRUE)


