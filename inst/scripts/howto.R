

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
cat("import(IRanges)\n",file="NAMESPACE",append=TRUE)
cat("import(viper)\n",file="NAMESPACE",append=TRUE)
cat("import(ChIPpeakAnno)\n",file="NAMESPACE",append=TRUE)
cat("import(DESeq)\n",file="NAMESPACE",append=TRUE)
cat("import(DiffBind)\n",file="NAMESPACE",append=TRUE)
cat("import(TxDb.Hsapiens.UCSC.hg19.knownGene)\n",file="NAMESPACE",append=TRUE)
cat("import(csaw)\n",file="NAMESPACE",append=TRUE)
cat("import(wordcloud)\n",file="NAMESPACE",append=TRUE)
cat("import(caTools)\n",file="NAMESPACE",append=TRUE)

