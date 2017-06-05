

# Building manuals and namespace using roxygen2
library(devtools)
document()

cat("\nimport(zoo)\nimport(IRanges)\nimport(viper)\n",file="NAMESPACE",append=TRUE)


