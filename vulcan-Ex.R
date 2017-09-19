pkgname <- "vulcan"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
options(pager = "console")
library('vulcan')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
cleanEx()
nameEx("average_fragment_length")
### * average_fragment_length

flush(stderr()); flush(stdout())

### Name: average_fragment_length
### Title: Define the average fragment length
### Aliases: average_fragment_length

### ** Examples

library(vulcandata)
sheetfile<-'deleteme.csv'
vulcandata::vulcansheet(sheetfile)
a<-read.csv(sheetfile,as.is=TRUE)
bams<-a$bamReads
unlink(sheetfile)
average_fragment_length(bams,plot=TRUE)



cleanEx()
nameEx("corr2p")
### * corr2p

flush(stderr()); flush(stdout())

### Name: corr2p
### Title: Convert correlation coefficient to p-value
### Aliases: corr2p

### ** Examples

set.seed(1)
a<-rnorm(1000)
b<-a+rnorm(1000,sd=10)
r<-cor(a,b,method='pearson')
corr2p(r,N=length(a))



cleanEx()
nameEx("densityauc")
### * densityauc

flush(stderr()); flush(stdout())

### Name: densityauc
### Title: densityauc - Calculate the AUC of a density object
### Aliases: densityauc

### ** Examples

set.seed(1)
a<-rnorm(1000)
d<-density(a)
window<-c(2,3)
da<-densityauc(d,window)

plot(d,main='')
abline(v=window,lty=2)
title(paste0('AUC between lines=',da))






cleanEx()
nameEx("dpareto")
### * dpareto

flush(stderr()); flush(stdout())

### Name: dpareto
### Title: Probability density of Pareto distributions
### Aliases: dpareto

### ** Examples

set.seed(1)
x<-abs(rnorm(1000))
n<-length(x)
exponent<-1+n/sum(log(x))
dp<-dpareto(x,exponent=exponent)
plot(dp)



cleanEx()
nameEx("fisherp")
### * fisherp

flush(stderr()); flush(stdout())

### Name: fisherp
### Title: Fisher integration of p-values
### Aliases: fisherp

### ** Examples

ps<-c(0.01,0.05,0.03,0.2)
fisherp(ps)



cleanEx()
nameEx("gsea")
### * gsea

flush(stderr()); flush(stdout())

### Name: gsea
### Title: GSEA
### Aliases: gsea

### ** Examples

reflist<-setNames(-sort(rnorm(1000)),paste0('gene',1:1000))
set<-paste0('gene',sample(1:200,50))
obj<-gsea(reflist,set,method='pareto',np=1000)
obj$p.value



cleanEx()
nameEx("kmgformat")
### * kmgformat

flush(stderr()); flush(stdout())

### Name: kmgformat
### Title: kmgformat - Nice Formatting of Numbers
### Aliases: kmgformat

### ** Examples

# Thousands
set.seed(1)
a<-runif(1000,0,1e4)
plot(a,yaxt='n')
kmg<-kmgformat(pretty(a))
axis(2,at=pretty(a),labels=kmg)

# Millions to Billions
set.seed(1)
a<-runif(1000,0,1e9)
plot(a,yaxt='n',pch=20,col=val2col(a))
kmg<-kmgformat(pretty(a))
axis(2,at=pretty(a),labels=kmg)



cleanEx()
nameEx("null_gsea")
### * null_gsea

flush(stderr()); flush(stdout())

### Name: null_gsea
### Title: Calculate Null Distribution for GSEA
### Aliases: null_gsea

### ** Examples

reflist<-setNames(-sort(rnorm(26)),LETTERS)
set<-c('A','B','D','F')
nulldist<-null_gsea(set,reflist)
nulldist[1:10]



cleanEx()
nameEx("p2corr")
### * p2corr

flush(stderr()); flush(stdout())

### Name: p2corr
### Title: Convert p-value to correlation coefficient
### Aliases: p2corr

### ** Examples

N<-100
p<-0.05
p2corr(p,N)



cleanEx()
nameEx("p2z")
### * p2z

flush(stderr()); flush(stdout())

### Name: p2z
### Title: p2z
### Aliases: p2z

### ** Examples

p<-0.05
p2z(p)



cleanEx()
nameEx("pareto.fit")
### * pareto.fit

flush(stderr()); flush(stdout())

### Name: pareto.fit
### Title: Estimate parameters of Pareto distribution
### Aliases: pareto.fit

### ** Examples

# Estimate the tail of a pospulation normally distributed
set.seed(1)
x<-rnorm(1000)
q95<-as.numeric(quantile(abs(x),0.95))
fit<-pareto.fit(abs(x),threshold=q95)
# We can infer the pvalue of a value very much right to the tail of the
# distribution
value<-5
pvalue<-ppareto(value, threshold=q95, exponent=fit$exponent,
lower.tail=FALSE)/20
plot(density(abs(x)),xlim=c(0,value+0.3),main='Pareto fit')
arrows(value,0.2,value,0)
text(value,0.2,labels=paste0('p=',signif(pvalue,2)))



cleanEx()
nameEx("plot_gsea")
### * plot_gsea

flush(stderr()); flush(stdout())

### Name: plot_gsea
### Title: Plot GSEA results
### Aliases: plot_gsea

### ** Examples

reflist<-setNames(-sort(rnorm(26)),LETTERS)
set<-c('A','B','D','F')
obj<-gsea(reflist,set,method='pareto')
plot_gsea(obj)



cleanEx()
nameEx("ppareto")
### * ppareto

flush(stderr()); flush(stdout())

### Name: ppareto
### Title: Cumulative distribution function of the Pareto distributions '
###   Gives NA on values < threshold
### Aliases: ppareto

### ** Examples

# Estimate the tail of a pospulation normally distributed
set.seed(1)
x<-rnorm(1000)
q95<-as.numeric(quantile(abs(x),0.95))
fit<-pareto.fit(abs(x),threshold=q95)
# We can infer the pvalue of a value very much right to the tail of the
# distribution
value<-5
pvalue<-ppareto(value, threshold=q95, exponent=fit$exponent,
lower.tail=FALSE)/20
plot(density(abs(x)),xlim=c(0,value+0.3),main='Pareto fit')
arrows(value,0.2,value,0)
text(value,0.2,labels=paste0('p=',signif(pvalue,2)))



cleanEx()
nameEx("rea")
### * rea

flush(stderr()); flush(stdout())

### Name: rea
### Title: REA: Rank EnrichmeNt Analysis
### Aliases: rea

### ** Examples

signatures<-setNames(-sort(rnorm(1000)),paste0('gene',1:1000))
set1<-paste0('gene',sample(1:200,50))
set2<-paste0('gene',sample(1:1000,50))
groups<-list(set1=set1,set2=set2)
obj<-rea(signatures=signatures,groups=groups)
obj



cleanEx()
nameEx("slice")
### * slice

flush(stderr()); flush(stdout())

### Name: slice
### Title: Slice
### Aliases: slice

### ** Examples

set.seed(1)
example<-matrix(rnorm(1000),nrow=100,ncol=10)
slice(example)



cleanEx()
nameEx("stouffer")
### * stouffer

flush(stderr()); flush(stdout())

### Name: stouffer
### Title: Stouffer integration of Z scores
### Aliases: stouffer

### ** Examples

zs<-c(1,3,5,2,3)
stouffer(zs)



cleanEx()
nameEx("textplot2")
### * textplot2

flush(stderr()); flush(stdout())

### Name: textplot2
### Title: textplot2 - An x y plot of non-overlapping text
### Aliases: textplot2

### ** Examples

obj_names<-apply(expand.grid(LETTERS,LETTERS),1,paste,collapse='')
a<-setNames(runif(26*26),obj_names)
b<-setNames(rnorm(26*26),obj_names)
plot(a,b,pch=20,col='grey')
top<-names(sort(-a))[1:50]
textplot2(a[top],b[top],words=top,new=FALSE,pointcolor='black')



cleanEx()
nameEx("val2col")
### * val2col

flush(stderr()); flush(stdout())

### Name: val2col
### Title: Convert a numeric vector into colors
### Aliases: val2col

### ** Examples

a<-rnorm(1000)
cols<-val2col(a)
plot(a,col=cols,pch=16)



cleanEx()
nameEx("vulcan")
### * vulcan

flush(stderr()); flush(stdout())

### Name: vulcan
### Title: VULCAN - VirtUaL Chipseq data Analysis using Networks
### Aliases: vulcan

### ** Examples

library(vulcandata)
# Generate an annotation file from the dummy ChIP-Seq dataset
vfile<-'deleteme.csv'
vulcandata::vulcansheet(vfile)
# Import BAM and BED information into a list object
vobj<-vulcan.import(vfile)
unlink(vfile)
# Annotate peaks to gene names
vobj<-vulcan.annotate(vobj,lborder=-10000,rborder=10000,method='sum')
# Normalize data for VULCAN analysis
vobj<-vulcan.normalize(vobj)
# Detect which conditions are present
names(vobj$samples)

# Load an ARACNe network
# This is a regulon object as specified in the VIPER package, named 'network'
load(system.file('extdata','network.rda',package='vulcandata',mustWork=TRUE))
# Run VULCAN
# We can reduce the minimum regulon size, since in this example only one
# chromosome
# was measured, and the networks would otherwise have too few hits
vobj_analysis<-vulcan(vobj,network=network,contrast=c('t90','t0'),minsize=5)
# Visualize output using the msviper plotting function
plot(vobj_analysis$msviper,mrs=7)




cleanEx()
nameEx("vulcan.annotate")
### * vulcan.annotate

flush(stderr()); flush(stdout())

### Name: vulcan.annotate
### Title: Function to annotate peaks for VULCAN analysis
### Aliases: vulcan.annotate

### ** Examples

library(vulcandata)
vfile<-'deleteme.csv'
vulcandata::vulcansheet(vfile)
vobj<-vulcan.import(vfile)
unlink(vfile)
vobj<-vulcan.annotate(vobj,lborder=-10000,rborder=10000,method='sum')



cleanEx()
nameEx("vulcan.import")
### * vulcan.import

flush(stderr()); flush(stdout())

### Name: vulcan.import
### Title: Function to import BAM files
### Aliases: vulcan.import

### ** Examples

library(vulcandata)
vfile<-'deleteme.csv'
vulcandata::vulcansheet(vfile)
vobj<-vulcan.import(vfile)
unlink(vfile)




cleanEx()
nameEx("vulcan.normalize")
### * vulcan.normalize

flush(stderr()); flush(stdout())

### Name: vulcan.normalize
### Title: Function to normalize promoter peak data
### Aliases: vulcan.normalize

### ** Examples

library(vulcandata)
vfile<-'deleteme.csv'
vulcandata::vulcansheet(vfile)
vobj<-vulcan.import(vfile)
unlink(vfile)
vobj<-vulcan.annotate(vobj,lborder=-10000,rborder=10000,method='sum')
vobj<-vulcan.normalize(vobj)



cleanEx()
nameEx("vulcan.pathways")
### * vulcan.pathways

flush(stderr()); flush(stdout())

### Name: vulcan.pathways
### Title: Function to calculate pathway enrichment over a ChIP-Seq profile
### Aliases: vulcan.pathways

### ** Examples

library(vulcandata)
vfile<-'deleteme.csv'
vulcandata::vulcansheet(vfile)
vobj<-vulcan.import(vfile)
unlink(vfile)
vobj<-vulcan.annotate(vobj,lborder=-10000,rborder=10000,method='sum')
vobj<-vulcan.normalize(vobj)
# Create a dummy pathway list (in entrez ids)
pathways<-list(
   pathwayA=rownames(vobj$normalized)[1:20],
   pathwayB=rownames(vobj$normalized)[21:50]
)
# Which contrast groups can be queried
names(vobj$samples)
results_gsea<-vulcan.pathways(vobj,pathways,contrast=c('t90','t0'),
method='GSEA')
results_rea<-vulcan.pathways(vobj,pathways,contrast=c('all'),method='REA')



cleanEx()
nameEx("wstouffer")
### * wstouffer

flush(stderr()); flush(stdout())

### Name: wstouffer
### Title: Weighted Stouffer integration of Z scores
### Aliases: wstouffer

### ** Examples

zs<-c(1,-3,5,2,3)
ws<-c(1,10,1,2,1)
wstouffer(zs,ws)



cleanEx()
nameEx("z2p")
### * z2p

flush(stderr()); flush(stdout())

### Name: z2p
### Title: z2p
### Aliases: z2p

### ** Examples

z<-1.96
z2p(z)



### * <FOOTER>
###
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
