# This small script will generate smaller version of the BAM files, keeping only
# one chromosome

# Input folder
cd /Users/giorgi01/Dropbox/rstudio/vulcan/chipseq/holding

# Output folder
cd /Users/giorgi01/Dropbox/rstudio/vulcanpackage/vulcan/inst/extdata


### Process beds
cd /Users/giorgi01/Dropbox/rstudio/vulcan/chipseq/holding
for f in *bed.gz
do
zcat $f | grep "chr22" | gzip --stdout > /Users/giorgi01/Dropbox/rstudio/vulcanpackage/vulcan/inst/extdata/$f
done

### Process bams
cd /Users/giorgi01/Dropbox/rstudio/vulcan/chipseq/holding
for f in *bam
do
samtools view -hb $f chr22 > /Users/giorgi01/Dropbox/rstudio/vulcanpackage/vulcan/inst/extdata/$f
done


### Generate new BAIs
cd /Users/giorgi01/Dropbox/rstudio/vulcanpackage/vulcan/inst/extdata
for f in *bam
do
samtools index $f
done


### Which genes is this chromosome covering?
# Generate an annotation file from the vulcandata reduced dataset
vfile<-"deleteme.csv"
vulcansheet(vfile)
# Import BAM and BED information into a list object
vobj<-vulcan.import(vfile)
unlink(vfile)
# Annotate peaks to gene names
vobj<-vulcan.annotate(vobj,lborder=-10000,rborder=10000,method="sum")
# Normalize data for VULCAN analysis
vobj<-vulcan.normalize(vobj)
# Detect which conditions are present
names(vobj$samples)

tt<-apply(vobj$normalized,1,function(x){
    tt<-t.test(x[1:4],x[9:12]);return(tt$stat)
})
genes<-names(sort(tt))


### Generate dummy network (network.rda)
set.seed(40)

network<-list()

for(tf in c("YBX1","CLOCK","RB1")){
    subgenes<-sample(genes[1:floor(length(genes)/2)],45)
    tfmode<-setNames(runif(length(genes),-1,0.5),subgenes)
    likelihood<-runif(length(genes),0,1)
    subnet<-list(tfmode=tfmode,likelihood=likelihood)
    network[[tf]]<-subnet
}
for(tf in c("ESR1","GATA3","RARA","FOXM1","RXRA")){
    subgenes<-sample(genes[ceiling(length(genes)/2):length(genes)],45)
    tfmode<-setNames(runif(length(genes),-0.5,1),subgenes)
    likelihood<-runif(length(genes),0,1)
    subnet<-list(tfmode=tfmode,likelihood=likelihood)
    network[[tf]]<-subnet
}


class(network)<-"regulon"
save(network,file="/Users/giorgi01/Dropbox/rstudio/vulcanpackage/vulcan/inst/extdata/network.rda")



# Run VULCAN
# We can reduce the minimum regulon size, since in this example only one chromosome
# was measured, and the networks would otherwise have too few hits
vobj_analysis<-vulcan(vobj,network=network,contrast=c("t90","t0"),minsize=5)
# Visualize output using the msviper plotting function
plot(vobj_analysis$msviper,mrs=7)







