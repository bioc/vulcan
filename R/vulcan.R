#' Function to import BAM files
#'
#' This function coalesces and annotates a set of BAM files into peak-centered
#' data
#'
#' @param sheetfile path to a csv annotation file containing sample information
#' and BAM location
#' @param intervals size of the peaks. If NULL (default) it is inferred from
#' the average fragment length observed in the dataset
#'
#' @return a list
#'
#' @examples
#' require(vulcandata)
#' vfile<-"deleteme.csv"
#' vulcansheet(vfile)
#' vobj<-vulcan.import(vfile)
#' unlink(vfile)
#'
#' @export
vulcan.import<-function(sheetfile,intervals=NULL){
    # Check the dataset
    sheet<-read.csv(sheetfile,as.is=TRUE)

    # Generate a DiffBind object
    dbobj<-dba(sampleSheet=sheetfile)
    message("Sheet loaded. You have ",nrow(sheet), " samples and ",
            length(unique(sheet$Condition))," conditions")

    # Select the interval size automatically (if not provided by the user)
    if(is.null(intervals)){
        # List of bam files
        bam.files <- sheet[,"bamReads"]
        intervals<-average_fragment_length(bam.files,plot=FALSE)*2
        message("Peak size automatically detected as ",intervals,"nt")
    }

    # Count reads in binding sites intervals
    dbcounts<-dba.count(dbobj,summits=intervals)

    # Extract counts from the dbacount object
    listcounts<-dbcounts$peaks
    names(listcounts)<-dbcounts$samples[,1]

    # Prepare RPKM matrix
    first<-listcounts[[1]]
    rawmat<-matrix(NA, nrow = nrow(first), ncol = length(listcounts)+3)
    colnames(rawmat)<-c("Chr","Start","End",names(listcounts))
    rownames(rawmat)<-1:nrow(rawmat)
    rawmat<-as.data.frame(rawmat)
    rawmat[,1]<-as.character(first[,1])
    rawmat[,2]<-as.integer(first[,2])
    rawmat[,3]<-as.integer(first[,3])
    for(i in 1:length(listcounts)){
        rawmat[,names(listcounts)[i]]<-as.numeric(listcounts[[i]]$RPKM)
    }
    peakrpkms<-rawmat
    rm(rawmat)

    # Prepare Count matrix
    first<-listcounts[[1]]
    rawmat<-matrix(NA, nrow = nrow(first), ncol = length(listcounts)+3)
    colnames(rawmat)<-c("Chr","Start","End",names(listcounts))
    rownames(rawmat)<-1:nrow(rawmat)
    rawmat<-as.data.frame(rawmat)
    rawmat[,1]<-as.character(first[,1])
    rawmat[,2]<-as.integer(first[,2])
    rawmat[,3]<-as.integer(first[,3])
    for(i in 1:length(listcounts)){
        rawmat[,names(listcounts)[i]]<-as.integer(listcounts[[i]]$Reads)
    }
    peakcounts<-rawmat
    rm(rawmat)


    # Create an annotation structure
    samples<-list()
    conditions<-unique(sheet$Condition)
    for(cond in conditions){
        heresamples<-sheet$SampleID[sheet$Condition==cond]
        samples[[cond]]<-heresamples
    }

    # Return output
    vobj<-list(peakcounts=peakcounts,samples=samples,peakrpkms=peakrpkms)
    return(vobj)
}

#' Function to annotate peaks for VULCAN analysis
#'
#' This function coalesces and annotates a set of BAM files into peak-centered
#' data
#'
#' @param vobj A list of peakcounts, samples and peakrpkms (i.e. the output of
#' the funcion vulcan.import)
#' @param method Method to deal with multiple peaks found within gene promoter
#' boundaries.
#' One of "closest","strongest","sum"(default),"topvar","farthest","lowvar"
#' @param lborder Boundary for peak annotation (in nucleotides) upstream of
#' the Transcription starting site (default: -10000)
#' @param lborder Boundary for peak annotation (in nucleotides) downstream of
#' the Transcription starting site (default: 10000)
#' @return a list
#' @examples
#' require(vulcandata)
#' vfile<-"deleteme.csv"
#' vulcansheet(vfile)
#' vobj<-vulcan.import(vfile)
#' unlink(vfile)
#' vobj<-vulcan.annotate(vobj,lborder=-10000,rborder=10000,method="sum")
#' @export
vulcan.annotate<-function(vobj,lborder=-10000,rborder=10000,
                          method=c("closest","strongest","sum","topvar",
                                   "farthest","lowvar")
){
    # Annotate (hg19)
    annotation<-toGRanges(TxDb.Hsapiens.UCSC.hg19.knownGene, feature="gene")

    ##### PROCESS RAW COUNTS
    gr<-GRanges(vobj$peakcounts)
    anno<-annotatePeakInBatch(gr,AnnotationData=annotation,output="overlapping",
                              FeatureLocForDistance="TSS",
                              bindingRegion=c(lborder,rborder))

    # Convert to a more handy data frame
    dfanno<-anno
    names(dfanno)<-1:length(dfanno)
    dfanno<-as.data.frame(dfanno)

    # Prepare the output table
    allsamples<-unique(unlist(vobj$samples))
    genes<-unique(dfanno$feature)
    peakspergene<-table(dfanno$feature)
    rawcounts<-matrix(NA,nrow=length(genes),ncol=length(allsamples))
    colnames(rawcounts)<-allsamples
    rownames(rawcounts)<-genes

    # All methods: if a gene has a single peak, you select that
    for(i in 1:length(genes)){
        gene<-genes[i]
        if(peakspergene[gene]==1){
            rawcounts[gene,allsamples]<-as.numeric(dfanno[dfanno$feature==gene,
                                                          allsamples])
        }
    }

    # Method closest: when multiple peaks are found, keep only the closest to
    # the TSS as the representative one
    if(method=="closest"){
        for(i in 1:length(genes)){
            gene<-genes[i]
            if(peakspergene[gene]>1){
                subanno<-dfanno[dfanno$feature==gene,]
                closest<-which.min(subanno$distanceToStart)
                # if(any(subanno$distanceToStart<0)){
                #   stop("Stop! Negative distances are not allowed")
                # }
                rawcounts[gene,allsamples]<-as.numeric(
                    subanno[closest,allsamples])
            }
        }
    }

    # Method farthest: when multiple peaks are found, keep only the closest
    # to the TSS as the representative one
    if(method=="farthest"){
        for(i in 1:length(genes)){
            gene<-genes[i]
            if(peakspergene[gene]>1){
                subanno<-dfanno[dfanno$feature==gene,]
                farthest<-which.max(subanno$distanceToStart)
                rawcounts[gene,allsamples]<-as.numeric(
                    subanno[farthest,allsamples])
            }
        }
    }


    # Method sum: when multiple peaks are found, keep the strongest as the
    # representative one
    if(method=="sum"){
        for(i in 1:length(genes)){
            gene<-genes[i]
            if(peakspergene[gene]>1){
                subanno<-dfanno[dfanno$feature==gene,]
                sums<-apply(subanno[,allsamples],2,sum)
                rawcounts[gene,allsamples]<-as.numeric(sums)
            }
        }
    }

    # Method strongest: when multiple peaks are found, keep the strongest as
    # the representative one
    if(method=="strongest"){
        for(i in 1:length(genes)){
            gene<-genes[i]
            if(peakspergene[gene]>1){
                subanno<-dfanno[dfanno$feature==gene,]
                sums<-apply(subanno[,allsamples],1,sum)
                top<-which.max(sums)
                rawcounts[gene,allsamples]<-as.numeric(subanno[top,allsamples])
            }
        }
    }

    # Method topvar: when multiple peaks are found, keep the
    # most varying as the representative one
    if(method=="topvar"){
        for(i in 1:length(genes)){
            gene<-genes[i]
            if(peakspergene[gene]>1){
                subanno<-dfanno[dfanno$feature==gene,]
                vars<-apply(subanno[,allsamples],1,var)
                top<-which.max(vars)
                rawcounts[gene,allsamples]<-as.numeric(subanno[top,allsamples])
            }
        }
    }

    # Method lowvar: when multiple peaks are found,
    # keep the least varying as the representative one
    if(method=="lowvar"){
        for(i in 1:length(genes)){
            gene<-genes[i]
            if(peakspergene[gene]>1){
                subanno<-dfanno[dfanno$feature==gene,]
                vars<-apply(subanno[,allsamples],1,var)
                top<-which.min(vars)
                rawcounts[gene,allsamples]<-as.numeric(subanno[top,allsamples])
            }
        }
    }

    ##### PROCESS RPKMS
    gr<-GRanges(vobj$peakrpkms)
    anno<-annotatePeakInBatch(gr,AnnotationData=annotation,output="overlapping",
                              FeatureLocForDistance="TSS",
                              bindingRegion=c(lborder,rborder))

    # Convert to a more handy data frame
    dfanno<-anno
    names(dfanno)<-1:length(dfanno)
    dfanno<-as.data.frame(dfanno)

    # Prepare the output table
    allsamples<-unique(unlist(vobj$samples))
    genes<-unique(dfanno$feature)
    peakspergene<-table(dfanno$feature)
    rpkms<-matrix(NA,nrow=length(genes),ncol=length(allsamples))
    colnames(rpkms)<-allsamples
    rownames(rpkms)<-genes

    # All methods: if a gene has a single peak, you select that
    for(i in 1:length(genes)){
        gene<-genes[i]
        if(peakspergene[gene]==1){
            rpkms[gene,allsamples]<-as.numeric(dfanno[dfanno$feature==gene,
                                                      allsamples])
        }
    }

    # Method closest: when multiple peaks are found, keep only the closest
    # to the TSS as the representative one
    if(method=="closest"){
        for(i in 1:length(genes)){
            gene<-genes[i]
            if(peakspergene[gene]>1){
                subanno<-dfanno[dfanno$feature==gene,]
                closest<-which.min(subanno$distanceToStart)
                # if(any(subanno$distanceToStart<0)){
                #   stop("Stop! Negative distances are not allowed")
                # }
                rpkms[gene,allsamples]<-as.numeric(subanno[closest,allsamples])
            }
        }
    }

    # Method farthest: when multiple peaks are found, keep only the closest
    # to the TSS as the representative one
    if(method=="farthest"){
        for(i in 1:length(genes)){
            gene<-genes[i]
            if(peakspergene[gene]>1){
                subanno<-dfanno[dfanno$feature==gene,]
                farthest<-which.max(subanno$distanceToStart)
                rpkms[gene,allsamples]<-as.numeric(subanno[farthest,allsamples])
            }
        }
    }


    # Method sum: when multiple peaks are found, keep the strongest
    # as the representative one
    if(method=="sum"){
        for(i in 1:length(genes)){
            gene<-genes[i]
            if(peakspergene[gene]>1){
                subanno<-dfanno[dfanno$feature==gene,]
                sums<-apply(subanno[,allsamples],2,sum)
                rpkms[gene,allsamples]<-as.numeric(sums)
            }
        }
    }

    # Method strongest: when multiple peaks are found, keep the strongest
    # as the representative one
    if(method=="strongest"){
        for(i in 1:length(genes)){
            gene<-genes[i]
            if(peakspergene[gene]>1){
                subanno<-dfanno[dfanno$feature==gene,]
                sums<-apply(subanno[,allsamples],1,sum)
                top<-which.max(sums)
                rpkms[gene,allsamples]<-as.numeric(subanno[top,allsamples])
            }
        }
    }

    # Method topvar: when multiple peaks are found, keep the most varying
    # as the representative one
    if(method=="topvar"){
        for(i in 1:length(genes)){
            gene<-genes[i]
            if(peakspergene[gene]>1){
                subanno<-dfanno[dfanno$feature==gene,]
                vars<-apply(subanno[,allsamples],1,var)
                top<-which.max(vars)
                rpkms[gene,allsamples]<-as.numeric(subanno[top,allsamples])
            }
        }
    }

    # Method lowvar: when multiple peaks are found, keep the least varying
    # as the representative one
    if(method=="lowvar"){
        for(i in 1:length(genes)){
            gene<-genes[i]
            if(peakspergene[gene]>1){
                subanno<-dfanno[dfanno$feature==gene,]
                vars<-apply(subanno[,allsamples],1,var)
                top<-which.min(vars)
                rpkms[gene,allsamples]<-as.numeric(subanno[top,allsamples])
            }
        }
    }

    ### Fix data types as needed
    for(j in 1:ncol(rawcounts)){
        rawcounts[,j]<-as.numeric(rawcounts[,j])
    }
    rawcounts<-as.matrix(rawcounts)

    for(j in 1:ncol(rpkms)){
        rpkms[,j]<-as.numeric(rpkms[,j])
    }
    rpkms<-as.matrix(rpkms)

    # Return object
    vobj$rawcounts<-rawcounts
    vobj$rpkms<-rpkms
    return(vobj)
}



# vulcan.signature<-function(vobj,contrast=1){
#   anno<-vobj[[contrast]]
#   ## As signature we use -log10(sign(fold)*p)
#   # Fold change is not recapitulating the replicate agreements
#   signature<-anno$p.value*sign(anno$Fold)
#   names(signature)<-anno$feature
#   # Uniform ultrasmall gaussian noise (if by chance no reflist genes get
#   # into the null GSEA function, we get a division by zero)
#   set.seed(1)
#   othergenes<-setdiff(allgenes,names(signature))
# gaussiannoise<-setNames(rnorm(length(othergenes),
#                               mean=0,sd=0.01),othergenes) # very small
#   signature<-c(signature,gaussiannoise)
# }



#' Function to normalize promoter peak data
#'
#' This function normalizes gene-centered ChIP-Seq data using VST
#'
#' @param vobj a list, the output of the \code{"vulcan.annotate"} function
#'
#' @return a list
#'
#' @examples
#' require(vulcandata)
#' vfile<-"deleteme.csv"
#' vulcansheet(vfile)
#' vobj<-vulcan.import(vfile)
#' unlink(vfile)
#' vobj<-vulcan.annotate(vobj,lborder=-10000,rborder=10000,method="sum")
#' vobj<-vulcan.normalize(vobj)
#' @export
vulcan.normalize<-function(vobj){
    # Extract raw counts from object
    samples<-vobj$samples
    rawcounts<-vobj$rawcounts
    allsamples<-unique(unlist(samples))
    allgenes<-rownames(rawcounts)

    # Generate a normalized abundance object
    conditions<-c()
    for(i in 1:length(samples)){
        conditions<-c(conditions,rep(names(samples)[i],length(samples[[i]])))
    }
    conditions<-factor(conditions)
    cds<-newCountDataSet(vobj$rawcounts,conditions)
    cds<-estimateSizeFactors(cds)
    cds<-estimateDispersions(cds,fitType="local")
    vsd<-varianceStabilizingTransformation(cds)
    normalized<-exprs(vsd)
    rownames(normalized)<-rownames(rawcounts)
    vobj$normalized<-normalized
    return(vobj)
}

#' VULCAN - VirtUaL Chipseq data Analysis using Networks
#'
#' This function calculates the enrichment of a gene regulatory network over a
#' ChIP-Seq derived signature
#'
#' @param vobj a list, the output of the \code{"vulcan.normalize"} function
#' @param network an object of class \code{"viper::regulon"}
#' @param contrast a vector of two fields, containing the condition names
#' to be compared (1 vs 2)
#' @param annotation an optional named vector to convert gene identifiers
#' (e.g. entrez ids to gene symbols)
#' Default (NULL) won't convert gene names.
#'
#' @return a list
#'
#' @examples
#' require(vulcandata)
#' # Generate an annotation file from the dummy ChIP-Seq dataset
#' vfile<-"deleteme.csv"
#' vulcansheet(vfile)
#' # Import BAM and BED information into a list object
#' vobj<-vulcan.import(vfile)
#' unlink(vfile)
#' # Annotate peaks to gene names
#' vobj<-vulcan.annotate(vobj,lborder=-10000,rborder=10000,method="sum")
#' # Normalize data for VULCAN analysis
#' vobj<-vulcan.normalize(vobj)
#' # Detect which conditions are present
#' names(vobj$samples)
#'
#' # Load an ARACNe network
#' # This is a regulon object as specified in the VIPER package, named "network"
#' load(system.file("extdata","network.rda",package="vulcandata",mustWork=TRUE))
#' # Run VULCAN
#' # We can reduce the minimum regulon size, since in this example only one
#' # chromosome
#' # was measured, and the networks would otherwise have too few hits
#' vobj_analysis<-vulcan(vobj,network=network,contrast=c("t90","t0"),minsize=5)
#' # Visualize output using the msviper plotting function
#' plot(vobj_analysis$msviper,mrs=7)
#'
#' @export
vulcan<-function(vobj,network,contrast,annotation=NULL,minsize=10){
    tfs<-names(network)
    samples<-vobj$samples
    normalized<-vobj$normalized

    # Prepare output objects
    msvipers<-matrix(NA,ncol=3,nrow=length(tfs))
    rownames(msvipers)<-tfs
    # Define contrast
    a<-samples[[contrast[1]]]
    b<-samples[[contrast[2]]]
    # Vulcan msviper implementation
    set.seed(1)
    signature<-rowTtest(normalized[,a],normalized[,b])$statistic
    dnull<-ttestNull(normalized[,a],normalized[,b],per=1000)
    msviper<-msviper(signature,network,dnull,minsize=minsize)
    # Annotate
    if(!is.null(annotation)){
        msviper<-msviperAnnot(msviper,annotation)
    }
    vobj$msviper<-msviper
    # Specific Master Regulators
    mrs<-cbind(msviper$es$nes,z2p(msviper$es$nes))
    colnames(mrs)<-c("NES","pvalue")
    vobj$mrs<-mrs

    return(vobj)
}


vulcan.pathways<-function(vobj,pathways,contrast=NULL,method=c("GSEA","REA")){
    normalized<-vobj$normalized
    samples<-vobj$samples
    allgenes<-unique(unlist(pathways))

    # Specific contrast
    if(!setequal(contrast,"all")){
        # Define contrast
        a<-samples[[contrast[1]]]
        b<-samples[[contrast[2]]]

        # Prepare signature
        set.seed(1)
        signature<-rowTtest(normalized[,a],normalized[,b])$statistic
        othergenes<-setdiff(allgenes,names(signature))
        gaussiannoise<-setNames(rnorm(length(othergenes),mean=0,sd=0.01),
                                othergenes) # very small
        signature<-c(signature,gaussiannoise)

        # GSEA
        if(method=="GSEA"){
            gsea.pathways<-setNames(rep(0,length(pathways)),names(pathways))
            message("Running GSEA for ",length(pathways)," pathways")
            pb<-txtProgressBar(0,length(pathways),style=3)
            i<-0
            for(pname in names(pathways)){
                p<-pathways[[pname]]
                obj<-gsea(reflist=signature,set=p,method="pareto",np=100)
                gsea.pathways[pname]<-obj$nes
                setTxtProgressBar(pb,i<-i+1)
            }
            return(gsea.pathways)
        }

        # REA
        if(method=="REA"){
            rea.pathways<-setNames(rep(0,length(pathways)),names(pathways))
            message("Running REA for ",length(pathways)," pathways")
            rea.pathways<-rea(signatures=signature,groups=pathways,minsize=1)
            return(rea.pathways)
        }
    } else {
        if(method!="REA"){
            stop("Multiple signatures supported only with method='REA'")
        }
        signatures<-t(scale(t(vobj$normalized)))
        rea.pathways<-rea(signatures=signatures,groups=pathways,minsize=1)
        return(rea.pathways)
    }



}











