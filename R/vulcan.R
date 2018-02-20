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
#' @return A list of components:
#' \describe{
#' \item{peakcounts}{A matrix of raw peak counts, peaks as rows, samples as
#' columns}
#' \item{peakrpkms}{A matrix of peak RPKMs, peaks as rows, samples as
#' columns}
#' \item{samples}{A vector of sample names and conditions}
#' }
#' @examples
#' library(vulcandata)
#' # Generate an annotation file from the dummy ChIP-Seq dataset
#' vfile<-tempfile()
#' vulcandata::vulcansheet(vfile)
#' # Import BAM and BED information into a list object
#' # vobj<-vulcan.import(vfile)
#' # This vobj is identical to the object returned by
#' # vulcandata::vulcanexample()
#' unlink(vfile)
#' @export
vulcan.import <- function(sheetfile, intervals = NULL) {
    # Check the dataset
    sheet <- read.csv(sheetfile, as.is = TRUE)

    # Generate a DiffBind object
    dbobj <- DiffBind::dba(sampleSheet = sheetfile)
    message("Sheet loaded. You have ", nrow(sheet),
            " samples and ", length(unique(sheet$Condition)),
            " conditions")

    # Select the interval size automatically
    # (if not provided by the user)
    if (is.null(intervals)) {
        # List of bam files
        bam.files <- sheet[, "bamReads"]
        intervals <- average_fragment_length(bam.files, plot = FALSE) * 2
        message("Peak size automatically detected as ",
                intervals, "nt")
    }

    # Count reads in binding sites intervals
    dbcounts <- DiffBind::dba.count(dbobj, summits = intervals)

    # Extract counts from the dbacount object
    listcounts <- dbcounts$peaks
    names(listcounts) <- dbcounts$samples[,1]

    # Prepare RPKM matrix
    first <- listcounts[[1]]
    rawmat <- matrix(NA, nrow = nrow(first), ncol = length(listcounts) + 3)
    colnames(rawmat) <- c("Chr", "Start", "End", names(listcounts))
    rownames(rawmat) <- 1:nrow(rawmat)
    rawmat <- as.data.frame(rawmat)
    rawmat[, 1] <- as.character(first[, 1])
    rawmat[, 2] <- as.integer(first[, 2])
    rawmat[, 3] <- as.integer(first[, 3])
    for (i in seq_len(length(listcounts))) {
        rawmat[, names(listcounts)[i]] <- as.numeric(listcounts[[i]]$RPKM)
    }
    peakrpkms <- rawmat

    # Prepare Count matrix
    first <- listcounts[[1]]
    rawmat <- matrix(NA, nrow = nrow(first), ncol = length(listcounts) + 3)
    colnames(rawmat) <- c("Chr", "Start", "End", names(listcounts))
    rownames(rawmat) <- 1:nrow(rawmat)
    rawmat <- as.data.frame(rawmat)
    rawmat[, 1] <- as.character(first[, 1])
    rawmat[, 2] <- as.integer(first[, 2])
    rawmat[, 3] <- as.integer(first[, 3])
    for (i in seq_len(length(listcounts))) {
        rawmat[, names(listcounts)[i]] <- as.integer(listcounts[[i]]$Reads)
    }
    peakcounts <- rawmat


    # Create an annotation structure
    samples <- list()
    conditions <- unique(sheet$Condition)
    for (cond in conditions) {
        heresamples <- sheet$SampleID[sheet$Condition == cond]
        samples[[cond]] <- heresamples
    }

    # Return output
    vobj <- list(peakcounts = peakcounts, samples = samples,
                peakrpkms = peakrpkms)
    return(vobj)
}

#' Function to annotate peaks for VULCAN analysis
#'
#' This function coalesces and annotates a set of BAM files into peak-centered
#' data. It implements the ChIPPeakANno methods, with specific choices dealing
#' with defining the genomic area around the promoter and which peaks to
#' include.
#'
#' @param vobj A list of peakcounts, samples and peakrpkms (i.e. the output of
#' the function vulcan.import)
#' @param method Method to deal with multiple peaks found within gene promoter
#' boundaries. One of sum (default), closest, strongest, topvar, farthest or
#' lowvar. This will affect only genes with multiple possible peaks. When a
#' single peak can be mapped to the promoter region of the gene, that peak
#' abundance will be considered as the gene promoter's occupancy.
#' \describe{
#' \item{sum}{when multiple peaks are found, sum their contributions}
#' \item{closest}{when multiple peaks are found, keep only the closest to the
#' TSS as the representative one}
#' \item{strongest}{when multiple peaks are found, keep the strongest as the
#' representative one}
#' \item{farthest}{when multiple peaks are found, keep only the closest to the
#' TSS as the representative one}
#' \item{topvar}{when multiple peaks are found, keep the most varying as the
#' representative one}
#' \item{lowvar}{when multiple peaks are found, keep the least varying as the
#' representative one}
#' }
#' @param lborder Boundary for peak annotation (in nucleotides) upstream of
#' the Transcription starting site (default: -10000)
#' @param rborder Boundary for peak annotation (in nucleotides) downstream of
#' the Transcription starting site (default: 10000)
#' @param TxDb TxDb annotation object containing the knownGene track. If NULL
#' (the default), TxDb.Hsapiens.UCSC.hg19.knownGene is loaded
#' @return A list of components:
#' \describe{
#' \item{peakcounts}{A matrix of raw peak counts, peaks as rows, samples as
#' columns}
#' \item{peakrpkms}{A matrix of peak RPKMs, peaks as rows, samples as
#' columns}
#' \item{rawcounts}{A matrix of raw gene counts, genes as rows, samples as
#' columns. The counts are associated to the promoter region of the gene}
#' \item{rpkms}{A matrix of RPKMs, genes as rows, samples as
#' columns. The RPKMs are associated to the promoter region of the gene}
#' \item{samples}{A vector of sample names and conditions}
#' }
#' @examples
#' library(vulcandata)
#' vobj<-vulcandata::vulcanexample()
#' vobj<-vulcan.annotate(vobj,lborder=-10000,rborder=10000,method='sum')
#' @export
vulcan.annotate <- function(vobj, lborder = -10000,
                            rborder = 10000, method = c("closest",
                                                        "strongest",
                                                        "sum",
                                                        "topvar",
                                                        "farthest",
                                                        "lowvar"),
                            TxDb=NULL) {

    if(!is.null(TxDb)){
        # Annotate (any genome)
        annotation <- toGRanges(TxDb,
                                feature = "gene")

    } else {
        # Annotate (hg19)
        annotation <- toGRanges(TxDb.Hsapiens.UCSC.hg19.knownGene,
                                feature = "gene")
    }


    ##### PROCESS RAW COUNTS
    gr <- GRanges(vobj$peakcounts)
    anno <- annotatePeakInBatch(gr, AnnotationData = annotation,
                                output = "overlapping",
                                FeatureLocForDistance = "TSS",
                                bindingRegion = c(lborder, rborder))

    # Convert to a more handy data frame
    dfanno <- anno
    names(dfanno) <- seq_len(length(dfanno))
    dfanno <- as.data.frame(dfanno)

    # Prepare the output table
    allsamples <- unique(unlist(vobj$samples))
    genes <- unique(dfanno$feature)
    peakspergene <- table(dfanno$feature)
    rawcounts <- matrix(NA, nrow = length(genes),
                        ncol = length(allsamples))
    colnames(rawcounts) <- allsamples
    rownames(rawcounts) <- genes

    # All methods: if a gene has a single
    # peak, you select that
    genesone <- names(peakspergene)[peakspergene ==
                                        1]
    for (gene in genesone) {
        rawcounts[gene, allsamples] <- as.numeric(dfanno[dfanno$feature ==
                                                            gene, allsamples])
    }

    # Other methods: they deal with cases
    # where multiple peaks are found
    genesmore <- names(peakspergene)[peakspergene > 1]
    rawcounts<-dist_calc(method,dfanno,rawcounts,genesmore,allsamples)

    ##### PROCESS RPKMS
    gr <- GRanges(vobj$peakrpkms)
    anno <- annotatePeakInBatch(gr, AnnotationData = annotation,
                                output = "overlapping",
                                FeatureLocForDistance = "TSS",
                                bindingRegion = c(lborder, rborder))

    # Convert to a more handy data frame
    dfanno <- anno
    names(dfanno) <- seq_len(length(dfanno))
    dfanno <- as.data.frame(dfanno)

    # Prepare the output table
    allsamples <- unique(unlist(vobj$samples))
    genes <- unique(dfanno$feature)
    peakspergene <- table(dfanno$feature)
    rpkms <- matrix(NA, nrow = length(genes),
                    ncol = length(allsamples))
    colnames(rpkms) <- allsamples
    rownames(rpkms) <- genes

    # All methods: if a gene has a single
    # peak, you select that
    genesone <- names(peakspergene)[peakspergene == 1]
    for (gene in genesone) {
        rpkms[gene, allsamples] <- as.numeric(
            dfanno[dfanno$feature == gene, allsamples]
        )
    }

    # Other methods: they deal with cases
    # where multiple peaks are found
    genesmore <- names(peakspergene)[peakspergene > 1]
    rpkms<-dist_calc(method,dfanno,rpkms,genesmore,allsamples)


    ### Fix data types as needed
    rawcounts<-matrix(as.numeric(rawcounts),nrow=nrow(rawcounts),
                    dimnames=dimnames(rawcounts))
    rpkms<-matrix(as.numeric(rpkms),nrow=nrow(rpkms),
                dimnames=dimnames(rpkms))

    # Return object
    vobj$rawcounts <- rawcounts
    vobj$rpkms <- rpkms
    return(vobj)
}


dist_calc<-function(method,dfanno,genematrix,genesmore,allsamples){
    # This function structure was strongly suggested
    # by the Bioconductor reviewer
    supportedMethods<-c(
        "closest",
        "strongest",
        "sum",
        "topvar",
        "farthest",
        "lowvar"
    )
    if(!method%in%supportedMethods){
        stop("unsupported method ", method)
    }

    # Method closest: when multiple peaks are
    # found, keep only the closest to the TSS
    # as the representative one

    # Method farthest: when multiple peaks
    # are found, keep only the closest to the
    # TSS as the representative one

    # Method sum: when multiple peaks are
    # found, sum their contributions

    # Method strongest: when multiple peaks
    # are found, keep the strongest as the
    # representative one

    # Method topvar: when multiple peaks are
    # found, keep the most varying as the
    # representative one

    # Method lowvar: when multiple peaks are
    # found, keep the least varying as the
    # representative one


    for (gene in genesmore) {
        subanno <- dfanno[dfanno$feature == gene, ]

        if (method == "closest") {
            closest <- which.min(subanno$distanceToStart)
            genematrix[gene, allsamples] <- as.numeric(subanno[closest,
                                                            allsamples])
        }

        if (method == "farthest") {
            farthest <- which.max(subanno$distanceToStart)
            genematrix[gene, allsamples] <- as.numeric(subanno[farthest,
                                                            allsamples])
        }

        if (method == "sum") {
            sums <- apply(subanno[, allsamples], 2, sum)
            genematrix[gene, allsamples] <- as.numeric(sums)
        }

        if (method == "strongest") {
            sums <- apply(subanno[, allsamples], 1, sum)
            top <- which.max(sums)
            genematrix[gene, allsamples] <- as.numeric(subanno[top,
                                                            allsamples])
        }

        if (method == "topvar") {
            vars <- apply(subanno[, allsamples], 1, var)
            top <- which.max(vars)
            genematrix[gene, allsamples] <- as.numeric(subanno[top,
                                                            allsamples])
        }

        if (method == "lowvar") {
            vars <- apply(subanno[, allsamples], 1, var)
            top <- which.min(vars)
            genematrix[gene, allsamples] <- as.numeric(subanno[top,
                                                            allsamples])
        }

    }
    return(genematrix)
}



#' Function to normalize promoter peak data
#'
#' This function normalizes gene-centered ChIP-Seq data using VST
#'
#' @param vobj a list, the output of the \code{'vulcan.annotate'} function
#'
#' @return A list of components:
#' \describe{
#' \item{peakcounts}{A matrix of raw peak counts, peaks as rows, samples as
#' columns}
#' \item{peakrpkms}{A matrix of peak RPKMs, peaks as rows, samples as
#' columns}
#' \item{rawcounts}{A matrix of raw gene counts, genes as rows, samples as
#' columns. The counts are associated to the promoter region of the gene}
#' \item{rpkms}{A matrix of RPKMs, genes as rows, samples as
#' columns. The RPKMs are associated to the promoter region of the gene}
#' \item{normalized}{A matrix of gene abundances normalized by
#' Variance-Stabilizing Transformation (VST), genes as rows, samples as
#' columns. The abundances are associated to the promoter
#' region of the gene}
#' \item{samples}{A vector of sample names and conditions}
#' }
#' @examples
#' \dontrun{
#' library(vulcandata)
#' vobj<-vulcandata::vulcanexample()
#' vobj<-vulcan.annotate(vobj,lborder=-10000,rborder=10000,method='sum')
#' vobj<-vulcan.normalize(vobj)
#' }
#' @export
vulcan.normalize <- function(vobj) {
    # Extract raw counts from object
    samples <- vobj$samples
    rawcounts <- vobj$rawcounts
    allsamples <- unique(unlist(samples))
    allgenes <- rownames(rawcounts)

    # Generate a normalized abundance object
    conditions <- c()
    for (i in seq_len(length(samples))) {
        conditions <- c(conditions, rep(names(samples)[i],
                                        length(samples[[i]])))
    }
    conditions <- factor(conditions)
    cds <- DESeq::newCountDataSet(vobj$rawcounts, conditions)
    cds <- DESeq::estimateSizeFactors(cds)
    cds <- DESeq::estimateDispersions(cds, fitType = "local")
    vsd <- DESeq::varianceStabilizingTransformation(cds)
    normalized <- Biobase::exprs(vsd)
    rownames(normalized) <- rownames(rawcounts)
    vobj$normalized <- normalized
    return(vobj)
}

#' VULCAN - VirtUaL Chipseq data Analysis using Networks
#'
#' This function calculates the enrichment of a gene regulatory network over a
#' ChIP-Seq derived signature
#'
#' @param vobj a list, the output of the \code{'vulcan.normalize'} function
#' @param network an object of class \code{'viper::regulon'}
#' @param contrast a vector of two fields, containing the condition names
#' to be compared (1 vs 2)
#' @param annotation an optional named vector to convert gene identifiers
#' (e.g. entrez ids to gene symbols)
#' Default (NULL) won't convert gene names.
#' @param minsize integer indicating the minimum regulon size for the analysis
#' to be run. Default: 10
#'
#' @return A list of components:
#' \describe{
#' \item{peakcounts}{A matrix of raw peak counts, peaks as rows, samples as
#' columns}
#' \item{peakrpkms}{A matrix of peak RPKMs, peaks as rows, samples as
#' columns}
#' \item{rawcounts}{A matrix of raw gene counts, genes as rows, samples as
#' columns. The counts are associated to the promoter region of the gene}
#' \item{rpkms}{A matrix of RPKMs, genes as rows, samples as
#' columns. The RPKMs are associated to the promoter region of the gene}
#' \item{normalized}{A matrix of gene abundances normalized by
#' Variance-Stabilizing Transformation (VST), genes as rows, samples as
#' columns. The abundances are associated to the promoter
#' region of the gene}
#' \item{samples}{A vector of sample names and conditions}
#' \item{msviper}{a multisample virtual proteomics object from the
#' viper package}
#' \item{mrs}{A table of master regulators for a specific signature, indicating
#' their Normalized Enrichment Score (NES) and p-value}
#' }
#' @examples
#' library(vulcandata)
#' # Get an example vulcan object (generated with vulcan.import() using the
#' # dummy dataset contained in the \textit{vulcandata} package)
#' vobj<-vulcandata::vulcanexample()
#' # Annotate peaks to gene names
#' vobj<-vulcan.annotate(vobj,lborder=-10000,rborder=10000,method='sum')
#' # Normalize data for VULCAN analysis
#' vobj<-vulcan.normalize(vobj)
#' # Detect which conditions are present
#' names(vobj$samples)
#'
#' # Load an ARACNe network
#' # This is a regulon object as specified in the VIPER package, named 'network'
#' load(system.file('extdata','network.rda',package='vulcandata',mustWork=TRUE))
#' # Run VULCAN
#' # We can reduce the minimum regulon size, since in this example only one
#' # chromosome
#' # was measured, and the networks would otherwise have too few hits
#' vobj_analysis<-vulcan(vobj,network=network,contrast=c('t90','t0'),minsize=5)
#' # Visualize output using the msviper plotting function
#' plot(vobj_analysis$msviper,mrs=7)
#'
#' @export
vulcan <- function(vobj, network, contrast, annotation = NULL,
                minsize = 10) {
    tfs <- names(network)
    samples <- vobj$samples
    normalized <- vobj$normalized

    # Prepare output objects
    msvipers <- matrix(NA, ncol = 3, nrow = length(tfs))
    rownames(msvipers) <- tfs
    # Define contrast
    a <- samples[[contrast[1]]]
    b <- samples[[contrast[2]]]
    # Vulcan msviper implementation
    set.seed(1)
    signature <- viper::rowTtest(normalized[, a], normalized[, b])$statistic
    dnull <- viper::ttestNull(normalized[, a], normalized[, b], per = 1000)
    msviper <- msviper(signature, network,
                    dnull, minsize = minsize)
    # Annotate
    if (!is.null(annotation)) {
        msviper <- viper::msviperAnnot(msviper, annotation)
    }
    vobj$msviper <- msviper
    # Specific Master Regulators
    mrs <- cbind(msviper$es$nes, z2p(msviper$es$nes))
    colnames(mrs) <- c("NES", "pvalue")
    vobj$mrs <- mrs

    return(vobj)
}


#' Function to calculate pathway enrichment over a ChIP-Seq profile
#'
#' This function applies Gene Set Enrichment Analysis or Rank Enrichment
#' Analysis over a ChIP-Seq signature contained in a vulcan package object
#'
#' @param vobj a list, the output of the \code{'vulcan.annotate'} function
#' @param contrast a vector with the name of the two conditions to
#' compare. If method=='REA', contrast can be set to 'all', and Rank Enrichment
#' Analysis will be performed for every sample independently, compared to the
#' mean of the dataset.
#' @param pathways a list of vectors, one vector of gene identifiers per
#' pathway
#' @param method either 'REA' for Rank Enrichment Analysis or 'GSEA' for Gene
#' Set Enrichment Analysis
#' @param np numeric, only for GSEA, the number of permutations to build the
#' null distribution. Default is 1000
#' @return if method=='GSEA', a named vector, with pathway names as names,
#' and the normalized enrichment score of either the GSEA as value.
#' If method=='REA', a matrix, with pathway names as rows and specific
#' contrasts as columns (the method 'REA' allows for multiple contrasts to
#' be calculated at the same time)
#' @examples
#' library(vulcandata)
#' vfile<-tempfile()
#' vulcandata::vulcansheet(vfile)
#' #vobj<-vulcan.import(vfile)
#' vobj<-vulcandata::vulcanexample()
#' unlink(vfile)
#' vobj<-vulcan.annotate(vobj,lborder=-10000,rborder=10000,method='sum')
#' vobj<-vulcan.normalize(vobj)
#' # Create a dummy pathway list (in entrez ids)
#' pathways<-list(
#'    pathwayA=rownames(vobj$normalized)[1:20],
#'    pathwayB=rownames(vobj$normalized)[21:50]
#' )
#' # Which contrast groups can be queried
#' names(vobj$samples)
#' results_gsea<-vulcan.pathways(vobj,pathways,contrast=c('t90','t0'),
#' method='GSEA')
#' results_rea<-vulcan.pathways(vobj,pathways,contrast=c('all'),method='REA')
#' @export
vulcan.pathways <- function(vobj, pathways,
                            contrast = NULL,
                            method = c("GSEA", "REA"), np=1000) {
    normalized <- vobj$normalized
    samples <- vobj$samples
    allgenes <- unique(unlist(pathways))

    # Specific contrast
    if (!setequal(contrast, "all")) {
        # Define contrast
        a <- samples[[contrast[1]]]
        b <- samples[[contrast[2]]]

        # Prepare signature
        signature <- viper::rowTtest(normalized[, a], normalized[, b])$statistic
        if (is.matrix(signature)) {
            signature <- signature[, 1]
        }
        othergenes <- setdiff(allgenes, names(signature))
        gaussiannoise <- setNames(rnorm(length(othergenes),
                                        mean = 0, sd = 0.01),
                                othergenes)  # very small
        signature <- c(signature, gaussiannoise)

        # GSEA
        if (method == "GSEA") {
            gsea.pathways <- setNames(rep(0, length(pathways)), names(pathways))
            message("Running GSEA for ",
                    length(pathways), " pathways")
            pb <- txtProgressBar(0, length(pathways),
                                style = 3)
            i <- 0
            for (pname in names(pathways)) {
                p <- pathways[[pname]]
                obj <- gsea(reflist = signature,
                            set = p, method = "pareto",
                            np = np)
                gsea.pathways[pname] <- obj$nes
                setTxtProgressBar(pb, i <- i + 1)
            }
            return(gsea.pathways)
        }

        # REA
        if (method == "REA") {
            rea.pathways <- setNames(rep(0, length(pathways)), names(pathways))
            message("Running REA for ", length(pathways), " pathways")
            rea.pathways <- rea(signatures = signature, groups = pathways,
                                minsize = 1)
            return(rea.pathways)
        }
    } else {
        if (method != "REA") {
            stop("Multiple signatures supported only with method='REA'")
        }
        signatures <- t(scale(t(vobj$normalized)))
        rea.pathways <- rea(signatures = signatures,
                            groups = pathways, minsize = 1)
        return(rea.pathways)
    }
}
