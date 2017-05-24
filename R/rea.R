#' REA: Rank EnrichmeNt Analysis
#'
#' REA Calculates enrichment of groups of objects over a vector of values
#' associated to a population of objects
#'
#' @param signatures a named vector, with values as signature values
#' (e.g. logFC) and names as object names (e.g. gene symbols)
#' @param groups a list of vectors of objects (e.g. pathways)
#' @param sweights weights associated to objects in the signature.
#' If NULL (default) all objects are treated according to the signature rank
#' @param gweights weights associated to association strength between each
#' object and each group. If NULL (default) all associations are treated equally
#' @param minsize integer. Minimum size of the groups to be analyzed. Default=1
#'
#' @examples
#' signatures<-setNames(-sort(rnorm(1000)),paste0("gene",1:1000))
#' set1<-paste0("gene",sample(1:200,50))
#' set2<-paste0("gene",sample(1:1000,50))
#' groups<-list(set1=set1,set2=set2)
#' obj<-rea(signatures=signatures,groups=groups)
#' obj
#' @export
rea<-function(
    signatures,
    groups,
    sweights=NULL,
    gweights=NULL,
    minsize=1
){
    ### Remove small groups
    groups<-groups[sapply(groups,length)>=minsize]

    ### Treat single "signature"
    if (is.null(nrow(signatures))){
        signatures <- matrix(signatures, length(signatures), 1,
                             dimnames=list(names(signatures), "sample1"))
    }


    ### Generate dummy signature weights
    if(is.null(sweights)){
        sweights<-matrix(1,nrow=nrow(signatures),ncol=ncol(signatures))
        dimnames(sweights)<-dimnames(signatures)
    }
    if(!identical(dim(signatures), dim(sweights))){
        stop("Signatures and Signature weights must be
             matrices of identical size")
    }

    ### Generate dummy group weights
    if(is.null(gweights)){
        gweights<-relist(rep(1,sum(sapply(groups,length))),skeleton=groups)
    }

    ### Apply weights to group belonging
    wgroups<-gweights
    for(i in 1:length(wgroups)){
        names(wgroups[[i]])<-groups[[i]]
    }

    ### Rank-transform columns
    ranks<-apply(signatures,2,rank,na.last="keep")
    ### Assign a 0 to signature weights where the signature was NA
    sweights[is.na(ranks)]<-0
    ### 0-1 bound ranks
    boundranks<-t(t(ranks)/(colSums(!is.na(signatures))+1))
    # Treat bound ranks as quantiles in a gaussian distribution (0=-Inf, 1=+Inf)
    gaussian <- qnorm(boundranks)
    ### Deal with NAs
    gaussian[is.na(gaussian)]<-0

    ### Apply signature weights to the normalized distribution
    gaussian<-gaussian*sweights


    ### Next, we see how each of the groups are behaving in these normalized
    ### signatures
    ### Create a boolean matrix with ngroup columns and signaturelength rows,
    ### indicating the matches
    matches <- sapply(wgroups, function(group, allElements) {
        hereMatches<-as.integer(allElements%in%names(group))
        names(hereMatches)<-allElements
        # Weigth by group belonging
        weightedMatches<-hereMatches
        weightedMatches[names(group)]<-weightedMatches[names(group)]*group
        return(weightedMatches)
    }, allElements=rownames(gaussian))
    # And then transpose it
    matches<-t(matches)

    # Number of matches per group
    groupmatches <- rowSums(matches)

    # Relative part of the signature that matches
    relativematches<-matches/groupmatches

    # This trick will overweight massively small groups with all their
    # components highly-ranked.
    # Extreme case is with a group with one gene at the top

    # The core linear algebra operation. The true magic of rea
    enrichmentScore <- relativematches %*% gaussian

    # Finally, every enrichment is square-rooted to respect the criterion
    # of normality
    normalizedEnrichmentScore<-enrichmentScore*sqrt(groupmatches)

    # Return output
    return(normalizedEnrichmentScore)
}




