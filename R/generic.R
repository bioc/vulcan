#' z2p
#'
#' This function gives a gaussian p-value corresponding to the provided Z-score
#'
#' @param z a Z score
#' @return a p-value
#' @examples
#' z<-1.96
#' z2p(z)
#' @export
z2p <- function(z) {
    pnorm(abs(z), lower.tail = FALSE) * 2
}

#' p2z
#'
#' This function gives a gaussian Z-score corresponding to the provided p-value
#' Careful: sign is not provided
#'
#' @param p a p-value
#' @return z a Z score
#' @examples
#' p<-0.05
#' p2z(p)
#' @export
p2z <- function(p) {
    qnorm(p/2, lower.tail = FALSE)
}

#' Stouffer integration of Z scores
#'
#' This function gives a gaussian Z-score corresponding to the provided p-value
#' Careful: sign is not provided
#'
#' @param x a vector of Z scores
#' @return Z an integrated Z score
#' @examples
#' zs<-c(1,3,5,2,3)
#' stouffer(zs)
#' @export
stouffer <- function(x) {
    Z <- sum(x)/sqrt(length(x))
    return(Z)
}

#' Weighted Stouffer integration of Z scores
#'
#' This function gives a gaussian Z-score corresponding to the provided p-value
#' Careful: sign is not provided
#'
#' @param x a vector of Z scores
#' @param w weight for each Z score
#' @return Z an integrated Z score
#' @examples
#' zs<-c(1,-3,5,2,3)
#' ws<-c(1,10,1,2,1)
#' wstouffer(zs,ws)
#' @export
wstouffer <- function(x, w) {
    Z <- sum(x * w)/sqrt(sum(w^2))
    return(Z)
}


#' Fisher integration of p-values
#'
#' This function applies the Fisher integration of pvalues
#'
#' @param ps a vector of p-values
#' @return p.val an integrated p-value
#' @examples
#' ps<-c(0.01,0.05,0.03,0.2)
#' fisherp(ps)
#' @export
fisherp <- function(ps) {
    Xsq <- -2 * sum(log(ps))
    p.val <- pchisq(Xsq, df = 2 * length(ps),
    lower.tail = FALSE)
    # p<-c(Xsq = Xsq, p.value = p.val)
    return(p.val)
}


#' Slice
#'
#' This function prints a slice of a matrix
#'
#' @param matrix A matrix
#' @return prints it
#' @examples
#' set.seed(1)
#' example<-matrix(rnorm(1000),nrow=100,ncol=10)
#' slice(example)
#' @export
slice <- function(matrix) {
    if (nrow(matrix) < 5) {
    stop("Input matrix has less than 5 rows")
    }
    if (ncol(matrix) < 5) {
    stop("Input matrix has less than 5 columns")
    }
    print(matrix[1:5, 1:5])
}



#' Convert correlation coefficient to p-value
#'
#' This functions converts an R value from a correlation calculation into a
#' p-value, using a T distribution with the provided number of samples N minus
#' 2 degrees of freedom
#' @param r a correlation coefficient
#' @param N a number of samples
#' @return p a p-value
#' @examples
#' set.seed(1)
#' a<-rnorm(1000)
#' b<-a+rnorm(1000,sd=10)
#' r<-cor(a,b,method='pearson')
#' corr2p(r,N=length(a))
#' @export
corr2p <- function(r, N) {
    # Get the t-distribution value
    t <- r/sqrt((1 - r^2)/(N - 2))
    # t follows the t distribution with
    # df=N?2
    p <- 1 - pt(abs(t), N - 2)
    return(p)
}

#' Convert p-value to correlation coefficient
#'
#' This functions converts an p-value into a the corresponding correlation
#' coefficient, using a T distribution with the provided number of samples N
#' minus 2 degrees of freedom
#' @param p a p-value
#' @param N a number of samples
#' @return r a correlation coefficient
#' @examples
#' N<-100
#' p<-0.05
#' p2corr(p,N)
#' @export
p2corr <- function(p, N) {
    # Get the t-value
    t <- abs(qt(p, N - 2))

    # Convert to correlation
    r <- sqrt(t^2/(N - 2 + t^2))

    return(r)
}


#' Convert a numeric vector into colors
#' @param z a vector of numbers
#' @param col1 a color name for the min value, default 'navy'
#' @param col2 a color name for the middle value, default 'white'
#' @param col3 a color name for the max value, default 'red3'
#' @param nbreaks Number of colors to be generated. Default is 30.
#' @param center boolean, should the data be centered? Default is TRUE
#' @param rank boolean, should the data be ranked? Default is FALSE
#' @return a vector of colors
#' @examples
#' a<-rnorm(1000)
#' cols<-val2col(a)
#' plot(a,col=cols,pch=16)
#' @export
val2col <- function(z, col1 = "navy", col2 = "white",
    col3 = "red3", nbreaks = 100, center = TRUE,
    rank = FALSE) {
    isMatrix <- FALSE
    if (is.matrix(z)) {
    isMatrix <- TRUE
    oriz <- z
    }
    if (is.character(z)) {
    z <- as.numeric(as.factor(z))
    }
    if (rank) {
    z <- rank(z)
    }
    if (center) {
    extreme = round(max(abs(z)))
    breaks <- seq(-extreme, extreme,
    length = nbreaks)
    z <- z - mean(z)
    } else {
    breaks <- seq(min(z), max(z), length = nbreaks)
    }
    ncol <- length(breaks) - 1
    col <- colorpanel(ncol, col1, col2, col3)
    CUT <- cut(z, breaks = breaks)
    # assign colors to heights for each point
    colorlevels <- col[match(CUT, levels(CUT))]
    names(colorlevels) <- names(z)
    if (isMatrix) {
    colormatrix <- matrix(colorlevels,
    ncol = ncol(oriz), nrow = nrow(oriz))
    dimnames(colormatrix) <- dimnames(oriz)
    return(colormatrix)
    }
    return(colorlevels)
}

#' kmgformat - Nice Formatting of Numbers
#'
#' This function will convert thousand numbers to K, millions to M, billions
#' to G, trillions to T, quadrillions to P
#'
#' @param input A vector of values
#' @param roundParam How many decimal digits you want
#' @return A character vector of formatted numebr names
#' @examples
#' # Thousands
#' set.seed(1)
#' a<-runif(1000,0,1e4)
#' plot(a,yaxt='n')
#' kmg<-kmgformat(pretty(a))
#' axis(2,at=pretty(a),labels=kmg)
#'
#' # Millions to Billions
#' set.seed(1)
#' a<-runif(1000,0,1e9)
#' plot(a,yaxt='n',pch=20,col=val2col(a))
#' kmg<-kmgformat(pretty(a))
#' axis(2,at=pretty(a),labels=kmg)
#' @export
kmgformat <- function(input, roundParam = 1) {
    signs <- sign(input)
    signs[signs == 1] <- ""
    signs[signs == -1] <- "-"
    absinput <- abs(input)
    output <- c()
    for (i in absinput) {
    if (i < 1000) {
    output <- c(output, i)
    } else if (i < 1e+06) {
    i <- round(i/1000, roundParam)
    i <- paste0(i, "K")
    output <- c(output, i)
    } else if (i < 1e+09) {
    i <- round(i/1e+06, roundParam)
    i <- paste0(i, "M")
    output <- c(output, i)
    } else if (i < 1e+12) {
    i <- round(i/1e+09, roundParam)
    i <- paste0(i, "G")
    output <- c(output, i)
    } else if (i < 1e+15) {
    i <- round(i/1e+12, roundParam)
    i <- paste0(i, "T")
    output <- c(output, i)
    } else if (i < 1e+18) {
    i <- round(i/1e+15, roundParam)
    i <- paste0(i, "P")
    output <- c(output, i)
    } else {
    output <- c(output, i)
    }
    }
    output <- paste0(signs, output)
    return(output)
}



#' densityauc - Calculate the AUC of a density object
#'
#' This function will calculate the AUC of a density object generated by the
#' \code{'density'} function.
#'
#' @param dens a density object
#' @param window a vector with two values, specifying the left and right borders
#' for the AUC to be calculated
#' @return a numeric value for the density AUC
#' @examples
#' set.seed(1)
#' a<-rnorm(1000)
#' d<-density(a)
#' window<-c(2,3)
#' da<-densityauc(d,window)
#'
#' plot(d,main='')
#' abline(v=window,lty=2)
#' title(paste0('AUC between lines=',da))
#'
#'
#'
#' @export
densityauc <- function(dens, window) {
    xt <- diff(dens$x[dens$x > window[1] &
    dens$x < window[2]])
    yt <- rollmean(dens$y[dens$x > window[1] &
    dens$x < window[2]], 2)
    sum(xt * yt)
}



#' textplot2 - An x y plot of non-overlapping text
#'
#' This function is an extension of the \code{'textplot'} function from the
#' \code{'wordcloud'} package, with the extra functionality of specifiying the
#' color of the points
#' @param x x coordinates
#' @param y y coordinates
#' @param words the text to plot
#' @param cex font size
#' @param new should a new plot be created
#' @param pointcolor a string specifying the color of the points
#' (default #FFFFFF00)
#' @param show.lines if true, then lines are plotted between x,y and the word,
#' for those words not covering their x,y coordinates
#' @param pch pch parameter for the plotted points
#' @param ...    Additional parameters to be passed to wordlayout and text.
#' @return nothing
#' @examples
#' obj_names<-apply(expand.grid(LETTERS,LETTERS),1,paste,collapse='')
#' a<-setNames(runif(26*26),obj_names)
#' b<-setNames(rnorm(26*26),obj_names)
#' plot(a,b,pch=20,col='grey')
#' top<-names(sort(-a))[1:50]
#' textplot2(a[top],b[top],words=top,new=FALSE,pointcolor='black')
#' @export
textplot2 <- function(x, y, words, cex = 1,
    pch = 16, pointcolor = "#FFFFFF00", new = TRUE,
    show.lines = TRUE, ...) {
    if (new) {
    plot(x, y, type = "n", ...)
    }
    lay <- wordlayout(x, y, words, cex, ...)
    if (show.lines) {
    for (i in seq_len(length(x))) {
    xl <- lay[i, 1]
    yl <- lay[i, 2]
    w <- lay[i, 3]
    h <- lay[i, 4]
    if (x[i] < xl || x[i] > xl +
    w || y[i] < yl || y[i] >
    yl + h) {
    points(x[i], y[i], pch = pch,
    col = pointcolor, cex = 0.5)
    nx <- xl + 0.5 * w
    ny <- yl + 0.5 * h
    lines(c(x[i], nx), c(y[i],
    ny), col = "grey")
    }
    }
    }
    text(lay[, 1] + 0.5 * lay[, 3], lay[,
    2] + 0.5 * lay[, 4], words, cex = cex,
    ...)
}

#' Define the average fragment length
#'
#' A function to get average fragment length from a ChIP-Seq experiment
#' @param bam.files a vector of BAM files locations
#' @param plot logical. Should a plot be generated?
#' @param max.dist numeric. Maximum fragment length accepted. Default=550
#' @return nothing
#' @examples
#' library(vulcandata)
#' sheetfile<-'deleteme.csv'
#' vulcandata::vulcansheet(sheetfile)
#' a<-read.csv(sheetfile,as.is=TRUE)
#' bams<-a$bamReads
#' unlink(sheetfile)
#' average_fragment_length(bams,plot=TRUE)
#' @export
average_fragment_length <- function(bam.files,
    plot = TRUE, max.dist = 550) {
    # obtain average fragment length from
    # cross-correlation of plus/minus strand
    # alignments
    x <- csaw::correlateReads(bam.files,
    max.dist = max.dist)
    # visualize (raw and smoothed)
    xs <- caTools::runmean(Rle(x), k = 101,
    endrule = "constant")
    frag.len <- which.max(xs)
    if (plot) {
    plot(0:max.dist, x, pch = "*", ylab = "CCF",
    xlab = "Delay (bp)")
    lines(0:max.dist, xs, col = "red",
    lwd = 3)
    abline(v = frag.len, col = "blue")
    legend("topright", paste0("fragment-length = ",
    frag.len, "bp"), bty = "n")
    }
    return(frag.len)
}

