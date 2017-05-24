#' FUnction to generate an import sheet CSV file for the vulcan test set
#'
#' This function generates a user-specific example sheet with appropriate paths
#' on a dummy ChIP-Seq dataset
#'
#' @param outfile Path and file name for the CSV sheet file
#' @return nothing
#' @examples
#' vfile<-"deleteme.csv"
#' vulcansheet(vfile)
#' tmp<-read.csv(vfile)
#' head(tmp)
#' unlink(vfile)
#' @export
vulcansheet<-function(outfile="sheetfile.csv"){
    header<-c("SampleID","Tissue","Factor","Condition","Replicate","bamReads",
              "Peaks","PeakCaller")
    row1<-c("T90_1","MCF7","ER","t90","1",
            system.file("extdata","A90_I.bam",package="vulcan",mustWork=TRUE),
            system.file("extdata","A90_peaks.bed.gz",package="vulcan",
                        mustWork=TRUE),"bed"
    )
    row2<-c("T90_2","MCF7","ER","t90","2",
            system.file("extdata","B90_I.bam",package="vulcan",mustWork=TRUE),
            system.file("extdata","B90_peaks.bed.gz",package="vulcan",
                        mustWork=TRUE),"bed"
    )
    row3<-c("T90_3","MCF7","ER","t90","3",
            system.file("extdata","C90_I.bam",package="vulcan",mustWork=TRUE),
            system.file("extdata","C90_peaks.bed.gz",package="vulcan",
                        mustWork=TRUE),"bed"
    )
    row4<-c("T90_4","MCF7","ER","t90","4",
            system.file("extdata","D90_I.bam",package="vulcan",mustWork=TRUE),
            system.file("extdata","D90_peaks.bed.gz",package="vulcan",
                        mustWork=TRUE),"bed"
    )
    row5<-c("T45_1","MCF7","ER","t45","1",
            system.file("extdata","A45_I.bam",package="vulcan",mustWork=TRUE),
            system.file("extdata","A45_peaks.bed.gz",package="vulcan",
                        mustWork=TRUE),"bed"
    )
    row6<-c("T45_2","MCF7","ER","t45","2",
            system.file("extdata","B45_I.bam",package="vulcan",mustWork=TRUE),
            system.file("extdata","B45_peaks.bed.gz",package="vulcan",
                        mustWork=TRUE),"bed"
    )
    row7<-c("T45_3","MCF7","ER","t45","3",
            system.file("extdata","C45_I.bam",package="vulcan",mustWork=TRUE),
            system.file("extdata","C45_peaks.bed.gz",package="vulcan",
                        mustWork=TRUE),"bed"
    )
    row8<-c("T45_4","MCF7","ER","t45","4",
            system.file("extdata","D45_I.bam",package="vulcan",mustWork=TRUE),
            system.file("extdata","D45_peaks.bed.gz",package="vulcan",
                        mustWork=TRUE),"bed"
    )
    row9<-c("T0_1","MCF7","ER","t0","1",
            system.file("extdata","A0_I.bam",package="vulcan",mustWork=TRUE),
            system.file("extdata","A0_peaks.bed.gz",package="vulcan",
                        mustWork=TRUE),"bed"
    )
    row10<-c("T0_2","MCF7","ER","t0","2",
             system.file("extdata","B0_I.bam",package="vulcan",mustWork=TRUE),
             system.file("extdata","B0_peaks.bed.gz",package="vulcan",
                         mustWork=TRUE),"bed"
    )
    row11<-c("T0_3","MCF7","ER","t0","3",
             system.file("extdata","C0_I.bam",package="vulcan",mustWork=TRUE),
             system.file("extdata","C0_peaks.bed.gz",package="vulcan",
                         mustWork=TRUE),"bed"
    )
    row12<-c("T0_4","MCF7","ER","t0","4",
             system.file("extdata","D0_I.bam",package="vulcan",mustWork=TRUE),
             system.file("extdata","D0_peaks.bed.gz",package="vulcan",
                         mustWork=TRUE),"bed"
    )
    sheetmat<-rbind(row1,row2,row3,row4,row5,row6,row7,row8,row9,
                    row10,row11,row12)
    colnames(sheetmat)<-header
    write.csv(sheetmat,file=outfile,row.names=FALSE)
}

