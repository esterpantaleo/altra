#*************************************************************************
#
#     get arguments
#
#************************************************************************
args               <- commandArgs(TRUE)

source_file        <- args[1]
locusStart         <- as.numeric(args[2])
locusEnd           <- as.numeric(args[3])
N                  <- as.numeric(args[4])
RL                 <- as.numeric(args[5])
M                  <- as.numeric(args[6])
C                  <- as.numeric(unlist(strsplit(args[7],"\\,")))
output_folder      <- args[8]
GenePredIn         <- args[9]
GenePredOut        <- args[10]
ExprFile           <- args[11]
LikelihoodFile     <- args[12]
JunctionFile       <- args[13]
breakpoints_file   <- args[14]
coverages_file     <- args[15]
source(source_file)

#*************************************************************************
#
#     set variables
#
#************************************************************************
individualExprOut  <- as.matrix(readTableNumbers(ExprFile), byrow=TRUE)
colExp             <- ncol(individualExprOut)
individualExprOut  <- as.matrix(individualExprOut[,seq(1, colExp, 2)],nrow=N,ncol=colExp/2)
NORM               <- sum(C) / N
individualExprOut  <- individualExprOut * C / NORM 
ifelse(N==1, ExprAverage <- individualExprOut, ExprAverage <- colMeans(individualExprOut))

#*************************************************************************
#
#     plot GenePred file and region identified as differentially expressed
#
#************************************************************************
pdf(file=paste(output_folder, "GenePredOut.pdf", sep=""))
#pdf(file=paste(output_folder, "GenePredOut.pdf", sep=""), height=50, width=150)
par(mar=c(3,3,3,3), mfrow=c(2,1), oma=c(3,1,1,1))

plotGenePred(GenePredOut=GenePredOut, GenePredIn=GenePredIn, plotStart=locusStart, plotEnd=locusEnd, Expr_v=ExprAverage)
plotGenePred(GenePredOut=GenePredIn, GenePredIn=GenePredIn, plotStart=locusStart, plotEnd=locusEnd)
#plotCoordinates(listCoordinate, plotStart=locusStart, plotEnd=locusEnd, cex.axis=2)
dev.off()


#***************************************************************************
#
#     plot log_likelihood
#
#***************************************************************************
pdf(file=paste(output_folder, "plot_log_likelihood.pdf", sep=""))
llk <- as.numeric(as.matrix(read.table(LikelihoodFile, 
                                       blank.lines.skip=TRUE, 
                                       comment.char="", 
                                       fill=TRUE, 
                                       header=FALSE))[,1])
m = min(max(llk), 100000000)
plot(1:length(llk), 
             llk, 
             col="black", 
             ylim=c(-1000 + m, m), 
             type='l')
dev.off()


#*****************************************************************************
#
#     plot reads, junctions and predictive
#
#*****************************************************************************
junctions <- NULL
strands   <- NULL
positions <- NULL
#I'm loading all the junctions for all individuals at once 
if (file.exists(JunctionFile) & file.info(JunctionFile)$size != 0){
    junctions_tmp <- read.table(JunctionFile)
    nc            <- ncol(junctions_tmp)
    nr            <- nrow(junctions_tmp)
    junctions     <- matrix(as.numeric(as.matrix(junctions_tmp[,5:nc])), nr, nc - 4)
    strands       <- junctions_tmp[,4]
    positions     <- matrix(as.numeric(as.matrix(junctions_tmp[,2:3])), nr, 2)
}
coverages <- data.frame(lapply(read.table(coverages_file, 
                                          header=FALSE, 
                                          comment.char=""), 
                               as.numeric), 
                        stringsAsFactors=FALSE)

pdf(file=paste(output_folder, "plot_reads.pdf", sep=""))
par(mar=c(1,1,1,1), mfrow=c(9,1), oma=c(1,1,1,1))
if (N==1){
    my_title <- ""
    plotJunctions(junctions, strands, positions, 1, locusStart, locusEnd)
    plotBreakpoints(breakpoints_file, locusStart, locusEnd)
    plotCoverage(GenePredOut, individualExprOut, locusStart, locusEnd, RL, M, my_title, coverages)
}else{
    for (i in 1:N){
        my_title <- i
        plotJunctions(junctions, strands, positions, i, locusStart, locusEnd)
        plotBreakpoints(breakpoints_file, locusStart, locusEnd)
        plotCoverage(GenePredOut, individualExprOut[i,], locusStart, locusEnd, RL, M, my_title, coverages[i,])
    }
}
dev.off()

