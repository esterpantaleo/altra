#*************************************************************************
#
#     get arguments
#
#************************************************************************
args               <- commandArgs(TRUE)

source_file        <- args[1]
readfile_labels_v  <- unlist(strsplit(args[2],"\\,"))
RL                 <- as.numeric(args[3])
M                  <- as.numeric(args[4])
C                  <- as.numeric(unlist(strsplit(args[5],"\\,")))
listCoordinate     <- as.numeric(unlist(strsplit(args[6],"\\,")))
output_folder      <- args[7]
GenePredIn         <- args[8]
GenePredOut        <- args[9]
ExprFile           <- args[10]
LikelihoodFile     <- args[11]
JunctionFile       <- args[12]
breakpoints_file   <- args[13]
coverages_file     <- args[14]
ListLociFile       <- args[15]
ListLociLine       <- args[16]

source(source_file)

#*************************************************************************
#
#     set variables
#
#************************************************************************
N                  <- length(readfile_labels_v)
locusStart         <- listCoordinate[1] + 0.5;
locusEnd           <- listCoordinate[length(listCoordinate)]-0.5
individualExprOut  <- as.matrix(readTableNumbers(ExprFile),byrow=TRUE)
colExp             <- ncol(individualExprOut)
individualExprOut  <- as.matrix(individualExprOut[,seq(1, colExp, 2)],nrow=N,ncol=colExp/2)
NORM               <- sum(C) / N
individualExprOut  <- individualExprOut * C / NORM 
if (N==1) ExprAverage <- individualExprOut else ExprAverage  <- colMeans(individualExprOut)


#*************************************************************************
#
#     plot GenePred file and region identified as differentially expressed
#
#************************************************************************
pdf(file=paste(output_folder, "GenePredOut.pdf", sep=""))
#pdf(file=paste(output_folder, "GenePredOut.pdf", sep=""), height=50, width=150)
par(mar=c(3,3,3,3), mfrow=c(2,1), oma=c(3,1,1,1))

plotGenePred(GenePredOut=GenePredOut, GenePredIn=GenePredIn, plotStart=locusStart, plotEnd=locusEnd, Expr_v=ExprAverage)
plotGenePred(GenePredOut=GenePredIn, GenePredIn=GenePredIn, plotStart=locusStart, plotEnd=locusEnd) # listCoordinate=listCoordinate, cex.axis=3)
if (file.exists(ListLociFile))
    plotListLociLine(locusStart, locusEnd, ListLociFile, ListLociLine)
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
	 positions     <- matrix(as.numeric(as.matrix(junctions_tmp[,2:3])), nr , 2)
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
         my_title <- readfile_labels_v[i]
	 plotJunctions(junctions, strands, positions, i, locusStart, locusEnd)
	 plotBreakpoints(breakpoints_file, locusStart, locusEnd)
	 plotCoverage(GenePredOut, individualExprOut[i,], locusStart, locusEnd, RL, M, my_title, coverages[i,])
   }
}
dev.off()

