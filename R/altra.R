#*************************************************************************
#
#     get arguments
#
#************************************************************************
#get basename
args               <- commandArgs(trailingOnly = FALSE)
scriptName         <- sub("--file=", "", args[grep("--file=", args)])
scriptBasename     <- dirname(scriptName)
sourceFile         <- file.path(scriptBasename,"utils.R")

args               <- commandArgs(trailingOnly = TRUE)  
locusStart         <- as.numeric(args[1])
locusEnd           <- as.numeric(args[2])
N                  <- as.numeric(args[3])
RL                 <- as.numeric(args[4])
M                  <- as.numeric(args[5])
C                  <- as.numeric(unlist(strsplit(args[6],"\\,")))
outputFolder       <- args[7]
VERBOSE            <- as.numeric(args[8])

GenePredIn         <- file.path(outputFolder,"GenePredInCut")
GenePredOut        <- file.path(outputFolder,"GenePredOut")
ExprFile           <- file.path(outputFolder,"ExprOut")
LikelihoodFile     <- file.path(outputFolder,"llkOut")
JunctionFile       <- file.path(outputFolder,"JunctionOut")
BreakpointFile     <- file.path(outputFolder,"BreakpointOut")
CoverageFile       <- file.path(outputFolder,"Coverage")
source(sourceFile)

#*************************************************************************
#
#     set variables
#
#************************************************************************
individualExprOut  <- as.matrix(readTableNumbers(ExprFile), byrow=TRUE)
colExp             <- ncol(individualExprOut)-2
individualExprOut  <- as.matrix(individualExprOut[,seq(1, colExp, 2)],nrow=N,ncol=colExp/2)
NORM               <- sum(C) / N
individualExprOut  <- individualExprOut * C / NORM 
ifelse(N==1, ExprAverage <- individualExprOut, ExprAverage <- colMeans(individualExprOut))

#*************************************************************************
#
#     plot GenePred file and region identified as differentially expressed
#
#************************************************************************
pdf(file=paste(outputFolder, "GenePredOut.pdf", sep=""))
par(mar=c(3,3,3,3), mfrow=c(2,1), oma=c(3,1,1,1))

plotGenePred(GenePredOut=GenePredOut, GenePredIn=GenePredIn, plotStart=locusStart, plotEnd=locusEnd, Expr_v=ExprAverage)
plotGenePred(GenePredOut=GenePredIn, GenePredIn=GenePredIn, plotStart=locusStart, plotEnd=locusEnd)

dev.off()


#***************************************************************************
#
#     plot log_likelihood
#
#***************************************************************************
if (VERBOSE){
    pdf(file=paste(outputFolder, "llkOut.pdf", sep=""))
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
}


#*****************************************************************************
#
#     plot reads, junctions and predictive
#
#*****************************************************************************
#loading all the junctions for all individuals at once 
if (file.exists(JunctionFile) & file.info(JunctionFile)$size != 0){
    junctions_tmp <- read.table(JunctionFile)
    nc            <- ncol(junctions_tmp)
    nr            <- nrow(junctions_tmp)
    junctions     <- matrix(as.numeric(as.matrix(junctions_tmp[,5:nc])), nr, nc - 4)
    strands       <- junctions_tmp[,4]
    positions     <- matrix(as.numeric(as.matrix(junctions_tmp[,2:3])), nr, 2)
}
#loading coverage file
coverages <- read.table(CoverageFile, header=FALSE, comment.char="")
       

pdf(file=paste0(outputFolder, "Coverage.pdf"))
par(mar=c(1,1,1,1), mfrow=c(9,1), oma=c(1,1,1,1))
if (N==1){
    #plotting coverage
    my_title <- ""
    plotJunctions(junctions, strands, positions, 1, locusStart, locusEnd)
    plotBreakpoints(BreakpointFile, locusStart, locusEnd)
    areaRelativeDifference <- plotCoverage(GenePredOut, individualExprOut, locusStart, locusEnd, RL, M, my_title, unlist(coverages[1], use.names=FALSE))
    write(areaRelativeDifference, append=TRUE, file=paste0(outputFolder, "areaRelativeDifference.txt"))
}else{
    for (i in 1:N){
        print(paste("plotting coverage for individual",i))
        my_title <- i
        plotJunctions(junctions, strands, positions, i, locusStart, locusEnd)
        plotBreakpoints(BreakpointFile, locusStart, locusEnd)
        areaRelativeDifference <- plotCoverage(GenePredOut, individualExprOut[i,], locusStart, locusEnd, RL, M, my_title, unlist(coverages[i,],use.names=FALSE))
        write(areaRelativeDifference, append=TRUE, file=paste0(outputFolder, "areaRelativeDifference_", i, ".txt"))
    }
}
dev.off()

