readTable <- function(GenePredIn){
    if (!file.exists(GenePredIn))
        print(paste("ERROR: file", GenePredIn, "does not exist"))
    return(data.frame(lapply(read.table(GenePredIn,
                                        fill=1,
                                        comment.char="",
                                        header=FALSE),
                             as.character),
                      stringsAsFactors=FALSE))
}

readTableNumbers <- function(ExprOut){
    if (!file.exists(ExprOut))
        print(paste("ERROR: file", GenePredIn, "does not exist"))
    return(as.matrix(data.frame(lapply(read.table(ExprOut,
                                                  fill=1,
                                                  comment.char="",
                                                  header=FALSE),
                                       as.numeric),
                                stringsAsFactors=FALSE)))
}

get_exons_start_end <- function(transcript){	
    exst <- as.numeric(strsplit(as.character(transcript[9]), ",")[[1]]) + 1
    exen <- as.numeric(strsplit(as.character(transcript[10]), ",")[[1]])
    return(list(exst=exst, exen=exen))
}

plotGenePred <- function(GenePredOut, GenePredIn=NULL, plotStart=NULL, plotEnd=NULL, is_xaxis=1, title=NULL, Expr_v=NULL, listCoordinate=NULL, cex.axis=1){
    # get data
    if (!is.null(GenePredIn)){
        if (!file.exists(GenePredIn))
            print(paste("ERROR: file", GenePredIn, "does not exist"))
        InTranscripts  <- readTable(GenePredIn)
        nrIn           <- nrow(InTranscripts)
    }
    if (!file.exists(GenePredOut))
        print(paste("ERROR: file", GenePredOut, "does not exist"))
    if (!file.info(GenePredOut)$size == 0){
        OutTranscripts <- readTable(GenePredOut)
        nrOut          <- nrow(OutTranscripts)
        chr            <- OutTranscripts[1,2]
        if (is.null(plotStart)==1)
            plotStart <- min(as.numeric(OutTranscripts[,4]))
        if (is.null(plotEnd)==1)
            plotEnd   <- max(as.numeric(OutTranscripts[,5]))
        plot("NA",
             axes="F",
             yaxt="n",
             main=title,
             xlim=c(plotStart,plotEnd),
             ylim=c(0,1),
             xlab=paste0("Position (Mb) on ", chr),
             ylab="",
             font.main=1)
        
        #plot axis
        if (is_xaxis==1){
            if (!is.null(listCoordinate)){
                axis(1,
                     at=listCoordinate,
                     lab=0:(length(listCoordinate)-1),
                     cex.axis=cex.axis)
                #par(tck=tck)
                abline(v=listCoordinate)
            }else{
                tck=axTicks(1)
                axis(1,
                     at=tck,
                     lab=format(tck/1000000))
            }
        }
        
        for (i in 1:nrOut){#set colors
            if (OutTranscripts[i,3] == "-") fill_col="red" else if (OutTranscripts[i,3] == "+") fill_col="blue" else  fill_col="black"
            border_col=fill_col;
            
            if (nrOut==1){y=0.5} else y=(nrOut+1-i)*(1/(nrOut+1))
            
            rect(OutTranscripts[i, 4],
                 y-0.001,
                 OutTranscripts[i, 5],
                 y+0.001,
                 col="black",
                 border="black")
            
            if (length(Expr_v) != 0) 
                text(x=OutTranscripts[i,5],
                     y=y+0.07,
                     labels=signif(Expr_v[i],1),
                     cex=0.7,
                     pos=3,
                     offset=0.1)
            trans <- get_exons_start_end(OutTranscripts[i,])	      
            ll    <- length(trans$exst)
            
            if (!is.null(GenePredIn)){
                for (j in 1:nrIn){
                    intrans <- get_exons_start_end(InTranscripts[j,])
                    mm      <- length(intrans$exst)
                    if (ll == mm){
                        if (all(trans$exst == intrans$exst) & all(trans$exen == intrans$exen)){
                            border_col="black";
                            break
                        }
                    }
                }
            }
            for (j in 1:ll) 
                rect(trans$exst[j],
                     y-0.02,
                     trans$exen[j],
                     y+0.02,
                     col=fill_col,
                     border=border_col)
        }
    }
}

#' Given a genePred file and a vector of expressions for each of the transcripts, plot the coverage
pred2coverage <- function(GenePredOut, indExprOut, locusStart, locusEnd, RL, M){
    h           <- rep(0, locusEnd - locusStart + 1)
    if (file.info(GenePredOut)$size == 0) return(h);
    Transcripts <- readTable(GenePredOut)
    nr          <- nrow(Transcripts)
    for (k in 1:nr){
        trans  <- get_exons_start_end(Transcripts[k,])
        E      <- length(trans$exst)
        
###################################################################
#  exst and exen are exon start and exon end in genomic coordinates (1-based)
#  ExSt is exon start in transcriptome coordinates (1-based)
#  cumInLength are cumulative intron lengths 
###################################################################
        cumInLength  <- (cumsum(c(trans$exst, 0) - c(-1, trans$exen) - 1) - trans$exst[1])[1:E]
        ExSt         <- cumsum(c(1, trans$exen - trans$exst + 1))
        
#################################################################
#  hh[x] is the coverage at position x in transcriptome coordinates
#################################################################
        hh          <- rep(0, ExSt[E + 1] - 1)
        Expression  <- as.numeric(indExprOut[k])
        if (Expression == 0) 
            next;
        #for each exon j
        for (j in 1:E){
###############coverage from reads with no junctions
            if (ExSt[j] < ExSt[j + 1] - RL)
                for (x in ExSt[j]:(ExSt[j + 1] - RL))
                    hh[x:(x+RL-1)] <- hh[x:(x+RL-1)] + Expression
###############coverage from reads with junctions
            if (j != E){
                if (ExSt[j + 1] - ExSt[j] > M){
                    for (x in max(ExSt[j], ExSt[j + 1] - RL + M):(ExSt[j + 1] - M - 1)){ 
                        read   <- ExSt[j + 1] - x
                        read_v <- read
                        rest   <- RL - read
                        for (jj in (j + 1):E){
                            read   <- min(rest, ExSt[jj + 1] - ExSt[jj])
                            read_v <- c(read_v, read)
                            rest   <- rest - read
                            if (rest <= 0)
                                break;
                        }
                        if (rest > 0 || any(read_v < M)){
                                        #print(rest)
                                        #print(read_v) 
                            break;
                        }
                        ##else
                        hh[x:(x+RL-1)] <- hh[x:(x+RL-1)] + Expression
                    }
                }
            }
        }
        
#################################################################
#  transform from transcriptome coordinates (x) (hh) to reference coordinates (y) (h)
#################################################################
        for (j in 1:E){
            x <- ExSt[j]:(ExSt[j+1]-1)
            y <- x + cumInLength[j] + trans$exst[1] - locusStart
            h[y] <- h[y] + hh[x]
        }
    } 
    return(h)
}

#' Plot coverage of one sample given number of reads that start at each base (individualCoverages)
plotCov <- function(individualCoverages, locusStart, locusEnd, my_title, ylim){
    l       <- length(individualCoverages)
    plot(locusStart:(locusStart + l - 1),
         individualCoverages,
         lwd=0.1,
         col="black",
         font.main=1,
         main=my_title,
         axes="F",
         type="h",
         xlab = "",
         ylab = "",
         xlim=range(c(locusStart, locusStart + l - 1)),
         ylim=ylim)
    axis(2, pos=locusStart)
}
    
plotCoverage <- function(GenePredOut, individualExprOut, locusStart, locusEnd, RL, M, my_title, individualCoverages){	
    l       <- length(individualCoverages)

    predCov <- pred2coverage(GenePredOut, individualExprOut, locusStart, locusEnd, RL, M)
    ylim    <- c(0,max(c(predCov,individualCoverages)))
    #plot coverage from bam files
    plotCov(individualCoverages, locusStart, locusEnd, my_title, ylim)
    par(new=TRUE)
    #plot predictive
    plot(locusStart:(locusStart + length(predCov) - 1), 
         predCov, 
         lwd=0.5,
         col="red",
         axes="F", 
         type="l", 
         xlab="", 
         ylab="", 
         xlim=range(c(locusStart, locusStart + l - 1)), 
         ylim=ylim)
    
    covDiff <- abs(individualCoverages-predCov[1:length(individualCoverages)])
       
    return(sum(covDiff)/sum(individualCoverages))
}

plotCoverageFromBam <- function(bamfile, chr, locusStart, locusEnd, my_title){
    individualCoverages <- system(paste0("cat ", chr, ":", locusStart, "-", locusEnd, " | coverageBed split -d -hist -abam ", bamfile, "-b stdin | awk '{printf(\"%s\",$5)} END{printf(\"\n\")}'"), intern=TRUE)
    plotCov(individualCoverages, locusStart, locusEnd, my_title)
}


#' plotJunctions
#'
#' Given the JunctionFile (a row in this file looks like "chr12 114455 114567 + 4 4 54 1 18") load counts (from column 5 to ncol) and strand (from col 4).
#' @param junction_counts
#' @param strands
#' @param positions
#' @param i: the individual (the column number in the matrix)
#' @param plotStart
#' @param plotEnd
plotJunctions <- function(junction_counts, strands, positions, i, plotStart, plotEnd){	 
    if (!is.null(junction_counts)){
        nr        <- nrow(junction_counts);
        max_count <- max(junction_counts)
        plot("NA", 
             xlim=c(plotStart, plotEnd), 
             ylim=c(0, 1), 
             axes=F, 
             ylab="", 
             xlab="")
        for (j in 1:nr){
            y <- (nr + 1 - j) / (nr+1)
            
            if (strands[j]=="+")
                clr="blue"
            if (strands[j]=="-")
                clr="red"
            line_width <- junction_counts[j, i]*10/max_count 
            lines(c(as.numeric(positions[j, 1]), as.numeric(positions[j, 2])), 
                  c(y, y),
                  col = clr, 
                  lwd = line_width,
                  lend=1) #line end type (0 is round)
            text(x = as.numeric(positions[j, 1]), 
                 y = y + 0.02, 
                 labels = junction_counts[j, i], 
                 cex=0.5, 
                 pos=2, 
                 offset=0)
        }
    }else plot(1, type="n", axes=F, xlab="", ylab="", xlim = c(plotStart, plotEnd))
}

plotBreakpoints <- function(breakpoints_file, plotStart, plotEnd){
    if (file.exists(breakpoints_file) & file.info(breakpoints_file)$size != 0){
        breakpoints <- readTable(breakpoints_file)
        breakpoints <- subset(breakpoints, !duplicated(breakpoints[,4])) 
        nr <- nrow(breakpoints)
        if (nr > 0){
            plot("NA", 
                 xlim=c(plotStart,plotEnd), 
                 ylim=c(0, 1), 
                 axes=F, 
                 ylab="", 
                 xlab="")
            for (j in 1:nr){
                y          <- (nr + 1 - j)/(nr+1)
                coordinate <- as.numeric(breakpoints[j, 4])
                if (breakpoints[j, 3] == "3SS")
                    arrows(coordinate - 10, y, coordinate, y, length = 0.05) 
                else if (breakpoints[j, 3] == "5SS")
                    arrows(coordinate + 10, y, coordinate, y, length = 0.05)
            }
        }
    }else plot(1, type="n", axes=F, xlab="", ylab="", xlim = c(plotStart, plotEnd))
}
    
