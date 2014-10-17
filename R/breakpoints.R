if (!"FLLat" %in% installed.packages())
   stop("R package FLLat is required; istall FLLat using install.packages(\"FLLat\")\n")

require("FLLat")

#===========================================================================
# get arguments
#===========================================================================
args                        <- commandArgs(TRUE)
MAX_NROW                    <- 15000


chr                         <- args[1]
locusStart                  <- as.numeric(args[2])
coverages_file              <- args[3]
is_pos                      <- as.numeric(args[4])
is_neg                      <- as.numeric(args[5])
MIN_EX_LEN                  <- as.numeric(args[6])
SSfile                      <- args[7]
SSfile_breakpoints          <- args[8]
SSfile_breakpoints_filtered <- args[9]
C                           <- as.numeric(unlist(strsplit(args[10],"\\,")))
DELTA_J                     <- as.numeric(args[11])

#===========================================================================
# load coverages_file
#==========================================================================
NORM            <- max(C)/C
h_matrix        <- t(as.matrix(data.frame(lapply(read.table(coverages_file, 
                                                           comment.char="", 
                                                           header=FALSE), 
                                                           as.numeric), 
                                                stringsAsFactors=FALSE))* NORM)

q               <- quantile(h_matrix[h_matrix!=0], probs = seq(0, 1, 0.05))[[20]]
l1              <- q/2
l2              <- 1000/q #1000 #100
nr              <- nrow(h_matrix)
start           <- 1
last_Beta       <- NULL
my_coordinate3  <- locusStart
my_coordinate5  <- locusStart

for (j in 1:ceiling(nr/MAX_NROW)){ 
    h           <- as.matrix(h_matrix[start:min(nr, start + MAX_NROW - 1),])
    trial       <- try(result  <- FLLat(h, J=1, B="pc", lam1=l1, lam2=l2), silent=TRUE)
    if (class(trial) != "try-error"){ #if there is no error
    	my_Beta <- c(last_Beta, result$Beta)
    	if (j==1){
       	   pdf(file=paste(SSfile,".pdf", sep=""));
       	   plot(result)
           dev.off()
    	}

    	if (length(my_Beta) > 1){
           for (i in 1:(length(my_Beta) - 1)){
       	       if (my_Beta[i] == 0 && my_Beta[i + 1] != 0)
               	  my_coordinate3 <- locusStart + start + i - 1
               else
		  if (my_Beta[i + 1] == 0 && my_Beta[i] != 0){
              	     my_coordinate5 <- locusStart + start + i - 1;	
	      	     if (my_coordinate5 - my_coordinate3 > MIN_EX_LEN){
	      	     	if (is_pos == 1){
                    	   write(paste(chr, 
                                       "+", 
                                       "3SS", 
                                       my_coordinate3,   
                                       "2", 
                                       sep=" "), 
                                 file=SSfile_breakpoints, 
                                 append=TRUE)
                           write(paste(chr, 
                                       "+", 
                                       "5SS", 
                                       my_coordinate5, 
                                       "2", 
                                       sep=" "), 
                                 file=SSfile_breakpoints, 
                                 append=TRUE)
                        }
                        if (is_neg == 1){
                           write(paste(chr, 
                                       "-", 
                                       "3SS", 
                                       my_coordinate3, 
                                       "2", 
                                       sep=" "), 
                                 file=SSfile_breakpoints, 
                                 append=TRUE)
  		           write(paste(chr, 
                                       "-", 
                                       "5SS", 
                                       my_coordinate5, 
                                       "2", 
                                       sep=" "), 
                                 file=SSfile_breakpoints, 
                                 append=TRUE)
                        }
                     }
                  }
           }
           last_Beta <- my_Beta[i + 1]
        }   
     }else{ #if FLLat returns error
        print(paste("FLLat returned an error; altra will not use breakpoints in region", j))
        last_Beta <- NULL
    }
     start     <- start + MAX_NROW
}


if (file.exists(SSfile_breakpoints) & file.info(SSfile_breakpoints)$size != 0){
   # load SS_file and SSfile_breakpoints
   ss_file             <- as.matrix(data.frame(lapply(read.table(SSfile, 
                                                                 comment.char="", 
                                                                 header=FALSE, 
                                                                 fill=TRUE), 
                                                      as.character), 
                                               stringsAsFactors=FALSE))
   ss_file_breakpoints <- as.matrix(data.frame(lapply(read.table(SSfile_breakpoints, 
                                                                 comment.char="", 
                                                                 header=FALSE, 
                                                                 fill=TRUE), 
                                                      as.character), 
                                               stringsAsFactors=FALSE))

   # write SSfile_breakpoints_filtered
   for (kk in 1:nrow(ss_file_breakpoints)){
      ss_type    <- ss_file_breakpoints[kk, 3]
      coordinate <- as.numeric(ss_file_breakpoints[kk, 4])
      is_already <- 0;
      for (j in 1:nrow(ss_file)){
         if (ss_file[j, 3] == ss_type){
            if (abs(coordinate - as.numeric(ss_file[j, 4])) < DELTA_J)
               is_already <- 1;
         }
      }
      if (is_already == 0)
          write(paste(chr, 
                      ss_file_breakpoints[kk, 2], 
                      ss_file_breakpoints[kk, 3], 
                      ss_file_breakpoints[kk, 4], 
                      "2", 
                      sep=" "), 
                 file=SSfile_breakpoints_filtered, 
                 append=TRUE)
   }
}

#Sys.sleep(120)
