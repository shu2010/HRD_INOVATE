# ------------------------------------------------------------------
#          Title: organise.segments_mean.R
#          Description: Computes seg.mean and BAF.mean
#          Input files: 
#          ASCAT segmentation profile output
#          
#          author contact: swetansu@gmail.com
#
# Software Dependency:
#    	   R version >= 3.6, ASCAT 2.4	 
# ------------------------------------------------------------------
##for seg.mean and BAF.mean calc

seg_mean_estimation <- function(input, ascat.output){
  
  seg_ref <- ascat.output$segments
  sampleName <- unique(seg_ref[,1])
  logr <- read.table(input,header=T,stringsAsFactors=F,sep="\t")
  
  logrPcfed <- read.table(paste0(sampleName, ".LogR.PCFed",".txt"),header=F,stringsAsFactors=F,sep="\t") 
  bafPcfed <- read.table(paste0(sampleName, ".BAF.PCFed",".txt"),header=F,stringsAsFactors=F,sep="\t")
  colnames(logrPcfed) <- c("probe","logR")
  colnames(bafPcfed) <- c("probe","BAF")
  pcfed <- merge(logrPcfed,bafPcfed,by.x=1,by.y=1,all.x=T,sort=F)
  pcfedAnn <- merge(logr[,1:3],pcfed,by.x=1,by.y=1,all.y=T,sort=F)
  pcfedAnn <- pcfedAnn[order(match(pcfedAnn[,1],logrPcfed[,1])),]
  
  #ploidy=ploidy,aberrantcellfraction=aberrantcellfraction
  
  
  #seg_ref <- ascat.output$segments
  
  #marker <- ascat.bc$SNPpos
  out.new <- data.frame()
  
  get_chr <- list()
  
  nProbes <- numeric()
  
  segMean <- numeric()
  
  bafMean <- numeric()
  
  chr_var <- unique(ascat.output$segments[,2])
  
  i <- integer()
  
  j <- integer()
  
  for(i in 1:length(chr_var)){
    
    get_chr  <- pcfedAnn[grep(paste0('^',chr_var[i],'$'), pcfedAnn[,2]), ]
    
    for(j in 1:length(seg_ref[,2])){
      
      if(seg_ref[j,2] == chr_var[i]){
        
        tmp <- get_chr[get_chr[,3] <= seg_ref[j,4] & get_chr[,3] >= seg_ref[j,3], ]
        
        nProbes[j] <- nrow(tmp)
        
        segMean[j] <- mean(tmp[,4], na.rm=TRUE)
        
        bafMean[j] <- mean(tmp[,5], na.rm=TRUE)
        
      }
      
    }
    
  }
  
  gc()
  
  out.new <- cbind(ascat.output$segments, "nProbes" = nProbes, "seg.Mean"=segMean, "baf.Mean"=bafMean)
  return(out.new)
}
