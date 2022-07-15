# ------------------------------------------------------------------
#          Title: calc_lst_new2.R
#          Description: Compute LST, Popova et al. (Cancer Res. 2012, 72: 5454-5462. 10.1158/0008-5472.CAN-12-1470.)
#          Input files: 
#          ASCAT output
#          
#          author contact: swetansu@gmail.com
#
# Software Dependency:
#    	   R version >= 3.6, ASCAT 2.4	 
# ------------------------------------------------------------------
calc_lst_new <- function(
  seg,                # Matrix of segmented output from ASCAT
  chrominfo,          # Matrix with information about the chromosomes
  nA = 7,             # The column index of 'seg' where copy number of A
  #   allele is found
  check.names = T,    # Rename any duplicated samples
  min.probes = 50,    # As described by Popova (50 for SNP6)  
  return.loc = FALSE, # Return location of LST sites
  chr.arm = 'no'      # Option to use chromosome arms defined during
  #   segmentation. The option must give a column
  #   that holds the chromosome arm information
  #   (or 'no' to not use this)
){
  seg[,1] <- as.character(seg[,1])
  if(check.names){
    seg <- check.names.fn(seg)
  }
  if(! all(seg[,8] <= seg[,7]) ){
    # In case ASCAT people change the algorithm  ### They did!!
    cat("Warning!! nB  not always <= nA!!  -- Correcting for internal use (only!)\n") 
    tmp <- seg
    seg[tmp[,8] > tmp[,7],7]  <- tmp[tmp[,8] > tmp[,7],8]
    seg[tmp[,8] > tmp[,7],8]  <- tmp[tmp[,8] > tmp[,7],7]
  }
  seg <- seg[seg[,5] >= min.probes,]
  nB <- nA+1
  samples <- unique(seg[,1])
  output <- setNames(rep(0,length(samples)), samples)
  if(return.loc) {
    out.seg <- matrix(0,0,10)
    colnames(out.seg) <- c(colnames(seg)[1:8],'LST breakpoint', 'chr. arm')
  }
  for(j in samples){
    sample.seg <- seg[seg[,1] %in% j,]
    sample.lst <- c()
    for(i in unique(sample.seg[,2])){
      sample.chrom.seg <- sample.seg[sample.seg[,2] %in% i,,drop=F]
      if(chr.arm !='no'){
        p.max <- if(any(sample.chrom.seg[,chr.arm] == 'p')){max(sample.chrom.seg[sample.chrom.seg[,chr.arm] == 'p',4])}
        q.min <- min(sample.chrom.seg[sample.chrom.seg[,chr.arm] == 'q',3])
      }
      if(nrow(sample.chrom.seg) < 2) {next}
      sample.chrom.seg.new <- sample.chrom.seg
      if(chr.arm == 'no'){
        # split into chromosome arms
        p.arm <- sample.chrom.seg.new[sample.chrom.seg.new[,3] <= chrominfo[i,4],] 
        q.arm <- sample.chrom.seg.new[sample.chrom.seg.new[,4] >= chrominfo[i,5],]
        q.arm<- tryCatch(shrink.seg.ai(q.arm), error=function(x) return(NULL))
        p.arm<- tryCatch(shrink.seg.ai(p.arm), error=function(x) return(NULL))
        if(!is.null(p.arm)){
          p.arm[nrow(p.arm),4] <- chrominfo[i,4]
          q.arm[1,3] <- chrominfo[i,5]
        }
      }
      if(chr.arm != 'no'){
        q.arm <- sample.chrom.seg.new[sample.chrom.seg.new[,chr.arm] == 'q',,drop=F]
        q.arm<- shrink.seg.ai(q.arm)
        q.arm[1,3] <- q.min
        if(any(sample.chrom.seg.new[,chr.arm] == 'p')){
          # split into chromosome arms
          p.arm <- sample.chrom.seg.new[sample.chrom.seg.new[,chr.arm] == 'p',,drop=F] 
          p.arm<- shrink.seg.ai(p.arm)
          p.arm[nrow(p.arm),4] <- p.max
        }
      }  
      if(!is.null(p.arm)){
        n.3mb <- which((p.arm[,4] - p.arm[,3]) < 3e6)
        while(length(n.3mb) > 0){
          p.arm <- p.arm[-(n.3mb[1]),]
          p.arm <- tryCatch(shrink.seg.ai(p.arm), error=function(x) return(NULL))
          #   if(!is.null(p.arm)){
          n.3mb <- which((p.arm[,4] - p.arm[,3]) < 3e6)
          #   }
        }
        if(!is.null(p.arm) && nrow(p.arm) >= 2){
          p.arm <- cbind(p.arm[,1:8], c(0,1)[match((p.arm[,4]-p.arm[,3]) >= 10e6, c('FALSE','TRUE'))])
          #        print("suspicious p.arm")
          for(k in 2:nrow(p.arm)){
            #          print("suspicious p.arm2")
            if(p.arm[k,9] == 1 & p.arm[(k-1),9]==1 & (p.arm[k,3]-p.arm[(k-1),4]) < 3e6){
              sample.lst <- c(sample.lst, 1)
              if(return.loc){
                print("suspicious p.arm3")
                ## Number indicates if start (1) or end (2) defines the breakpoint
                a<- cbind(p.arm[(k-1),1:8], 2,'p-arm') 
                b <- cbind(p.arm[k,1:8], 1,'p-arm')
                colnames(a)[9:10]<- colnames(b)[9:10]<- c('LST breakpoint', 'chr. arm')
                out.seg <- rbind(out.seg, a,b)
              }
            }
          }
        }
      }
      
      if(!is.null(q.arm)){
        n.3mb <- which((q.arm[,4] - q.arm[,3]) < 3e6)
        while(length(n.3mb) > 0){
          q.arm <- q.arm[-(n.3mb[1]),]
          q.arm <- tryCatch(shrink.seg.ai(q.arm), error=function(x) return(NULL))
          n.3mb <- which((q.arm[,4] - q.arm[,3]) < 3e6)
        }
        if(!is.null(p.arm) && nrow(q.arm) >= 2){
          q.arm <- cbind(q.arm[,1:8], c(0,1)[match((q.arm[,4]-q.arm[,3]) >= 10e6, c('FALSE','TRUE'))])
          for(k in 2:nrow(q.arm)){
            if(q.arm[k,9] == 1 & q.arm[(k-1),9]==1 & (q.arm[k,3]-q.arm[(k-1),4]) < 3e6){
              sample.lst <- c(sample.lst, 1)
              if(return.loc){
                ## Number indicates if start (1) or end (2) defines the breakpoint
                a<- cbind(q.arm[(k-1),1:8], 2,'q-arm') 
                b <- cbind(q.arm[k,1:8], 1,'q-arm')
                colnames(a)[9:10]<- colnames(b)[9:10]<- c('LST breakpoint', 'chr. arm')
                out.seg <- rbind(out.seg, a,b)
              }
            }
          }
        }
      }
    }
    output[j] <- sum(sample.lst)
  }
  if(return.loc){
    return(out.seg)
  } else {
    return(output)
  }
}


