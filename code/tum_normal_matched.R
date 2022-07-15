# ------------------------------------------------------------------
#          Title: tum_normal_matched.R
#          Description: Computes HRD scores using ASCAT output (Timms et al. Breast Cancer Research (2014) 16:475)
#          Input files: 
#          SNP array data (700K array)
#          
#          author contact: swetansu@gmail.com
#
# Software Dependency:
#    	   R version >= 3.6, ASCAT 2.4	 
# ------------------------------------------------------------------


.libPaths(c( "~/R/x86_64-pc-linux-gnu-library/3.6", .libPaths() ) )
args <- commandArgs(trailingOnly = TRUE)
arg_file <- args[1]

#run R in separate directories
#arg_path <- getwd()
#new_path <- paste(arg_path,arg_file, sep="/")
#setwd(new_path)
#formatting file names for I/O (3 column input "id", "tum_id", "nor_id")

f1 <- readChar(arg_file, file.info(arg_file)$size)
f1 <- gsub("\t", " ", f1)
f1 <- gsub("\n", "", f1)
files <- strsplit(f1, " ")
aocs_id <- files[[1]][1]
out.file <- paste("Results_", files[[1]][1], sep="" )
files <- files[[1]][-1]

sam1 <- list()
for(i in 1:length(files)){
  sam1[[i]] <- read.delim(files[i] , header = T, sep = "\t", stringsAsFactors = F)
}

#make ascat BAF and LRR


##checks if arrays of two different sizes are used, matches the SNP IDs and retain only the common ones
library(dplyr)
library(purrr)
ini_file <- sam1 %>% reduce(full_join, by = "SNP.Name")
ini_file <- ini_file[complete.cases(ini_file), ]


#run R in separate directories
arg_path <- getwd()
mod_dir <- paste(arg_file,"tab",sep=".")
new_path <- paste(arg_path,mod_dir, sep="/")
setwd(new_path)

ascat_run <- function(ini_file){
#Write dataframes to output files
##process BAF

BAF_df <- ini_file[,c(1:4,8)]
colnames(BAF_df) <- c("SNP.Name", "Chr", "Position", "X1", "X2")
BAF_file <- paste("BAF_",aocs_id,sep = "")
write.table(BAF_df[,-5], file=paste(BAF_file,"tum.txt",sep= ".") , quote = FALSE, row.names = FALSE, sep = "\t")
write.table(BAF_df[,-4], file=paste(BAF_file,"nor.txt",sep= ".") , quote = FALSE, row.names = FALSE, sep = "\t")

##process LRR

#LRR_df <- data.frame(matrix(unlist(ini_LRR), ncol = length(files) ,byrow = FALSE),stringsAsFactors=FALSE)
LRR_df <- ini_file[,c(1:3,5,9)]
colnames(LRR_df) <- c("SNP.Name", "Chr", "Position", "X1", "X2")
LRR_file <- paste("LRR_",aocs_id,sep = "")
write.table(LRR_df[,-5], file=paste(LRR_file,"tum.txt",sep= ".") , quote = FALSE, row.names = FALSE, sep = "\t")
write.table(LRR_df[,-4], file=paste(LRR_file,"nor.txt",sep= ".") , quote = FALSE, row.names = FALSE, sep = "\t")


#Run ASCAT2.4
#for matched germline data 
#the path to aspcf function in ascat.R function may be hardcoded..warning!!
library(RColorBrewer)
library(ASCAT)
ascat.bc <- ascat.loadData(Tumor_LogR_file=paste(LRR_file,"tum.txt",sep= "."), Tumor_BAF_file=paste(BAF_file,"tum.txt",sep= "."), 
                           Germline_LogR_file=paste(LRR_file,"nor.txt",sep= "."), Germline_BAF_file=paste(BAF_file,"nor.txt",sep= "."))

ascat.plotRawData(ascat.bc)
ascat.bc <- ascat.aspcf(ascat.bc)
ascat.plotSegmentedData(ascat.bc)
ascat.output = ascat.runAscat(ascat.bc)

##segmentation profile
#convert ascat output to tangible input for scoring the lesions
source("organise.segments_mean.R")
ascat.segments <- seg_mean_estimation(paste(LRR_file,"tum.txt",sep= "."), ascat.output)
cols <- colnames(ascat.segments)
cols <- gsub("nMajor", "nA", cols)
cols <- gsub("nMinor", "nB", cols)
colnames(ascat.segments) <- cols
#add the allele specific copy numbers to get total copy number (tcn)
ascat.segments.tcn <- cbind(ascat.segments[,c(1:4,7)], ascat.segments[,5]+ascat.segments[,6], ascat.segments[,5:6])
#save(ascat.bc, ascat.output, ascat.segments, file = "ascat_files.RData") ## uncomment to save segmentation profile
colnames(ascat.segments.tcn)[6] <- "tcn"

#Extract ploidy calls
ploidy <- ascat.output$ploidy
names(ploidy) <- colnames(ascat.output$nA)

#Add ploidy info
ascat.segments.tcn <- cbind(ascat.segments.tcn[,1:8], ploidy[match(ascat.segments.tcn[,1], names(ploidy))])
colnames(ascat.segments.tcn)[9] <- "ploidy"

#Add aberrant cell fraction info
cellularity <- ascat.output$aberrantcellfraction
names(cellularity) <- colnames(ascat.output$nA)
ascat.segments.tcn <- cbind(ascat.segments.tcn[,1:9], cellularity[match(ascat.segments.tcn[,1], names(cellularity))])
colnames(ascat.segments.tcn)[10] <- "cellularity"
#save(ascat.bc, ascat.output, ascat.segments, ascat.segments.tcn, file = "ascat_files.RData")

#remove chromosome 17 due to ubiquitous LOH for HRD computation
#ascat.segments.tcn.hrd <- ascat.segments.tcn[ ascat.segments.tcn[,2] != 17, ] #commented since chr17 may not be ubiquitously deleted and may hold unique information
ascat.segments.tcn.hrd <- ascat.segments.tcn
ascat.segments.tcn.hrd[,2] <- gsub("X", "23", ascat.segments.tcn.hrd[,2])
ascat.segments.tcn.hrd[,2] <- gsub("Y", "24", ascat.segments.tcn.hrd[,2])
class(ascat.segments.tcn.hrd[,2]) <- "numeric"
rm(ascat.segments.tcn)
ascat.segments.tcn.hrd[,1] <- gsub("Log.R.Ratio",aocs_id,ascat.segments.tcn.hrd[,1])
save(ascat.segments.tcn.hrd, file = "ascat_seg_fin.RData")
return(ascat.segments.tcn.hrd)
}

system.time({ seg_ascat_data <- ascat_run(ini_file) })

run_HRD_est_fun <- function(ascat_seg_dat){
source("NTAI_LST_LOH_scores.R")
source("calc_lst_new2.R")
chrominfo <-GetChrominfo()
ntai <- calc.ai(ascat_seg_dat, chrominfo, check.names = FALSE)
#lst <- calc.lst(ascat_seg_dat, chrominfo, check.names = FALSE)
lst <- calc_lst_new(ascat_seg_dat, chrominfo, check.names = FALSE)
hrd <- calc.hrd(ascat_seg_dat,check.names = FALSE)
#hrd_loc <- calc.hrd(ascat.segments.tcn.hrd,check.names = FALSE,return.loc= TRUE)
#samples <- gsub(".ascat", "", files)
comb_res <- cbind(ntai, "LST"= lst, "HRD" = hrd, "sample_id" = aocs_id)
#save(ascat.segments.tcn,ascat.segments.tcn.hrd, ntai,lst,hrd,comb_res, file = "ascat_check.RData") ##uncomment for QC
write.table(comb_res, file = out.file , col.names = TRUE, row.names = F, sep = "\t", quote = F)
}
system.time({ run_HRD_est_fun(ascat_seg_dat = seg_ascat_data) })
