

#capture file name for parallel processing:modified for testing on Liz samples -> modification is generic -> it has new GetChrominfo function 
args <- commandArgs(trailingOnly = TRUE)
arg_file <- args[1]
print(arg_file)
#run R in separate directories
#arg_path <- getwd()
#new_path <- paste(arg_path,arg_file, sep="/")
#setwd(new_path)

#formatting file names for I/O
#out.file <- paste("Results_", args[1], sep="")
#f1 <- readChar(arg_file, file.info(arg_file)$size)
#f1 <- gsub("\n","",f1)
#f2 <- readChar(f1, file.info(f1)$size)
f2 <- readChar(arg_file, file.info(arg_file)$size)
f2 <- gsub("\t", " ", f2)
f2 <- gsub("\n","",f2)
files <- strsplit(f2, " ")
sam_id <- files[[1]][1]
f2 <- files[[1]][-1]
#f2 <- gsub("\n"," ",f2)
#file_list <- gsub("^G","G",file_list, perl = TRUE)
#files <-strsplit(f2, " ")
#files <- list()
#files[[1]] <- f2
#print(files)
print(class(f2))
str(f2)
#sam1 <- list()
#for(i in 1:length(files[[1]])){
 # sam1[[i]] <- read.delim(files[[1]][i] , header = T, sep = "\t")
#}
#sam1 <- read.delim(files[[1]][1], header = T, sep = "\t")
sam1 <- read.delim(f2, header = T, sep = "\t")
#make ascat Tum_BAF and Tum_LRR
#ini_BAF <- list()
#ini_LRR <- list()
#for (k in 1: length(sam1)){
#  ini_BAF[[k]] <- (sam1[[k]][,4]) 
#  ini_LRR[[k]] <- sam1[[k]][,5]
#}

ini_BAF <- sam1[,1:4]
ini_LRR <- sam1[,-(4)]


#run R in separate directories
arg_path <- getwd()
mod_dir <- paste(arg_file,"tab",sep=".")
new_path <- paste(arg_path,mod_dir, sep="/")
setwd(new_path)

#Write dataframes to output files
#BAF_df <- data.frame(matrix(unlist(ini_BAF), ncol = length(files[[1]]), byrow = FALSE),stringsAsFactors=FALSE)
#dim1 <- dim(BAF_df)[2]
#lim1 <- dim(BAF_df)[2] + 1
#lim2 <- dim(BAF_df)[2] + 3
#BAF_df[,lim1:lim2] <- c(sam1[[k]][,1:3])
#BAF_df <- BAF_df[,c(lim1:lim2,1:dim1)]
#BAF_df[,1] <- gsub("-","_", BAF_df[,1])
BAF_file <- paste("BAF_",sam_id,sep = "")
write.table(ini_BAF, file=paste(BAF_file,"txt",sep= ".") , quote = FALSE, row.names = FALSE, sep = "\t")
#make ascat Tum_LRR
#ini_LRR <- list()
#for (k in 1: length(sam1)){
#  ini_LRR[[k]] <- sam1[[k]][,5] 
#}

#LRR_df <- data.frame(matrix(unlist(ini_LRR), ncol = length(files[[1]]) ,byrow = FALSE),stringsAsFactors=FALSE)
#dim1 <- dim(LRR_df)[2]
#lim1 <- dim(LRR_df)[2] + 1
#lim2 <- dim(LRR_df)[2] + 3
#LRR_df[,lim1:lim2] <- c(sam1[[k]][,1:3])
#LRR_df <- LRR_df[,c(lim1:lim2,1:dim1)]
#LRR_df[,1] <- gsub("-","_", LRR_df[,1])
LRR_file <- paste("LRR_",sam_id,sep = "")
write.table(ini_LRR, file=paste(LRR_file,"txt",sep= "."), quote = FALSE, row.names = FALSE, sep = "\t")

#for NO germline data (unmatched)
#the path to aspcf function in ascat.R function may be hardcoded..warning!!
#source("/share/ClusterShare/thingamajigs/swepat/HRD signature/ASCAT2.1/ascat.R")
#Using ASCAT2.4.. installed locally
library(ASCAT)
ascat.bc = ascat.loadData(paste(LRR_file,"txt",sep= "."), paste(BAF_file,"txt",sep= "."))
#save(ascat.bc, file = "ascat_bc1.RData")
#ascat.plotRawData(ascat.bc)
#source("/share/ClusterShare/thingamajigs/swepat/HRD signature/ASCAT2.1/predictGG.R")
#platform = "Illumina2.5M"
platform = "Illumina700k"
ascat.gg = ascat.predictGermlineGenotypes(ascat.bc, platform)
#save(ascat.gg, file = "ascat_gg.RData")
#source("/share/ClusterShare/thingamajigs/swepat/HRD signature/ASCAT2.1/aspcf.R")
ascat.bc = ascat.aspcf(ascat.bc,ascat.gg=ascat.gg)
#ascat.plotSegmentedData(ascat.bc)
ascat.output = ascat.runAscat(ascat.bc)
#convert ascat output to tangible input for scoring the lesions
#source("/share/ClusterShare/thingamajigs/swepat/HRD signature/ASCAT2.1/organize.ascat.segments.R")
#ascat.segments <- organize.ascat.segments(ascat.output, ascat.bc$SNPpos)

#ASCAT2.4 output segmentation
source("/share/ClusterShare/thingamajigs/swepat/Liz/Batch3_ILOE24-11837/ILOE24-11837/Analysis/ASCAT/ASCAT_24/organise.segments_mean.R")
ascat.segments <- seg_mean_estimation(paste(LRR_file,"txt",sep= "."), ascat.output)
cols <- colnames(ascat.segments)
cols <- gsub("nMajor", "nA", cols)
cols <- gsub("nMinor", "nB", cols)
colnames(ascat.segments) <- cols
#add the allele specific copy numbers to get total copy number (tcn)
ascat.segments.tcn <- cbind(ascat.segments[,c(1:4,7)], ascat.segments[,5]+ascat.segments[,6], ascat.segments[,5:6], "seg.Mean"= ascat.segments[,8],"baf.Mean"= ascat.segments[,9])
#save(ascat.bc, ascat.output, ascat.segments, file = "ascat_files.RData")
colnames(ascat.segments.tcn)[6] <- "tcn"
#Extract ploidy calls
ploidy <- ascat.output$ploidy
names(ploidy) <- colnames(ascat.output$nA)
#Add ploidy info
ascat.segments.tcn <- cbind(ascat.segments.tcn[,1:10], ploidy[match(ascat.segments.tcn[,1], names(ploidy))])
colnames(ascat.segments.tcn)[11] <- "ploidy"
#head(ascat.segments.tcn)
#Add aberrant cell fraction info
cellularity <- ascat.output$aberrantcellfraction
names(cellularity) <- colnames(ascat.output$nA)
ascat.segments.tcn <- cbind(ascat.segments.tcn[,1:11], cellularity[match(ascat.segments.tcn[,1], names(cellularity))])
colnames(ascat.segments.tcn)[12] <- "cellularity"
ascat.segments.tcn[,1] <- gsub("Log.R.Ratio",sam_id,ascat.segments.tcn[,1])
#Write segments for breakpoint phylogeny
#save(ascat.bc, ascat.output, ascat.segments, ascat.segments.tcn, file = "ascat_files.RData")
#write.table(ascat.segments.tcn,file=paste0(sam_id,"_segmented.csv"),col.names=T,row.names=F,quote=F,sep="\t")

#For genomic lesion scores

#source("/share/ClusterShare/thingamajigs/swepat/Liz/Batch3_ILOE24-11837/ILOE24-11837/Analysis/ASCAT/ASCAT_24/organise.segments_mean.R")
ascat.gs <- seg_mean_estimation(paste(LRR_file,"txt",sep= "."), ascat.output)
cols <- colnames(ascat.gs)
cols <- gsub("nMajor", "nA", cols)
cols <- gsub("nMinor", "nB", cols)
colnames(ascat.gs) <- cols
#add the allele specific copy numbers to get total copy number (tcn)
#add the allele specific copy numbers to get total copy number (tcn)
ascat.gs.tcn <- cbind(ascat.gs[,c(1:4,7)], ascat.gs[,5]+ascat.gs[,6], ascat.gs[,5:6])
#save(ascat.bc, ascat.output, ascat.segments, file = "ascat_files.RData")
colnames(ascat.gs.tcn)[6] <- "tcn"
#Extract ploidy calls
ploidy <- ascat.output$ploidy
names(ploidy) <- colnames(ascat.output$nA)
#Add ploidy info
ascat.gs.tcn <- cbind(ascat.gs.tcn[,1:8], ploidy[match(ascat.gs.tcn[,1], names(ploidy))])
colnames(ascat.gs.tcn)[9] <- "ploidy"
#head(ascat.segments.tcn)
#Add aberrant cell fraction info
cellularity <- ascat.output$aberrantcellfraction
names(cellularity) <- colnames(ascat.output$nA)
ascat.gs.tcn <- cbind(ascat.gs.tcn[,1:9], cellularity[match(ascat.gs.tcn[,1], names(cellularity))])
colnames(ascat.gs.tcn)[10] <- "cellularity"

#add the allele specific copy numbers to get total copy number (tcn)
##ascat.segments.tcn <- cbind(ascat.segments[,1:5], ascat.segments[,6]+ascat.segments[,7], ascat.segments[,6:7])
#save(ascat.bc, ascat.output, ascat.segments, file = "ascat_files.RData")

#remove chromosome 17 due to ubiquitous LOH for HRD computation
#ascat.gs.tcn.hrd <- ascat.gs.tcn[ ascat.gs.tcn[,2] != 17, ]
#ascat.gs.tcn.hrd[,2] <- gsub("X", "23", ascat.gs.tcn.hrd[,2])
#ascat.gs.tcn.hrd[,2] <- gsub("Y", "24", ascat.gs.tcn.hrd[,2])
#class(ascat.gs.tcn.hrd[,2]) <- "numeric"

#ascat.gs.tcn.new <- ascat.gs.tcn[ ascat.gs.tcn[,2] != 17, ]
ascat.gs.tcn[,2] <- gsub("X", "23", ascat.gs.tcn[,2])
ascat.gs.tcn[,2] <- gsub("Y", "24", ascat.gs.tcn[,2])
class(ascat.gs.tcn[,2]) <- "numeric"

#ascat.gs.tcn.hrd[,1] <- gsub("Log.R.Ratio",sam_id,ascat.gs.tcn.hrd[,1])
ascat.gs.tcn[,1] <- gsub("Log.R.Ratio",sam_id,ascat.gs.tcn[,1])
#save(ascat.gs.tcn,ascat.gs.tcn.hrd, file = "ascat_gs_fin.RData")
save(ascat.gs.tcn,file = "ascat_gs_fin.RData")

source("/share/ClusterShare/thingamajigs/swepat/HRD signature/NTAI_LST_LOH_scores.R")
source("/share/ClusterShare/thingamajigs/swepat/Cristina/proj_2021/INOVATE_HRD_July_2021/Analysis/paired_Tumour/calc_lst_new2.R")
chrominfo <-GetChrominfo()
#ntai <- calc.ai(ascat.gs.tcn.hrd, chrominfo, check.names = FALSE)
ntai <- calc.ai(ascat.gs.tcn, chrominfo, check.names = FALSE)
#lst <- calc.lst(ascat.gs.tcn, chrominfo, check.names = FALSE)
lst <- calc_lst_new(ascat.gs.tcn, chrominfo, check.names = FALSE)
hrd <- calc.hrd(ascat.gs.tcn,check.names = FALSE)
#hrd_loc <- calc.hrd(ascat.segments.tcn.hrd,check.names = FALSE,return.loc= TRUE)
#samples <- gsub(".ascat", "", files)
comb_res <- cbind(ntai, "LST"= lst, "HRD" = hrd, "sample_id" = sam_id)
#save(ascat.segments.tcn,ascat.segments.tcn.hrd, ntai,lst,hrd,comb_res, file = "ascat_check.RData")
write.table(comb_res, file = paste("Result",sam_id,"txt",sep=".") , col.names = TRUE, row.names = F, sep = "\t", quote = F)
