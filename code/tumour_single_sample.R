# ------------------------------------------------------------------
#          Title: tumour_single_sample.R
#          Description: Computes HRD scores using ASCAT output (Timms et al. Breast Cancer Research (2014) 16:475) using only tumour samples.
#          Matched normals are not needed for baseline computation.
#          Input files: 
#          SNP array data (700K array)
#          
#          author contact: swetansu@gmail.com
#
# Software Dependency:
#    	   R version >= 3.6, ASCAT 2.4	 
# ------------------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)
arg_file <- args[1]
print(arg_file)

#formatting file names for I/O

f2 <- readChar(arg_file, file.info(arg_file)$size)
f2 <- gsub("\t", " ", f2)
f2 <- gsub("\n","",f2)
files <- strsplit(f2, " ")
sam_id <- files[[1]][1]
f2 <- files[[1]][-1]

print(class(f2))
str(f2)
sam1 <- read.delim(f2, header = T, sep = "\t")
#make ascat Tum_BAF and Tum_LRR
ini_BAF <- sam1[,1:4]
ini_LRR <- sam1[,-(4)]

#run R in separate directories
arg_path <- getwd()
mod_dir <- paste(arg_file,"tab",sep=".")
new_path <- paste(arg_path,mod_dir, sep="/")
setwd(new_path)

#Write dataframes to output files

BAF_file <- paste("BAF_",sam_id,sep = "")
write.table(ini_BAF, file=paste(BAF_file,"txt",sep= ".") , quote = FALSE, row.names = FALSE, sep = "\t")
#make ascat Tum_LRR

LRR_file <- paste("LRR_",sam_id,sep = "")
write.table(ini_LRR, file=paste(LRR_file,"txt",sep= "."), quote = FALSE, row.names = FALSE, sep = "\t")

#for NO germline data (unmatched)
#the path to aspcf function in ascat.R function may be hardcoded..warning!!
#Using ASCAT2.4.. installed locally
library(ASCAT)
ascat.bc = ascat.loadData(paste(LRR_file,"txt",sep= "."), paste(BAF_file,"txt",sep= "."))
##platform needs to be specified for ASCAT to run efficiently. ASCAT provides support for most of Illumina arrays
#platform = "Illumina2.5M"
platform = "Illumina700k"
ascat.gg = ascat.predictGermlineGenotypes(ascat.bc, platform)
ascat.bc = ascat.aspcf(ascat.bc,ascat.gg=ascat.gg)
#ascat.plotSegmentedData(ascat.bc)
ascat.output = ascat.runAscat(ascat.bc)
#convert ascat output to tangible input for scoring the lesions

#ASCAT2.4 output segmentation
source("organise.segments_mean.R")
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
ascat.gs <- seg_mean_estimation(paste(LRR_file,"txt",sep= "."), ascat.output)
cols <- colnames(ascat.gs)
cols <- gsub("nMajor", "nA", cols)
cols <- gsub("nMinor", "nB", cols)
colnames(ascat.gs) <- cols

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

ascat.gs.tcn[,2] <- gsub("X", "23", ascat.gs.tcn[,2])
ascat.gs.tcn[,2] <- gsub("Y", "24", ascat.gs.tcn[,2])
class(ascat.gs.tcn[,2]) <- "numeric"
ascat.gs.tcn[,1] <- gsub("Log.R.Ratio",sam_id,ascat.gs.tcn[,1])
#save(ascat.gs.tcn,ascat.gs.tcn.hrd, file = "ascat_gs_fin.RData")
save(ascat.gs.tcn,file = "ascat_gs_fin.RData")

source("NTAI_LST_LOH_scores.R")
source("calc_lst_new2.R")
chrominfo <-GetChrominfo()
ntai <- calc.ai(ascat.gs.tcn, chrominfo, check.names = FALSE)
#lst <- calc.lst(ascat.gs.tcn, chrominfo, check.names = FALSE)
lst <- calc_lst_new(ascat.gs.tcn, chrominfo, check.names = FALSE)
hrd <- calc.hrd(ascat.gs.tcn,check.names = FALSE)
#hrd_loc <- calc.hrd(ascat.segments.tcn.hrd,check.names = FALSE,return.loc= TRUE)
#samples <- gsub(".ascat", "", files)
comb_res <- cbind(ntai, "LST"= lst, "HRD" = hrd, "sample_id" = sam_id)
#save(ascat.segments.tcn,ascat.segments.tcn.hrd, ntai,lst,hrd,comb_res, file = "ascat_check.RData")
write.table(comb_res, file = paste("Result",sam_id,"txt",sep=".") , col.names = TRUE, row.names = F, sep = "\t", quote = F)
