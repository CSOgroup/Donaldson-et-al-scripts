#! /usr/bin/Rscript



###############################################################################################################################################################
calculate_MoC_with_domainTypes <- function(set1, set2, chr_len, gaps_as_clusters = FALSE, file_as_input = FALSE) {
	#******************************************************************************************** FUNCTION DEFINITION
    getOverlap <- function(a, b){
        # smallest end - biggest start    
        return(max(0, min(a[2], b[2]) - max(a[1], b[1]) +1))
    }
    prepareDT <- function(dt, chr_len) {
        stopifnot(is.numeric(chr_len))
        stopifnot(ncol(dt) == 3)
        # do not take the 1st column with the "chr6"    
        dt <- dt[,2:3]
        colnames(dt) <- c("Start", "End")
        # add a column filled with 0
        dt$is_gap <- 0
        # test that nrow dt is bigger than 1 otherwise the 2:... will create NA
        if(nrow(dt) > 1) {
          # create data frame that will hold the gap for the 1st dataset
          # start of the gap = end of the end of the TAD + 1 (do not take the last row)    
          # end of the gap = start of the TAD - 1 (do not take the first row)
          dt_gaps <- data.frame( Start = (dt$End[1:(nrow(dt)-1)] + 1),
                                  End = (dt$Start[2:nrow(dt)] -1))
          stopifnot(is.numeric(dt_gaps$Start[1]))
          stopifnot(is.numeric(dt_gaps$End[1]))    
          # select only the row with end > start
          dt_gaps <- dt_gaps[dt_gaps$Start < dt_gaps$End,]
        } else{
          dt_gaps <- data.frame(Start = numeric(), End=numeric())
        }
        # add gaps at the beginning until first TAD and at the end until end of chromo size
        # FIRST DOMAIN APPENDED SHOULD START WITH 1 NOT WITH 0
        # if needed, add gap before 1st TAD
        if(dt$Start[1] > 1) {
            tmpDT <- data.frame(Start = 1, End = dt$Start[1]-1)
            dt_gaps <- rbind(tmpDT, dt_gaps)
        }
        # if needed, add gap until chromosome end    
        if(dt$End[nrow(dt)] < chr_len) {
            tmpDT <- data.frame(Start = dt$End[nrow(dt)] + 1, End = chr_len)
            dt_gaps <- rbind(dt_gaps, tmpDT)        
        }
        # add a column to indicate there are gaps
        if(nrow(dt_gaps) > 0) {
            dt_gaps$is_gap <- 1
            dt_final <- rbind(dt, dt_gaps)
        } else{
            dt_final <- dt
        }
        dt_final <- dt_final[order(dt_final$Start),]
        return(dt_final)
    }
	#********************************************************************************************
    if(file_as_input) {
        if(file.info(set1)$size == 0 & file.info(set2)$size == 0)
            return(NA)
        if(file.info(set1)$size == 0 & file.info(set2)$size > 0)
            return(0)            
        set1DT <- read.delim(set1, header=F, col.names=c("chromo", "start", "end"))
        set2DT <- read.delim(set2, header=F, col.names=c("chromo", "start", "end"))
    } else {
        set1DT <- set1
        set2DT <- set2
        colnames(set1DT) <- c("chromo", "start", "end")
        colnames(set2DT) <- c("chromo", "start", "end")
    }
    stopifnot(is.numeric(set1DT$start[1]))
    stopifnot(is.numeric(set2DT$start[1]))
    stopifnot(is.numeric(set1DT$end[1]))        
    stopifnot(is.numeric(set2DT$end[1]))            
    
    # prepare data for 1st caller
    ptot1 = prepareDT(set1DT, chr_len)

    # prepare data for 2nd caller
    ptot2 = prepareDT(set2DT, chr_len)

    # number of clusters for each TAD caller (# of rows)
    Nclust1 = nrow(ptot1)
    Nclust2 = nrow(ptot2)

    # in MoC definition: if I = J = 1 => MoC = 1
    if(Nclust1==1 & Nclust2==1)
        return(1)
        
    # added for the problem of negative value if compared 1 single whole-chromosome domain to more than 1 domain
    # in this case return NA
    if( (Nclust1 == 1 & Nclust2 > 1) | (Nclust1 > 1 & Nclust2 == 1) )
        return(NA)

    # interTAD regions are considered as domains
    if(gaps_as_clusters){
        crossum = 0
        # iterate over dt1 -> for each domain retrieve the overlapping domains ("fragments") and compute mutual concordance
        for(i in 1:Nclust1) {
            # index of the 1st end in dt2 that has end > than current start
            firstj = min(which(ptot2$End > ptot1$Start[i]))
            # index of the 1st start in dt2 with start < current end
            lastj = max(which(ptot2$Start < ptot1$End[i]))
            
            if(is.infinite(firstj) | is.infinite(lastj)) next
            
            for(j in (firstj:lastj)) {
                crossum = crossum + (getOverlap(ptot1[i,1:2], ptot2[j,1:2]))**2/
                                        ( (ptot1[i,2] - ptot1[i,1] + 1)*(ptot2[j,2] - ptot2[j,1] + 1) )
            }            
        }            
        MoC = 1/(sqrt(Nclust1*Nclust2)-1)*(crossum - 1)
        return(MoC)
    } else{          
    # if not interTAD_as_domains => concordance term for <TAD <-> interTAD regions> = 0
        crossum = 0
        for(i in 1:Nclust1 ){
            firstj = min(which(ptot2$End > ptot1$Start[i]))
            lastj = max(which(ptot2$Start < ptot1$End[i] ))            
            for(j in (firstj:lastj)){
                crossum = crossum + ( as.numeric(ptot1$is_gap[i] == ptot2$is_gap[j]) )*(getOverlap(ptot1[i,1:2], ptot2[j,1:2]))**2/
                               ( (ptot1[i,2] - ptot1[i,1] + 1)*(ptot2[j,2] - ptot2[j,1] + 1) )              
            }             
        }
        MoC = 1/(sqrt(Nclust1*Nclust2)-1)*(crossum - 1)
        return(MoC)
    }
}    
###############################################################################################################################################################


SSHFS <- F

### !!! WARNING
# - hard-coded tolerance threshold = 2*binSize
# - the resolution of the input data is inferred from the name of input folder _(#)kb_
# (so for the moment only implemented comparison between of the data coming from the same resolution)

cat(paste0("> Start compareTADs.R\n"))
startTime <- Sys.time()

options(scipen = 100)

suppressPackageStartupMessages(library(foreach, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
suppressPackageStartupMessages(library(doMC, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
suppressPackageStartupMessages(library(reshape2, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
suppressPackageStartupMessages(library(optparse, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
suppressPackageStartupMessages(library(IRanges, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
suppressPackageStartupMessages(library(dplyr, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))

option_list = list(
  make_option(c("-f", "--folder"), type="character", default=NULL,
              help="input folder", metavar="character"),
  
  make_option(c("-p", "--pattern"), type="character", default=NULL,
              help="matching pattern to retrieve files", metavar="character"),
  
  make_option(c("-c", "--ncpu"), type="integer", default=NULL,
              help="matching pattern to retrieve files", metavar="integer"),
  
  make_option(c("-b", "--caller1"), type="character", default=NULL,
              help="caller 1", metavar="character"),
  
  make_option(c("-d", "--cl1"), type="character", default=NULL,
              help="cell line 1", metavar="character"),
  
  make_option(c("-e", "--treat1"), type="character", default=NULL,
              help="treatment 1", metavar="character"),

  make_option(c("-B", "--caller2"), type="character", default=NULL,
              help="caller 2", metavar="character"),

  make_option(c("-D", "--cl2"), type="character", default=NULL,
              help="cell line 2", metavar="character"),
  
  make_option(c("-E", "--treat2"), type="character", default=NULL,
              help="treatment 2", metavar="character"),
  
  make_option(c("-o", "--outFold"), type="character", default=NULL,
              help="output folder", metavar="character")
)

opt_parser <- OptionParser(option_list=option_list);
opt <- parse_args(opt_parser);
if(is.null(opt$folder)| is.null(opt$pattern)  | is.null(opt$outFold) | is.null(opt$caller1) | is.null(opt$caller2) |
   is.null(opt$cl1) | is.null(opt$cl2) |is.null(opt$treat1) | is.null(opt$treat2)) {
  stop("Error - missing input argument !")
}

matchPat <- opt$pattern

inFolder <- opt$folder
resol <- as.numeric(gsub(".+_(.+?)kb_.+", "\\1", inFolder))
stopifnot(!is.na(resol))

treat1 <- opt$treat1
treat2 <- opt$treat2

cell_line1 <- opt$cl1
cell_line2 <- opt$cl2

caller1 <- opt$caller1
caller2 <- opt$caller2

outFold <- opt$outFold
system(paste0("mkdir -p ", outFold))

ncpu <- ifelse(is.null(opt$ncpu), 1, opt$ncpu)
registerDoMC(ncpu)




cat("!!! WARNING: hard-coded tolerance threshold = 2*bin size !!!\n")
tolThresh <- 2*resol*10^3
cat(paste0("... found tolerance threshold: ", tolThresh, "\n"))




#cond1, cond2, cond1_cond2
my_col <- c("skyblue", "royalblue2", "darkolivegreen3") # for the unique1, unique2, shared
my_col2 <- c("dodgerblue1", "darkolivegreen3") # for not conserved, conserved
names(my_col2) <- c("not cons.", "conserved")

cond1 <- paste0(cell_line1, "_", treat1)
cond2 <- paste0(cell_line2, "_", treat2)

inFolder1 <- file.path(inFolder, cell_line1, treat1, "BEDfiles", caller1)
inFolder2 <- file.path(inFolder, cell_line2, treat2, "BEDfiles", caller2)

all_files1 <- list.files(inFolder1, full.names = T)
all_files1 <- all_files1[grep(matchPat, all_files1)]
stopifnot(length(all_files1) > 0)
all_files2 <- list.files(inFolder2, full.names = T)
all_files2 <- all_files2[grep(matchPat, all_files2)]
stopifnot(length(all_files2) > 0)

all_chr1 <- gsub(".+_(chr.+?)_.+", "\\1", basename(all_files1))
all_chr2 <- gsub(".+_(chr.+?)_.+", "\\1", basename(all_files2))

all_chr <- intersect(all_chr1, all_chr2)

cat(paste0("... from folder1: ", length(all_chr), "/", length(all_chr1), " matching chromosomes\n"))
cat(paste0("... from folder2: ", length(all_chr), "/", length(all_chr2), " matching chromosomes\n"))


all_chr_overlapBD <- foreach(chromo = all_chr, .combine='rbind') %dopar% {
  inFile1 <- all_files1[grep(paste0("_", chromo, "_"), all_files1)]
  stopifnot(length(inFile1) == 1)
  stopifnot(file.exists(inFile1))
  
  inFile2 <- all_files2[grep(paste0("_", chromo, "_"), all_files2)]
  stopifnot(length(inFile2) == 1)
  stopifnot(file.exists(inFile2))
  
  set1DT <- read.delim(inFile1, header=F, col.names=c("chromo", "start", "end"), stringsAsFactors = F)
  set2DT <- read.delim(inFile2, header=F, col.names=c("chromo", "start", "end"), stringsAsFactors = F)
  
  set1DT <- set1DT[set1DT$chromo == chromo,]
  if(nrow(set1DT) < 1)
    stop("NO DATA IN SET1 DT\n")
  set2DT <- set2DT[set2DT$chromo == chromo,]
  if(nrow(set2DT) < 1)
    stop("NO DATA IN SET2 DT\n")
  
  ###====================================== ASYMMETRIC: MATCHING STARTS AND ENDS
  #== A) strict: conservation of TAD as it
  set1DT_matching_strict <- unlist(foreach(i=1:nrow(set1DT), .combine='c') %dopar% {
    any( (abs(set2DT$start - set1DT$start[i]) <= tolThresh) &
           (abs(set2DT$end - set1DT$end[i]) <= tolThresh) )
  })
  set1DT_nStrictMatch <- sum(set1DT_matching_strict)
  stopifnot(set1DT_nStrictMatch <= nrow(set1DT))
  cat(paste0(cond1, " TADs strictly matching in ", cond2, ": ", set1DT_nStrictMatch, "/", nrow(set1DT), "\n" ))
  
  set2DT_matching_strict <- unlist(foreach(i=1:nrow(set2DT), .combine='c') %dopar% {
    any( (abs(set1DT$start - set2DT$start[i]) <= tolThresh) &
           (abs(set1DT$end - set2DT$end[i]) <= tolThresh) )
  })
  set2DT_nStrictMatch <- sum(set2DT_matching_strict)
  stopifnot(set2DT_nStrictMatch <= nrow(set2DT))
  cat(paste0(cond2, " TADs strictly matching in ", cond1, ": ", set2DT_nStrictMatch, "/", nrow(set2DT), "\n" ))
  
  
  #== B) less strict: conservation of starts and ends
  set1DT_matching <- unlist(foreach(i=1:nrow(set1DT), .combine='c') %dopar% {
    any(abs(set2DT$start - set1DT$start[i]) <= tolThresh) &
      any(abs(set2DT$end - set1DT$end[i]) <= tolThresh)   
  })
  set1DT_nLooseMatch <- sum(set1DT_matching)
  stopifnot(set1DT_nLooseMatch <= nrow(set1DT))
  stopifnot(set1DT_nLooseMatch >= set1DT_nStrictMatch)
  cat(paste0(cond1, " TADs loosely matching in ", cond2, ": ", set1DT_nLooseMatch, "/", nrow(set1DT), "\n" ))
  
  set2DT_matching <- unlist(foreach(i=1:nrow(set2DT), .combine='c') %dopar% {
    any(abs(set1DT$start - set2DT$start[i]) <= tolThresh) &
      any(abs(set1DT$end - set2DT$end[i]) <= tolThresh)   
  })
  set2DT_nLooseMatch <- sum(set2DT_matching)
  stopifnot(set2DT_nLooseMatch <= nrow(set2DT))
  stopifnot(set2DT_nLooseMatch >= set2DT_nStrictMatch)
  cat(paste0(cond2, " TADs loosely matching in ", cond1, ": ", set2DT_nLooseMatch, "/", nrow(set2DT), "\n" ))
  
  
  ###====================================== SYMMETRIC: BOUNDARY REGIONS CONSERVATION
  cat(paste0("> Calculate boundary overlap - ", chromo, "\n"))
  overlapVect <- get_boundariesOverlap_pipeline(set1DT, set2DT, tolRad=tolThresh, verbose = F,
                                                myTit = paste0("Boundary overlap - ", chromo),
                                                set1_name=cond1, set2_name=cond2,
                                                set1_col = "darkgreen", set2_col = "darkblue",
                                                gridNew = F, plotVenn=F, plotAtCall =F, returnPlot = F)
  
  
  ###====================================== get MoC to annotate the plots
  sizefile1 <- paste0(gsub("(.+)/BEDfiles/.+", "\\1", inFile1),"/NORM/", chromo, ".size")
  sizefile2 <- paste0(gsub("(.+)/BEDfiles/.+", "\\1", inFile2),"/NORM/", chromo, ".size") 
  stopifnot(file.exists(sizefile1))
  stopifnot(file.exists(sizefile2))
  size1 <- read.delim(sizefile1, header=F)[1,2]
  size2 <- read.delim(sizefile2, header=F)[1,2]
  stopifnot(is.numeric(size1))
  stopifnot(is.numeric(size2))
  chrSize <- max(c(size1, size2))
  
  r_moc <- calculate_MoC_with_domainTypes(set1=inFile1, set2=inFile2, 
                                          chr_len=chrSize, gaps_as_clusters = TRUE, file_as_input = TRUE)
  
  data.frame(chromo = chromo, 
             
            uniqueBR1 = overlapVect[1],
            uniqueBR2 = overlapVect[2],
            sharedBR = overlapVect[3],
            
            looseMatchTAD1 = set1DT_nLooseMatch,
            looseUniqueTAD1 = nrow(set1DT) - set1DT_nLooseMatch,
            
            strictMatchTAD1 = set1DT_nStrictMatch,
            strictUniqueTAD1 = nrow(set1DT) - set1DT_nStrictMatch,
            
            looseMatchTAD2 = set2DT_nLooseMatch,
            looseUniqueTAD2 = nrow(set2DT) - set2DT_nLooseMatch,
            
            strictMatchTAD2 = set2DT_nStrictMatch,
            strictUniqueTAD2 = nrow(set2DT) - set2DT_nStrictMatch,
            
            MoC = r_moc,
            
            stringsAsFactors =F
             )
}

