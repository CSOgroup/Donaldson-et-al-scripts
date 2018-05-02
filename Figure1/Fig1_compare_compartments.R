#! /usr/bin/Rscript

### !!! WARNING
# - the resolution of the input data is inferred from the name of input folder _(#)kb_
# (so for the moment only implemented comparison between of the data coming from the same resolution)

SSHFS <- F
setDir <- ifelse(SSHFS, "/media/electron", "")
startTime <- Sys.time()

cat(paste0("> Start compare_compartments.R\n"))

options(scipen = 100)

suppressPackageStartupMessages(library(dplyr, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
suppressPackageStartupMessages(library(foreach, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
suppressPackageStartupMessages(library(doMC, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
suppressPackageStartupMessages(library(reshape2, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
suppressPackageStartupMessages(library(optparse, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))

option_list = list(
  make_option(c("-f", "--folder1"), type="character", default=NULL,
              help="input folder1", metavar="character"),
  
  make_option(c("-F", "--folder2"), type="character", default=NULL,
              help="input folder2", metavar="character"),
  
  make_option(c("-p", "--pattern"), type="character", default=NULL,
              help="matching pattern to retrieve files", metavar="character"),
  
  make_option(c("-c", "--ncpu"), type="integer", default=NULL,
              help="matching pattern to retrieve files", metavar="integer"),
  
  make_option(c("-d", "--cl1"), type="character", default=NULL,
              help="cell line 1", metavar="character"),
  
  make_option(c("-e", "--treat1"), type="character", default=NULL,
              help="treatment 1", metavar="character"),
  
  make_option(c("-D", "--cl2"), type="character", default=NULL,
              help="cell line 2", metavar="character"),
  
  make_option(c("-E", "--treat2"), type="character", default=NULL,
              help="treatment 2", metavar="character"),
  
  make_option(c("-o", "--outFold"), type="character", default=NULL,
              help="output folder", metavar="character")
)

# to be compatible with both MAPQFILTER AND GSE63525
# => add possible folder2
# => treat2 is not mandatory anymore

opt_parser <- OptionParser(option_list=option_list);
opt <- parse_args(opt_parser);
if(is.null(opt$folder1)| is.null(opt$outFold) |
   is.null(opt$cl1) | is.null(opt$cl2) |is.null(opt$treat1)) {
  stop("Error - missing input argument !")
}

matchPat <- ifelse(is.null(opt$pattern), "_binPCscore.txt$|_cptmtCoord.txt$", opt$pattern)

inFolder1 <- opt$folder1
if(is.null(opt$folder2)) {
  # this will be the case for MAPQFILTER AND REPLICATES
  inFolder2 <- inFolder1
} else {
  inFolder2 <- opt$folder2
}

treat1 <- opt$treat1
treat2 <- opt$treat2  # may be NULL for GSE63525

cell_line1 <- opt$cl1
cell_line2 <- opt$cl2

outFold <- opt$outFold
system(paste0("mkdir -p ", outFold))

registerDoMC(ifelse(is.null(opt$ncpu), 1, opt$ncpu))

plot_cols <- c(even = "darkgrey", odd = "gray29", "0" = "white", "-1" = "blue", "1" = "red")

chromoDelimCol <- "gray41"

# plotType <- "png"
plotType <- "svg"


resol1 <- gsub(".+_(.+?)kb_.+", "\\1", inFolder1)
resol2 <- gsub(".+_(.+?)kb_.+", "\\1", inFolder2)
stopifnot(resol1 == resol2)
resol <- resol1
stopifnot(!is.na(as.numeric(resol)))

if(as.numeric(resol) < 1000 ) {
  plotType <- "png"
}
myWidth <- 12
myHeight <- 7

myWidth2 <- ifelse(plotType == "png", 600,12)
myHeight2 <- ifelse(plotType == "png", 480,7)



inFolder1 <- file.path(inFolder1, cell_line1, treat1, "CPTMTS")

if(is.null(treat2)) {
  # this will be the case for GSE63525
  inFolder2 <- file.path(inFolder2, cell_line2, "CPTMTS")
} else{
  # this will be the case for MAQFILTER or REPLICATES
  inFolder2 <- file.path(inFolder2, cell_line2, treat2, "CPTMTS")
}

cat(paste0("inFolder1: ", inFolder1, "\n"))
cat(paste0("inFolder2: ", inFolder2, "\n"))

all_files1 <- list.files(inFolder1, full.names = T)
all_files1 <- all_files1[grep(matchPat, all_files1)]
stopifnot(length(all_files1) > 0)
all_files2 <- list.files(inFolder2, full.names = T)
all_files2 <- all_files2[grep(matchPat, all_files2)]
stopifnot(length(all_files2) > 0)

all_chr1 <- unique(gsub(".+_(chr.+?)_.+", "\\1", basename(all_files1)))
all_chr2 <- unique(gsub(".+_(chr.+?)_.+", "\\1", basename(all_files2)))

if(all(regexpr("_", all_chr2) > 0)) {
  all_chr2 <- unique(gsub("^(chr.+?)_.+", "\\1", basename(all_files2)))
}
if(all(regexpr("_", all_chr1) > 0)) {
  all_chr1 <- unique(gsub("^(chr.+?)_.+", "\\1", basename(all_files1)))
}
stopifnot(!any(regexpr("_", all_chr1)>0))
stopifnot(!any(regexpr("_", all_chr2)>0))

all_chr <- intersect(all_chr1, all_chr2)

stopifnot(length(all_chr) > 0)

chrOrder <- paste0("chr", as.character(sort(as.numeric(as.character(gsub("chr", "", all_chr))))))

evenChromo <- chrOrder[c(2*(1:(length(chrOrder)/2)))]

cat(paste0("... from folder1: ", length(all_chr), "/", length(all_chr1), " matching chromosomes\n"))
cat(paste0("... from folder2: ", length(all_chr), "/", length(all_chr2), " matching chromosomes\n"))

cond1 <- paste0(cell_line1, "_", treat1)
if(is.null(treat2)) {
  cond2 <- paste0(cell_line2)
}else {
  cond2 <- paste0(cell_line2, "_", treat2) 
}

all_chr_cptmtScores <- foreach(chromo = all_chr, .combine='rbind') %do% {
  inFile1 <- all_files1[grep(paste0("_", chromo, "_"), basename(all_files1))]
  inFile1 <- inFile1[grep("_binPCscore.txt$", inFile1)]
  if(length(inFile1) == 0) {
    inFile1 <- all_files1[grep(paste0("^", chromo, "_"), basename(all_files1))]
    inFile1 <- inFile1[grep("_binPCscore.txt$", inFile1)]
  }
  stopifnot(length(inFile1) == 1)
  stopifnot(file.exists(inFile1))
  inFile2 <- all_files2[grep(paste0("_", chromo, "_"), basename(all_files2))]
  inFile2 <- inFile2[grep("_binPCscore.txt$", inFile2)]
  if(length(inFile2) == 0) {
    inFile2 <- all_files2[regexpr(paste0("^", chromo, "_"), basename(all_files2))>0]
    inFile2 <- inFile2[grep("_binPCscore.txt$", inFile2)]
  }
  stopifnot(length(inFile2) == 1)
  stopifnot(file.exists(inFile2))
  
  pcDT1 <- read.delim(inFile1, header=T,stringsAsFactors = F)
  pcDT2 <- read.delim(inFile2, header=T,stringsAsFactors = F)
  stopifnot(all(colnames(pcDT1) == colnames(pcDT2)))
  stopifnot( (pcDT1$binEnd[1]-pcDT1$binStart[1]) == (pcDT2$binEnd[1]-pcDT2$binStart[1]))
  
  minCoord <- min(c(nrow(pcDT1), nrow(pcDT2)))
  pcDT1$condition <- cond1
  pcDT2$condition <- cond2
  pcDT1 <- pcDT1[1:minCoord,]
  pcDT2 <- pcDT2[1:minCoord,]
  
  stopifnot(all(pcDT1$start == pcDT2$start))
  stopifnot(all(pcDT1$end == pcDT2$end))
  rawCor <- cor(pcDT1$PCscore[!is.na(pcDT1$PCscore) & !is.na(pcDT2$PCscore)],pcDT2$PCscore[!is.na(pcDT1$PCscore) & !is.na(pcDT2$PCscore)])
  changeSignCor <- cor(pcDT1$PCscore[!is.na(pcDT1$PCscore) & !is.na(pcDT2$PCscore)],-pcDT2$PCscore[!is.na(pcDT1$PCscore) & !is.na(pcDT2$PCscore)])
  if(changeSignCor > rawCor) {
    pcDT2$PCscore <- -1*pcDT2$PCscore
    pcDT2$sign <- "changed"
  } else {
    pcDT2$sign <- "unchanged"
  }
  pcDT1$sign <- "unchanged"
  rbind(pcDT1, pcDT2)
}

all_chr_cptmtCoord <- foreach(chromo = all_chr, .combine='rbind') %dopar% {
  inFile1 <- all_files1[grep(paste0("_", chromo, "_"), all_files1)]
  inFile1 <- inFile1[grep("_cptmtCoord.txt$", inFile1)]
  if(length(inFile1) == 0) {
    inFile1 <- all_files1[regexpr(paste0("^", chromo, "_"), basename(all_files1))>0]
    inFile1 <- inFile1[grep("_cptmtCoord.txt$", inFile1)]
  }
  stopifnot(length(inFile1) == 1)
  stopifnot(file.exists(inFile1))
  
  inFile2 <- all_files2[grep(paste0("_", chromo, "_"), all_files2)]
  inFile2 <- inFile2[grep("_cptmtCoord.txt$", inFile2)]
  if(length(inFile2) == 0) {
    inFile2 <- all_files2[regexpr(paste0("^", chromo, "_"), basename(all_files2))>0]
    inFile2 <- inFile2[grep("_cptmtCoord.txt$", inFile2)]
  }
  stopifnot(length(inFile2) == 1)
  stopifnot(file.exists(inFile2))
  
  pcDT1 <- read.delim(inFile1, header=T,stringsAsFactors = F)
  pcDT2 <- read.delim(inFile2, header=T,stringsAsFactors = F)
  
  stopifnot(all(colnames(pcDT1) == colnames(pcDT2)))
  stopifnot( (pcDT1$binEnd[1]-pcDT1$binStart[1]) == (pcDT2$binEnd[1]-pcDT2$binStart[1]))
  # maxCoord <- max(c(nrow(pcDT1), nrow(pcDT2)))
  
  # update the max end with the one retrieved
  maxCoord <- max(all_chr_cptmtScores$binEnd[all_chr_cptmtScores$chromo == chromo], na.rm=T)
  
  pcDT1 <- pcDT1[pcDT1$start < maxCoord,]
  pcDT2 <- pcDT2[pcDT2$start < maxCoord,]
  
  pcDT1$end[nrow(pcDT1)] <- min(c(pcDT1$end[nrow(pcDT1)], maxCoord))
  pcDT2$end[nrow(pcDT2)] <- min(c(pcDT2$end[nrow(pcDT2)], maxCoord))
  
  pcDT1$condition <- cond1
  pcDT2$condition <- cond2
  changeSign <- unique(all_chr_cptmtScores$sign[all_chr_cptmtScores$chromo == chromo & all_chr_cptmtScores$condition== cond2])
  stopifnot(length(changeSign) == 1)
  if(changeSign == "changed") {
    pcDT2$cptmtLab <- ifelse(regexpr("neg", pcDT2$cptmtLab) > 0, gsub("neg", "pos", pcDT2$cptmtLab),gsub("pos", "neg", pcDT2$cptmtLab))
    pcDT2$region <- ifelse(regexpr("neg", pcDT2$region) > 0, gsub("neg", "pos", pcDT2$region),gsub("pos", "neg", pcDT2$region))
  }
  rbind(pcDT1, pcDT2)
}

cond1_txt <- gsub("_", " ",cond1)
cond2_txt <- gsub("_", " ",cond2)
binSize <- unique(all_chr_cptmtScores$binEnd-all_chr_cptmtScores$binStart+1)
stopifnot(length(binSize) == 1)
binSize <- binSize/1000

maxChromoPositions <- c(by(all_chr_cptmtScores,all_chr_cptmtScores$chromo, function(x) max(x$binEnd,na.rm=T)))
chromPos <- data.frame(chromo = chrOrder, pos = maxChromoPositions[chrOrder])

chromPos_cumsum <- data.frame(
                  chromo = chrOrder,
                  posShift = c(0,cumsum(as.numeric(chromPos$pos))[-nrow(chromPos)]),
                  stringsAsFactors = F
                    )

###### ADJUST START AND END GENOME WIDE FOR PC SCORES
all_chr_cptmtScores$chromo <- factor(as.character(all_chr_cptmtScores$chromo),levels = chrOrder)
all_chr_cptmtScores$adjStart <- unlist(sapply(1:nrow(all_chr_cptmtScores), function(i) {
  all_chr_cptmtScores$binStart[i] + chromPos_cumsum$posShift[as.character(chromPos_cumsum$chromo) == as.character(all_chr_cptmtScores$chromo[i])]
}))
all_chr_cptmtScores$adjEnd <- all_chr_cptmtScores$adjStart + all_chr_cptmtScores$binEnd - all_chr_cptmtScores$binStart
stopifnot(all(all_chr_cptmtScores$adjEnd[all_chr_cptmtScores$chromo == "chr1"] == all_chr_cptmtScores$binEnd[all_chr_cptmtScores$chromo == "chr1"] ))
stopifnot(all(all_chr_cptmtScores$adjStart[all_chr_cptmtScores$chromo == "chr1"] == all_chr_cptmtScores$binStart[all_chr_cptmtScores$chromo == "chr1"] ))

chrlabelpos <- max(all_chr_cptmtScores$PCscore,na.rm=T)+0.01
yaxmax <- chrlabelpos + 0.02

###### ADJUST START AND END GENOME WIDE FOR COMPARTMENT COORD
all_chr_cptmtCoord$chromo <- factor(as.character(all_chr_cptmtCoord$chromo),levels = chrOrder)
all_chr_cptmtCoord$adjStart <- unlist(sapply(1:nrow(all_chr_cptmtCoord), function(i) {
  all_chr_cptmtCoord$start[i] + chromPos_cumsum$posShift[as.character(chromPos_cumsum$chromo) == as.character(all_chr_cptmtCoord$chromo[i])]
}))
all_chr_cptmtCoord$adjEnd <- all_chr_cptmtCoord$adjStart + all_chr_cptmtCoord$end - all_chr_cptmtCoord$start
stopifnot(all(all_chr_cptmtCoord$adjEnd[all_chr_cptmtCoord$chromo == "chr1"] == all_chr_cptmtCoord$end[all_chr_cptmtCoord$chromo == "chr1"] ))
stopifnot(all(all_chr_cptmtCoord$adjStart[all_chr_cptmtCoord$chromo == "chr1"] == all_chr_cptmtCoord$start[all_chr_cptmtCoord$chromo == "chr1"] ))

xaxisRange <- c(0, max(c(all_chr_cptmtCoord$start, all_chr_cptmtCoord$adjEnd, all_chr_cptmtCoord$adjStart, all_chr_cptmtCoord$adjEnd), na.rm=T))

#################################################################################################### PLOT THE PC SCORES FOR COND1 -- version continuous X
score1DT <- all_chr_cptmtScores[all_chr_cptmtScores$condition == cond1,]
score1DT <- score1DT[order(as.numeric(score1DT$chromo), score1DT$binStart,score1DT$binEnd),]
score1DT$binNbr <- paste0("bin",1:nrow(score1DT))
score1DT$binNbr <- factor(score1DT$binNbr, levels=score1DT$binNbr)
score1DT$chromoN <- ifelse(as.character(score1DT$chromo) %in% evenChromo,"even","odd")
score1DT$midPos <- (score1DT$adjStart+score1DT$adjEnd)/2

#################################################################################################### PLOT THE PC SCORES FOR COND2 -- version continuous X

score2DT <- all_chr_cptmtScores[all_chr_cptmtScores$condition == cond2,]
score2DT <- score2DT[order(as.numeric(score2DT$chromo), score2DT$binStart,score2DT$binEnd),]
score2DT$chromoN <- ifelse(as.character(score2DT$chromo) %in% evenChromo,"even","odd")
score2DT$midPos <- (score2DT$adjStart+score2DT$adjEnd)/2

#################################################################################################### PLOT COR DT1 VS DT2
score1DT_b <- score1DT
score1DT_b <- score1DT_b[,c("chromo", "binStart", "binEnd", "PCscore")]
colnames(score1DT_b) <- c("chromo", "binStart", "binEnd", "PCscore1")
score2DT_b <- score2DT
score2DT_b <- score2DT_b[,c("chromo", "binStart", "binEnd", "PCscore")]
colnames(score2DT_b) <- c("chromo", "binStart", "binEnd", "PCscore2")
corrScoresDT <- inner_join(score1DT_b, score2DT_b, by=c("chromo", "binStart", "binEnd"))

# test cor by chromo
corr_by_chr_DT <- do.call('rbind', by(corrScoresDT, corrScoresDT$chromo, function(x) {
  tmp <- cor.test(x$PCscore2,x$PCscore1)
  setNames(c(round(tmp$estimate,4), round(tmp$p.value, 4)), c("PC_corr_coeff", "PC_corr_pval"))
}))
corr_by_chr_DT <- as.data.frame(corr_by_chr_DT)
corr_by_chr_DT$chromo <- rownames(corr_by_chr_DT)

#################################################################################################### PUT TOGETHER PC SCORES PLOT FOR COND1 AND COND2
############################ comput how many change
score1DT_tmp <- score1DT
score2DT_tmp <- score2DT
score1DT_tmp$bin_id <- paste0(score1DT_tmp$binStart,"_",score1DT_tmp$binEnd)
score2DT_tmp$bin_id <- paste0(score2DT_tmp$binStart,"_",score2DT_tmp$binEnd)

score1DT_tmp$sign_cond1 <- sign(score1DT_tmp$PCscore)
score2DT_tmp$sign_cond2 <- sign(score2DT_tmp$PCscore)

score1DT_tmp <- na.omit(score1DT_tmp)
score2DT_tmp <- na.omit(score2DT_tmp)

interDT <- inner_join(score1DT_tmp[,c("chromo", "bin_id", "sign_cond1")],score2DT_tmp[,c("chromo", "bin_id", "sign_cond2")], by =c("chromo", "bin_id"))

binsChangingCptmt <- sum(interDT$sign_cond1 != interDT$sign_cond2)
interDT_tmp <- interDT
interDT_tmp$changing <- as.numeric(interDT_tmp$sign_cond1 != interDT_tmp$sign_cond2)

binsChangingCptmt_chr <- inner_join( aggregate(changing ~chromo, data=interDT_tmp, FUN=length),
            aggregate(changing ~chromo, data=interDT_tmp, FUN=sum), by="chromo")
colnames(binsChangingCptmt_chr) <- c("chromo", "nbr_bins", "nbr_changing")
stopifnot(all(binsChangingCptmt_chr$nbr_bins >= binsChangingCptmt_chr$nbr_changing))
stopifnot(sum(binsChangingCptmt_chr$nbr_changing) == binsChangingCptmt)

# write in file the number of bins changing compartments
outDT <- data.frame(
          cell_line1 = cell_line1,
          treat1 = treat1,
          cell_line2 = cell_line2,
          treat2 = ifelse(is.null(treat2), NA, treat2),
          nbr_bins=nrow(interDT),
					nbr_changing = binsChangingCptmt,
					PC_corr_coeff = round(testcor$estimate,4),
					PC_corr_pval = round(testcor$p.value, 4)
					)

# write in file the number of bins changing compartments
outFile <- paste0(outFold, "/", "gw_bins_changing_pcscore_corr.txt")
write.table(outDT, file = outFile, col.names=T, row.names=F, quote=F, sep="\t")
cat(paste0("... written: ", outFile, "\n"))

stopifnot(all(corr_by_chr_DT$chromo %in% binsChangingCptmt_chr$chromo ))
stopifnot(all(binsChangingCptmt_chr$chromo %in% corr_by_chr_DT$chromo ))

outDT_chr <- inner_join(corr_by_chr_DT, binsChangingCptmt_chr, by = c("chromo"))
outDT_chr$cell_line1 <- cell_line1
outDT_chr$treat1 <- treat1
outDT_chr$cell_line2 <- cell_line2
if(is.null(treat2)){
  outDT_chr$treat2 <- NA
} else {
  outDT_chr$treat2 <- treat2
}

outDT_chr <- outDT_chr[,c("chromo", "cell_line1", "treat1", "cell_line2", "treat2", 
                          "nbr_bins", "nbr_changing", "PC_corr_coeff", "PC_corr_pval")]

# write in file the number of bins changing compartments
outFile <- paste0(outFold, "/", "all_chr_bins_changing_pcscore_corr.txt")
write.table(outDT_chr, file = outFile, col.names=T, row.names=F, quote=F, sep="\t")
cat(paste0("... written: ", outFile, "\n"))


