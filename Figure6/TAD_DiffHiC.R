#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
c1 = args[1] # Cell line / condition 1
c2 = args[2] # Cell line / condition 2
Interactome_c1 = args[3] # Interactome tables for c1 - HiC-DC results with columns as c("binI", "binJ", "counts", "D", "pvalue", "qvalue")
Interactome_c2 = args[4] # Interactome tables for c2 - HiC-DC results with columns as c("binI", "binJ", "counts", "D", "pvalue", "qvalue")
chr = args[5] # Chromosome
TAD_start = as.numeric(args[6]) # 1-based coordinate of TAD start (e.g. 20001)
TAD_end = as.numeric(args[7]) # 0-based coordinate of TAD end (e.g. 80000)
binsize = as.numeric(args[8]) # bin size used for HiC-DC (fixed bin mode)
pval_adj_thresh = as.numeric(args[9]) # FDR cutoff to call significant changes
OutFile = args[10] # File name for output table

#### TAD_DiffHiC.R
# Script to call significantly different HiC interactions occurring in a given genomic domain (e.g. a TAD) between two cell lines / conditions.
# Input: see command line arguments.
# Output: A table containing, for each pair of Hi-C bins (rows) the following information (columns)
# - chr, binI, binJ, binI_start, binI_end, binJ_start, binJ_end: location of the two bins
# - minus_Log10_qval_HiCDC_c1, minus_Log10_qval_HiCDC_c2: interactome scores (computed as HiC-DC's -log10(q-values)) for the two conditions
# - c1_minus_c2: difference between the interactome scores 
# - pval_difference: p-value of the difference 
# - qval_difference: q-value of the difference
# - stronger_in: '1' iff (c1_minus_c2>0) & (qval_difference<pval_adj_thresh) & (qval_HiCDC_c1<pval_adj_thresh); vice-versa, '2'; '0' if there was no significant difference 
# Refer to Donaldson, Sungalee, Zufferey, Tavernari et al. for details.
# 
# Contact daniele.tavernari@unil.ch for enquiries

SestrinTad_start = floor(TAD_start/binsize)+1 # in 20kb: first bin is 108840001-108860000
SestrinTad_end = floor(TAD_end/binsize) # in 20 kb: last bin is 109420001-109440000

###################

############ INTERACTOME MATRICES ############

Interactome1 <- read.table(Interactome_c1, quote = '', stringsAsFactors=FALSE, sep = "\t", header = TRUE)

Interactome2 <- read.table(Interactome_c2, quote = '', stringsAsFactors=FALSE, sep = "\t", header = TRUE)

interacs = Interactome1
colnames(interacs) <- c("binI", "binJ", "counts", "D", "pvalue", "qvalue")
interacs[is.na(interacs[,"qvalue"]),"qvalue"] = 1

interacs$qvalue_original = interacs$qvalue
interacs$qvalue[interacs$qvalue==0] <- min(interacs$qvalue[interacs$qvalue!=0]) # Clipping to avoid singularities

interacs_f1 <- data.frame(qval_score = -log10(interacs$qvalue), binI = NA, binJ = NA)

splitted_binI <- strsplit(interacs$binI,"-")
splitted_binJ <- strsplit(interacs$binJ,"-")
interacs_f1$binI <- sapply(1:length(splitted_binI), function(k) floor(as.numeric(splitted_binI[[k]][2])/binsize)+1)
interacs_f1$binJ <- sapply(1:length(splitted_binJ), function(k) floor(as.numeric(splitted_binJ[[k]][2])/binsize)+1)
interacs_f1[,"binI_binJ"] = paste0(interacs_f1[,"binI"], "_", interacs_f1[,"binJ"])
colnames(interacs_f1) = c("qvalS_f1", "binI", "binJ", "binI_binJ" )


interacs <- Interactome2
colnames(interacs) <- c("binI", "binJ", "counts", "D", "pvalue", "qvalue")
interacs[is.na(interacs[,"qvalue"]),"qvalue"] = 1

interacs$qvalue_original = interacs$qvalue
interacs$qvalue[interacs$qvalue==0] <- min(interacs$qvalue[interacs$qvalue!=0]) # Clipping to avoid singularities

interacs_f2 <- data.frame(qval_score = -log10(interacs$qvalue), binI = NA, binJ = NA)

splitted_binI <- strsplit(interacs$binI,"-")
splitted_binJ <- strsplit(interacs$binJ,"-")
interacs_f2$binI <- sapply(1:length(splitted_binI), function(k) floor(as.numeric(splitted_binI[[k]][2])/binsize)+1)
interacs_f2$binJ <- sapply(1:length(splitted_binJ), function(k) floor(as.numeric(splitted_binJ[[k]][2])/binsize)+1)
interacs_f2[,"binI_binJ"] = paste0(interacs_f2[,"binI"], "_", interacs_f2[,"binJ"])
colnames(interacs_f2) = c("qvalS_f2", "binI", "binJ", "binI_binJ" )

interacs_f12 = merge(interacs_f1, interacs_f2, by = "binI_binJ")
interacs_f12$binI.y = NULL
interacs_f12$binJ.y = NULL
colnames(interacs_f12) = c( "binI_binJ", "qvalS_f1","binI", "binJ","qvalS_f2" )
interacs_f12$qvalS_f1 = as.numeric(interacs_f12$qvalS_f1)
interacs_f12$qvalS_f2 = as.numeric(interacs_f12$qvalS_f2)

interacs_f12[,"qvalS_diff"] = (interacs_f12[,"qvalS_f1"]-interacs_f12[,"qvalS_f2"])

interacs_sestrin_f12 = interacs_f12[( (interacs_f12$binI>=SestrinTad_start) & (interacs_f12$binI<=SestrinTad_end) & (interacs_f12$binJ>=SestrinTad_start) & (interacs_f12$binJ<=SestrinTad_end) & (interacs_f12$binI!=interacs_f12$binJ) ),]

##############################################

############## STATISTICAL TESTS #############
pval_emp <- function(query, bg, alternative){
   if (alternative=="greater")
   {
      return( sum(1*(bg>=query))/length(bg) )
   }
   if (alternative=="less")
   {
      return( sum(1*(bg<=query))/length(bg) )
   }
}

# 1 minus 2
interacs_sestrin_f12[,"pval"] = NA
interacs_sestrin_f12[,"pval_adj"] = NA
interacs_sestrin_f12[,"has_adjSignif_diff_and_interacs"] = 0

for (diag in 1:(SestrinTad_end-SestrinTad_start) )
{
   interacs_f12_diag = interacs_f12[(interacs_f12[,"binJ"]-interacs_f12[,"binI"])==diag,]
   bg = interacs_f12_diag[,"qvalS_diff"] # 
   # print(length(bg))
   interacs_sestrin_f12_diag = interacs_sestrin_f12[(interacs_sestrin_f12[,"binJ"]-interacs_sestrin_f12[,"binI"])==diag,]
   for (roww in rownames(interacs_sestrin_f12_diag) )
   {
      query = interacs_sestrin_f12_diag[roww,"qvalS_diff"]
      # print(query)
      pos_query = interacs_sestrin_f12_diag[roww,"binI_binJ"]
      if (query>=0)
      {
         interacs_sestrin_f12[interacs_sestrin_f12[,"binI_binJ"]==pos_query,"pval"] = pval_emp(query, bg, alternative = "greater")
      }
      else {
         interacs_sestrin_f12[interacs_sestrin_f12[,"binI_binJ"]==pos_query,"pval"] = pval_emp(query, bg, alternative = "less")
      }
      
      # print(sum(is.na(bg)))
   }
}

# interacs_sestrin_f12 = interacs_sestrin_f12[!is.na(interacs_sestrin_f12$pval),]

interacs_sestrin_f12[,"pval_adj"] = p.adjust( interacs_sestrin_f12[,"pval"], method = "BH" )
interacs_sestrin_f12[ (interacs_sestrin_f12[,"pval_adj"]<pval_adj_thresh) & ( (interacs_sestrin_f12[,"qvalS_diff"]>0) & (interacs_sestrin_f12[,"qvalS_f1"]>(-log10(pval_adj_thresh))) ) ,"has_adjSignif_diff_and_interacs"] = 1
interacs_sestrin_f12[ (interacs_sestrin_f12[,"pval_adj"]<pval_adj_thresh) & ( (interacs_sestrin_f12[,"qvalS_diff"]<0) & (interacs_sestrin_f12[,"qvalS_f2"]>(-log10(pval_adj_thresh))) ) ,"has_adjSignif_diff_and_interacs"] = 2

tabtosave = interacs_sestrin_f12
tabtosave$binI_binJ=NULL
colnames(tabtosave) = c("minus_Log10_qval_HiCDC_c1", "binI", "binJ", "minus_Log10_qval_HiCDC_c2", "c1_minus_c2", "pval_difference", "qval_difference", "stronger_in" )
tabtosave$chr = chr
tabtosave$binI_end = tabtosave$binI*binsize
tabtosave$binI_start = tabtosave$binI_end-binsize+1
tabtosave$binJ_end = tabtosave$binJ*binsize
tabtosave$binJ_start = tabtosave$binJ_end-binsize+1
tabtosave$binI = tabtosave$binI-SestrinTad_start+1
tabtosave$binJ = tabtosave$binJ-SestrinTad_start+1
tabtosave$c1 = c1
tabtosave$c2 = c2
tabtosave$binI_end = format(tabtosave$binI_end, scientific = F)
tabtosave$binI_start = format(tabtosave$binI_start, scientific = F)
tabtosave$binJ_end = format(tabtosave$binJ_end, scientific = F)
tabtosave$binJ_start = format(tabtosave$binJ_start, scientific = F)

tabtosave2 = tabtosave[,c( "c1", "c2", "chr", "binI", "binI_start", "binI_end", "binJ", "binJ_start", "binJ_end", "minus_Log10_qval_HiCDC_c1", "minus_Log10_qval_HiCDC_c2", "c1_minus_c2", "pval_difference", "qval_difference", "stronger_in" )]
write.table(tabtosave2, file = OutFile, row.names = F, col.names = T, quote = F, sep = "\t" )

##############################################
