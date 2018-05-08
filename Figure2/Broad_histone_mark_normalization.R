#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
LY19WT_bdg = args[1]
LY19Y646F_bdg = args[2]
LY19WT_peaks = args[3]
LY19Y646F_peaks = args[4]
cutoff_top = as.numeric(args[5])
cutoff_pval = as.numeric(args[6])
chr = as.numeric(args[7])
OutFolder = args[8]

######### Functions #########
# compute read-enrichment for each bin (useful with large windows, e.g. 20kb)
winEnrich = function(bed, win, q = 0.9){
	
	ranked = sort(bed$reads)
	cutoff = ranked[ceiling(q*length(ranked))]
	peaks = sum(bed$reads > cutoff)
	notpeaks = nrow(bed) - peaks
	
	win$p = rep(1,nrow(win))
	for(i in 1:nrow(win)){
		bins = bed[bed$start >= win$start[i] & bed$end <= win$end[i],]
		win.peaks = sum(bins$reads > cutoff)
		win.notpeaks = nrow(bins) - win.peaks
		out.peaks = peaks - win.peaks
		out.notpeaks = notpeaks - win.notpeaks
		m = matrix(c(win.peaks, out.peaks, win.notpeaks, out.notpeaks),nrow=2,byrow=T)
		p = fisher.test(m,alternative="greater")$p.value
		win$p[i] = p
		if(i %% 100 == 0)
			cat(".")
	}

	return(win)
}

sharedBins = function(bed.list){
	bed = bed.list[[1]]
	regions = bed$start
	for(i in 2:length(bed.list)){
		bed = bed.list[[i]]
		r = bed$start
		regions = regions[regions %in% r]
	}
	return(regions)
}

scalingFactorBin = function(bed, regions, q = 0.9){

	ranked = sort(bed$reads)
	cutoff = ranked[ceiling(q*length(ranked))]

	shared.bed = bed[bed$start %in% regions,]
	
	reads = 0
	nbins = 0
	
	bins = shared.bed[shared.bed$reads > cutoff,]
	scaling.factor = sum(bins$reads)/nrow(bins)
	
	return(scaling.factor)
}

# extract bins that overlap with a given set of peaks
peakToBin = function(peaks, bed){
	bed$sig = rep(FALSE,nrow(bed))
	for(i in 1:nrow(peaks)){
		s = bed$start >= peaks$start[i]
		e = bed$e <= peaks$end[i]
		bed$sig[s&e] = TRUE
	}
	bed = bed[bed$sig,]
	return(bed)
}
#############################

Q = 1-cutoff_top

NormFactors_df = data.frame(row.names = paste0("chr", chr), chr = chr, LY19WT_DMSO = NA, LY19Y646F_DMSO = NA )



kh.bin.100 = read.table( LY19WT_bdg, header = FALSE ) 
colnames(kh.bin.100) = c("chr","start","end","reads")
k.peaks.chr = read.table( LY19WT_peaks, header = FALSE )
colnames(k.peaks.chr) = c("chr","start","end","reads", "c1", "c2", "c3", "c4")
k.peaks.chr = k.peaks.chr[k.peaks.chr$chr == paste0("chr",chr),]
k.peaks.chr = k.peaks.chr[,c("chr","start","end","reads")]
kh.bin = winEnrich(bed = kh.bin.100, win = k.peaks.chr, q=Q)
T = -log10(cutoff_pval)
kh.peaks.sig = kh.bin[-log10(kh.bin$p) > T,]
LY19WT_DMSO.bin.sig = peakToBin(peaks = kh.peaks.sig, bed = kh.bin.100)

kh.bin.100 = read.table( LY19Y646F_bdg, header = FALSE ) 
colnames(kh.bin.100) = c("chr","start","end","reads")
# kh.bin.100 = kh.bin.100[kh.bin.100$reads>0,]
k.peaks.chr = read.table( LY19Y646F_peaks, header = FALSE )
colnames(k.peaks.chr) = c("chr","start","end","reads", "c1", "c2", "c3", "c4")
k.peaks.chr = k.peaks.chr[k.peaks.chr$chr == paste0("chr",chr),]
k.peaks.chr = k.peaks.chr[,c("chr","start","end","reads")]
kh.bin = winEnrich(bed = kh.bin.100, win = k.peaks.chr, q=Q)
T = -log10(cutoff_pval)
kh.peaks.sig = kh.bin[-log10(kh.bin$p) > T,]
LY19Y646F_DMSO.bin.sig = peakToBin(peaks = kh.peaks.sig, bed = kh.bin.100)




regions = sharedBins(list(LY19WT_DMSO.bin.sig, LY19Y646F_DMSO.bin.sig))

LY19WT_DMSO.bin.sig_factor = scalingFactorBin(bed = LY19WT_DMSO.bin.sig, regions = regions, q=Q)
LY19Y646F_DMSO.bin.sig_factor = scalingFactorBin(bed = LY19Y646F_DMSO.bin.sig, regions = regions, q=Q)

minF = min(LY19WT_DMSO.bin.sig_factor, LY19Y646F_DMSO.bin.sig_factor)

NormFactors_df[paste0("chr",chr),"LY19WT_DMSO"] = minF/LY19WT_DMSO.bin.sig_factor
NormFactors_df[paste0("chr",chr),"LY19Y646F_DMSO"] = minF/LY19Y646F_DMSO.bin.sig_factor

write.table(NormFactors_df, file = paste0(OutFolder,"NormFactors_chr",chr,"_BHMN.txt") , row.names = F )

Cellline = "LY19WT"
Treatment = "DMSO"
kh.bin.100 = read.table( LY19WT_bdg, header = FALSE )
colnames(kh.bin.100) = c("chr","start","end","reads")
kh.bin.100$norm = kh.bin.100$reads*NormFactors_df[paste0("chr",chr),paste0(Cellline,"_",Treatment) ]
norm.kh.bin.100 = kh.bin.100[,c("chr","start","end","norm")]
write.table(norm.kh.bin.100, paste0(OutFolder,Cellline, "_chr",chr,"_BHM_normalized.bdg"),row.names=F,quote=F,sep="\t",col.names=F)

Cellline = "LY19Y646F"
Treatment = "DMSO"
kh.bin.100 = read.table( LY19Y646F_bdg, header = FALSE )
colnames(kh.bin.100) = c("chr","start","end","reads")
kh.bin.100$norm = kh.bin.100$reads*NormFactors_df[paste0("chr",chr),paste0(Cellline,"_",Treatment) ]
norm.kh.bin.100 = kh.bin.100[,c("chr","start","end","norm")]
write.table(norm.kh.bin.100, paste0(OutFolder,Cellline, "_chr",chr,"_BHM_normalized.bdg"),row.names=F,quote=F,sep="\t",col.names=F)


