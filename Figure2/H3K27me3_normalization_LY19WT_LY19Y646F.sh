#!/bin/bash

#### H3K27me3_normalization_LY19WT_LY19Y646F.sh
# Pipeline for jointly-normalizing broad histone mark ChIP-seq data across two datasets (LY19WT and LY19Y646F).
# Input: see "INPUT" section.
# Output: normalized ChIP-seq tracks for LY19WT and LY19Y646F
# Refer to Donaldson, Sungalee, Zufferey, Tavernari et al. for details.
# 
# Requires: R, Python (numpy, pandas), bedtools
# 
# Contact daniele.tavernari@unil.ch for enquiries

############ INPUT ############
LY19WT_alignments_bed="Datadir/LY19WT_DMSO_H3K27me3-rep0_alignments_sorted_mapq1.bed" # Alignment of ChIP-seq reads in .bed format for LY19WT
LY19Y646F_alignments_bed="Datadir/LY19Y646F_DMSO_H3K27me3-rep0_alignments_sorted_mapq1.bed" # Alignment of ChIP-seq reads in .bed format for LY19Y646F
LY19WT_sicer_peaks_bed="Datadir/LY19WT_DMSO_rep0_SICERpeaks_FDR0.001.bed" # SICER peaks for LY19WT
LY19Y646F_sicer_peaks_bed="Datadir/LY19Y646F_DMSO_rep0_SICERpeaks_FDR0.001.bed" # SICER peaks for LY19Y646F
binsize=100 # interval size for binning of ChIP-seq values 
cutoff_top=0.1 # cutoff_top*100 is the top percentile of most populated bins used to test enrichment of peaks 
cutoff_pval=0.00001 # p-value cutoff to retain peaks enriched in most populated bins
firstChr=1 # chromosomes on which performing the normalization
lastChr=22
Outdir="Outdir/" # folder in which results should be saved
Scriptsdir="Scriptsdir/" # folder where this scripts and the other functions are stored
###############################

echo "Broad histone mark normalization pipeline"

echo "--- Same-depth normalization"

for chr in $(seq ${firstChr} ${lastChr})
do

	echo "--------- Processing chr"${chr}
	
	binnedFile=${Outdir}/"chr"${chr}_"binsize"${binsize}".bed"
	${Scriptsdir}/BinChr.py -o ${binnedFile} -c ${chr} -b ${binsize}

	echo "------------ Binning ChIP data for LY19WT"

	Cellline="LY19WT"

	# Remove duplicates
	filein=${LY19WT_alignments_bed}
	fileout_dup=${Outdir}/${Cellline}"_withdup_chr"${chr}".bed"
	fileout_nodup_wt=${Outdir}/${Cellline}"_nodup_chr"${chr}".bed"
	awk -v chr="chr${chr}" '{if ($1==chr) {print $1"\t"$2"\t"$3"\t"$6}}' ${filein} > ${fileout_dup}
	sort -u -k2,2n -k3,3n -k4,4 ${fileout_dup} > ${fileout_nodup_wt}
	rm ${fileout_dup}

	# Bin ChIP data
	binned_wt=${Outdir}/${Cellline}"_nodup_chr"${chr}"_binsize"${binsize}
	bedtools coverage -a ${binnedFile} -b ${fileout_nodup_wt} > ${binned_wt}

	echo "------------ Binning ChIP data for LY19Y646F"

	Cellline="LY19Y646F"

	# Remove duplicates
	filein=${LY19Y646F_alignments_bed}
	fileout_dup=${Outdir}/${Cellline}"_withdup_chr"${chr}".bed"
	fileout_nodup_mut=${Outdir}/${Cellline}"_nodup_chr"${chr}".bed"
	awk -v chr="chr${chr}" '{if ($1==chr) {print $1"\t"$2"\t"$3"\t"$6}}' ${filein} > ${fileout_dup}
	sort -u -k2,2n -k3,3n -k4,4 ${fileout_dup} > ${fileout_nodup_mut}
	rm ${fileout_dup}

	# Bin ChIP data
	binned_mut=${Outdir}/${Cellline}"_nodup_chr"${chr}"_binsize"${binsize}
	bedtools coverage -a ${binnedFile} -b ${fileout_nodup_mut} > ${binned_mut}

	
	echo "------------ Same-depth normalizing ..."

	# Same-depth normalization
	${Scriptsdir}/SameDepthNormalize.py -c ${chr} -w ${binned_wt} -m ${binned_mut} -W ${fileout_nodup_wt} -M ${fileout_nodup_mut}

	echo "------------ Same-depth normalization done"

	echo "------------ Broad histone mark normalization"

	LY19WT_bdg=${binned_wt}"sameDepth.bdg"
	LY19Y646F_bdg=${binned_mut}"sameDepth.bdg"

	Rscript ${Scriptsdir}/Broad_histone_mark_normalization.R ${LY19WT_bdg} ${LY19Y646F_bdg} ${LY19WT_sicer_peaks_bed} ${LY19Y646F_sicer_peaks_bed} ${cutoff_top} ${cutoff_pval} ${chr} ${Outdir}

	rm ${binnedFile} ${fileout_nodup_wt} ${binned_wt} ${binned_wt}"sameDepth.bdg" ${fileout_nodup_mut} ${binned_mut} ${binned_mut}"sameDepth.bdg"

done


echo "Done"
