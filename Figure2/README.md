**H3K27me3_normalization_LY19WT_LY19Y646F.sh** 
This script jointly normalizes the H3K27me3 ChIP-seq tracks for LY19WT and LY19Y646F cell lines.
See methods section in Donaldson, Sungalee, Zufferey, Tavernari et al. for details.

_Input_
LY19WT_alignments_bed: alignment of ChIP-seq reads in .bed format for LY19WT
LY19Y646F_alignments_bed: alignment of ChIP-seq reads in .bed format for LY19Y646F
LY19WT_sicer_peaks_bed: SICER peaks for LY19WT
LY19Y646F_sicer_peaks_bed: SICER peaks for LY19Y646F
binsize: interval size for binning of ChIP-seq values 
cutoff_top: cutoff_top*100 is the top percentile of most populated bins used to test enrichment of peaks 
cutoff_pval: p-value cutoff to retain peaks enriched in most populated bins
firstChr: first chromosome on which performing the normalization
lastChr: last chromosome on which performing the normalization
Outdir: folder in which results should be saved
Scriptsdir: folder where this scripts and the other functions are stored

_Output_
Chromosome-wise normalized bedGraph (.bdg) H3K27me3 ChIP-seq tracks for LY19WT and LY19Y646F.  

_Requires_
R, Python (numpy, pandas), bedtools

_Calls_
**BinChr.py**, **SameDepthNormalize.py**, **Broad_histone_mark_normalization.R**
