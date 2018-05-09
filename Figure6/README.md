**TAD_DiffHiC.R** 

This script computes differences between Hi-C interactomes (computed by HiC-DC [[Carty et al. (2017)]](https://www.nature.com/articles/ncomms15454) as significant Hi-C interactions) and assesses their statistical significance.  
See methods section in Donaldson, Sungalee, Zufferey, Tavernari et al. for details.

_Input_

* c1: Cell line / condition 1
* c2: Cell line / condition 2
* Interactome_c1: Interactome tables for c1 - HiC-DC results with columns as c("binI", "binJ", "counts", "D", "pvalue", "qvalue")
* Interactome_c2: Interactome tables for c2 - HiC-DC results with columns as c("binI", "binJ", "counts", "D", "pvalue", "qvalue")
* chr: Chromosome
* TAD_start: 1-based coordinate of TAD start (e.g. 20001)
* TAD_end: 0-based coordinate of TAD end (e.g. 80000)
* binsize: bin size used for HiC-DC (fixed bin mode)
* pval_adj_thresh: FDR cutoff to call significant changes
* OutFile: File name for output table

_Output_

A table containing, for each pair of Hi-C bins (rows) the following information (columns)
* chr, binI, binJ, binI_start, binI_end, binJ_start, binJ_end: location of the two bins
* minus_Log10_qval_HiCDC_c1, minus_Log10_qval_HiCDC_c2: interactome scores (computed as HiC-DC's -log10(q-values)) for the two conditions
* c1_minus_c2: difference between the interactome scores 
* pval_difference: p-value of the difference 
* qval_difference: q-value of the difference
* stronger_in: '1' iff (c1_minus_c2>0) & (qval_difference<pval_adj_thresh) & (qval_HiCDC_c1<pval_adj_thresh); vice-versa, '2'; '0' if there was no significant difference 
