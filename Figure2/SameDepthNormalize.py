#!/usr/bin/env python

from __future__ import division
import time
import pandas as pd
import subprocess
import sys, getopt


def main(argv):
	start_time1 = time.time()

	try:
		opts, args = getopt.getopt(argv,"hc:w:m:W:M:",["chrom=","wt=","mut=","wtBed=","mutBed="])
	except getopt.GetoptError:
		print 'Usage: SameDepthNormalize.py ... etc'
		sys.exit(2)
	for opt, arg in opts:
		if opt == '-h':
			print 'Usage: SameDepthNormalize.py ... etc'
			sys.exit()
		elif opt in ("-c", "--chrom"):
		    chrom = int(arg)
		elif opt in ("-w", "--wt"):
			wt = arg
		elif opt in ("-m", "--mut"):
			mut = arg
		elif opt in ("-W", "--wtBed"):
			wtBed = arg
		elif opt in ("-M", "--mutBed"):
			mutBed = arg

	# nl_WtDMSO
	inFile = wtBed
	# inFile_df = pd.read_csv(inFile, sep='\t', header=None, names=['chr','start', 'end', 'c4', 'c5', 'c6'], index_col=False)
	nl_WtDMSO = subprocess.check_output(['wc', '-l', inFile]).split(' ')
	nl_WtDMSO = int((nl_WtDMSO[0]))

	# nl_MutDMSO
	inFile = mutBed
	# inFile_df = pd.read_csv(inFile, sep='\t', header=None, names=['chr','start', 'end', 'c4', 'c5', 'c6'], index_col=False)
	nl_MutDMSO = subprocess.check_output(['wc', '-l', inFile]).split(' ')
	nl_MutDMSO = int((nl_MutDMSO[0]))

	mode = 'min'

	if mode == 'min':
		factor_WtDMSO = min([ nl_WtDMSO, nl_MutDMSO ])/nl_WtDMSO
		factor_MutDMSO = min([ nl_WtDMSO, nl_MutDMSO ])/nl_MutDMSO

	if mode == 'max':
		factor_WtDMSO = max([ nl_WtDMSO, nl_MutDMSO ])/nl_WtDMSO
		factor_MutDMSO = max([ nl_WtDMSO, nl_MutDMSO ])/nl_MutDMSO

	# nl_WtDMSO
	binFile = wt
	outFile = wt+"sameDepth.bdg"
	binFile_df = pd.read_csv(binFile, sep='\t', header=None, names=['chr','start', 'end', 'strand', 'nreads', 'c6', 'c7', 'c8'], index_col=False)
	outFile_df = pd.DataFrame(columns=['chr','start', 'end','nreads'])
	outFile_df.loc[:,'chr'] = binFile_df.loc[:,'chr']
	outFile_df.loc[:,'start'] = binFile_df.loc[:,'start']
	outFile_df.loc[:,'end'] = binFile_df.loc[:,'end']
	outFile_df.loc[:,'nreads'] = binFile_df.loc[:,'nreads']*float(factor_WtDMSO)
	outFile_df.to_csv(outFile, sep='\t', header=False, index=False)
	

	# nl_MutDMSO
	binFile = mut
	outFile = mut+"sameDepth.bdg"
	binFile_df = pd.read_csv(binFile, sep='\t', header=None, names=['chr','start', 'end', 'strand', 'nreads', 'c6', 'c7', 'c8'], index_col=False)
	outFile_df = pd.DataFrame(columns=['chr','start', 'end','nreads'])
	outFile_df.loc[:,'chr'] = binFile_df.loc[:,'chr']
	outFile_df.loc[:,'start'] = binFile_df.loc[:,'start']
	outFile_df.loc[:,'end'] = binFile_df.loc[:,'end']
	outFile_df.loc[:,'nreads'] = binFile_df.loc[:,'nreads']*float(factor_MutDMSO)
	outFile_df.to_csv(outFile, sep='\t', header=False, index=False)

	elapsed_time1 = time.time() - start_time1

if __name__ == "__main__":
	import warnings
	warnings.filterwarnings("ignore")
	main(sys.argv[1:])
