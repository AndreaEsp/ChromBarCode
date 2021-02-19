# -*- coding: utf-8 -*-

# The following script computes Pearson correlations among the binding domains and a the of chromatin
# marks and also develops a control model to assess the statistical significance of such correlations.

# Chromatin marks data have been downloaded from ENCODE. See the Data folder for the list of ENCODE 
# accession numbers.

import numpy as np
import sys,os
from scipy.stats.mstats import pearsonr as masked_pearsonr

def build_random_model(binding_domain):
	'''Randomize a binding domain by bootstrapping'''
	tot = np.sum(binding_domain)
	random_model = np.zeros(binding_domain.shape)
	while tot != 0:
		random_model[np.random.randint(len(random_model))] += 1
		tot -= 1
	return random_model

def calculate_masked_correlation(epigenetic_signal, random_model):
	'''Compute the Pearson correlation masking the bin with zero binding sites'''
	random_model_masked =  np.ma.masked_where(random_model==0, random_model, copy=True)
	epigenetic_signal_masked = np.ma.masked_where(random_model==0, epigenetic_signal, copy=True)
	r = masked_pearsonr(random_model_masked,epigenetic_signal_masked)[0]
	return r

replicas = 100              # Number of random models for each binding domain
perc_inf,perc_sup = 5,95    # Significance threshold 

chrs = [i for i in range(2,23,2)]
F_name = ["H3K4me1","H3K4me2","H3K4me3","H3K9ac","H3K27ac","H3K36me3","H3K79me2","H4K20me1","H3K9me3","H3K27me3","H2AFZ","CTCF","POLR2AphosphoS2","POLR2AphosphoS5","EZH2","STAT1","ZNF687","RUNX3","CEBPZ","ZNF217","TCF12","ZNF592","USF1","HCFC1","MXI1","HSF1","ZEB1","MYB","EP300","KLF5","CEBPB","TAF1","RAD21","ZNF384","CUX1","ELF1","RB1","ZBTB33","TBP","NFYB","LARP7","WRNIP1","JUNB","EBF1","BHLHE40","STAT3","SIN3A","MEF2C","NKRF","PBX3","MTA3","ESRRA","NR2F1","ZSCAN29","NR2C1","NFXL1","IKZF1","IRF3","TBL1XR1","PAX5","HDAC2","KDM1A","YY1","CHD4","RXRA","DPF2","SMAD1","MAZ","NR2C2","SIX5","CHD1","ZFP36","GATAD2B","BCL3","SRF","RBBP5","BCLAF1","BRCA1","ZNF24","SMARCA5","ZNF207","NFATC1","REST","UBTF","RELB","ZNF622","SMC3","NFIC","CBX5","FOXK2","CBFB","CBX3","HDAC6","MEF2B","PAX8","MTA2","ATF7","TBX21","E2F8","ETS1","IRF5","ELK1","NBN","IRF4","MAX","RFX5","CHD2","BATF","STAT5A","JUND","BCL11A","ARID3A","ZNF143","BACH1","NRF1","SMAD5","TARDBP","SPI1","RCOR1","TRIM22","NFATC3","RAD51","GABPA","BMI1","PKNOX1","YBX1","ASH2L","EED","SKIL","E2F4","ZBTB40","SUZ12","ATF2","MAFK","NFYA","ETV6","EGR1","USF2","E4F1","MLLT1","IKZF2","TCF7","CREM","ZBED1","POLR2A","ARNT","HDGF","MEF2A","ZZZ3","DNase","RepG1","RepG2","RepS1","RepS2","RepS3","RepS4","RNAseq_total"]

# Reading histone marks data
feature = []
for chr in chrs:
	vals = []
	for k in range(len(F_name)):
		vals.append(np.genfromtxt("chr%d/%s.bed"%(chr,F_name[k]),usecols=4))
	feature.append(np.array(vals))

# Control model
rand_corr = []
for chr,m in zip(chrs,np.arange(len(chrs))):
	print("Chomosome %d\n"%chr)
	data_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))+os.sep+'Data'
	pol = np.genfromtxt(data_path+os.sep+'sbs_polymers'+os.sep+'best_polymer_chr%d.txt'%(chr))
	for i in range(1,len(pol[0])):
		for j in range(replicas):
			rand_model = build_random_model(pol[:,i])
			r = []
			for k in range(len(F_name)):
				r.append( calculate_masked_correlation(feature[m][k],rand_model) )
			rand_corr.append( np.array(r) )
rand_corr = np.array(rand_corr)

# Uncomment the following line to save the data
np.savetxt("control_sample.txt",rand_corr,fmt='%.8f')

# Actual correlations
actual_corr = []
for chr,m in zip(chrs,np.arange(len(chrs))):
	print("Chomosome %d\n"%chr)
	data_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))+os.sep+'Data'
	pol = np.genfromtxt(data_path+os.sep+'sbs_polymers'+os.sep+'best_polymer_chr%d.txt'%(chr))
	for i in range(1,len(pol[0])):
		r = []
		for k in range(len(F_name)):
			r.append( calculate_masked_correlation(feature[m][k],pol[:,i]) )
		actual_corr.append( np.array(r) )
actual_corr = np.array(actual_corr)

# Uncomment the following line to save the data
#np.savetxt("actual_corr.txt",actual_corr,fmt='%.8f')

actual_corr[np.where(np.isnan(actual_corr)==True)] = 0
rand_corr[np.where(np.isnan(rand_corr)==True)] = 0
significant_sample = np.zeros(actual_corr.shape)

# Retain only the significant correlation values
for i in np.arange(0,rand_corr.shape[1]):
    inf=scipy.percentile(rand_corr[:,i],perc_inf)
    sup=scipy.percentile(rand_corr[:,i],perc_sup)
    (significant_sample[:,i])[np.where((actual_corr[:,i]<=inf))] = (actual_corr[:,i])[np.where((actual_corr[:,i]<=inf))]
    (significant_sample[:,i])[np.where((actual_corr[:,i]>=sup))] = (actual_corr[:,i])[np.where((actual_corr[:,i]>=sup))]

# Uncomment the following line to save the data
np.savetxt("significant_sample.txt",significant_sample,fmt='%.8f')
