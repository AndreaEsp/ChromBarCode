# -*- coding: utf-8 -*-

# The following script computes Pearson correlations among the binding domains and a set of chromatin
# marks and also develops a control model to assess the statistical significance of such correlations.

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
F_name = ["H3K4me3","H3K4me1","H3K36me3","H3K9me3","H3K27me3"]

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
