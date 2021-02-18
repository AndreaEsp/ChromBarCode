# -*- coding: utf-8 -*-

# The following script predicts the location of the binding sites of odd-numbered
# chromosomes on the basis of their histone mark profiles.

import numpy as np
import sys,os
from scipy import stats 
import matplotlib.pyplot as plt
import pandas as pd

data_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))+os.sep+'Data'

ncol   = 30         # The number of binding domains per chromosome
nclass = 9          # The number of epigenetic classes
exp_resol = 500 	# The resolution (in bp) of the considered histone marks
pol_resol = 5000	# The resolution (in bp) of the sbs polymers
gray_bead = 5		# Number of gray sites per bin

F_name = ["H3K4me3","H3K4me1","H3K36me3","H3K9me3","H3K27me3"]
odd_chrs  = [i for i in range(1,22,2)]
class_order = [1,2,0,3,4,5,6,7,8]

r = int(pol_resol/exp_resol)

centroid = np.genfromtxt("9-centroid_matrix.txt")
centroid = centroid[class_order]

for chr in odd_chrs:
	pol = np.genfromtxt(data_path+os.sep+'sbs_polymers'+os.sep+'best_polymer_chr%d.txt'%(chr))
	null_bins = np.where(np.sum(pol,1)==0)[0]
	pol_new = np.zeros((len(pol),nclass+1))

	feature = []
	for k in range(len(F_name)):
		feature.append(np.genfromtxt('chr%d/%s.bed'%(chr,F_name[k]),usecols=4))
	feature = np.array(feature).T

	null_idx = []
	cors = []
	for i in range(len(feature)):
		if np.isnan(np.max(np.corrcoef(feature[i],centroid)[0][1:])) == False:
			cors.append(np.max(np.corrcoef(feature[i],centroid)[0][1:]))
			k = np.argmax(np.corrcoef(feature[i],centroid)[0][1:]) + 1
			pol_new[int(i/r),k] += 1
		else:
			null_idx.append(i)
			pol_new[int(i/r),0] += 1

	pol_new[:,0] += gray_bead
	
	#Exclude the non-mappable bins
	pol_new[null_bins] = 0

	# Uncomment the following line to save the data
	np.savetxt("pol_predicted_chr%d.txt"%(chr), pol_new, fmt='%d')
