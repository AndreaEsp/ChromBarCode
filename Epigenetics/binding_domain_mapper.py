
# -*- coding: utf-8 -*-

# The following script maps a 30-domain SBS polymer to a 9-domain SBS-polymer 
# by summing up all the binding sites belonging to the same epigenetic class 

import numpy as np
import sys,os

ncol   = 30     # The number of binding domains per chromosome
nclass = 9      # The number of epigenetic classes

data_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))+os.sep+'Data'

chrs = np.array([i for i in range(2,23,2)])
labels = np.genfromtxt(data_path+os.sep+'epigenetic_classes'+os.sep+'labels.txt')
labels_per_chr = {'%d'%i:labels[j*ncol:(j+1)*ncol] for i,j in zip(range(2,23,2),range(11))}

for chr in chrs:
	pol_tmp = np.genfromtxt(data_path+os.sep+'sbs_polymers'+os.sep+'best_polymer_chr%d.txt'%(chr))
	cols = len(np.unique(labels_per_chr['%s'%chr]))
	pol = np.zeros((len(pol_tmp),nclass+1))
	pol[:,0] = pol_tmp[:,0]
	for i,j in zip(range(1,cols+1),np.unique(labels_per_chr['%s'%chr])):
		idx = np.where(labels_per_chr['%s'%chr]==j)[0] + 1
		pol[:,i] = np.sum(pol_tmp[:,idx],1)
	# Uncomment the following line to save the data
	np.savetxt("polymer_chr%d_9-domain.txt"%(chr), pol, fmt='%d')
