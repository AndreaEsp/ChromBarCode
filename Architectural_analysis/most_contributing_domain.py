# -*- coding: utf-8 -*-

# The following script computes the most contributing binding domain to each pair interaction

import numpy as np
import sys,os

data_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))+os.sep+'Data'

# Select an even-numbered chromosome
chr = 20

ncol   = 30     # The number of binding domains per chromosome
nclass = 9      # The number of epigenetic classes

labels = np.genfromtxt(data_path+os.sep+'epigenetic_classes'+os.sep+'labels.txt')
labels_per_chr = {'%d'%i:labels[j*ncol:(j+1)*ncol] for i,j in zip(range(2,23,2),range(11))}
case = labels_per_chr['%d'%chr]

pol = np.genfromtxt(data_path+os.sep+'sbs_polymers'+os.sep+'best_polymer_chr%d.txt'%(chr))
N   = len(pol)

# Excluding the first column (gray beads)
pol = pol[:,1:].astype(int)

count_couples = np.zeros(ncol)
count_couples_per_class = np.zeros(nclass)
matrix_first = np.zeros((N,N))
rank_first_contrib = np.zeros(ncol)
tot=0
loc=[]

# Counting couples and map on epigenetic classes
for j in np.arange(0,N):
	for i in np.arange(j,N):
		count_couples=pol[j]*pol[i]
		for k in np.arange(nclass):
			count_couples_per_class[k]=count_couples[np.where(case==k+1)].sum()
		count_couples_sorted_indexes = np.argsort(count_couples_per_class,kind='mergesort')
		
		#First contributing
		if count_couples_per_class[count_couples_sorted_indexes[-1]] != 0:
			if count_couples_per_class[count_couples_sorted_indexes[-1]] == count_couples_per_class[count_couples_sorted_indexes[-2]]:
				loc.append((j,i))
			matrix_first[i,j]= 1+count_couples_sorted_indexes[-1]
		else:
			matrix_first[i,j]= -1
		tot+=1

# Uncomment the following line to save the data
# np.savetxt('1st_contribution.txt',matrix_first,fmt="%d")
