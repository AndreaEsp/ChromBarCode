# -*- coding: utf-8 -*-

# The following script computes the first and second most abundant binding site type per bin for a specific chromosome

import numpy as np
import sys,os
import matplotlib
import matplotlib.pyplot as plt
import scipy

#Add path
analysis_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))+os.sep+'Modules'
sys.path.append(analysis_path)

import Chroms_info as c_info

case = {'%s'%c_info.chrs[i]:{'c_start':c_info.c_start[i],'c_end':c_info.c_end[i],'bin':c_info.N_bin[i]} for i in range(len(c_info.chrs))}

ncol   = 30     # The number of binding domains per chromosome
nclass = 9      # The number of epigenetic classes

data_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))+os.sep+'Data'
labels = np.genfromtxt(data_path+os.sep+'epigenetic_classes'+os.sep+'labels.txt')
labels_per_chr = {'%d'%i:labels[j*ncol:(j+1)*ncol] for i,j in zip(range(2,23,2),range(11))}

# Select an even-numbered chromosome
chr = 20

# Map the binding domains on the 9 epigenetic classes
pol_tmp = np.genfromtxt(data_path+os.sep+'sbs_polymers'+os.sep+'best_polymer_chr%d.txt'%(chr))
cols = len(np.unique(labels_per_chr['%s'%chr]))
pol = np.zeros((len(pol_tmp),nclass+1))
pol[:,0] = pol_tmp[:,0]
for i,j in zip(range(1,cols+1),np.unique(labels_per_chr['%s'%chr])):
	idx = np.where(labels_per_chr['%s'%chr]==j)[0] + 1
	pol[:,i] = np.sum(pol_tmp[:,idx],1)

# Select a genomic region (in Mb)
start = case['%s'%chr]['c_start']
end   = case['%s'%chr]['c_end']
N = len(pol)

# Excluding the first column (gray beads)
pol = pol[:,1:].astype(int)

# Read the matrix of first contribution
matrix_first  = np.genfromtxt('1st_contribution.txt')
matrix_first  = np.tril(matrix_first) + np.triu(matrix_first.T, 1)

# First most present class
cols1=[]
cols1_indices=[]
for i in range(N):
	top = pol[i,:].max()
	deg = np.where(pol[i,:]==top)[0]
	if pol[i,deg[0]] == 0:
		cols1.append(0)
		cols1_indices.append(0)
	elif len(deg) == 1:
		cols1.append(deg[0]+1)
		cols1_indices.append(deg[0])
	else:
		multyp=[]
		for j in deg:
			multyp.append(np.where(matrix_first[i]==j+1)[0].size)
		multyp=np.array(multyp)
		cols1.append(deg[multyp.argmax()]+1)
		cols1_indices.append(deg[multyp.argmax()])
cols1 = np.array(cols1)
cols1_indices = np.array(cols1_indices)

# Second most present class
cols2=[]
for i in range(N):
	row = np.ma.array(pol[i,:])
	row[cols1_indices[i]] = np.ma.masked
	top = row.max()
	deg = np.where(row==top)[0]
	if  row[deg[0]] == 0:
		cols2.append(0)
	elif len(deg) == 1:
		cols2.append(deg[0]+1)
	else:
		multyp=[]
		for j in deg:
			multyp.append(np.where(matrix_first[i]==j+1)[0].size)
		multyp=np.array(multyp)
		cols2.append(deg[multyp.argmax()]+1)
cols2 = np.array(cols2)


# Plot the results as linear colored bars
binning = np.linspace(start,end,int(round((end - start)/c_info.resol))+1)
palette = ['#ffffff','#29ed02','#808000','#1EB001','#c71f0c','#00A2FF','#0000FF','#FF00FF','#B47638','#2F2F2F']

# First most present
fig1, ax1 = plt.subplots(figsize=(6, 0.4))
ax1.set_xlim((start,end))
for i in range(len(binning)-1):
	k = int(round((binning[i]-case[chr]['c_start'])/c_info.resol))
	ax1.axvspan(binning[i], binning[i+1], color=palette[cols1[k]])
ax1.xaxis.set_visible(False)
ax1.yaxis.set_visible(False)
ax1.spines['top'].set_visible(False)
ax1.spines['bottom'].set_visible(False)
ax1.spines['left'].set_visible(False)
ax1.spines['right'].set_visible(False)
# Uncomment the following line to save the plot
#plt.savefig("1st_most_present_chr%d[%.3f-%.3f].pdf"%(chr,start,end))

# Second most present
fig2, ax2 = plt.subplots(figsize=(6, 0.4))
ax2.set_xlim((start,end))
for i in range(len(binning)-1):
	k = int(round((binning[i]-case[chr]['c_start'])/c_info.resol))
	ax2.axvspan(binning[i], binning[i+1], color=palette[cols2[k]])
ax2.xaxis.set_visible(False)
ax2.yaxis.set_visible(False)
ax2.spines['top'].set_visible(False)
ax2.spines['bottom'].set_visible(False)
ax2.spines['left'].set_visible(False)
ax2.spines['right'].set_visible(False)
# Uncomment the following line to save the plot
#plt.savefig("2nd_most_present_chr%d[%.3f-%.3f].pdf"%(chr,start,end))
