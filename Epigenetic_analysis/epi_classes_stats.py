# -*- coding: utf-8 -*-

# Statistical analyses of the 9 epigenetic classes of binding domains

import numpy as np
import sys,os
from scipy import stats 
import matplotlib.pyplot as plt
import matplotlib.colors as col

import seaborn as sns

ncol   = 30     # The number of binding domains per chromosome
nclass = 9      # The number of epigenetic classes

data_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))+os.sep+'Data'

chrs = np.array([i for i in range(2,23,2)])

labels = np.genfromtxt(data_path+os.sep+'epigenetic_classes'+os.sep+'labels.txt')
labels_per_chr = {'%d'%i:labels[j*ncol:(j+1)*ncol] for i,j in zip(range(2,23,2),range(11))}
pol_per_chr = {'%d'%i:np.genfromtxt(data_path+os.sep+'sbs_polymers'+os.sep+'best_polymer_chr%d.txt'%i) for i in chrs}
class_order = [2,3,1,4,5,6,7,8,9]

# Color palette
palette = ['#808000','#1EB001','#29ed02','#c71f0c','#00A2FF','#0000FF','#FF00FF','#B47638','#2F2F2F']


# ---- Genomic Coverage ----
tot_mass = 0
for i in pol_per_chr:
	tot_mass += np.sum(pol_per_chr[i])

coverage = []
for k in class_order:
	count = 0
	for i in pol_per_chr:
		idx = np.where(labels_per_chr[i]==k)[0] + 1
		count += np.sum(pol_per_chr[i][:,idx])
	coverage.append(count)
coverage = np.array(coverage)/tot_mass

# Barplot
fig1, ax1 = plt.subplots(figsize=(7,4))
ax1.bar(range(nclass),coverage*100,align='center',width=.8,linewidth=2.3,edgecolor='black',color=palette)
ax1.set_ylim((0,15))
ax1.set_yticks([0,5,10,15])
ax1.set_xticks([])
ax1.tick_params(axis='both',which='both',bottom=False,top=False,right=False,left=True,direction='in',width=1.5,labelsize=25)
# Uncomment the following line to save the plot
plt.savefig("coverage.pdf")


# ---- Distribution of the different classes over chromosomes ----
class_dist = np.zeros((nclass,len(chrs)))
for i in range(nclass):
	for j in range(len(chrs)):
		class_dist[i,j] = len(np.where( labels_per_chr['%d'%chrs[j]]==class_order[i] )[0])

# Barplots
fig2, ax2 = plt.subplots(figsize=(5,5), nrows=nclass, ncols=1, sharex=False)
fig2.subplots_adjust(hspace=0.17)
for ax in range(nclass):
	ax2[ax].bar(chrs, class_dist[ax], align='center', width=1.6, linewidth=1, edgecolor='black', color=palette[ax])
	ax2[ax].axhline(y=np.mean(class_dist[ax]), ls='--', lw=1, color='black')
	ax2[ax].set_ylim((0,class_dist[ax].max()))
	ax2[ax].set_xlim((1,23))
	ax2[ax].set_xticks([])
	ax2[ax].set_yticks([])
	ax2[ax].tick_params(axis='both',which='both', labelsize=6, pad=1, length=2.3, width=1.3)
ax2[-1].set_xticks(chrs)
ax2[-1].set_xticklabels([])
# Uncomment the following line to save the plot
plt.savefig("classes_distribution.pdf")

# ---- Correlation of the genomic location of the different classes over chromosomes ----
correl_mat = np.zeros((nclass,nclass))
for m in range(nclass):
	for k in range(m,nclass):
		print("Pair %d-%d"%(class_order[m],class_order[k]))
		out = []
		for i in range(len(chrs)):
			idx1 = np.where( labels_per_chr['%d'%chrs[i]] == class_order[m])[0] + 1
			tmp1 = pol_per_chr['%d'%chrs[i]][:,idx1]
			idx2 = np.where( labels_per_chr['%d'%chrs[i]] == class_order[k])[0] + 1
			tmp2 = pol_per_chr['%d'%chrs[i]][:,idx2]
			if tmp1.shape[1] > 0 and tmp2.shape[1] > 0:
				for t in np.transpose(tmp1):
					for s in np.transpose(tmp2):
						out.append(stats.pearsonr(t,s)[0])
		correl_mat[m,k] = np.mean( np.array(out) )
		correl_mat[k,m] = np.mean( np.array(out) )


fig3,ax3 = plt.subplots(figsize=(8,6))
cmap = col.ListedColormap(['#000074','#9797FF','#FFAFAF','#8E0000'])
cmap.set_bad('white')
bounds = [-0.3,-0.15,0,0.15,0.3]
norm = col.BoundaryNorm(bounds,cmap.N)
sns.heatmap(correl_mat,square=True,cmap=cmap,linewidths=1,edgecolor='black',linecolor='w',xticklabels=False,
	yticklabels=False,norm=norm,ax=ax3,cbar_kws={'ticks':[-0.3,-0.15,0,0.15,0.3]})
cax=plt.gcf().axes[-1]
cax.tick_params(labelsize=25)

# Uncomment the following line to save the plot
plt.savefig("correlation_matrix.pdf")

plt.show()
