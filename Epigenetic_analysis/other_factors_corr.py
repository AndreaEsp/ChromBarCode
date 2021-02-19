# -*- coding: utf-8 -*-

# The following script analyzes the Pearson correlation coefficient among the domains of each
# epigenetic class and the epigenetic data not-used for the hierarchical clustering.

import numpy as np
import sys,os
import scipy
import matplotlib.colors as col
import matplotlib.pyplot as plt 
from matplotlib.colors import ListedColormap
from scipy.cluster.hierarchy import dendrogram, linkage, cophenet, fcluster, set_link_color_palette
from scipy.spatial.distance import pdist
from scipy import stats
import seaborn as sns
import pandas as pd

def noise_centroids(X_rand, N_clust, N_step=100, method='ward', chrs=22, ncol=30):
    """Run the noise analysis for the cluster centroids"""
    sample=[]
    h=0
    for i in np.arange(0,N_step):
        R = X_rand[(chrs*ncol)*i:(chrs*ncol)*(i+1),:]    
        try:
            Z = linkage(R, method)
            c, coph_dists = cophenet(Z, pdist(R))
            if((scipy.cluster.hierarchy.is_valid_linkage(Z, warning=False, throw=False, name=None)) is False):
                print("ERROR::LINKAGE FAILED")
                sys.exit(0)
            partition = fcluster(Z,N_clust,criterion='maxclust')
            Centroids=[]
            for j in np.arange(1,N_clust+1):
                Centroids.append(np.mean(R[np.where(partition==j)],axis=0))
            sample.append(np.array(Centroids))
            h+=1
        except ValueError:
            print('STEP %d: Clustering failed, keep going to the next step'%(i+1))
    return np.reshape(np.array(sample),(N_clust*h,X_rand.shape[1]))

def cut_off(noisy_sample, rand_sample, perc_inf=2.5, perc_sup=97.5):
    """Retain only the significant correlation values"""
    clean_sample = np.zeros(noisy_sample.shape)
    for i in np.arange(0,rand_sample.shape[1]):
        inf = scipy.percentile(rand_sample[:,i],perc_inf)
        sup = scipy.percentile(rand_sample[:,i],perc_sup)
        (clean_sample[:,i])[np.where((noisy_sample[:,i]<=inf))] = (noisy_sample[:,i])[np.where((noisy_sample[:,i]<=inf))]
        (clean_sample[:,i])[np.where((noisy_sample[:,i]>=sup))] = (noisy_sample[:,i])[np.where((noisy_sample[:,i]>=sup))]
    return clean_sample

ncol   = 30     # The number of binding domains per chromosome
nclass = 9      # The number of epigenetic classes

data_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))+os.sep+'Data'
labels = np.genfromtxt(data_path+os.sep+'epigenetic_classes'+os.sep+'labels.txt')
chrs = np.array([i for i in range(2,23,2)])

class_order = [2,3,1,4,5,6,7,8,9]

# Select the list of columns/marks to show
marks = ["H3K4me2","H3K9ac","H3K27ac","H3K79me2","H4K20me1","H2AFZ"]
#marks = ["DNase","RepG1","RepS1","RepS2","RepS3","RepS4","RepG2","RNAseq_total"]
#marks = ["CTCF","RAD21","SMC3","POLR2A"]
#marks = ["POLR2AphosphoS2","POLR2AphosphoS5","EZH2","STAT1","ZNF687","RUNX3","CEBPZ","ZNF217","TCF12","ZNF592","USF1","HCFC1","MXI1","HSF1","ZEB1","MYB","EP300","KLF5","CEBPB","TAF1","ZNF384","CUX1","ELF1","RB1","ZBTB33","TBP","NFYB","LARP7","WRNIP1","JUNB","EBF1","BHLHE40","STAT3","SIN3A","MEF2C","NKRF","PBX3","MTA3","ESRRA","NR2F1","ZSCAN29","NR2C1","NFXL1","IKZF1","IRF3","TBL1XR1","PAX5","HDAC2","KDM1A","YY1","CHD4","RXRA","DPF2","SMAD1","MAZ","NR2C2","SIX5","CHD1","ZFP36","GATAD2B","BCL3","SRF","RBBP5","BCLAF1","BRCA1","ZNF24","SMARCA5","ZNF207","NFATC1","REST","UBTF","RELB","ZNF622","NFIC","CBX5","FOXK2","CBFB","CBX3","HDAC6","MEF2B","PAX8","MTA2","ATF7","TBX21","E2F8","ETS1","IRF5","ELK1","NBN","IRF4","MAX","RFX5","CHD2","BATF","STAT5A","JUND","BCL11A","ARID3A","ZNF143","BACH1","NRF1","SMAD5","TARDBP","SPI1","RCOR1","TRIM22","NFATC3","RAD51","GABPA","BMI1","PKNOX1","YBX1","ASH2L","EED","SKIL","E2F4","ZBTB40","SUZ12","ATF2","MAFK","NFYA","ETV6","EGR1","USF2","E4F1","MLLT1","IKZF2","TCF7","CREM","ZBED1","ARNT","HDGF","MEF2A","ZZZ3"]

# Import data
x = pd.read_csv(data_path+os.sep+'epigenetic_correlations'+os.sep+'significant_correlations.csv', usecols=marks)[marks]
X = x.to_numpy()

corr_all = []
corr_avg = []
for k in class_order:
	corr_all.append(X[np.where(labels==k)[0]])
	corr_avg.append(np.mean(X[np.where(labels==k)[0]],0))
corr_all = np.vstack(corr_all)
corr_avg = np.vstack(corr_avg)

# Noise analysis of the centroid matrix
x_rand = pd.read_csv(data_path+os.sep+'epigenetic_correlations'+os.sep+'control_correlations.csv.gz', usecols=marks)[marks]
X_rand = x_rand.to_numpy()
rnd_com = noise_centroids(X_rand, N_clust=nclass, N_step=100, method='ward', chrs=22, ncol=30)
corr_avg = cut_off(corr_avg, rnd_com, perc_inf=10, perc_sup=90)

# Custom colorbar
cmap = col.ListedColormap(['#d675d6','#8e4e8e','#472747','#387e73','#33c199','#87ff9d'])
cmap.set_bad('black')
norm = col.BoundaryNorm([-.3,-.2,-.1,0,.1,.2,.3], cmap.N)

# Custom colorbar (centroids)
cmap_centroid = col.ListedColormap(['#d675d6','#472747','#387e73','#33c199','#87ff9d'])
cmap_centroid.set_bad('black')
norm_centroid = col.BoundaryNorm([-.2,-.1,0,.1,.2,.3], cmap.N)

fig1, ax1 = plt.subplots()
mask = np.zeros(corr_all.shape)
mask[np.where(corr_all==0)]=1
percsup , percinf = np.percentile(corr_all,98) , np.percentile(corr_all,2)
sns.heatmap(corr_all, cmap=cmap, norm=norm, vmax=.3, vmin=-.3, square=False, linewidths=0, annot=False, linecolor='white',
    edgecolor='black', cbar=False, xticklabels=marks, yticklabels=False, mask=mask, ax=ax1)
ax1.tick_params(axis='both',labelsize=8,labelbottom=True,labeltop=False,labelleft=False,length=0)
plt.xticks(rotation=90, horizontalalignment='center',fontweight='bold')
plt.tight_layout()
# Uncomment the following line to save the plot
# plt.savefig("other_factors_correlations.pdf")

fig2, ax2 = plt.subplots()
mask = np.zeros_like(corr_avg)
mask[np.where(corr_avg==0)]=1
sns.heatmap(corr_avg, cmap=cmap_centroid, norm=norm_centroid, vmax=.3, vmin=-.2, square=True, linewidths=1, annot=False,linecolor='white',
    edgecolor='black', cbar=False, xticklabels=marks, yticklabels=False, mask=mask, ax=ax2)#, cbar_kws={"shrink":.5, "ticks":[]})
ax2.tick_params(axis='both',labelsize=8,labelbottom=True,labeltop=False,labelleft=False,length=0)
plt.xticks(rotation=90, horizontalalignment='center',fontweight='bold')
plt.tight_layout()

# Uncomment the following line to save the plot
# plt.savefig("other_factors_correlations_centroids.pdf")

plt.show()
