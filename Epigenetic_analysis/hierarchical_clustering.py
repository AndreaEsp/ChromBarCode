# -*- coding: utf-8 -*-

# The following script performs a hierarchical clustering of the matrix of correlations
# among binding domains and chromatin marks.

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

F_name = ["H3K4me3","H3K4me1","H3K36me3","H3K9me3","H3K27me3"]
data_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))+os.sep+'Data'

# Import data
x = pd.read_csv(data_path+os.sep+'epigenetic_correlations'+os.sep+'significant_correlations.csv', usecols=F_name)[F_name]
X = x.to_numpy()
X_norm = stats.zscore(X)

clust_method='ward'        # Clustering method

# Clustering procedure
Z = linkage(X_norm,clust_method)
L = scipy.cluster.hierarchy.leaves_list(Z)
c, coph_dists = cophenet(Z, pdist(X_norm))

# Dendrogram
fig,ax = plt.subplots(figsize=(5,5))
ax.set_frame_on(False),ax.set_yticks([])
#scipy.cluster.hierarchy.set_link_color_palette(['r','g','b','magenta'])
scipy.cluster.hierarchy.set_link_color_palette(['#29ed02','#808000','#1EB001','#c71f0c','#00A2FF','#0000FF','#FF00FF','#B47638','#2F2F2F'])

D = dendrogram(
    Z,
    truncate_mode='none',           # show only the last p merged clusters ('none' or 'lastp')
    p=40,                           # show only the last p merged clusters
    show_leaf_counts=False,         # otherwise numbers in brackets are counts
    #show_contracted=True,          # to get a distribution impression in truncated branches
    leaf_rotation=90.,              # rotates the x axis labels
    leaf_font_size=7.,              # font size for the x axis labels
    orientation='top',
    #count_sort='descending',
    no_labels=True,
    color_threshold=7.2,
    above_threshold_color='black',
    #annotate_above=5,
    #max_d=0.75,
)
plt.show()
# Uncomment the following line to save the plot
plt.savefig("dendrogram.pdf")

# Custom colorbar
cmap = col.ListedColormap(['#d675d6','#8e4e8e','#472747','#387e73','#33c199','#87ff9d'])
cmap.set_bad('black')
norm = col.BoundaryNorm([-.3,-.2,-.1,0,.1,.2,.3], cmap.N)

# Sorting and plotting the matrix of observation
XX = np.copy(X)
X_sorted = XX[L,:]
fig, ax = plt.subplots()
mask = np.zeros(X.shape)
mask[np.where(X_sorted==0)]=1
percsup , percinf = np.percentile(X_sorted,98) , np.percentile(X_sorted,2)
sns.heatmap(X_sorted, cmap=cmap, norm=norm, vmax=.3, vmin=-.3, square=False, linewidths=0, annot=False,linecolor='white',
    edgecolor='black', cbar=False, xticklabels=F_name, yticklabels=False, mask=mask, ax=ax)
ax.tick_params(axis='both',labelsize=8,labelbottom=True,labeltop=False,labelleft=False,length=0)
plt.xticks(rotation=45, horizontalalignment='center',fontweight='bold')
plt.tight_layout()
# Uncomment the following line to save the plot
# plt.savefig("correlation_matrix.pdf")

# Compute centroid matrix
N_clusters = 9
final_clustering = fcluster(Z,N_clusters,criterion='maxclust')
centroid=[]
for i in np.arange(1,N_clusters+1):
    centroid.append(np.mean(X[np.where(final_clustering==i)],axis=0))
centroid=np.array(centroid)

# Noise analysis of the centroid matrix
x_rand = pd.read_csv(data_path+os.sep+'epigenetic_correlations'+os.sep+'control_correlations.csv.gz', usecols=F_name)[F_name]
X_rand = x_rand.to_numpy()
rnd_com = noise_centroids(X_rand, N_clust=N_clusters, N_step=100, method=clust_method, chrs=11, ncol=30)
final_centroid = cut_off(centroid, rnd_com, perc_inf=1, perc_sup=99)
# Uncomment the following line to save the data
# np.savetxt("%d-centroid_matrix.txt"%(N_clusters),final_centroid,fmt='%.8f')

# Sorting rows
order = (1,2,0,3,4,5,6,7,8)
final_centroid = final_centroid[order,:]

# Custom colorbar for centroids
cmap_centroid = col.ListedColormap(['#d675d6','#472747','#387e73','#33c199','#87ff9d'])
cmap_centroid.set_bad('black')
norm_centroid = col.BoundaryNorm([-.2,-.1,0,.1,.2,.3], cmap.N)

fig, ax = plt.subplots()
mask = np.zeros_like(final_centroid)
mask[np.where(final_centroid==0)]=1
sns.heatmap(final_centroid, cmap=cmap_centroid, norm=norm_centroid, vmax=.3, vmin=-.2, square=True, linewidths=1, annot=False,linecolor='white',
    edgecolor='black', cbar=False, xticklabels=F_name, yticklabels=False, mask=mask, ax=ax)
ax.tick_params(axis='both',labelsize=8,labelbottom=True,labeltop=False,labelleft=False,length=0)
plt.xticks(rotation=60, horizontalalignment='center',fontweight='bold')
plt.tight_layout()
# Uncomment the following line to save the plot
plt.savefig("%d-centroids_matrix.pdf"%(N_clusters))

plt.show()
