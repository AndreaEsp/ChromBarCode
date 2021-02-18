# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as col

# Plot the 1st contribution matrix as a heatmap

data = np.genfromtxt('1st_contribution.txt') 
nclass = 9.

# Domain contribution is colored according to the epigetic class it belongs to
cmap = col.ListedColormap(['#ffffff','#00A2FF','#0000FF','#9437FF','#FF7F50','#B47638','#1EB001','#2F2F2F','#808000','#FF00FF','#EE230C'])

uniq_vals = set(np.reshape(data,np.size(data)))
bounds = [i for i in uniq_vals]
bounds.append(nclass+1)
bounds.sort()
bounds=list(np.array(bounds)-0.5)
norm = col.BoundaryNorm(bounds, cmap.N)

fig, ax = plt.subplots()
B = plt.imshow(data, cmap=cmap, alpha=1, interpolation='nearest', norm=norm)
ax.set_xticklabels([], minor=False, fontsize=10)
ax.set_yticklabels([], minor=False, fontsize=10)
ax.tick_params(axis='x',labelsize=20,labelbottom='off',labeltop='off')

plt.show()

# Uncomment the following line to save the plot
#plt.savefig("1st_contribution.pdf")