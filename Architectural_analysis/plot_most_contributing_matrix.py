# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as col

# Plot the 1st contribution matrix as a heatmap

data = np.genfromtxt('1st_contribution.txt') 
nclass = 9.

# Domain contribution is colored according to the epigetic class it belongs to
cmap = col.ListedColormap(['#ffffff','#ffffff','#29ed02','#808000','#1EB001','#c71f0c','#00A2FF','#0000FF','#FF00FF','#B47638','#2F2F2F'])

uniq_vals = set([-1.0, 0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0])
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