# -*- coding: utf-8 -*-

# The following script plots the binding domains in such a way that:
# - their signal is smoothed
# - each row will have a fixed number of domains
# - the domains will be ordered from the first to the last one

import numpy as np
import sys,os
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from scipy import signal

#Add path
analysis_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))+os.sep+'Modules'
sys.path.append(analysis_path)

import Polymer as poly
import Chroms_info as c_info

case = {'%s'%c_info.chrs[i]:{'c_start':c_info.c_start[i],'c_end':c_info.c_end[i],'bin':c_info.N_bin[i]} for i in range(len(c_info.chrs))}

ncol   = 30         # The number of binding domains per chromosome
n_rows = 3			# The number of rows to plot
n_dom_per_row = 10	# The number of domains per row

chrs = np.array([i for i in range(2,23,2)])

chrom = '20'        # Choose a SBS polymer to plot

start = case[chrom]['c_start']
end   = case[chrom]['c_end']
resol = c_info.resol

data_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))+os.sep+'Data'
pol = poly.polymer(data_path+os.sep+'sbs_polymers'+os.sep+'best_polymer_chr%s.txt'%(chrom),coords=[start, end])

# Choose the region of the selected chr to plot
z_start = start
z_end   = end

# Color code
palette = ['#000085','#610a00','#8c0e00','#c9c902','#570057','#808000','#ad00ad','#3838ff','#1EB001','#ff00ff',
		   '#404000','#545454','#c71f0c','#ff8fff','#B47638','#29ed02','#ff1900','#f59f49','#ffbfb8','#005c91',
		   '#84bd79','#5c5c9c','#fa7a6b','#00A2FF','#137300','#ffd9ff','#adade0','#85d2ff','#d9f1ff','#b39ab3']

if palette is None:
	cmap = plt.cm.brg
	colors = [cmap(i) for i in np.linspace(0,1,n_dom_per_row)]
	#colors.insert(0, 'gray')
else:
	colors = palette.copy()
	#colors.insert(0, 'gray')
		
bins = np.linspace(z_start+resol/2, z_end+resol/2, int( round((z_end - z_start)/resol) ), endpoint=False)
idx_start = int( round((z_start - start)/resol) ) 
idx_end = int( round((z_end - start)/resol) ) 

# Binding domains plot
fig, ax = plt.subplots(figsize=(7,2), nrows=n_rows, ncols=1, sharex=False, sharey=False)
fig.subplots_adjust(hspace=0)
col=0
for i in range(n_rows):
	maxes=[]
	for k,c in zip(range(n_dom_per_row), colors):
		y = pol.body[:,i*n_dom_per_row+k+1]
		y_smooth = signal.savgol_filter(y, 201, 3, mode='constant', cval=0)
		ax[i].fill_between(bins, y_smooth, 0, facecolor=colors[col], alpha=1, lw=.4, edgecolors='black')
		maxes.append(np.max(y_smooth))
		ax[i].spines['right'].set_visible(False)
		ax[i].spines['left'].set_visible(False)
		ax[i].spines['top'].set_visible(False)
		ax[i].spines['bottom'].set_linewidth(1.5)
		ax[i].xaxis.set_visible(False)
		ax[i].yaxis.set_visible(False)
		col+=1
	ax[i].axis([z_start, z_end, 0, np.max(maxes)])

# Uncomment the following line to save the plot
#plt.savefig('domains_chr%s_[%.3f-%.3f]Mb.pdf'%(chrom,z_start,z_end))

plt.show()
