# -*- coding: utf-8 -*-

import numpy as np
import sys,os
import matplotlib
import matplotlib.pyplot as plt
import scipy
from scipy import signal

#Add path
analysis_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))+os.sep+'Modules'
sys.path.append(analysis_path)

import Polymer as poly
import Tools as tools
import Chroms_info as c_info

case = {'%s'%c_info.chrs[i]:{'c_start':c_info.c_start[i],'c_end':c_info.c_end[i],'bin':c_info.N_bin[i]} for i in range(len(c_info.chrs))}

overlaps     = []
rnd_overlaps = []
iteration    = 100

for chr in ['%d'%i for i in c_info.chrs]:
	print('Chr%s'%chr)
	data_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))+os.sep+'Data'
	pol = poly.polymer(data_path+os.sep+'sbs_polymers'+os.sep+'best_polymer_chr%s.txt'%(chr),coords=[case[chr]['c_start'],case[chr]['c_end']])

	# Compute overlap among the domains
	overlaps.append(tools.intra_overlaps(pol))

	# Compute overlap among bootstrapped domains (random model)
	rnd_model = pol.copy()
	for k in range(iteration):
		print("Iteration %d"%k)
		for i in np.arange(1,1+pol.n):
			rnd_model.body[:,i] = tools.build_random_model(pol.body[:,i])
		rnd_overlaps.append(tools.intra_overlaps(rnd_model))

overlaps = np.array(overlaps).flatten()
rnd_overlaps = np.array(rnd_overlaps).flatten()

print('Overlaps --> mean: %f\t std: %f'%(np.mean(overlaps), np.std(overlaps)))
print('Random overlaps --> mean: %f\t std: %f'%(np.mean(rnd_overlaps), np.std(rnd_overlaps)))

# Histogram of results
fig1, ax1 = plt.subplots()
bins = np.linspace(0,0.8,30)
ax1.hist(rnd_overlaps, bins=bins, weights=np.ones_like(rnd_overlaps)/float(len(rnd_overlaps)), histtype='stepfilled', edgecolor='black', color='#1E78B4', linewidth=2)
ax1.hist(overlaps, bins=bins, weights=np.ones_like(overlaps)/float(len(overlaps)), histtype='stepfilled', edgecolor='black', color='#FF7F7E', linewidth=2, alpha=0.9)
plt.yscale('linear')
plt.xscale('linear')
ax1.tick_params(labelsize=20)
ax1.set_xlim(0,0.7)
ax1.set_ylim(0,0.2)
ax1.set_xticks([0,0.25,0.50,0.75])
ax1.set_yticks([0.0,0.10,0.20])
plt.setp(ax1.spines.values(), linewidth=1)
plt.show()

# Uncomment the following line to save the plot
# plt.savefig("bd_overlaps.pdf")
