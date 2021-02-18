# -*- coding: utf-8 -*-

import numpy as np
import sys,os
import matplotlib
import matplotlib.pyplot as plt
import scipy
from scipy import signal,optimize

#Add path
analysis_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))+os.sep+'Modules'
sys.path.append(analysis_path)

import Polymer as poly
import Tools as tools
import Chroms_info as c_info

case = {'%s'%c_info.chrs[i]:{'c_start':c_info.c_start[i],'c_end':c_info.c_end[i],'bin':c_info.N_bin[i]} for i in range(len(c_info.chrs))}

mean_length = np.mean(c_info.c_end - c_info.c_start)
sizes   = []
int_len = []
iteration = 100

for chr in ['%d'%i for i in c_info.chrs]:
	print('Chr%s'%chr)
	data_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))+os.sep+'Data'
	pol = poly.polymer(data_path+os.sep+'sbs_polymers'+os.sep+'best_polymer_chr%s.txt'%(chr),coords=[case[chr]['c_start'],case[chr]['c_end']])

	# Sort the binding domains by (ascending) order of their barycenter
	pol.sort()

	# Compute domain size and range of interaction
	sizes.append( pol.bd_size()[1:] )
	int_len.append( pol.bd_interaction_length()[1:] )

sizes = np.array(sizes).flatten()
int_len = np.array(int_len).flatten()

# Do the same on bootstrapped domains (random model)
print("\nRandom controls")
rnd_int_len_all = []
for chr,start,end in zip(['%d'%i for i in c_info.chrs],c_info.c_start,c_info.c_end):
	print('Chr%s rand'%chr)
	pol = poly.polymer(data_path+os.sep+'sbs_polymers'+os.sep+'best_polymer_chr%s.txt'%(chr),coords=[case[chr]['c_start'],case[chr]['c_end']])
	rnd_int_len = []
	for k in range(iteration):
		rnd_model = pol.copy()
		rnd_model.sort()
		for i in np.arange(1,1+pol.n):
			rnd_model.body[:,i] = tools.build_random_model(pol.body[:,i])
		rnd_int_len.append( rnd_model.bd_interaction_length()[1:] )
	rnd_int_len = np.array(rnd_int_len).flatten()
	rnd_int_len = rnd_int_len/(end-start)
	rnd_int_len_all.append(rnd_int_len)
rnd_int_len_all = (np.array(rnd_int_len_all).flatten())*mean_length

print('Size --> mean: %f\t std: %f'%(np.mean(sizes), np.std(sizes)))
print('Int length --> mean: %f\t std: %f'%(np.mean(int_len), np.std(int_len)))
print('Random Int length --> mean: %f\t std: %f'%(np.mean(rnd_int_len_all), np.std(rnd_int_len_all)))

# Histogram of results
fig, ax = plt.subplots()
bins = np.linspace(5,175,15)
ax.hist(int_len, bins=bins, weights=np.ones_like(int_len)/float(len(int_len)), histtype='stepfilled', edgecolor='black', color='#FF7F7E',linewidth=2, zorder=1, alpha=1)
ax.hist(rnd_int_len_all, bins=bins, weights=np.ones_like(rnd_int_len_all)/float(len(rnd_int_len_all)), histtype='stepfilled', edgecolor='black', color='#1E78B4',linewidth=2, alpha=0.7, zorder=10)

# Power-law fit of the interaction length
p_law = lambda x,a,b,c: (a/(x-b)) + c
y,x = np.histogram(int_len, bins=bins, weights=np.ones_like(int_len)/float(len(int_len)))
x_smooth = np.zeros(len(x)-1)
for i in range(len(x)-1):
	x_smooth[i] = (x[i+1] + x[i])/2
x_fit = x_smooth[5:]
y_fit = y[5:]
popt1, pcov1 = scipy.optimize.curve_fit(f=p_law, xdata=x_fit, ydata=y_fit, p0=[20,0.1,-0.1])

xval = np.linspace(60, 170,100)
fit  = p_law(xval,*popt1)
ax.plot(xval, fit, color='darkred', linestyle='--',linewidth=5, zorder=20)

ax.set_xlim(0,180)
ax.set_ylim(0,0.2)
ax.set_xticks([0,50,100,150])
ax.tick_params(labelsize=20)
plt.setp(ax.spines.values(), linewidth=1)
plt.yscale('linear')
plt.xscale('linear')
plt.show()

# Uncomment the following line to save the plot
#plt.savefig("bd_interaction_length.pdf")
