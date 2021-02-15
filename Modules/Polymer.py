# -*- coding: utf-8 -*-
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from scipy import signal
import os,sys

class polymer:

	"""Class to read and manage a SBS polymer. Coordinates must be set in Mb"""
	
	def __init__(self, path, coords=[0,0]):
		
		self.path = path
		self.coords = coords 
		self.body = np.genfromtxt(self.path)
		self.N = self.body.shape[0]
		self.n = self.body.shape[1]-1
		self.r = int( np.max(np.sum(self.body,1)) )
		self.start = self.coords[0]
		self.end = self.coords[1]
		self.nmb = len(np.where(np.sum(self.body,1)==0)[0])

		if self.start==0 and self.end==0:
			self.resol = 0
		else:
			self.resol = round( (self.end - self.start)/self.N, 3)

	def info(self):
		"""Print out SBS polymer's features"""
		print('\nNumber of bins: %d \nNumber of domains: %d \nModel resolution: %d \nNumber of beads: %d (excluding possible non-mappable regions)'%(self.N, self.n, self.r, np.sum(self.body)))

	def save(self, path):
		"""Save the SBS polymer to file"""
		np.savetxt(path, self.body, fmt='%d', delimiter=' ')

	def copy(self):
		"""Copy the SBS polymer"""
		return polymer(self.path, [self.start, self.end])

	def sort(self):
		"""Sort the binding domains by (ascending) order of their barycenter"""
		if self.resol == 0:
			print('Coordinates must be set to use this function')
		else:
			self.binning = np.linspace(self.start,self.end,self.N)
			self.masses = []
			for i in range(0,self.n):
				self.masses.append(np.sum(self.binning*self.body[:,i+1])/np.sum(self.body[:,i+1]))
			self.body = np.c_[self.body[:,0], self.body[:,1:][:,np.argsort(self.masses)]]
			self.com = (np.array(self.masses)[np.argsort(self.masses)])
			return self.body

	def draw(self, z_start, z_end, ylim, palette=None, save=False):
		"""Plot the binding domains of a specific region

		Parameters
		----------
		z_start: Start coordinate of the region (must be in Mb).

		z_end:   End coordinate of the region (must be in Mb).

		ylim:    Top y-axis limit for each domain.

		palette: The list of colors (one per domain) to use. If None (default), set an automatic palette.

		save:    If True, save the plot. (default is False)"""
		
		if palette is None:
			cmap = plt.cm.brg
			self.colors = [cmap(i) for i in np.linspace(0,1,self.n)]
			self.colors.insert(0, 'gray')
		else:
			self.colors = palette.copy()
			self.colors.insert(0, 'gray')
		
		self.bins = np.linspace(z_start, z_end, int( round((z_end - z_start)/self.resol) ), endpoint=False)
		self.idx_start = int( round((z_start - self.start)/self.resol) ) 
		self.idx_end = int( round((z_end - self.start)/self.resol) ) 

		fig, ax = plt.subplots(nrows=self.n+1, ncols=1, sharex=True, sharey=True)
		fig.subplots_adjust(hspace=.1)		
		for i,c in zip(range(self.n+1), self.colors):
			ax[i].bar(self.bins, self.body[self.idx_start:self.idx_end,i], width=self.resol, align='edge', color=c)
			ax[i].axis([z_start, z_end, 0, ylim])
			ax[i].spines['right'].set_visible(False)
			ax[i].spines['left'].set_visible(False)
			ax[i].spines['top'].set_visible(False)
			ax[i].xaxis.set_visible(False)
			ax[i].yaxis.set_visible(False)
		ax[-1].xaxis.set_visible(True)
		ax[-1].xaxis.set_major_locator(ticker.AutoLocator())
		ax[-1].tick_params(axis='x',labelsize=12)

		if save==True: plt.savefig('domains_[%.3f-%.3f]Mb.pdf'%(z_start,z_end), transparent=False, bbox_inches='tight')

	def draw_smooth(self, z_start, z_end, ylim, w_length, polyorder, hspace=.1, lw=.5, x_lw=1, lim='equal', palette=None, save=False, save_as='png'):
		"""Plot a smooth version of the binding domains of a specific region

		Parameters
		----------
		z_start: Start coordinate of the region (must be in Mb).

		z_end:   End coordinate of the region (must be in Mb).

		ylim:    Top y-axis limit for each domain.

		w_length: The length of the filter window. Must be a positive odd integer.

		polyorder: The order of the polynomial used to fit the samples. Must be less than window_length.

		hspace: The separation among the subplots.

		lw: The thickness of the curve border.

		x_lw: The thickness of the x axes.

		lim: If 'equal' (default), each plot will have the same y limit (ylim). If 'each', the y limit of each plot will be (0,max).

		palette: ThThe list of colors (one per domain) to use. If None (default), set an automatic palette.

		save:    If True, save the plot. (defalut is False)"""
		
		if palette is None:
			cmap = plt.cm.brg
			self.colors = [cmap(i) for i in np.linspace(0,1,self.n)]
			self.colors.insert(0, 'gray')
		else:
			self.colors = palette.copy()
			self.colors.insert(0, 'gray')
		
		self.bins = np.linspace(z_start+self.resol/2, z_end+self.resol/2, int( round((z_end - z_start)/self.resol) ), endpoint=False)
		self.idx_start = int( round((z_start - self.start)/self.resol) ) 
		self.idx_end = int( round((z_end - self.start)/self.resol) ) 

		fig, ax = plt.subplots(nrows=self.n+1, ncols=1, sharex=True, sharey=False)
		fig.subplots_adjust(hspace=hspace)		
		for i,c in zip(range(self.n+1), self.colors):
			y = self.body[self.idx_start:self.idx_end,i]
			y_smooth = signal.savgol_filter(y, w_length, polyorder, mode='constant', cval=0)
			ax[i].fill_between(self.bins, y_smooth, 0, facecolor=c, alpha=1, lw=lw, edgecolors='black')
			if lim == 'equal': ax[i].axis([z_start, z_end, 0, ylim])
			if lim == 'each' : ax[i].axis([z_start, z_end, 0, np.max(y_smooth)])
			ax[i].spines['right'].set_visible(False)
			ax[i].spines['left'].set_visible(False)
			ax[i].spines['top'].set_visible(False)
			ax[i].spines['bottom'].set_linewidth(x_lw)
			ax[i].xaxis.set_visible(False)
			ax[i].yaxis.set_visible(False)
		ax[-1].xaxis.set_visible(True)
		ax[-1].xaxis.set_major_locator(ticker.AutoLocator())
		ax[-1].tick_params(axis='x',labelsize=12)

		if save==True: plt.savefig('domains_smooth_w%d_p%d_[%.3f-%.3f]Mb.%s'%(w_length,polyorder,z_start,z_end,save_as), transparent=False, bbox_inches='tight')


	def draw_subset(self, z_start, z_end, ylim, margin, save=False):
		"""Plot the binding domains of a specific region. This function plots only the domains
		whose barycenter falls in the specified region ± the given margin.

		Parameters
		----------
		z_start: Start coordinate of the region (must be in Mb).

		z_end:   End coordinate of the region (must be in Mb).

		ylim:    Top y-axis limit for each domain.

		margin:  Left and right margin for the considered region (must be in Mb). 

		save:    If True, save the plot. (defalut is False)"""

		self.domais_subset = np.where( (self.com >= z_start-margin) & (self.com <= z_end+margin) )[0] + 1
		self.domais_subset = np.insert(self.domais_subset,0,0)

		cmap = plt.cm.brg
		self.colors = [cmap(i) for i in np.linspace(0,1,len(self.domais_subset)-1)]
		self.colors.insert(0, 'gray')
		
		self.bins = np.linspace(z_start, z_end, int( round((z_end - z_start)/self.resol) ), endpoint=False)
		self.idx_start = int( round((z_start - self.start)/self.resol) ) 
		self.idx_end = int( round((z_end - self.start)/self.resol) )

		if len(self.domais_subset)==1: print("No barycenter between %.3f-%.3f Mb\n"%(z_start-margin,z_end+margin))

		else:
			fig, ax = plt.subplots(nrows=len(self.domais_subset), ncols=1, sharex=True, sharey=True)
			for i,c in zip(range(len(self.domais_subset)), self.colors):
				ax[i].bar(self.bins, self.body[self.idx_start:self.idx_end,self.domais_subset[i]], width=self.resol, align='edge', color=c)
				ax[i].axis([z_start, z_end, 0, ylim])
				ax[i].spines['right'].set_visible(False)
				ax[i].spines['left'].set_visible(False)
				ax[i].spines['top'].set_visible(False)
				ax[i].xaxis.set_visible(False)
				ax[i].yaxis.set_visible(False)
			ax[-1].xaxis.set_visible(True)
			ax[-1].xaxis.set_major_locator(ticker.AutoLocator())
			ax[-1].tick_params(axis='x',labelsize=12)

			if save==True: plt.savefig('subset_of_domains_[%.3f-%.3f]Mb.png'%(z_start,z_end), transparent=False, bbox_inches='tight', dpi=100)

	def bd_size(self):
		"""Compute the size of the binding domains. 
		For each domain the size is defined as the fraction of its binding sites multiplied by that polymer length"""
		self.length = self.end - self.start
		return (np.sum(self.body,0) / np.sum(self.body)) * self.length

	def bd_interaction_length(self):
		"""Compute the interaction length of the binding domains
		For each domain the interaction length is defined as two times the standard deviation of the barycenter of that domain"""
		std = []
		for i in range(self.n+1):
			avg = np.average(self.binning, weights=self.body[:,i] )
			std.append( np.sqrt(np.average((self.binning-avg)**2, weights=self.body[:,i])) )
		return 2*np.array(std)