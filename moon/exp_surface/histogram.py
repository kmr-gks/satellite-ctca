import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import pandas as pd
import numpy as np
import os

def save_hist2d(df, x, y, bins, title, xlabel, ylabel, file_name):
	plt.clf()
	hist, xedges, yedges, im = plt.hist2d(df[x], df[y], bins=bins, norm=LogNorm(vmin=max(vmin,1), vmax=vmax))
	#colorbar for number of real particles
	hist_real_par=hist*real_par_num_per_sup_par
	plt.clf()
	plt.pcolormesh(xedges, yedges, hist_real_par.T, norm=LogNorm(vmin=max(vmin * real_par_num_per_sup_par, 1), vmax=vmax * real_par_num_per_sup_par))
	plt.title(title)
	plt.yscale('log')
	plt.colorbar(label='number of actual particles')
	plt.xlabel(xlabel)
	plt.ylabel(ylabel)
	#get average
	x_bin_indices = np.digitize(df[x], xedges) - 1
	valid_indices = (x_bin_indices >= 0) & (x_bin_indices < len(xedges) - 1)
	x_bin_indices = x_bin_indices[valid_indices]
	y_valid = df[y][valid_indices]
	y_means = np.zeros(len(xedges) - 1)
	np.add.at(y_means, x_bin_indices, y_valid)
	bin_counts = np.bincount(x_bin_indices, minlength=len(xedges) - 1)
	y_means = np.divide(y_means, bin_counts, where=bin_counts > 0)
	x_bin_centers = (xedges[:-1] + xedges[1:]) / 2
	plt.plot(x_bin_centers, y_means, color='red', label='average energy')
	plt.legend()
	plt.savefig(file_name, dpi=300)

#file_name = os.environ["OUTPUT_FILE_NAME"]
file_name="output,sy=16,sz=256,nt=10_12.csv"
real_par_num_per_sup_par=13020

#load data from file
df=pd.read_csv(file_name)

log_bins=[len(df['time'].value_counts()), np.logspace(np.log10(df['energy'].min()),np.log10(df['energy'].max()), num=50)]
hist,xedges,yedges =np.histogram2d(df['time'], df['energy'], bins=log_bins)
vmin, vmax = hist.min(), hist.max()

#plot histogram of electron energy
save_hist2d(df[df['species']==1], 'time', 'energy', log_bins, "Energy distribution (electron)", "time[sec]", "energy [eV]", file_name.replace('.csv','')+"Energy distribution (electron)")

#plot histogram of ion energy
save_hist2d(df[df['species']==2], 'time', 'energy', log_bins, "Energy distribution (ion)", "time[sec]", "energy [eV]", file_name.replace('.csv','')+"Energy distribution (ion)")

#plot histogram of all particles energy
save_hist2d(df, 'time', 'energy', log_bins, "Energy distribution (all particles)", "time[sec]", "energy [eV]", file_name.replace('.csv','')+"Energy distribution (all particles)")
