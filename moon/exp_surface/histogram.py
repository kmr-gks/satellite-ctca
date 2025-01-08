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

file_name = os.environ["OUTPUT_FILE_NAME"]
#file_name="output,sy=16,sz=256,nt=10_12.csv"
real_par_num_per_sup_par=13020

#load data from file
df=pd.read_csv(file_name)

#csv data is binned by time-step and energy, so we can use pivot_table to create a 2D histogram
hist_data = df.pivot_table(index='energy(10*log10eV)', columns='time-step', values='sup-par-count', aggfunc='sum', fill_value=0)

plt.imshow(hist_data, aspect='auto', origin='lower', extent=[1, 100, -100, 100], norm=LogNorm())

plt.colorbar(label='Particle Count')
plt.xlabel('Time-step')
plt.ylabel('Energy (10 * log10 eV)')
plt.savefig("out.png", dpi=300)
