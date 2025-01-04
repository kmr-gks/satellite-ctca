import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import pandas as pd
import numpy as np
import os

def save_hist2d(df, x, y, bins, title, xlabel, ylabel, file_name):
	plt.clf()
	plt.hist2d(df[x], df[y], bins=bins, norm=LogNorm())
	plt.title(title)
	plt.yscale('log')
	plt.colorbar(label='number of super particles')
	plt.xlabel(xlabel)
	plt.ylabel(ylabel)
	plt.savefig(file_name, dpi=300)

file_name = os.environ["OUTPUT_FILE_NAME"]

#load data from file
df=pd.read_csv(file_name)

log_bins=[len(df['step'].value_counts()), np.logspace(np.log10(df['energy'].min()),np.log10(df['energy'].max()), num=50)]

#plot histogram of electron energy
save_hist2d(df[df['species']==1], 'step', 'energy', log_bins, "Energy distribution (electron)", "x[m]", "energy density [eV/m^3]", file_name.replace('.csv','')+"Energy distribution (electron)")

#plot histogram of ion energy
save_hist2d(df[df['species']==2], 'step', 'energy', log_bins, "Energy distribution (ion)", "x[m]", "energy density [eV/m^3]", file_name.replace('.csv','')+"Energy distribution (ion)")

#plot histogram of all particles energy
save_hist2d(df, 'step', 'energy', log_bins, "Energy distribution (all particles)", "x[m]", "energy density [eV/m^3]", file_name.replace('.csv','')+"Energy distribution (all particles)")
