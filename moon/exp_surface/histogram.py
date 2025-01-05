import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import pandas as pd
import numpy as np
import os

def save_hist2d(df, x, y, bins, title, xlabel, ylabel, file_name):
	plt.clf()
	plt.hist2d(df[x], df[y], bins=bins, norm=LogNorm(vmin=max(vmin,1), vmax=vmax))
	plt.title(title)
	plt.yscale('log')
	plt.colorbar(label='number of super particles')
	plt.xlabel(xlabel)
	plt.ylabel(ylabel)
	plt.savefig(file_name, dpi=300)

file_name = os.environ["OUTPUT_FILE_NAME"]
#file_name="output,sy=16,sz=256,nt=10_11.csv"

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
