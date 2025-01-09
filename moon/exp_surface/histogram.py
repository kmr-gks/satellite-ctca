import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib.ticker import FuncFormatter
import pandas as pd
import numpy as np
import os
import re

# custom formatter for y-axis
def log_formatter(value, tick_number):
    return f"$10^{{{int(value/10)}}}$"

def save_hist2d(df, title, file_name):
	plt.clf()
	#csv data is binned by time-step and energy, so we can use pivot_table to create a 2D histogram
	hist_data = df.pivot_table(index='energy(10*log10eV)', columns='time-step', values='par-count', aggfunc=np.sum, fill_value=0)
	#plot histogram
	plt.imshow(hist_data, aspect='auto', origin='lower', extent=[time_min, time_max, energy_min, energy_max], norm=LogNorm(vmin=count_min, vmax=count_max))
	plt.colorbar(label=colorbar_label)
	plt.title(title)
	plt.xlabel(xlabel)
	plt.ylabel(ylabel)
	#set log10 scale for energy
	plt.gca().yaxis.set_major_formatter(FuncFormatter(log_formatter))
	plt.savefig(file_name, dpi=300)

file_name = os.environ["OUTPUT_FILE_NAME"]
#file_name="output,sy=16,sz=256,nt=10_9.csv"
colorbar_label='number of actual particles'
xlabel='Time-step'
ylabel='Energy [eV]'
output_file_name = file_name.replace('.csv','')+"Energy distribution"

#load data from file
df_all=pd.read_csv(file_name)
#get nonzero data
df_non0=df_all[df_all['par-count'] > 0]
#get extent of data
time_min = df_non0['time-step'].min()
time_max = df_non0['time-step'].max()
energy_min = df_non0['energy(10*log10eV)'].min()
energy_max = df_non0['energy(10*log10eV)'].max()
count_min = df_non0['par-count'].min()
count_max = df_non0['par-count'].max()
#filter data to only include the range of interest
df_all=df_all[(df_all['time-step'] >= time_min) & (df_all['time-step'] <= time_max) & (df_all['energy(10*log10eV)'] >= energy_min) & (df_all['energy(10*log10eV)'] <= energy_max)]
#split data by species
df_ele=df_all[df_all['species'] == 1]
df_ion=df_all[df_all['species'] == 2]

#save histograms
save_hist2d(df_all, "Energy distribution (all particles)", output_file_name + " (all particles)")
save_hist2d(df_ele, "Energy distribution (electron)", output_file_name + " (electron)")
save_hist2d(df_ion, "Energy distribution (ion)", output_file_name + " (ion)")
