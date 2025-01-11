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

def save_hist2d(fig, axe, df, column, title):
	#csv data is binned by time and energy, so we can use pivot_table to create a 2D histogram
	hist_data = df.pivot_table(index='energy(10*log10eV)', columns='time', values=column, aggfunc=np.sum, fill_value=0)
	#plot histogram
	im= axe.imshow(hist_data, aspect='auto', origin='lower', extent=[time_min, time_max, energy_min, energy_max], norm=LogNorm(vmin=count_min, vmax=count_max))
	fig.colorbar(im, label=colorbar_label)
	axe.set_title(title)
	axe.set_xlabel(xlabel)
	axe.set_ylabel(ylabel)
	#set scale for energy(10*log10eV -> linear scale)
	axe.yaxis.set_major_formatter(FuncFormatter(log_formatter))

if 'OUTPUT_FILE_NAME' in os.environ: #called from job script
	file_name = os.environ["OUTPUT_FILE_NAME"]
else: #default file name if not called from job script
	file_name="output.csv"
colorbar_label='number of actual particles'
xlabel='time [sec]'
ylabel='Energy [eV]'

#load data from file
df_all=pd.read_csv(file_name)

#save histograms for energy
column, description='par-count', 'Energy'
xlabel, ylabel = 'time [sec]', 'Energy [eV]'
#get nonzero data
df_non0 = df_all[(df_all[column] > 0)]
#get extent of data
time_min, time_max = df_non0['time'].agg(['min', 'max'])
energy_min, energy_max = df_non0['energy(10*log10eV)'].agg(['min', 'max'])
count_min, count_max = df_non0[column].agg(['min', 'max'])
print(f"extent of histogram of {description}: time: {time_min} to {time_max}, energy: {energy_min} to {energy_max}, count: {count_min} to {count_max}")
df_allpar = df_all[
	df_all['time'].between(time_min, time_max) &
	df_all['energy(10*log10eV)'].between(energy_min, energy_max)
]
#split data by species and save histograms
df_ele, df_ion = df_allpar[df_allpar['species'] == 1], df_allpar[df_allpar['species'] == 2]

#save histograms for energy
fig, axes = plt.subplots(1, 2, figsize=(14, 6))
save_hist2d(fig, axes[0], df_ele, column, description+"(electron)")
save_hist2d(fig, axes[1], df_ion, column, description+"(ion)")
plt.tight_layout()
plt.savefig(description+' of all particles.png', dpi=300)

#save histograms for velocity
columns = ['vx-count', 'vy-count', 'vz-count', 'vx-p', 'vy-p', 'vz-p', 'vx-n', 'vy-n', 'vz-n']
descriptions = ['x-comp. of vel.', 'y-comp. of vel.', 'z-comp. of vel.', 'x-comp. of vel. (+)', 'y-comp. of vel. (+)', 'z-comp. of vel. (+)', 'x-comp. of vel. (-)', 'y-comp. of vel. (-)', 'z-comp. of vel. (-)']
xlabels = ['time [sec]'] * 9
ylabels = ['x-comp. of vel. [m/s]', 'y-comp. of vel. [m/s]', 'z-comp. of vel. [m/s]', 'x-comp. of vel. (+)[m/s]', 'y-comp. of vel. (+)[m/s]', 'z-comp. of vel. (+)[m/s]', 'x-comp. of vel. (-)[m/s]', 'y-comp. of vel. (-)[m/s]', 'z-comp. of vel. (-)[m/s]']
#get nonzero data
df_non0 = df_all[(df_all[columns] > 0).any(axis=1)].replace(0, np.nan)
#get extent of data
time_min, time_max = df_non0['time'].agg(['min', 'max'])
energy_min, energy_max = df_non0['energy(10*log10eV)'].agg(['min', 'max'])
count_min, count_max = df_non0[columns].min().min(), df_non0[columns].max().max()
print(f"extent of histogram of {description}: time: {time_min} to {time_max}, energy: {energy_min} to {energy_max}, count: {count_min} to {count_max}")
df_allpar = df_all[
	df_all['time'].between(time_min, time_max) &
	df_all['energy(10*log10eV)'].between(energy_min, energy_max)
]

#save histograms
for column, description,xlabel,ylabel in zip(columns, descriptions,xlabels,ylabels):
	#split data by species
	df_ele, df_ion = df_allpar[df_allpar['species'] == 1], df_allpar[df_allpar['species'] == 2]
	fig, axes = plt.subplots(1, 2, figsize=(14, 6))
	save_hist2d(fig, axes[0], df_ele, column, description+"(electron)")
	save_hist2d(fig, axes[1], df_ion, column, description+"(ion)")
	plt.tight_layout()
	plt.savefig(description+' of all particles.png', dpi=300)
