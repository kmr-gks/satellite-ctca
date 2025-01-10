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

def save_hist2d(df, column, title):
	plt.clf()
	#csv data is binned by time and energy, so we can use pivot_table to create a 2D histogram
	hist_data = df.pivot_table(index='energy(10*log10eV)', columns='time', values=column, aggfunc=np.sum, fill_value=0)
	#plot histogram
	plt.imshow(hist_data, aspect='auto', origin='lower', extent=[time_min, time_max, energy_min, energy_max], norm=LogNorm(vmin=count_min, vmax=count_max))
	plt.colorbar(label=colorbar_label)
	plt.title(title)
	plt.xlabel(xlabel)
	plt.ylabel(ylabel)
	#set scale for energy(10*log10eV -> linear scale)
	plt.gca().yaxis.set_major_formatter(FuncFormatter(log_formatter))
	plt.savefig(title+'.png', dpi=300)

if 'OUTPUT_FILE_NAME' in os.environ: #called from job script
	file_name = os.environ["OUTPUT_FILE_NAME"]
else: #default file name if not called from job script
	file_name="output.csv"
colorbar_label='number of actual particles'
xlabel='time [sec]'
ylabel='Energy [eV]'

#load data from file
df_all=pd.read_csv(file_name)
#get extent of data
time_min, time_max = df_all['time'].agg(['min', 'max'])
#get nonzero data
df_non0 = df_all[(df_all[['par-count', 'vx-count', 'vy-count', 'vz-count']] > 0).any(axis=1)].replace(0, np.nan)
#get extent of data
energy_min, energy_max = df_non0['energy(10*log10eV)'].agg(['min', 'max'])
count_min, count_max = df_non0[['par-count', 'vx-count', 'vy-count', 'vz-count']].min().min(), df_non0[['par-count', 'vx-count', 'vy-count', 'vz-count']].max().max()
#filter data to only include the range of interest
df_allpar = df_all[
	df_all['time'].between(time_min, time_max) &
	df_all['energy(10*log10eV)'].between(energy_min, energy_max)
]

columns = ['par-count', 'vx-count', 'vy-count', 'vz-count']
descriptions = ['Energy', 'x-comp. of vel.', 'y-comp. of vel.', 'z-comp. of vel.']
xlabels = ['time [sec]', 'time [sec]', 'time [sec]', 'time [sec]']
ylabels = ['Energy [eV]', 'x-comp. of vel. [m/s]', 'y-comp. of vel. [m/s]', 'z-comp. of vel. [m/s]']

#save histograms
for column, description,xlabel,ylabel in zip(columns, descriptions,xlabels,ylabels):
	#split data by species
	df_ele, df_ion = df_allpar[df_allpar['species'] == 1], df_allpar[df_allpar['species'] == 2]
	save_hist2d(df_allpar, column, description+"(all particles)")
	save_hist2d(df_ele, column, description+"(electron)")
	save_hist2d(df_ion, column, description+"(ion)")
