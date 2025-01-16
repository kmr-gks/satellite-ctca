import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.colors import LogNorm
from matplotlib.ticker import FuncFormatter
import pandas as pd
import numpy as np
import os

# custom formatter for y-axis
def log_formatter(value, tick_number):
    return f"$10^{{{value/10}}}$"

def save_hist2d(axe, df, column, title):
	#csv data is binned by time and energy, so we can use pivot_table to create a 2D histogram
	hist_data = df.pivot_table(index='energy(10*log10eV)', columns='time', values=column, aggfunc=np.sum, fill_value=0)
	#plot histogram
	im= axe.imshow(hist_data, aspect='auto', origin='lower', extent=[time_min, time_max, energy_min, energy_max], norm=LogNorm(vmin=count_min, vmax=count_max))
	#axe.set_title(title)
	axe.text(0.5, 0.05, title, transform=axe.transAxes, fontsize=12, va='bottom', ha='center', bbox=dict(facecolor='white', alpha=0.8, edgecolor='none'))
	#set scale for energy(10*log10eV -> linear scale)
	axe.yaxis.set_major_formatter(FuncFormatter(log_formatter))
	return im

if 'OUTPUT_FILE_NAME' in os.environ: #called from job script
	file_name = os.environ["OUTPUT_FILE_NAME"]
else: #default file name if not called from job script
	file_name="output.csv"

colorbar_label='number of actual particles'
dpi=600

#load data from file
df_all=pd.read_csv(file_name)

#save histograms for energy
column, description='par-count', 'Energy'
xlabel, ylabel = 'time [sec]', 'Energy [eV]'

#get nonzero data
df_non0 = df_all[(df_all[column] > 0)]
#get extent of data
if df_non0.empty:
	df_non0 = df_all
	count_min, count_max = 1,10
else:
	count_min, count_max = df_non0[column].agg(['min', 'max'])
time_min, time_max = df_all['time'].agg(['min', 'max'])
energy_min, energy_max = df_non0['energy(10*log10eV)'].agg(['min', 'max'])

df_allpar = df_all[
	df_all['energy(10*log10eV)'].between(energy_min, energy_max)
]
#split data by species and save histograms
df_ele, df_ion = df_allpar[df_allpar['species'] == 1], df_allpar[df_allpar['species'] == 2]

#save histograms for energy
fig, axes = plt.subplots(1, 2, figsize=(6, 2.5), sharex="col", sharey="row", constrained_layout=True)
norm = mpl.colors.LogNorm(vmin=count_min, vmax=count_max)
cmap = plt.get_cmap()

save_hist2d(axes[0], df_ion, column, description+"(ion)")
save_hist2d(axes[1], df_ele, column, description+"(electron)")
fig.supxlabel(xlabel)
fig.supylabel(ylabel)
cbar=fig.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=cmap), ax=axes.ravel().tolist(), aspect=25, location="right", pad=0.1)
cbar.set_label(colorbar_label)
plt.savefig(description, dpi=dpi)

#save histograms for velocity
columns = ['vx-p', 'vy-p', 'vz-p', 'vx-n', 'vy-n', 'vz-n']
descriptions = ['x-comp. of vel. (+)', 'y-comp. of vel. (+)', 'z-comp. of vel. (+)', 'x-comp. of vel. (-)', 'y-comp. of vel. (-)', 'z-comp. of vel. (-)']
file_name=['x-comp. of vel.', 'y-comp. of vel.', 'z-comp. of vel.']
xlabel = 'time [sec]'
ylabel = "velocity [m/s]"
#get nonzero data
df_non0 = df_all[(df_all[columns] > 0).any(axis=1)].replace(0, np.nan)
if df_non0.empty:
	df_non0 = df_all
	count_min, count_max = 1,10
else:
	count_min, count_max = df_non0[columns].min().min(), df_non0[columns].max().max()
#get extent of data
energy_min, energy_max = df_non0['energy(10*log10eV)'].agg(['min', 'max'])

df_allpar = df_all[
	df_all['energy(10*log10eV)'].between(energy_min, energy_max)
]

df_ele, df_ion = df_allpar[df_allpar['species'] == 1], df_allpar[df_allpar['species'] == 2]

#save histograms for velocity
for i in range(3):
	fig, axes = plt.subplots(2, 2, figsize=(6, 5), sharex="col", sharey="row", constrained_layout=True)
	norm = mpl.colors.LogNorm(vmin=count_min, vmax=count_max)
	cmap = plt.get_cmap()
	save_hist2d(axes[0,0], df_ion, columns[i], descriptions[i]+"(ion)")
	save_hist2d(axes[0,1], df_ele, columns[i], descriptions[i]+"(electron)")
	save_hist2d(axes[1,0], df_ion, columns[i+3], descriptions[i+3]+"(ion)")
	save_hist2d(axes[1,1], df_ele, columns[i+3], descriptions[i+3]+"(electron)")
	fig.supxlabel(xlabel)
	fig.supylabel(ylabel)
	cbar=fig.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=cmap), ax=axes.ravel().tolist(), aspect=25, location="right", pad=0.1)
	cbar.set_label(colorbar_label)
	plt.savefig(file_name[i], dpi=dpi)
