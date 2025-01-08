import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import pandas as pd
import numpy as np
import os

def save_hist2d(df, title, file_name):
	plt.clf()
	#csv data is binned by time-step and energy, so we can use pivot_table to create a 2D histogram
	hist_data = df.pivot_table(index='energy(10*log10eV)', columns='time-step', values='sup-par-count', aggfunc=np.sum, fill_value=0)
	#plot histogram
	plt.imshow(hist_data, aspect='auto', origin='lower', extent=[1, 100, -100, 100], norm=LogNorm())
	#set labels
	plt.colorbar(label=colorbar_label)
	plt.title(title)
	plt.xlabel(xlabel)
	plt.ylabel(ylabel)
	plt.savefig(file_name, dpi=300)

#file_name = os.environ["OUTPUT_FILE_NAME"]
file_name="output,sy=16,sz=256,nt=10_9.csv"
real_par_num_per_sup_par=13020
colorbar_label='number of actual particles'
xlabel='Time-step'
ylabel='Energy (10 * log10 eV)'
output_file_name = file_name.replace('.csv','')+"Energy distribution"

#load data from file
df_all=pd.read_csv(file_name)
df_ele=df_all[df_all['species'] == 1]
df_ion=df_all[df_all['species'] == 2]

#save histograms
save_hist2d(df_all, "Energy distribution (all particles)", output_file_name + " (all particles)")
save_hist2d(df_ele, "Energy distribution (electron)", output_file_name + " (electron)")
save_hist2d(df_ion, "Energy distribution (ion)", output_file_name + " (ion)")
