import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import pandas as pd
import numpy as np
import os

#file_name = os.environ["OUTPUT_FILE_NAME"]
file_name="output,sy=16,sz=256,nt=10_11.csv"

#load data from file
df_original=pd.read_csv(file_name)

df=df_original[df_original['species']==1]
log_bins=[9, np.logspace(np.log10(df['energy'].min()),np.log10(df['energy'].max()), num=50)]
#plot histogram of electron energy
plt.clf()
plt.hist2d(df['step'], df['energy'], bins=log_bins, norm=LogNorm())
title="Energy distribution (electron)"
plt.title(title)
plt.yscale('log')
plt.colorbar(label='number of super particles')
plt.xlabel("x[m]")
plt.ylabel("energy density [eV/m^3]")
plt.savefig(file_name.replace('.csv','')+title, dpi=300)

df=df_original[df_original['species']==2]
log_bins=[9, np.logspace(np.log10(df['energy'].min()),np.log10(df['energy'].max()), num=50)]
#plot histogram of ion energy
plt.clf()
plt.hist2d(df['step'], df['energy'], bins=log_bins, norm=LogNorm())
title="Energy distribution (ion)"
plt.title(title)
plt.yscale('log')
plt.colorbar(label='number of super particles')
plt.xlabel("x[m]")
plt.ylabel("energy density [eV/m^3]")
plt.savefig(file_name.replace('.csv','')+title, dpi=300)
