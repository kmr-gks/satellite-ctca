import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import pandas as pd
import numpy as np
import os

file_name = os.environ["OUTPUT_FILE_NAME"]

#load data from file
df=pd.read_csv(file_name)
log_bins=[9, np.logspace(np.log10(df['energy'].min()),np.log10(df['energy'].max()), num=50)]

#plot histogram
plt.hist2d(df['step'], df['energy'], bins=log_bins, norm=LogNorm())
title="Energy distribution"
plt.title(title)
plt.yscale('log')
plt.colorbar(label='number of particles')
plt.xlabel("x")
plt.ylabel("energy")

file_name=file_name.replace('.csv','')
plt.savefig(file_name+title, dpi=300)
