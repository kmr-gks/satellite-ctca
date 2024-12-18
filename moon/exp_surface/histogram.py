import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import pandas as pd
import os

file_name = os.environ["OUTPUT_FILE_NAME"]

#load data from file
df=pd.read_csv(file_name)

#plot histogram
plt.hist2d(df['step'], df['energy'], bins=(10,10), norm=LogNorm())
plt.colorbar(label='number of particles (pbuf)')
plt.xlabel("step/satellite x[pbuf%x]")
plt.ylabel("energy [pbuf%x^2]")
#plt.show()

file_name=file_name.replace('.csv','')
plt.savefig(file_name, dpi=300)
