import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import pandas as pd

#load data from file
df=pd.read_csv("output.csv")

#plot histogram
plt.hist2d(df['step'], df['energy'], bins=(6,10), norm=LogNorm())
plt.colorbar(label='number of particles (pbuf)')
plt.xlabel("step/satellite x[pbuf%x]")
plt.ylabel("energy [pbuf%x^2]")
#plt.show()
plt.savefig('histogram.high.png', dpi=300)
