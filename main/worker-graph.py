import matplotlib.pyplot as plt
from matplotlib import animation
import numpy as np
import pandas as pd
import io


#load data from file
file_name="dshield1.sh.2928515.out"

fig=plt.figure()
ims=[]

#read line one by one
with open(file_name,encoding='utf-8') as f:
	lines = f.readlines()
	start_str="113: ,   "
	#if line starts with start_str, then it is the line we want
	for i in range(len(lines)):
		if lines[i].startswith(start_str):
			#delete the first 7 characters
			lines[i]=lines[i][len(start_str):]
			#delete comma and space
			#lines[i]=lines[i].replace(",","")
			lines[i]=lines[i].replace(" ","")
			#read as csv
			df = pd.read_csv(io.StringIO(lines[i]), header=None)
			#to numpy array
			data = df.to_numpy()
			data*=1000000
			#print(data[0])
			#plot
			x=np.arange(0,len(data[0]),1)
			print(x,data[0])
			plt.title("worker phi")
			#set color red
			
			plt.xlabel("x")
			plt.ylabel("phi (EMSES unit x10^-6)")

			img=plt.plot(x, data[0],label='worker',color='black')
			
			ims.append(img)

ani=animation.ArtistAnimation(fig, ims, interval=100)
plt.show()
ani.save('worker.gif', writer='imagemagick', fps=10)
