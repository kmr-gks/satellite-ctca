import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm


#load data from file
file_name="job.sh.3218323.out"

with open(file_name,encoding='utf-8') as f:
	lines = f.readlines()
	start_str="129:  step"
	data_str="129:  energy  "
	step=0
	particle_x=[]
	particle_energy=[]
	#if line starts with start_str, then it is the line we want
	for i in range(len(lines)):
		if lines[i].startswith(start_str):
			step+=1
		if lines[i].startswith(data_str):
			lines[i]=lines[i][len(data_str):-1]
			particle_x.append(step)
			particle_energy.append(float(lines[i]))
	plt.hist2d(particle_x, particle_energy, norm=LogNorm())
	plt.colorbar(label='number of particles (pbuf)')
	plt.xlabel("step/satellite x[pbuf%x]")
	plt.ylabel("energy [pbuf%x^2]")
	#plt.show()
	plt.savefig('histogram.png')
