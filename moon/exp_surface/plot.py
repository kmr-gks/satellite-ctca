import emout
import matplotlib.pyplot as plt
import os
import sys

job = sys.argv[1]
dirname=job+'_plot/'

data = emout.Emout('output/' + job)
print("read:", job)
os.makedirs(job + '_plot',exist_ok=True)

print(data.phisp[:].val_si.shape)

if sys.argv[2]=='x':
	x=int(sys.argv[3])
	data.phisp[-1, :, :, x].val_si.plot()
	data.j1yz[-1, :, :, x].val_si.plot(mode='stream',title=f'j1yz[A/m^2], x={x}')
	plt.savefig(dirname+f'j1yz_x{x}.png')
	plt.clf()

	data.phisp[-1, :, :, x].val_si.plot()
	data.j2yz[-1, :, :, x].val_si.plot(mode='stream',title=f'j2yz[A/m^2], x={x}')
	plt.savefig(dirname+f'j2yz_x{x}.png')
	plt.clf()

	data.phisp[-1, :, :, x].val_si.plot()
	data.jyz[-1, :, :, x].val_si.plot(mode='stream',title=f'jyz[A/m^2], x={x}')
	plt.savefig(dirname+f'jyz_x{x}.png')
	plt.clf()

	data.phisp[-1, :, :, x].val_si.plot(savefilename=dirname+f'phisp_x{x}.png',title=f'Potential[V] at x={x} [m]')

	data.phisp[:, :, :, x].val_si.gifplot(action='save', filename=dirname+f'phisp_x{x}.gif',title=f'Potential[V] at x={x} [m]')
	plt.clf()

	data.nd1p[-1, :, :, x].val_si.plot(savefilename=dirname+f'nd1p_x{x}.png',title=f'Electron density[/cc] at x={x} [m]')

	data.nd1p[:, :, :, x].val_si.gifplot(action='save', filename=dirname+f'nd1p_x{x}.gif',title=f'Electron density[/cc] at x={x} [m]')
	plt.clf()

	data.nd2p[-1, :, :, x].val_si.plot(savefilename=dirname+f'nd2p_x{x}.png',title=f'Ion density[/cc] at x={x} [m]')

	data.nd2p[:, :, :, x].val_si.gifplot(action='save', filename=dirname+f'nd2p_x{x}.gif',title=f'Ion density[/cc] at x={x} [m]')
	plt.clf()
elif sys.argv[2]=='y':
	y=int(sys.argv[3])
	data.phisp[-1, :, y, :].val_si.plot()
	data.j1xz[-1, :, y, :].val_si.plot(mode='stream',title=f'j1xz[A/m^2], y={y}')
	plt.savefig(dirname+f'j1xz_y{y}.png')
	plt.clf()

	data.phisp[-1, :, y, :].val_si.plot()
	data.j2xz[-1, :, y, :].val_si.plot(mode='stream',title=f'j2xz[A/m^2], y={y}')
	plt.savefig(dirname+f'j2xz_y{y}.png')
	plt.clf()

	data.phisp[-1, :, y, :].val_si.plot()
	data.jxz[-1, :, y, :].val_si.plot(mode='stream',title=f'jxz[A/m^2], y={y}')
	plt.savefig(dirname+f'jxz_y{y}.png')
	plt.clf()

	data.phisp[-1, :, y, :].val_si.plot(savefilename=dirname+f'phisp_y{y}.png',title=f'Potential[V] at y={y} [m]')

	data.phisp[:, :, y, :].val_si.gifplot(action='save', filename=dirname+f'phisp_y{y}.gif',title=f'Potential[V] at y={y} [m]')
	plt.clf()

	data.nd1p[-1, :, y, :].val_si.plot(savefilename=dirname+f'nd1p_y{y}.png',title=f'Electron density[/cc] at y={y} [m]')

	data.nd1p[:, :, y, :].val_si.gifplot(action='save', filename=dirname+f'nd1p_y{y}.gif',title=f'Electron density[/cc] at y={y} [m]')
	plt.clf()

	data.nd2p[-1, :, y, :].val_si.plot(savefilename=dirname+f'nd2p_y{y}.png',title=f'Ion density[/cc] at y={y} [m]')

	data.nd2p[:, :, y, :].val_si.gifplot(action='save', filename=dirname+f'nd2p_y{y}.gif',title=f'Ion density[/cc] at y={y} [m]')
	plt.clf()
elif sys.argv[2]=='z':
	z=int(sys.argv[3])
	data.phisp[-1, z, :, :].val_si.plot()
	data.j1xy[-1, z, :, :].val_si.plot(mode='stream',title=f'j1xy[A/m^2], z={z}')
	plt.savefig(dirname+f'j1xy_z{z}.png')
	plt.clf()

	data.phisp[-1, z, :, :].val_si.plot()
	data.j2xy[-1, z, :, :].val_si.plot(mode='stream',title=f'j2xy[A/m^2], z={z}')
	plt.savefig(dirname+f'j2xy_z{z}.png')
	plt.clf()

	data.phisp[-1, z, :, :].val_si.plot()
	data.jxy[-1, z, :, :].val_si.plot(mode='stream',title=f'jxy[A/m^2], z={z}')
	plt.savefig(dirname+f'jxy_z{z}.png')
	plt.clf()

	data.phisp[-1, z, :, :].val_si.plot(savefilename=dirname+f'phisp_z{z}.png',title=f'Potential[V] at z={z} [m]')

	data.phisp[:, z, :, :].val_si.gifplot(action='save', filename=dirname+f'phisp_z{z}.gif',title=f'Potential[V] at z={z} [m]')
	plt.clf()

	data.nd1p[-1, z, :, :].val_si.plot(savefilename=dirname+f'nd1p_z{z}.png',title=f'Electron density[/cc] at z={z} [m]')

	data.nd1p[:, z, :, :].val_si.gifplot(action='save', filename=dirname+f'nd1p_z{z}.gif',title=f'Electron density[/cc] at z={z} [m]')
	plt.clf()

	data.nd2p[-1, z, :, :].val_si.plot(savefilename=dirname+f'nd2p_z{z}.png',title=f'Ion density[/cc] at z={z} [m]')

	data.nd2p[:, z, :, :].val_si.gifplot(action='save', filename=dirname+f'nd2p_z{z}.gif',title=f'Ion density[/cc] at z={z} [m]')
	plt.clf()
