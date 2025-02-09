import emout
import matplotlib.pyplot as plt
import os
import sys

job = sys.argv[1]
dirname = job+'_plot/'
grid_length = 0.5

data = emout.Emout('output/' + job)
print("read:", job)
os.makedirs(job + '_plot',exist_ok=True)

print(data.phisp[:].val_si.shape)

if sys.argv[2]=='x':
	x=int(sys.argv[3])
	data.phisp[-1, :, :, x].val_si.plot()
	data.j1yz[-1, :, :, x].val_si.plot(mode='stream',title=f'Electron current density [A/$m^2$] at x={x*grid_length} [m]')
	plt.savefig(dirname+f'j1yz_x{x*grid_length}.png')
	plt.clf()

	data.phisp[-1, :, :, x].val_si.plot()
	data.j2yz[-1, :, :, x].val_si.plot(mode='stream',title=f'Ion current density [A/$m^2$] at x={x*grid_length} [m]')
	plt.savefig(dirname+f'j2yz_x{x*grid_length}.png')
	plt.clf()

	data.phisp[-1, :, :, x].val_si.plot()
	data.jyz[-1, :, :, x].val_si.plot(mode='stream',title=f'Current density [A/$m^2$] at x={x*grid_length} [m]')
	plt.savefig(dirname+f'jyz_x{x*grid_length}.png')
	plt.clf()

	data.phisp[-1, :, :, x].val_si.plot(savefilename=dirname+f'phisp_x{x*grid_length}.png',title=f'Potential [V] at x={x*grid_length} [m]')
	data.phisp[:, :, :, x].val_si.gifplot(action='save', filename=dirname+f'phisp_x{x*grid_length}.gif',title=f'Potential [V] at x={x*grid_length} [m]')
	plt.clf()

	data.nd1p[-1, :, :, x].val_si.plot(savefilename=dirname+f'nd1p_x{x*grid_length}.png',title=f'Electron density [/cc] at x={x*grid_length} [m]')
	data.nd1p[:, :, :, x].val_si.gifplot(action='save', filename=dirname+f'nd1p_x{x*grid_length}.gif',title=f'Electron density [/cc] at x={x*grid_length} [m]')
	plt.clf()

	data.nd2p[-1, :, :, x].val_si.plot(savefilename=dirname+f'nd2p_x{x*grid_length}.png',title=f'Ion density [/cc] at x={x*grid_length} [m]')
	data.nd2p[:, :, :, x].val_si.gifplot(action='save', filename=dirname+f'nd2p_x{x*grid_length}.gif',title=f'Ion density [/cc] at x={x*grid_length} [m]')
	plt.clf()

	data.bx[-1, :, :, x].val_si.plot(savefilename=dirname+f'bx_x{x*grid_length}.png',title=f'X comp. of Magnetic field [A/m] at x={x*grid_length} [m]')
	data.by[-1, :, :, x].val_si.plot(savefilename=dirname+f'by_x{x*grid_length}.png',title=f'Y comp. of Magnetic field [A/m] at x={x*grid_length} [m]')
	data.bz[-1, :, :, x].val_si.plot(savefilename=dirname+f'bz_x{x*grid_length}.png',title=f'Z comp. of Magnetic field [A/m] at x={x*grid_length} [m]')

	data.ex[-1, :, :, x].val_si.plot(savefilename=dirname+f'ex_x{x*grid_length}.png',title=f'X comp. of Electric field [V/m] at x={x*grid_length} [m]')
	data.ey[-1, :, :, x].val_si.plot(savefilename=dirname+f'ey_x{x*grid_length}.png',title=f'Y comp. of Electric field [V/m] at x={x*grid_length} [m]')
	data.ez[-1, :, :, x].val_si.plot(savefilename=dirname+f'ez_x{x*grid_length}.png',title=f'Z comp. of Electric field [V/m] at x={x*grid_length} [m]')
elif sys.argv[2]=='y':
	y=int(sys.argv[3])
	data.phisp[-1, :, y, :].val_si.plot()
	data.j1xz[-1, :, y, :].val_si.plot(mode='stream',title=f'Electron current density [A/$m^2$] at y={y*grid_length} [m]')
	plt.savefig(dirname+f'j1xz_y{y*grid_length}.png')
	plt.clf()

	data.phisp[-1, :, y, :].val_si.plot()
	data.j2xz[-1, :, y, :].val_si.plot(mode='stream',title=f'Ion current density [A/$m^2$] at y={y*grid_length} [m]')
	plt.savefig(dirname+f'j2xz_y{y*grid_length}.png')
	plt.clf()

	data.phisp[-1, :, y, :].val_si.plot()
	data.jxz[-1, :, y, :].val_si.plot(mode='stream',title=f'Current density [A/$m^2$] at y={y*grid_length} [m]')
	plt.savefig(dirname+f'jxz_y{y*grid_length}.png')
	plt.clf()

	data.phisp[-1, :, y, :].val_si.plot(savefilename=dirname+f'phisp_y{y*grid_length}.png',title=f'Potential [V] at y={y*grid_length} [m]')
	data.phisp[:, :, y, :].val_si.gifplot(action='save', filename=dirname+f'phisp_y{y*grid_length}.gif',title=f'Potential [V] at y={y*grid_length} [m]')
	plt.clf()

	data.nd1p[-1, :, y, :].val_si.plot(savefilename=dirname+f'nd1p_y{y*grid_length}.png',title=f'Electron density [/cc] at y={y*grid_length} [m]')
	data.nd1p[:, :, y, :].val_si.gifplot(action='save', filename=dirname+f'nd1p_y{y*grid_length}.gif',title=f'Electron density [/cc] at y={y*grid_length} [m]')
	plt.clf()

	data.nd2p[-1, :, y, :].val_si.plot(savefilename=dirname+f'nd2p_y{y*grid_length}.png',title=f'Ion density [/cc] at y={y*grid_length} [m]')
	data.nd2p[:, :, y, :].val_si.gifplot(action='save', filename=dirname+f'nd2p_y{y*grid_length}.gif',title=f'Ion density [/cc] at y={y*grid_length} [m]')
	plt.clf()

	data.bx[-1, :, y, :].val_si.plot(savefilename=dirname+f'bx_y{y*grid_length}.png',title=f'X comp. of Magnetic field [A/m] at y={y*grid_length} [m]')
	data.by[-1, :, y, :].val_si.plot(savefilename=dirname+f'by_y{y*grid_length}.png',title=f'Y comp. of Magnetic field [A/m] at y={y*grid_length} [m]')
	data.bz[-1, :, y, :].val_si.plot(savefilename=dirname+f'bz_y{y*grid_length}.png',title=f'Z comp. of Magnetic field [A/m] at y={y*grid_length} [m]')

	data.ex[-1, :, y, :].val_si.plot(savefilename=dirname+f'ex_y{y*grid_length}.png',title=f'X comp. of Electric field [V/m] at y={y*grid_length} [m]')
	data.ey[-1, :, y, :].val_si.plot(savefilename=dirname+f'ey_y{y*grid_length}.png',title=f'Y comp. of Electric field [V/m] at y={y*grid_length} [m]')
	data.ez[-1, :, y, :].val_si.plot(savefilename=dirname+f'ez_y{y*grid_length}.png',title=f'Z comp. of Electric field [V/m] at y={y*grid_length} [m]')
elif sys.argv[2]=='z':
	z=int(sys.argv[3])
	data.phisp[-1, z, :, :].val_si.plot()
	data.j1xy[-1, z, :, :].val_si.plot(mode='stream',title=f'Electron current density [A/$m^2$] at z={z*grid_length} [m]')
	plt.savefig(dirname+f'j1xy_z{z*grid_length}.png')
	plt.clf()

	data.phisp[-1, z, :, :].val_si.plot()
	data.j2xy[-1, z, :, :].val_si.plot(mode='stream',title=f'Ion current density [A/$m^2$] at z={z*grid_length} [m]')
	plt.savefig(dirname+f'j2xy_z{z*grid_length}.png')
	plt.clf()

	data.phisp[-1, z, :, :].val_si.plot()
	data.jxy[-1, z, :, :].val_si.plot(mode='stream',title=f'Current density [A/$m^2$] at z={z*grid_length} [m]')
	plt.savefig(dirname+f'jxy_z{z*grid_length}.png')
	plt.clf()

	data.phisp[-1, z, :, :].val_si.plot(savefilename=dirname+f'phisp_z{z*grid_length}.png',title=f'Potential [V] at z={z*grid_length} [m]')
	data.phisp[:, z, :, :].val_si.gifplot(action='save', filename=dirname+f'phisp_z{z*grid_length}.gif',title=f'Potential [V] at z={z*grid_length} [m]')
	plt.clf()

	data.nd1p[-1, z, :, :].val_si.plot(savefilename=dirname+f'nd1p_z{z*grid_length}.png',title=f'Electron density [/cc] at z={z*grid_length} [m]')
	data.nd1p[:, z, :, :].val_si.gifplot(action='save', filename=dirname+f'nd1p_z{z*grid_length}.gif',title=f'Electron density [/cc] at z={z*grid_length} [m]')
	plt.clf()

	data.nd2p[-1, z, :, :].val_si.plot(savefilename=dirname+f'nd2p_z{z*grid_length}.png',title=f'Ion density [/cc] at z={z*grid_length} [m]')
	data.nd2p[:, z, :, :].val_si.gifplot(action='save', filename=dirname+f'nd2p_z{z*grid_length}.gif',title=f'Ion density [/cc] at z={z*grid_length} [m]')
	plt.clf()

	data.bx[-1, z, :, :].val_si.plot(savefilename=dirname+f'bx_z{z*grid_length}.png',title=f'X comp. of Magnetic field [A/m] at z={z*grid_length} [m]')
	data.by[-1, z, :, :].val_si.plot(savefilename=dirname+f'by_z{z*grid_length}.png',title=f'Y comp. of Magnetic field [A/m] at z={z*grid_length} [m]')
	data.bz[-1, z, :, :].val_si.plot(savefilename=dirname+f'bz_z{z*grid_length}.png',title=f'Z comp. of Magnetic field [A/m] at z={z*grid_length} [m]')

	data.ex[-1, z, :, :].val_si.plot(savefilename=dirname+f'ex_z{z*grid_length}.png',title=f'X comp. of Electric field [V/m] at z={z*grid_length} [m]')
	data.ey[-1, z, :, :].val_si.plot(savefilename=dirname+f'ey_z{z*grid_length}.png',title=f'Y comp. of Electric field [V/m] at z={z*grid_length} [m]')
	data.ez[-1, z, :, :].val_si.plot(savefilename=dirname+f'ez_z{z*grid_length}.png',title=f'Z comp. of Electric field [V/m] at z={z*grid_length} [m]')
