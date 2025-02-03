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
x, y, z = int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4])

data.phisp[-1, :, :, x].val_si.plot()
data.j1yz[-1, :, :, x].val_si.plot(mode='stream')
plt.savefig(dirname+f'j1yz_x{x}.png')
plt.clf()

data.phisp[-1, :, y, :].val_si.plot()
data.j1xz[-1, :, y, :].val_si.plot(mode='stream')
plt.savefig(dirname+f'j1xz_y{y}.png')
plt.clf()

data.phisp[-1, z, :, :].val_si.plot()
data.j1xy[-1, z, :, :].val_si.plot(mode='stream')
plt.savefig(dirname+f'j1xy_z{z}.png')
plt.clf()

data.phisp[-1, :, :, x].val_si.plot()
data.j2yz[-1, :, :, x].val_si.plot(mode='stream')
plt.savefig(dirname+f'j2yz_x{x}.png')
plt.clf()

data.phisp[-1, :, y, :].val_si.plot()
data.j2xz[-1, :, y, :].val_si.plot(mode='stream')
plt.savefig(dirname+f'j2xz_y{y}.png')
plt.clf()

data.phisp[-1, z, :, :].val_si.plot()
data.j2xy[-1, z, :, :].val_si.plot(mode='stream')
plt.savefig(dirname+f'j2xy_z{z}.png')
plt.clf()

data.phisp[-1, :, :, x].val_si.plot()
data.jyz[-1, :, :, x].val_si.plot(mode='stream')
plt.savefig(dirname+f'jyz_x{x}.png')
plt.clf()

data.phisp[-1, :, y, :].val_si.plot()
data.jxz[-1, :, y, :].val_si.plot(mode='stream')
plt.savefig(dirname+f'jxz_y{y}.png')
plt.clf()

data.phisp[-1, z, :, :].val_si.plot()
data.jxy[-1, z, :, :].val_si.plot(mode='stream')
plt.savefig(dirname+f'jxy_z{z}.png')
plt.clf()

data.phisp[-1, :, :, x].val_si.plot(savefilename=dirname+f'phisp_x{x}.png')
data.phisp[-1, :, y, :].val_si.plot(savefilename=dirname+f'phisp_y{y}.png')
data.phisp[-1, z, :, :].val_si.plot(savefilename=dirname+f'phisp_z{z}.png')

data.phisp[:, :, :, x].val_si.gifplot(action='save', filename=dirname+f'phisp_x{x}.gif')
data.phisp[:, :, y, :].val_si.gifplot(action='save', filename=dirname+f'phisp_y{y}.gif')
data.phisp[:, z, :, :].val_si.gifplot(action='save', filename=dirname+f'phisp_z{z}.gif')
plt.clf()

data.nd1p[-1, :, :, x].val_si.plot(savefilename=dirname+f'nd1p_x{x}.png')
data.nd1p[-1, :, y, :].val_si.plot(savefilename=dirname+f'nd1p_y{y}.png')
data.nd1p[-1, z, :, :].val_si.plot(savefilename=dirname+f'nd1p_z{z}.png')

data.nd1p[:, :, :, x].val_si.gifplot(action='save', filename=dirname+f'nd1p_x{x}.gif')
data.nd1p[:, :, y, :].val_si.gifplot(action='save', filename=dirname+f'nd1p_y{y}.gif')
data.nd1p[:, z, :, :].val_si.gifplot(action='save', filename=dirname+f'nd1p_z{z}.gif')
plt.clf()

data.nd2p[-1, :, :, x].val_si.plot(savefilename=dirname+f'nd2p_x{x}.png')
data.nd2p[-1, :, y, :].val_si.plot(savefilename=dirname+f'nd2p_y{y}.png')
data.nd2p[-1, z, :, :].val_si.plot(savefilename=dirname+f'nd2p_z{z}.png')

data.nd2p[:, :, :, x].val_si.gifplot(action='save', filename=dirname+f'nd2p_x{x}.gif')
data.nd2p[:, :, y, :].val_si.gifplot(action='save', filename=dirname+f'nd2p_y{y}.gif')
data.nd2p[:, z, :, :].val_si.gifplot(action='save', filename=dirname+f'nd2p_z{z}.gif')
plt.clf()
