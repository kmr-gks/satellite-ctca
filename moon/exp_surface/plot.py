import emout
import matplotlib.pyplot as plt
import os

job = '3982757'
dirname=job+'_plot/'

data = emout.Emout('output/' + job)
print("read:", job)
os.makedirs(job + '_plot',exist_ok=True)

print(data.phisp[:].val_si.shape)
x, y, z = 128, 128, 60 #grid point

data.phisp[-1, :, :, x].val_si.plot()
data.j1yz[-1, :, :, x].val_si.plot(mode='stream')
plt.savefig(dirname+'j1yz_x.png')
plt.clf()

data.phisp[-1, :, y, :].val_si.plot()
data.j1xz[-1, :, y, :].val_si.plot(mode='stream')
plt.savefig(dirname+'j1xz_y.png')
plt.clf()

data.phisp[-1, z, :, :].val_si.plot()
data.j1xy[-1, z, :, :].val_si.plot(mode='stream')
plt.savefig(dirname+'j1xy_z.png')
plt.clf()

data.phisp[-1, :, :, x].val_si.plot()
data.j2yz[-1, :, :, x].val_si.plot(mode='stream')
plt.savefig(dirname+'j2yz_x.png')
plt.clf()

data.phisp[-1, :, y, :].val_si.plot()
data.j2xz[-1, :, y, :].val_si.plot(mode='stream')
plt.savefig(dirname+'j2xz_y.png')
plt.clf()

data.phisp[-1, z, :, :].val_si.plot()
data.j2xy[-1, z, :, :].val_si.plot(mode='stream')
plt.savefig(dirname+'j2xy_z.png')
plt.clf()

data.phisp[-1, :, :, x].val_si.plot()
data.jyz[-1, :, :, x].val_si.plot(mode='stream')
plt.savefig(dirname+'jyz_x.png')
plt.clf()

data.phisp[-1, :, y, :].val_si.plot()
data.jxz[-1, :, y, :].val_si.plot(mode='stream')
plt.savefig(dirname+'jxz_y.png')
plt.clf()

data.phisp[-1, z, :, :].val_si.plot()
data.jxy[-1, z, :, :].val_si.plot(mode='stream')
plt.savefig(dirname+'jxy_z.png')
plt.clf()

data.phisp[-1, :, :, x].val_si.plot(savefilename=dirname+"phisp_x.png")
data.phisp[-1, :, y, :].val_si.plot(savefilename=dirname+"phisp_y.png")
data.phisp[-1, z, :, :].val_si.plot(savefilename=dirname+"phisp_z.png")

data.phisp[:, :, :, x].val_si.gifplot(action='save', filename=dirname+'phisp_x.gif')
data.phisp[:, :, y, :].val_si.gifplot(action='save', filename=dirname+'phisp_y.gif')
data.phisp[:, z, :, :].val_si.gifplot(action='save', filename=dirname+'phisp_z.gif')

data.nd1p[-1, :, :, x].val_si.plot(savefilename=dirname+"nd1p_x.png")
data.nd1p[-1, :, y, :].val_si.plot(savefilename=dirname+"nd1p_y.png")
data.nd1p[-1, z, :, :].val_si.plot(savefilename=dirname+"nd1p_z.png")

data.nd1p[:, :, :, x].val_si.gifplot(action='save', filename=dirname+'nd1p_x.gif')
data.nd1p[:, :, y, :].val_si.gifplot(action='save', filename=dirname+'nd1p_y.gif')
data.nd1p[:, z, :, :].val_si.gifplot(action='save', filename=dirname+'nd1p_z.gif')

data.nd2p[-1, :, :, x].val_si.plot(savefilename=dirname+"nd2p_x.png")
data.nd2p[-1, :, y, :].val_si.plot(savefilename=dirname+"nd2p_y.png")
data.nd2p[-1, z, :, :].val_si.plot(savefilename=dirname+"nd2p_z.png")

data.nd2p[:, :, :, x].val_si.gifplot(action='save', filename=dirname+'nd2p_x.gif')
data.nd2p[:, :, y, :].val_si.gifplot(action='save', filename=dirname+'nd2p_y.gif')
data.nd2p[:, z, :, :].val_si.gifplot(action='save', filename=dirname+'nd2p_z.gif')
