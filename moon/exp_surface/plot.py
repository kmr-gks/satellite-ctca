import emout
import matplotlib.pyplot as plt
import os

job = '3822856'
dirname=job+'_plot/'

data = emout.Emout('output/' + job)
os.makedirs(job + '_plot',exist_ok=True)

#data.phisp[-1, :, :, 256].val_si.plot()
print(data.phisp[:].val_si.shape)
x, y, z = 128, 128, 60

data.phisp[:, z, :, :].val_si.gifplot(action='save', filename=dirname+'phisp_z.gif') # for save on a file


data.j1yz[-1, :, :, x].val_si.plot(mode='stream', savefilename=dirname+"j1yz_x.png")
data.j1xz[-1, :, y, :].val_si.plot(mode='stream', savefilename=dirname+"j1xz_y.png")
data.j1xy[-1, z, :, :].val_si.plot(mode='stream', savefilename=dirname+"j1xy_z.png")

data.j2yz[-1, :, :, x].val_si.plot(mode='stream', savefilename=dirname+"j2yz_x.png")
data.j2xz[-1, :, y, :].val_si.plot(mode='stream', savefilename=dirname+"j2xz_y.png")
data.j2xy[-1, z, :, :].val_si.plot(mode='stream', savefilename=dirname+"j2xy_z.png")


data.jyz[-1, :, :, x].val_si.plot(mode='stream', savefilename=dirname+"jyz_x.png")
data.jxz[-1, :, y, :].val_si.plot(mode='stream', savefilename=dirname+"jxz_y.png")
data.jxy[-1, z, :, :].val_si.plot(mode='stream', savefilename=dirname+"jxy_z.png")

data.phisp[-1, :, :, x].val_si.plot(savefilename=dirname+"phisp_x.png")
data.phisp[-1, :, y, :].val_si.plot(savefilename=dirname+"phisp_y.png")
data.phisp[-1, z, :, :].val_si.plot(savefilename=dirname+"phisp_z.png")

data.nd1p[-1, :, :, x].val_si.plot(savefilename=dirname+"nd1p_x.png")
data.nd1p[-1, z, :, :].val_si.plot(savefilename=dirname+"nd1p_z.png")
data.nd1p[-1, :, y, :].val_si.plot(savefilename=dirname+"nd1p_y.png")

data.nd2p[-1, :, :, x].val_si.plot(savefilename=dirname+"nd2p_x.png")
data.nd2p[-1, :, y, :].val_si.plot(savefilename=dirname+"nd2p_y.png")
data.nd2p[-1, z, :, :].val_si.plot(savefilename=dirname+"nd2p_z.png")

