import emout
import matplotlib.pyplot as plt

data=emout.Emout("output/3688582")

y=int(data.inp.ny//2)
y_real=data.unit.length.reverse(y)

#show phi real unit
plt.clf()
#0 or -1
ax=data.phisp[-1, :, int(data.inp.ny//2), :].val_si.plot()
plt.title("Phi_sp at y={}[m]".format(y_real))
plt.savefig("phisp,y={}m.png".format(y_real))

#show Ex, Ey, Ez real unit
plt.clf()
ax=data.ex[-1, :, int(data.inp.ny//2), :].val_si.plot()
plt.title("Ex at y={}[m]".format(y_real))
plt.savefig("ex,y={}m.png".format(y_real))

plt.clf()
ax=data.ey[-1, :, int(data.inp.ny//2), :].val_si.plot()
plt.title("Ey_sp at y={}[m]".format(y_real))
plt.savefig("ey,y={}m.png".format(y_real))

plt.clf()
ax=data.ez[-1, :, int(data.inp.ny//2), :].val_si.plot()
plt.title("Ez at y={}[m]".format(y_real))
plt.savefig("ez,y={}m.png".format(y_real))
