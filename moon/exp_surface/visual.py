import emout
import matplotlib.pyplot as plt

data=emout.Emout("output/3688582")

y=int(data.inp.ny//2)
y_real=data.unit.length.reverse(y)

#show phi real unit
plt.clf()
ax=data.phisp[-1, :, int(data.inp.ny//2), :].val_si.plot()
plt.title("Phi_sp at y={}[m]".format(y_real))
plt.savefig("phisp,y={}m.png".format(y_real))
