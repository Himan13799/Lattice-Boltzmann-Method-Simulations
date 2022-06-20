import numpy as np
import matplotlib.pyplot as plt

x=65
r=10
y=45
u_time_series=[]
vts=[]
t=[]
for i in range(2000):
	U=np.loadtxt("ux"+str(i*100)+".txt").T
	V=np.loadtxt("uy"+str(i*100)+".txt").T
	u_time_series.append(U[x+r+2][y])
	vts.append(V[x+r+2][y])
	t.append(i*100)
plt.title("u at rare stagnation point(Re=60) vs t")
plt.xlabel("t")
plt.ylabel("velocity time series")
plt.plot(t,u_time_series)
plt.plot(t,vts)
plt.legend(["Axial velocity","y velocity"])
plt.show()