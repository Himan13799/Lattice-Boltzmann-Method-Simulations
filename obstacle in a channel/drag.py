import numpy as np
import matplotlib.pyplot as plt
L=10
U=0.04
avg=[]
t=[i for i in range(200000)]
avg=np.loadtxt("avgrho.txt")
drag = np.loadtxt("drag.txt").T
lift = np.loadtxt("lift.txt").T
dc = (drag[1]/(avg))/(U*U*L)
lc = (lift[1]/(avg))/(U*U*L)
plt.plot(t,dc)
plt.xlabel("t")
plt.ylabel("C_d")
plt.title("C_d vs t (Re=60)")
plt.show()
plt.plot(t,lc)
plt.xlabel("t")
plt.ylabel("C_l")
plt.title("C_l vs t (Re=60)")
plt.show()