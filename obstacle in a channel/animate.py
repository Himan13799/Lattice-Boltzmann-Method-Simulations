import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

w = 1
Y, X = np.mgrid[0:w:90j, 0:w:260j]
fig = plt.figure()
ax=fig.add_axes([0,0,1,1])
ims = []

"""
for i in range(72):
	print(i+1)
	namex="ux"+str((i+1)*1000)+".txt"
	namey="uy"+str((i+1)*1000)+".txt"
	U = np.loadtxt(namex, dtype=float)
	V = np.loadtxt(namey, dtype=float)
	im = plt.streamplot(X, Y, U, V, density = 2.7)
	ims.append([im])
"""
U = np.loadtxt("ux190000.txt", dtype=float)
V = np.loadtxt("uy190000.txt", dtype=float)
color = U**2+V**2
fig, at = plt.subplots()
stream = at.streamplot(X,Y,U,V, color=color, density=4, cmap='jet',arrowsize=1)

def animate(iter):
    at.collections = [] # clear lines streamplot
    at.patches = [] # clear arrowheads streamplot
    namex="ux"+str((iter+1901)*100)+".txt"
    namey="uy"+str((iter+1901)*100)+".txt"
    U = np.loadtxt(namex, dtype=float)
    V = np.loadtxt(namey, dtype=float)
    color = U**2+V**2
    plt.axis('off')
    fig.subplots_adjust(left=0, bottom=0, right=1, top=1, wspace=None, hspace=None)
    stream = plt.streamplot(X,Y,U,V,color=color,density=4,cmap='jet',arrowsize=1)
    print(iter)
    return stream

anim =   animation.FuncAnimation(fig, animate, frames=99, interval=200, blit=False, repeat=False)
anim.save('./Streamlines.gif', writer='imagemagick', fps=7)