import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

w = 1
Y, X = np.mgrid[0:w:256j, 0:w:256j]
fig = plt.figure()
ax=fig.add_axes([0,0,1,1])
ims = []

U = np.loadtxt("ux1000.txt", dtype=float)
V = np.loadtxt("uy1000.txt", dtype=float)
color = U**2+V**2
fig, at = plt.subplots()
stream = at.streamplot(X,Y,U,V, color=color, density=2.7, cmap='jet',arrowsize=1)

def animate(iter):
    at.collections = [] # clear lines streamplot
    at.patches = [] # clear arrowheads streamplot
    namex="ux"+str((iter+1)*1000)+".txt"
    namey="uy"+str((iter+1)*1000)+".txt"
    U = np.loadtxt(namex, dtype=float)
    V = np.loadtxt(namey, dtype=float)
    color = U**2+V**2
    plt.axis('off')
    fig.subplots_adjust(left=0, bottom=0, right=1, top=1, wspace=None, hspace=None)
    stream = plt.streamplot(X,Y,U,V,color=color,density=2.7,cmap='jet',arrowsize=1)
    print(iter)
    return stream

anim =   animation.FuncAnimation(fig, animate, frames=72, interval=200, blit=False, repeat=False)
anim.save('./animStreamlines.gif', writer='imagemagick', fps=4)
