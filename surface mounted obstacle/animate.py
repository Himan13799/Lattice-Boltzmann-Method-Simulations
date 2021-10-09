import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

w = 1
Y, X = np.mgrid[0:w:100j, 0:w:300j]
fig = plt.figure()
ax=fig.add_axes([0,0,1,1])
ims = []

U = np.loadtxt("ux0.txt", dtype=float)
V = np.loadtxt("uy0.txt", dtype=float)
color = U**2+V**2
fig, at = plt.subplots()
stream = ax.streamplot(X,Y,U,V, color=color, density=2.7, cmap='jet',arrowsize=1)

def animate(iter):
    at.collections = [] # clear lines streamplot
    at.patches = [] # clear arrowheads streamplot
    namex="ux"+str((iter+1)*50)+".txt"
    namey="uy"+str((iter+1)*50)+".txt"
    U = np.loadtxt(namex, dtype=float)
    V = np.loadtxt(namey, dtype=float)
    color = U**2+V**2
    plt.axis('off')
    fig.subplots_adjust(left=0, bottom=0, right=1, top=1, wspace=None, hspace=None)
    stream = plt.streamplot(X,Y,U,V,color=color,density=2.7,cmap='jet',arrowsize=1)
    print(iter)
    return stream

anim =   animation.FuncAnimation(fig, animate, frames=200, interval=400, blit=False, repeat=False)
anim.save('./animation.gif', writer='imagemagick', fps=2)
