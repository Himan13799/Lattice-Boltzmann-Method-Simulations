import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
import matplotlib.animation as animation

w = 1
Y, X = np.mgrid[0:w:100j, 0:w:300j]
fig = plt.figure(figsize = (6, 6))
ax=fig.add_axes([0,0,1,1])

U = np.loadtxt("ux195000.txt", dtype=float)
V = np.loadtxt("uy195000.txt", dtype=float)
fig, at = plt.subplots()
im=plt.imshow(np.sqrt(U**2+V**2),cmap=cm.jet)

def animate(iter):
    at.collections = [] # clear lines streamplot
    at.patches = [] # clear arrowheads streamplot
    namex="ux"+str((iter+1951)*100)+".txt"
    namey="uy"+str((iter+1951)*100)+".txt"
    U = np.loadtxt(namex, dtype=float)
    V = np.loadtxt(namey, dtype=float)
    plt.axis('off')
    fig.subplots_adjust(left=0, bottom=0, right=1, top=1, wspace=None, hspace=None)
    im=plt.imshow(np.sqrt(U**2+V**2),cmap=cm.jet,origin='lower')
    print(iter)
    return im

anim =   animation.FuncAnimation(fig, animate, frames=49, interval=400, blit=False, repeat=False)
anim.save('./animation.gif', writer='imagemagick', fps=7)