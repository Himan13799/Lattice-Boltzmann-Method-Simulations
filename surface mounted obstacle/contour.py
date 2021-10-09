import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
import matplotlib.animation as animation

w = 1
Y, X = np.mgrid[0:w:100j, 0:w:300j]
fig = plt.figure(figsize = (6, 6))
ax=fig.add_axes([0,0,1,1])

U = np.loadtxt("ux0.txt", dtype=float)
V = np.loadtxt("uy0.txt", dtype=float)
color = U**2+V**2
fig, at = plt.subplots()
fig.colorbar(cm.ScalarMappable(cmap=cm.jet))
im=plt.imshow(np.sqrt(U**2+V**2),cmap=cm.jet)

def animate(iter):
    at.collections = [] # clear lines streamplot
    at.patches = [] # clear arrowheads streamplot
    namex="ux"+str((iter+1)*50)+".txt"
    namey="uy"+str((iter+1)*50)+".txt"
    U = np.loadtxt(namex, dtype=float)
    V = np.loadtxt(namey, dtype=float)
    plt.axis('off')
    fig.subplots_adjust(left=0, bottom=0, right=1, top=1, wspace=None, hspace=None)
    im=plt.imshow(np.sqrt(U**2+V**2),cmap=cm.jet,origin='lower')
    print(iter)
    return im

anim =   animation.FuncAnimation(fig, animate, frames=200, interval=400, blit=False, repeat=False)
anim.save('./contour.gif', writer='imagemagick', fps=2)