import numpy as np
import matplotlib.pyplot as plt

w = 1
Y, X = np.mgrid[0:w:90j, 0:w:260j]
l = [166900,167800,168300,169500,170200,170700,172700,173700,174500,175700,176400,177500,178700,179700,181100]
for i in l:
	t=str(i)
	U = np.loadtxt("ux"+t+".txt",dtype=float)
	V = np.loadtxt("uy"+t+".txt",dtype=float)
	color = np.sqrt(U**2+V**2)
	plt.clf()
	stream = plt.streamplot(X,Y,U,V, color=color, density=5, cmap='jet',arrowsize=1)
	plt.title("streamlines @ t="+t)
	plt.show()