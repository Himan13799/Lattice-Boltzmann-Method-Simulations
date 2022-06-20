from scipy.fftpack import fft,fftfreq
import numpy as np
# Number of samplepoints
N = 1000
x=65
r=10
h=45
y=[]
for i in range(N):
	U = np.loadtxt("ux"+str((i+2000-N)*100)+".txt").T
	y.append(U[x+r+2][h])
yf = fft(y)
yf=yf[1:]
xf = fftfreq(len(yf))
i = xf>0
import matplotlib.pyplot as plt
print(np.argmax(np.abs(yf[i])))
plt.plot(xf[i], np.abs(yf[i]))
plt.grid()
plt.show()