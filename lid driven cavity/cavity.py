import numpy as np
import math
import matplotlib.pyplot as plt
from matplotlib import cm

#Definition of lattice constants
rho0=2.7 #Lattice density
U=0.1 #Lattice velocity
NJ=256 #Number of lattices in x-direction
NI=256 #Number of lattices in y-direction
NSTEPS=10000 #Total number of time-steps : change if needed(to ensure convergence)
RE=2000 #Reynolds Number
WI=50 #Write Interval - Number of time-steps after which plots and files are prepared
DL_COUNT=0

#Calculation of LBM parameters
kinematic_viscosity=(U*NI)/RE #Viscosity 
OMEGA=1/(3.*kinematic_viscosity + 0.5) #Time-relaxation parameter

C=np.array([[1,1],[1,0],[1,-1],[0,1],[0,0],[0,-1],[-1,1],[-1,0],[-1,-1]])
"""
6 3 0
 \|/
7-4-1   velocity set D2Q9
 /|\ 
8 5 2
"""
W=np.array([1/36, 1/9, 1/36, 1/9, 4/9, 1/9, 1/36, 1/9, 1/36])
COL1 = np.array([0,1,2])  #Left wall of Lattice
COL2 = np.array([3,4,5])  #Interior of Lattice
COL3 = np.array([6,7,8])  #Right wall of Lattice
ROW1 = np.array([0,3,6])  #Bottom wall of Lattice
ROW2 = np.array([1,4,7])  #Interior of Lattice
ROW3 = np.array([2,5,8])  #Top wall of Lattice

#Initialization of Arrays
FIN=np.zeros((9,NI,NJ))
tmp=np.zeros((9,NI,NJ))
FEQ=np.zeros((9,NI,NJ))
FOUT=np.zeros((9,NI,NJ))
RHO=np.zeros((NI,NJ))
vel_field=np.zeros((2,NI,NJ))
P=np.zeros((NI,NJ))
VORTICITY=np.zeros((NI,NJ))
DRAG=np.zeros(NSTEPS+1)
LIFT=np.zeros(NSTEPS+1)

def CALC_QUANTITIES(): #Calculation of Macroscopic Quantities
    global vel_field,P,RHO
    RHO=np.sum(FIN, axis=0) #Density
    P=RHO/3  #Pressure
    vel_field[:,:,:]=0
    for k in range(9):
        vel_field[0,:,:] += C[k,0] * FIN[k, :, :]
        vel_field[1,:,:] += C[k,1] * FIN[k, :, :]
    vel_field /= RHO #Velocity
    RHO[:,NJ-2]=rho0
    vel_field[0,:,NJ-2]=U
    vel_field[1,:,NJ-2]=0

def CALC_EQUILIBRIUM(): #Calculation of Equilibrium Function
    USQR = 3./2. * (vel_field[0]**2 + vel_field[1]**2)
    for k in range(9):
        VU = 3.0 * (C[k,0]*vel_field[0,:,:] + C[k,1]*vel_field[1,:,:])
        FEQ[k,:,:] = RHO*W[k] * (1 + VU + 0.5*VU**2 - USQR)

def CALC_VORTICITY(): #Vorticity Calculation
    for i in range(NI):
        for j in range(NJ):
            if i==0 or j==0 or i==(NI-1) or j==(NJ-1):
                VORTICITY[i,j]=0
            else:
                VORTICITY[i,j] = (vel_field[1,i+1,j] - vel_field[1,i-1,j] - vel_field[0,i,j+1] + vel_field[0,i,j-1])/2

def COLLISION():
    global FOUT,tmp
    tmp=FOUT
    FOUT = FIN - OMEGA*(FIN-FEQ)

def STREAMING():
    for i in range(NI-2):
        for j in range(NJ-2):
            for k in range(9):
                x=i+1-C[k,0]
                y=j+1-C[k,1]
                if x==0 or x==NI-1 or y==0 or y==NJ-1: FIN[k,i+1,j+1]=FIN[8-k,i+1,j+1]
                else: FIN[k,i+1,j+1] = FOUT[k,x,y]
        
def WRITE_FILES(): #Create post-processing files
    f1=open("ux"+str(t)+".txt","w+")
    f2=open("uy"+str(t)+".txt","w+")
    f3=open("vort"+str(t)+".txt","w+")
    for j in range(NJ):
        for i in range(NI):
            f1.write(str(vel_field[0,i,j])+" ")
            f2.write(str(vel_field[1,i,j])+" ")
            #f2.write(P[i,j])
            f3.write(str(VORTICITY[i,j])+" ")
        f1.write('\n')
        f2.write('\n')
        f3.write('\n')
    f1.close()
    f2.close()
    f3.close()

def SAVE_PLOTS(): #Create post-processing plots
    plt.clf() 
    plt.axis('off')
    plt.imshow(np.sqrt(vel_field[0]**2+vel_field[1]**2).transpose(),cmap=cm.jet)
    plt.savefig("vel_"+str(t)+".png")
    


#Initialization of solution
RHO[:,:] = rho0
vel_field[0,:,NJ-2]=U
vel_field[1,:,NJ-2]=0
CALC_EQUILIBRIUM()
FIN=FEQ.copy()
FOUT=FIN.copy()

def error(U0,vel_field):
  num=np.sum((U0[0]-vel_field[0])**2+(U0[1]-vel_field[1])**2)
  denom=np.sum(vel_field[0]*vel_field[0]+vel_field[1]*vel_field[1])+1e-30
  return num/denom


color = vel_field[0]**2+vel_field[1]**2
fig, ax = plt.subplots()
im=plt.imshow(np.sqrt(vel_field[0]**2+vel_field[1]**2),cmap=cm.jet)
for t in range(NSTEPS+1):
    print(t)
    STREAMING()
    U0=vel_field
    CALC_QUANTITIES()
    CALC_EQUILIBRIUM()
    COLLISION()
    FIN=tmp
    if t%1000==0:
      e=error(U0,vel_field)
      print(e)
    if t % WI ==0:
        CALC_VORTICITY()
        WRITE_FILES()
        SAVE_PLOTS()