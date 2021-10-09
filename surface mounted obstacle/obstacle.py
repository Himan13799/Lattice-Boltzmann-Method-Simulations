import numpy as np
import math
import matplotlib.pyplot as plt
from matplotlib import cm

#Definition of lattice constants
rho0=1 #Lattice density
U=0.04 #Lattice velocity
NJ=100 #Number of lattices in x-direction
NI=300 #Number of lattices in y-direction
NSTEPS=10000 #Total number of time-steps
RE=200 #Reynolds Number
WI=50 #Write Interval - Number of time-steps after which plots and files are prepared
WI_DL=50 #Write Interval for Drag and Lift forces
DL_COUNT=0 #Counter

#Obstacle parameters
R=20 #Radius of circle or Edge-length of Square
CX=NI // 2 #X-coordinate of Center of obstacle
CY=R/2 #Y-coordinate of Center of obstacle

OBSTACLE=np.zeros((NI,NJ),dtype=bool)

#Major and minor and minor axes of Ellipse
B=5
A=1.5*B

#Airfoil
AIRFOIL=30 #Last two digits of airfoil (NACA00XX)
CHORD=30 #Chord length in terms of lattices
AOA=-5 #Angle of Attack (Positive->Anti-clockwise)

#Calculation of LBM parameters
kinematic_viscosity=(U*2*R)/RE #Viscosity 
OMEGA=1/(3.*kinematic_viscosity + 0.5) #Time-relaxation parameter

C=np.array([[1,1],[1,0],[1,-1],[0,1],[0,0],[0,-1],[-1,1],[-1,0],[-1,-1]])
W=np.array([ 1/36, 1/9, 1/36, 1/9, 4/9, 1/9, 1/36, 1/9, 1/36])
COL1 = np.array([0,1,2])  #Left wall of Lattice
COL2 = np.array([3,4,5])  #Interior of Lattice
COL3 = np.array([6,7,8])  #Right wall of Lattice
ROW1 = np.array([0,3,6])  #Bottom wall of Lattice
ROW2 = np.array([1,4,7])  #Interior of Lattice
ROW3 = np.array([2,5,8])  #Top wall of Lattice

#Initialization of Arrays
FIN=np.zeros((9,NI,NJ))
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
    vel_field[0,int(CX-R/2):int(CX+R/2),int(CY-R/2):int(CY+R/2)]=0

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
                
def CALC_DRAG_LIFT(): #Calculation of Drag and Lift Forces
    global DL_COUNT
    DRAG[DL_COUNT]=0
    LIFT[DL_COUNT]=0
    for i in range(1,NI-1):
        for j in range(1,NJ-1):
            if (OBSTACLE[i-1,j-1] or OBSTACLE[i-1,j] or OBSTACLE[i-1,j+1] or OBSTACLE[i,j-1] or OBSTACLE[i,j+1] or OBSTACLE[i+1,j-1] or OBSTACLE[i+1,j] or OBSTACLE[i+1,j+1]):
                for k in range(9):
                    DRAG[DL_COUNT]=DRAG[DL_COUNT] + C[k,0]*(FIN[k,i,j]+FOUT[8-k,i,j])
                    LIFT[DL_COUNT]=LIFT[DL_COUNT] + C[k,1]*(FIN[k,i,j]+FOUT[8-k,i,j]) 
    DL_COUNT +=1

def OBSTACLE_SQUARE():
    for i in range(R):
      for j in range(R):
        OBSTACLE[int(CX-R/2)+i][int(CY-R/2)+j]=True

"""
def OBSTACLE_CIRCLE(X,Y):  
    return (X-CX)**2 + (Y-CY)**2 < (R)**2  #Equation of Circle

def OBSTACLE_ROTATED_SQUARE(X,Y):
    return abs(X-CX) + abs(Y-CY) < (R/math.sqrt(2))  #Equation of Rotated Square

def OBSTACLE_ELLIPSE(X,Y):
    return ((X-CX)/A)**2 +((Y-CY)/B)**2 < 1  #Equation of Ellipse
    
def OBSTACLE_AIRFOIL(X,Y):
    AOA_RAD=(AOA*np.pi/180)
    y_u=(5*0.3*CHORD*(0.2969*(((X-CX)/CHORD)**0.5) - 0.126*((X-CX)/CHORD) - 0.3516*(((X-CX)/CHORD)**2) + 0.2843*(((X-CX)/CHORD)**3) - 0.1036*(((X-CX)/CHORD)**4)))
    y_l=-y_u
    Y_U=y_u*np.cos(AOA_RAD) - (X-CX)*np.sin(AOA_RAD)
    Y_L=y_l*np.cos(AOA_RAD) - (X-CX)*np.sin(AOA_RAD)
    return ((X-CX)>0) & ((X-CX)<(CHORD*np.cos(AOA_RAD))) & ((Y-CY)>Y_L) & ((Y-CY)<Y_U) #Equation of NACA00XX Airfoil
"""

def BC_OUTLET(): #Outflow Boundary Condition 
    FIN[COL3,-1,:] = FIN[COL3,-2,:]

def BC_INLET(): #Velocity at Inlet
    vel_field[0,0,:] = U
    vel_field[1,0,:] = 0
    RHO[0,:] = 1./(1.-vel_field[0,0,:]) * (np.sum(FIN[COL2,0,:],axis=0) + 2*np.sum(FIN[COL3,0,:],axis=0))
    CALC_EQUILIBRIUM()
    FIN[COL1,0,:] = FEQ[COL1,0,:] + FIN[[8,7,6],0,:] - FEQ[[8,7,6],0,:]
   
def BC_FREESLIP(): #Free-slip Boundary Condition at the top and bottom boundaries
    FIN[ROW1,:,0]=FIN[[8,5,2],:,0]
    FIN[ROW3,:,-1]=FIN[[6,3,0],:,-1]

def BC_OBSTACLE(): #Bounce-back
    for k in range(9):
        FOUT[k, OBSTACLE]=FIN[8-k, OBSTACLE]
        
def COLLISION():
    global FOUT
    FOUT = FIN - OMEGA*(FIN-FEQ)

def STREAMING():
    for k in range(9):
        FIN[k,:,:] = np.roll(np.roll(FOUT[k,:,:],C[k,0],axis=0),C[k,1],axis=1)
        
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
    
"""
def WRITE_DRAG_LIFT(): #Create forces file
    f4=open("drag_lift.txt","w+")
    f4.write('VARIABLES = "T", "DRAG", "LIFT" \n')
    for n in range(DL_COUNT):
        f4.write(str(n*WI_DL)+" "+str(DRAG[n])+" "+str(LIFT[n])+"\n")
    f4.close()    
"""    
#VEL = np.fromfunction(CALC_INIVEL, (2,NI,NJ))
OBSTACLE_SQUARE()
#OBSTACLE = np.fromfunction(OBSTACLE_SQUARE, (NI,NJ)) #Changable obstacle definition 

#Initialization of solution
RHO[:,:] = rho0
vel_field[0,0,:] = U
CALC_EQUILIBRIUM()
FIN=FEQ.copy()

def error(U0,vel_field):
  num=np.sum((U0[0]-vel_field[0])**2+(U0[1]-vel_field[1])**2)
  denom=np.sum(vel_field[0]*vel_field[0]+vel_field[1]*vel_field[1])+1e-30
  return num/denom


color = vel_field[0]**2+vel_field[1]**2
fig, ax = plt.subplots()
im=plt.imshow(np.sqrt(vel_field[0]**2+vel_field[1]**2),cmap=cm.jet)
for t in range(NSTEPS+1):
    if t%1000==0:
      print(t)
    BC_OUTLET()
    BC_FREESLIP()
    U0=vel_field
    CALC_QUANTITIES()
    BC_INLET()
    COLLISION()
    BC_OBSTACLE()
    STREAMING()
    if t%1000==0:
      e=error(U0,vel_field)
      print(e)
    if t % WI_DL==0:
        CALC_DRAG_LIFT()
    if t % WI ==0:
        CALC_VORTICITY()
        WRITE_FILES()
        SAVE_PLOTS()