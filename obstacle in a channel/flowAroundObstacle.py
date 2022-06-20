
#
# 2D flow around a cylinder
#

from numpy import *; from numpy.linalg import *
import matplotlib.pyplot as plt; from matplotlib import cm
###### Flow definition #########################################################
maxIter = 200000 # Total number of time iterations.
Re      = 220.0  # Reynolds number.
nx = 260; ny = 90; ly=ny-1.0; q = 9 # Lattice dimensions and populations.
cx = nx/4; cy=ny/2; r=ny/9;          # Coordinates of the cylinder.
uLB     = 0.04                       # Velocity in lattice units.
nulb    = uLB*r/Re; omega = 1.0 / (3.*nulb+0.5); # Relaxation parameter.
DRAG=zeros(maxIter+1)
LIFT=zeros(maxIter+1)
DL_COUNT=0
###### Lattice Constants #######################################################
c = array([(x,y) for x in [0,-1,1] for y in [0,-1,1]]) # Lattice velocities.
t = 1./36. * ones(q)                                   # Lattice weights.
t[asarray([norm(ci)<1.1 for ci in c])] = 1./9.; t[0] = 4./9.
noslip = [c.tolist().index((-c[i]).tolist()) for i in range(q)] 
i1 = arange(q)[asarray([ci[0]<0  for ci in c])] # Unknown on right wall.
i2 = arange(q)[asarray([ci[0]==0 for ci in c])] # Vertical middle.
i3 = arange(q)[asarray([ci[0]>0  for ci in c])] # Unknown on left wall.

###### Function Definitions ####################################################
sumpop = lambda fin: sum(fin,axis=0) # Helper function for density computation.
def equilibrium(rho,u):              # Equilibrium distribution function.
    cu   = 3.0 * dot(c,u.transpose(1,0,2))
    usqr = 3./2.*(u[0]**2+u[1]**2)
    feq = zeros((q,nx,ny))
    for i in range(q): feq[i,:,:] = rho*t[i]*(1.+cu[i]+0.5*cu[i]**2-usqr)
    return feq

def WRITE_FILES(vel_field,fin): #Create post-processing files
    f=open("avgrho.txt","w+")
    rho = sum(fin,axis=0)
    r=average(rho)
    f.write(str(r)+" ")
    f.close()
    f1=open("ux"+str(time)+".txt","w+")
    f2=open("uy"+str(time)+".txt","w+")
    #f3=open("vort"+str(t)+".txt","w+")
    #savetxt("fin"+str(time)+".txt",fin.reshape(fin.shape[0],-1))
    for j in range(ny):
        for i in range(nx):
            f1.write(str(vel_field[0,i,j])+" ")
            f2.write(str(vel_field[1,i,j])+" ")
        f1.write('\n')
        f2.write('\n')
    f1.close()
    f2.close()
def WRITE_DRAG_LIFT(): #Create forces file
    f4=open("drag.txt","w+")
    f5=open("lift.txt","w+")
    for n in range(DL_COUNT):
        f4.write(str(n)+" "+str(DRAG[n])+"\n")
        f5.write(str(n)+" "+str(LIFT[n])+"\n")
    f4.close()
    f5.close()

def CALC_DRAG_LIFT(FIN,FOUT): #Calculation of Drag and Lift Forces
    global DL_COUNT
    DRAG[DL_COUNT]=0
    LIFT[DL_COUNT]=0
    for i in range(1,nx-1):
        for j in range(1,ny-1):
            if (obstacle[i-1,j-1] or obstacle[i-1,j] or obstacle[i-1,j+1] or obstacle[i,j-1] or obstacle[i,j+1] or obstacle[i+1,j-1] or obstacle[i+1,j] or obstacle[i+1,j+1]):
                for k in range(9):
                    DRAG[DL_COUNT]=DRAG[DL_COUNT] + c[k,0]*(FIN[k,i,j]+FOUT[8-k,i,j])
                    LIFT[DL_COUNT]=LIFT[DL_COUNT] + c[k,1]*(FIN[k,i,j]+FOUT[8-k,i,j]) 
    DL_COUNT +=1
###### Setup: cylindrical obstacle and velocity inlet with perturbation ########
obstacle = fromfunction(lambda x,y: (abs((x-cx)+(y-cy))+abs((x-cx)-(y-cy)))<2*r, (nx,ny))
y=array([(i)*(ny-i-1) for i in range(ny)])
vel=zeros((2,nx,ny))
vel[0,0,:] = 4*y*uLB/((ny-1)**2)
vel[1,0,:] = 0
feq = equilibrium(1.0,vel); fin = feq.copy()

###### Main time loop ##########################################################
for time in range(maxIter):
    fin[i1,-1,:] = fin[i1,-2,:] # Right wall: outflow condition.
    rho = sumpop(fin)           # Calculate macroscopic density and velocity
    u = dot(c.transpose(), fin.transpose((1,0,2)))/rho

    u[:,0,:] =vel[:,0,:] # Left wall: compute density from known populations.
    rho[0,:] = 1./(1.-u[0,0,:]) * (sumpop(fin[i2,0,:])+2.*sumpop(fin[i1,0,:]))

    feq = equilibrium(rho,u) # Left wall: Zou/He boundary condition.
    fin[i3,0,:] = fin[i1,0,:] + feq[i3,0,:] - fin[i1,0,:]
    fout = fin - omega * (fin - feq)  # Collision step.
    for i in range(q): fout[i,obstacle] = fin[noslip[i],obstacle]
    for i in range(q): # Streaming step.
        fin[i,:,:] = roll(roll(fout[i,:,:],c[i,0],axis=0),c[i,1],axis=1)

    CALC_DRAG_LIFT(fin,fout)
    if (time%100==0): # Visualization
        print(time)
        #plt.clf(); plt.imshow(sqrt(u[0]**2+u[1]**2).transpose(),cmap=cm.Reds)
        #plt.savefig("vel."+str(time/100).zfill(4)+".png")
        #WRITE_FILES(u,fin)

WRITE_DRAG_LIFT()

