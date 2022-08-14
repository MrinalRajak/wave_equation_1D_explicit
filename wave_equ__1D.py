



# wave equation in 1D (Explicit) ,constant number C=vel*(dt/dx)
import numpy as np
import matplotlib.pyplot as plt
L=1.0;T=1.0;dt=0.001;Nt=int(round(T/dt))
t=np.linspace(0,Nt*dt,Nt+1) # Mesh points in time
vel=10. # velocity of wave
c=0.9 # (constant number should be kept less than 1 for numerical stability)
dx=(dt*vel)/float(c)
Nx=int(round(L/dx))
x=np.linspace(0,L,Nx+1) # Mesh points in space 
c2=c**2 # Help variable in the scheme
s='t%f'
u=np.zeros(Nx+1) # solution array at new time level
u_n=np.zeros(Nx+1) # solution at 1 time level back
u_nml=np.zeros(Nx+1) # solution at 2 time level back
# Initial distribution at t=0
sig=0.1

def I(x):
    return (1/np.sqrt(2*np.pi*sig))*np.exp((-(x-L/2)**2)/(2*sig**2))
'''
def I(x):
    return x*(L-x)
'''
# set initial condition u(x,0)=I(x)
for i in range(0,Nx+1):
    u_n[i]=I(x[i])
# Apply special formula for the first step
for i in range(1,Nx):
    u[i]=u_n[i]+0.5*c2*(u_n[i+1]-2*u_n[i]+u_n[i-1])
# Incorporating u=0 ,boundary conditions4
u[0]=0;u[Nx]=0 # at both boundaries
# switch variables before next step
u_nml[:],u_n[:]=u_n,u
for t in range(1,Nt):
    # update all inner mesh points
    for i in range(1,Nx):
        u[i]=2*u_n[i]-u_nml[i]+c2*(u_n[i+1]-2*u_n[i]+u_n[i-1])
    # insert boundary conditions
    u[0]=0;u[Nx]=0
    # switch variables before next time step
    u_nml[:],u_n[:]=u_n,u
    if(t%10==0):
        plt.plot(x,u,label=s%t)
        plt.legend()
print(u[i])
plt.show()
    































