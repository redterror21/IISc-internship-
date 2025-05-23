import numpy as np 
import matplotlib.pyplot as plt
from matplotlib import cm  #for color maps
import matplotlib.animation as animation
import time

# grid and domain properties while using staggered grid 
# defining no. of points inside boundary

x,y=16,16
Lx,Ly=1,1
dx,dy=Lx/x,Ly/y

visc=0.01

# Problem: Top wall moving, rest are fixed,
Ut=2 
Ub=0
Vl=0
Vr=0

# check our Reynolds number
print('Reynolds number', Ut*Lx/visc) #ensuring laminar flow only with low numbers
# choose timestep size based on advection-diffusion stability criteria # USING Stability conditions
# rule 1:
dt1 = 0.5/(visc*(1/(dx*dx) + 1/(dy*dy)))# For Diffusive part von Neumann stability analysis (done on explicit time integration) to check numerical instability 
# to ensure diffusion doesn't propagate too quicly over timesteps
# rule 2:
C=2 #safety factor between 1 and 2
dt2 = C*visc/(Ut*Ut) #For advection part, wheere we let the viscous part of time to dominate  
dt = min(dt1, dt2) # so satisafies both advection and diffusion rule
print(dt)

# poission solver SOR
# Build the coefficient arrays 
Ap= np.zeros([y+2, x+2])
Ae = 1/(dx*dx)*np.ones([y+2, x+2])
Aw= 1/(dx*dx)*np.ones([y+2, x+2]) 
An=1/(dy*dy)*np.ones([y+2, x+2]) 
As=1/(dy*dy)*np.ones([y+2, x+2])

# left wall coefs: Aw = 0 
Aw[1:-1, 1] = 0
# right wall coefs: Ae = 0 
Ae[1:-1,-2] = 0
#top wall coefs: An = 0 
An[-2,1:-1] = 0
# bottom wall: As = 0 
As[1,1:-1] = 0
Ap=-(Ae+ Aw + An + As)


def sor_solver (p, S, dx, dy):
    pk = np.empty_like(p)
    it = 0
    err = 1e10
    tol = 1e-8
    maxit=10000
    B = 1.25
    while err > tol and it < maxit:
        pk = p.copy()
        for i in range(1, x+1):
            for j in range(1, y+1): 
                ap = Ap[j,i]
                an = An[j,i] 
                aso= As[j,i] 
                ae = Ae[j,i] 
                aw = Aw[j,i]
                pe = p[j,i+1]
                pw = p[j,i-1]
                pn = p[j+1,i]
                ps = p[j-1,i]
                result = S[j,i] - (ae*pe + aw*pw + an*pn + aso*ps)
                p[j,i] = B*result/ap + (1-B)*pk [j,i]
        err = np. linalg.norm(p.ravel()- pk.ravel(),2)
        it += 1
    print('Poisson solver error:', err, 'iterations', it)
    return p,err

#time advance code
p=np.zeros([y+2,x+2])
u=np.zeros([y+2,x+2]) # same mesh just displaced by half cell
v=np.zeros([y+2,x+2])
# from Solution Algorithm 
ut=np.zeros_like(u) # u tilda velocities values 
vt=np.zeros_like(v) #v tilda velocities values
divut=np.zeros([y+2,x+2])

t=0
tsteps=5

for n in range(0,tsteps):
    # setting boundary conditions:  done inside loop for varying veloctiy at boundary
    u[:,1]=0 #left wall impermeable so all j values at first i point = 0
    u[:,-1]=0 #right wall impermeable so all j values at last i point = 0
    # top wall: ughost = 2*Ut-uinterior
    u[-1,:] = 2*Ut-u[-2,:] 
    # bottom wall:
    u[0,:] = 2*Ub-u[1,:] # bottom wall
    
    v[-1,:] = 0 # top wall, impermeable 
    v[1,:] = 0 # bottom wall, impermeable 
    # left wall: vghost = 2V - vinterior 
    v[:,0] = 2*Vl-v[:,1]
    # right wall:
    v[:,-1]= 2*Vr-v[:,-2] 


    # x momentum equation
    for i in range(2,x +1): 
        for j in range(1,y +1): 
            ue = 0.5*(u[j,i+1] + u[j,i]) 
            uw = 0.5*(u[j,i] + u[j, i-1]) 
            un = 0.5*(u[j+1,i] + u[j,i]) 
            us = 0.5*(u[j,i] + u[j-1,i])
            vn = 0.5*(v[j+1, i-1] + v[j+1,i]) 
            vs=0.5*(v[j,i-1]+ v[j, i])
            convection=-(ue*ue-uw*uw)/dx-(un*vn-us*vs)/dy
            diffusion = visc* (u[j,i-1] - 2*u[j,i] + u[j,i+1]/(dx*dx) + (u[j-1,i] - 2*u[j,i] + u[j+1,i])/(dy*dy))
            ut[j,i]= u[j,i] + dt*(convection+diffusion)

    # y momentum equation
    for i in range(1,x +1): 
        for j in range(2,y +1): 
            ve = 0.5*(v[j,i+1]+v[j,i]) 
            vw = 0.5*(v[j,i] + v[j, i-1]) 
            ue = 0.5*(u[j,i+1]+u[j-1,i+1]) 
            uw = 0.5*(u[j,i] + u[j-1,i])
            vn = 0.5*(v[j+1,i]+v[j,i]) 
            vs = 0.5*(v[j,i] + v[j-1,i])
            convection= -(ue*ve-uw*vw)/dx-(vn*vn-vs*vs)/dy
            diffusion = visc* (v[j,i+1] - 2*v[j,i] + v[j,i-1]/(dx*dx) + (v[j+1,i] - 2*v[j,i] + v[j-1,i])/(dy*dy))
            vt[j,i]= v[j,i] + dt*(convection+diffusion)

    divut[1:-1,1:-1]=(ut[1:-1,2:] - ut[1:-1,1:-1])/dx + (vt[2:,1:-1] - vt[1:-1,1:-1])/dy
    
    p_rhs=divut/dt

    p,err=  sor_solver(p,p_rhs,dx,dy)

    
    u[1:-1,2:-1] = ut[1:-1, 2:-1] - dt*(p[1:-1, 2:-1]- p[1:-1,1:-2])/dx 
    v[2:-1, 1:-1] = vt[2:-1, 1:-1] - dt*(p[2:-1, 1:-1] - p[1:-2,1:-1])/dy
    t+=dt
plt.quiver(u,v,scale=5)
divu=np. zeros_like(p)
for i in range(1, x+1):
    for j in range(1, y+1):
        divu[j, i]=(u[j,i+1] - u[j,i])/dx + (v[j+1,i] -v[j,i])/dy
plt.imshow(divu, origin="lower")
plt.colorbar()
plt.show()

    





                                                        # for animations
# ani = animation.FuncAnimation(fig, update, frames=100, interval=50)

# Save as MP4 (requires ffmpeg)
# ani.save('animation.mp4', writer='ffmpeg')

# OR as a GIF (requires Pillow)
# ani.save('animation.gif', writer='pillow')

# plt.show()
                                                        # ----------
