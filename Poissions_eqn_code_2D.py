import numpy as np
import matplotlib.pyplot as plt

from timeit import default_timer as timer

start = timer()

Lx=1
Ly=1
n=51 #grid points
m=51
dx= Lx/n-1 # assuming uniformly spaced nodes
dy= Ly/m-1
x = np.linspace(0, Lx, n)
y = np.linspace(0, Ly, m)
xx,yy=np.meshgrid(x,y)

p=np.zeros([n,m])

# # source term function (=0 for now making poisson's eqn )
S=np.zeros([n,m])

tolerance = 1e-4
error=1e10
iter_count=0
max_iter_count=6000
while error > tolerance and iter_count<max_iter_count:
    p_k=p.copy()
    
    # boundary conditions applied
    ## already dirchilet BCs applied by default
    #P_right=0
    p[:,-1]=0
    #p_left=0
    p[:,0]=0
    p[10:-10,0]=1 #flow in a box
    #p_bottom=0
    p[0,:]=0
    #p_top=0
    p[-1,:]=0

    for j in range(1,n-1):
        for i in range(1,m-1):
            #gauss siedel iterative solver
            p[i,j]=(S[i,j])*((dx*dx*dy*dy)/(2*(dx*dx+dy*dy))) + ((dy*dy)/(2*(dx*dx+dy*dy)))*(p_k[i,j+1]+p[i,j-1]) + ((dx*dx)/(2*(dx*dx+dy*dy)))*(p_k[i+1,j]+p[i-1,j])
            #jacobi iterative solver
            # p[i,j]=(S[i,j])*((dx*dx*dy*dy)/(2*(dx*dx+dy*dy))) + ((dy*dy)/(2*(dx*dx+dy*dy)))*(p_k[i,j+1]+p_k[i,j-1]) + ((dx*dx)/(2*(dx*dx+dy*dy)))*(p_k[i+1,j]+p_k[i-1,j])
    error= np.linalg.norm(p-p_k,2) #L2 norm of difference between value of p at k (p_k) and k+1 (p) iternation
    iter_count+=1

if iter_count==max_iter_count:
    print('diverging solution with: ',iter_count,' many iterations')
    print('error is ',error)
else:
    print ('solution converged in ',iter_count,' iterations')

end = timer()

print('time taken: ', end - start)

plt.contourf(xx,yy,p,cmap="inferno")
plt.colorbar()
plt.show()
