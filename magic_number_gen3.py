import numpy as np
s=31
ss=s-(s&1)
s2=s//2
l=np.linspace(-ss,ss,s)
X,Y=np.meshgrid(l,l)
Z=(X*X+Y*Y)**.5

a=np.full_like(Z,np.inf)
a[0:2]=np.log(Z[0:2])
def rot(a):
    a=np.minimum(a,a.T)
    a=np.minimum(a,a[::-1,::-1])
    
    return a
for i in range(2,s//2+1):
    #a=rot(a)
    a[i,i:-i]=4*a[i-1,i:-i]-a[i-1,i-1:-i-1]-a[i-1,i+1:s-i+1]-a[i-2,i:-i]
import matplotlib.pyplot as plt
a=rot(a)
for i in np.linspace(0.001,5,30):
    a=np.log(Z+i)
    a[1:-1,1:-1]-=0.25*(a[2:,2:]+a[2:,:-2]+a[:-2,2:]+a[:-2,:-2])
    a[s2,s2]=0
    print(i,np.abs(a[1:-1,1:-1]).sum())
#a[1:-1,1:-1]-=0.25*(a[2:,2:]+a[2:,:-2]+a[:-2,2:]+a[:-2,:-2])
plt.imshow(a)
plt.show()
