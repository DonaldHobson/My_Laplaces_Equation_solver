import numpy as np
def f(x):
    return np.sin(x)
    if x==0:
        return 0
    return -0.25j*x*x*(2*np.log(x)-3)
a=f(.5+.5j)-f(.5)-f(.5j)+f(0)
a=4*a.real+0j
#print(a*4)
def g(x,y):
    z=x+1j*y
    return f(z+.5+.5j)+f(z-.5-.5j)-f(z+.5-.5j)-f(z-.5+.5j)
    #return sum((-1)**(i!=j)*f(x+.5*i+y*1j+.5j*j) for i,j in ((1,1),(1,-1),(-1,1),(-1,-1)))
def disp(x):
    print("Complex{%s,%s},"%(x.real,x.imag),end="")
u_input=1#int(input("Max jump (default=2)"))
disp(a)
for i in range(1,u_input):
    for j in range(i+1):
        disp(g(i,j))
l=np.arange(20)
X,Y=np.meshgrid(l,l)
import matplotlib.pyplot as plt
Q=np.log(X+1j*Y)
plt.imshow(Q.real)
plt.show()
Z=np.zeros_like(X,complex)
for i in range(20):
    for j in range(20):
        Z[i,j]=g(j-9,i-9)
Z[9,9]=a
w=Z[2:,1:-1]+Z[:-2,1:-1]+Z[1:-1,2:]+Z[1:-1,:-2]-4*Z[1:-1,1:-1]
plt.imshow(w.real)
plt.show()
plt.imshow(w.imag)
plt.show()
plt.imshow(Z.real)
plt.show()
##X,Y=np.meshgrid(np.linspace(1.5,2.5,1000),np.linspace(2.5,3.5,1000))
##Z=np.log(X+1j*Y).mean()
##print(Z)
plt.imshow((Z-1j*np.arctan2(Y,X)).imag)
plt.show()
plt.imshow((Z-0.5*np.log(X*X+Y*Y)).real)
plt.show()
plt.imshow((Q-Z).imag)
plt.show()
