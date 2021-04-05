import subprocess
import numpy as np#"-np", "1",
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
#def f(logs):
logs=6
s=2**logs
##sm=s-1
##xy=np.zeros([sm*4,2],int)
##xy[:sm,0]=np.arange(sm)
##xy[sm:sm*2,0]=sm
##xy[sm*2:sm*3,0]=np.arange(sm,0,-1)
##def cycle(a,b,c):
##    if c==0:
##        b[:]=a
##        return
##    b[:c]=a[-c:]
##    b[c:]=a[:-c]
##cycle(xy[:,0],xy[:,1],sm)
##q=np.zeros_like(xy)
##
##
###proc = subprocess.Popen(["mpirun" ,"-np", "1", "./p3",str(logs),"y"],stdin=subprocess.PIPE,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
##
###import select
###q = select.poll()
###q.register(proc.stdout,select.POLLIN)
###print(q.poll(0))
##
##
##def get_l(xy,v):
##    proc = subprocess.Popen(["mpirun" ,"-np", "1", "./p3",str(logs),"n"],stdin=subprocess.PIPE,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
##    #p=proc.stderr.read()
##    #print(p)
##    #proc = subprocess.Popen(["./p3",str(logs),"n"],stdin=subprocess.PIPE,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
##    def write(x,y,v):
##        #print(x,y,v)
##        proc.stdin.write(int(x).to_bytes(4,"little",signed=True))
##        proc.stdin.write(int(y).to_bytes(4,"little",signed=True))
##        proc.stdin.write(np.float64(v).data.tobytes())
##        proc.stdin.flush()
###while True:
###                     
###    print(proc.stdout.read(50))
###    input()
##    for (x,y),v in zip(xy,v):
##        write(x,y,v)
##    write(0,0,np.NaN)
##    a=proc.stdout.read(s**2*8)
##    p=proc.stderr.read()
##    #print(p)
##    if len(a)<s**2*8:
##        print(a)
##    dt=np.dtype(np.float64)
##    dt=dt.newbyteorder('<')
##    a=np.frombuffer(a,dt,s**2)
##
##    return a.reshape(s,s).astype(np.float32)
##
###plt.imshow(x_c)
###plt.show()
##import cv2
##src=np.zeros([s,s],np.uint8)
##g=8
##for i in range(0,s,g):
##    for j in range(i%(g*2),s,g*2):
##        src[i:i+g,j:j+g]=255
##
##fig, ax = plt.subplots()
##
##def get_w(t):
##    
##    cycle(xy,q,t)
##    x_c=get_l(xy,q[:,0])
##    y_c=get_l(xy,q[:,1])
##    return cv2.remap(src,x_c,y_c,cv2.INTER_LINEAR,borderMode=cv2.BORDER_REPLICATE)
##img=0
###ax.imshow(src,"gray")
###plt.show()
##def init():
##    global img
##    img=ax.imshow(get_w(0),"gray")
##    return img,
##
##def update(frame):
##    img.set_data(get_w(frame))
##    return img,
##
###plt.imshow(,"gray")
##
##ani = FuncAnimation(fig, update, frames=range(s*4),init_func=init, blit=True,interval=30)
##plt.show()


proc = subprocess.Popen(["mpirun" ,"-np", "1", "./p3",str(logs),"n","20"],stdin=subprocess.PIPE,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
#p=proc.stderr.read()
#print(p)
#proc = subprocess.Popen(["./p3",str(logs),"n"],stdin=subprocess.PIPE,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
def write(x,y,v):
    #print(x,y,v)
    proc.stdin.write(int(x).to_bytes(4,"little",signed=True))
    proc.stdin.write(int(y).to_bytes(4,"little",signed=True))
    proc.stdin.write(np.float64(v).data.tobytes())
    proc.stdin.flush()

for i in range(s):    
    write(i,0,0)
    write(i,s-1,0)
for i in range(1,s-1):    
    write(0,i,0)
    write(s-1,i,0)
tl=[]
bl=[]
from math import *
for y in range(s//8,(7*s)//8):    
    x=s/10*(-cos(8/6*2*pi*(y-s//8)/s)+1)
    print(x,y)
    x1=int(s//4+x+1+y/2)
    write(x1,y,1)
    x2=int(s//4-x+y/2)
    write(x2,y,1)
    tl.append((x1,y))
    bl.append((x2,y))
write(0,0,np.NaN)
a=proc.stdout.read(s**2*8)
p=proc.stderr.read()
print(p)
if len(a)<s**2*8:
    print(a)
dt=np.dtype(np.float64)
dt=dt.newbyteorder('<')
a=np.frombuffer(a,dt,s**2)

a=a.reshape(s,s).astype(np.float32)
plt.imshow(a)
plt.plot(*list(zip(*tl))[::-1],"+k")
plt.plot(*list(zip(*bl))[::-1],"+k")
plt.show()
##p=1
###while True:
##p=proc.stderr.read()                 
##print(p)#input()
##
##
##
##plt.imshow(a)
##plt.show()
##import cProfile as cp
#cp.run("f(7)")
    #print(proc.stdout.read(5))
    #print(q.poll(0))
#write(10,10,-1.)
#write(10,11,-1.)
#write(0,0,np.NAN)

#for i in range(8):
#    value[i] = ord(proc.stdout.read(1))
#    print "value i -> " + str(value[i])

#proc.stdin.write('q')
#proc.wait()
