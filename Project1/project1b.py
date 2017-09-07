import numpy as np 
from matplotlib import pyplot as plt 
import sys

#n = int(sys.argv[0])
n = 100

a = np.zeros(n)-1.0
b = np.zeros(n)+2.0
c = np.zeros(n)-1.0

h = 1./(n+1)
v = np.zeros(n)
x = np.zeros(n)
f = np.zeros(n)
u = np.zeros(n)
b_ = np.zeros(n)
temp = np.zeros(n)
#a = np.array([-1,-1,-1,-1])
#b = np.array([2,2,2,2])
#c = np.array([-1,-1,-1,-1])

A = np.zeros((n,n))
# Boundaries
A[0,0] = b[0]
A[0,1] = c[0]
A[-1,-2] = a[-2]
A[-1,-1] = b[-1]

# Loop
for i in range(1, n-1):
    A[i, i-1] = a[i-1]
    A[i,i] = b[i]
    A[i,i+1] = c[i]

for i in range(n):
    x[i] = i*h
    f[i] = 100*np.exp(-10*x[i])
    u[i] = 1 - (1-np.exp(-10))*x[i]-np.exp(-10*x[i])
    b_[i]= h**2 * f[i]

btemp = b[0]
v[0] = b_[0]/btemp

for i in range(1, n): 
    temp[i] = (c[i-1])/btemp
    btemp = b[i] - a[i]*temp[i]
    print btemp
    v[i] = (b_[i] - a[i]*v[i-1])/btemp

for i in range(n-1, 0):
    v[i] = v[i] - temp[i+1]*v[i+1]

plt.plot(x, v, 'r', x, u, 'b')
plt.show()

