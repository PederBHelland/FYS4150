import numpy as np 
import scipy

n = 100
a = np.zeros(n)-1.0
b = np.zeros(n)+2.0
c = np.zeros(n)-1.0
A = np.zeros((n,n))
h = 1./(n+1)
v = np.zeros(n)
x = np.zeros(n)
f = np.zeros(n)
u = np.zeros(n)
b_ = np.zeros(n)

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

x = np.linalg.solve(A,b_)      # solve the equation
P, L, U = scipy.linalg.lu(A)   # LU decomposition
print (A-P*L*U)                # check if zero
