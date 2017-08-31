import numpy as np 
from matplotlib import pyplot as plt 

n = 100

A = np.zeros((n,n))

for i in range(n):
	for j in range(n):
		if i==j:
			A[i,j] = 2
		elif i==j-1 or i == j+1:
			A[i,j] = -1

h = 1./(n+1)
v = np.zeros(n)
x = np.zeros(n); x[-1] = (n-1)*h

for i in range(1,n-1):
	x[i] = i*h
	v[i] = 1 - (1-np.exp(-10))*x[i]-np.exp(-10*x[i])

b = A.dot(v)
print b

plt.plot(x, b/h**2, 'r', x, 100*np.exp(-10*x), 'b')
plt.show()
#f = np.zeros(n)

#for i in range(1,n):
#	f[i-1] = -(v[i+1] + v[i-1]-2*v[i])/h**2
