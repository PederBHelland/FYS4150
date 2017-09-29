from matplotlib import pyplot as plt
import numpy as np
import sys
import pylab
import sys

input = sys.argv[1]
file = open(input, "r")
"""
A = np.loadtxt(input)

plt.figure()
for i in xrange(A.shape[0]) :
    
    plt.plot(np.insert(A[:,i],[0, A.shape[0]], [0, 0]),label='i=%d'%i)
    plt.hold('on')
plt.legend()
plt.show()
sys.exit(1)
"""
i = 0
for line in file: 
    line_ = line.split()
    if i==0:
        n=int(line_[0])
        rho = np.zeros(n)
        u1 = np.zeros(n)
        u2 = np.zeros(n)
        u3 = np.zeros(n)
    elif 0 < i <= n:
        rho[i-1] = float(line_[0])
        u1[i-1] = float(line_[1])**2   
    elif n < i <= 2*n:        
        u2[i-n-1] = float(line_[1])**2
    elif 2*n < i <= 3*n:
        u3[i-2*n-1] = float(line_[1])**2
    i += 1

print u2

plt.plot(rho, u1, 'r', rho, u2, 'b', rho, u3, 'g')
plt.title('n=%s' % n, fontsize = 22)
plt.xlabel('$\omega_r$', fontsize = 22)
plt.ylabel('The wave function', fontsize = 22)
pylab.xticks(fontsize=16)
pylab.yticks(fontsize=16)
plt.legend(['Different values of $\omega_r$'], fontsize = 16) 
plt.show()


