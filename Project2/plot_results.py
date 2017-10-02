from matplotlib import pyplot as plt
import numpy as np
import sys
import pylab
import sys

input = sys.argv[1]
file = open(input, "r")

i = 0
for line in file: 
    line_ = line.split()
    if i==0:
        n=int(line_[0])
	w_r = float(line_[1])
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

plt.plot(rho, u1, 'r', rho, u2, 'b', rho, u3, 'g')
plt.title('Energy levels for the three first eigenvalues, $w_r$ = %.2f' %w_r , fontsize = 20)
plt.xlabel('rho', fontsize = 18)
plt.ylabel('The probability distribution', fontsize = 18)
pylab.xticks(fontsize=16)
pylab.yticks(fontsize=16)
plt.legend(['Ground state', 'First excited state', 'Second excited state'], fontsize = 16) 
plt.show()


