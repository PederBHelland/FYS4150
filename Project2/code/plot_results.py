from matplotlib import pyplot as plt
import numpy as np
import sys
import pylab

input = sys.argv[1]
file = open(input, "r")

i = 0
for line in file: 
    line_ = line.split()
    if i==0:
        n=int(line_[0])
        r = np.zeros(n)
        omega = np.zeros(n)
    else:
        r[i-1] = line_[0]
        omega[i-1] = line_[1]
    i += 1

plt.plot(r, omega)
plt.title('n=%s' % n, fontsize = 22)
plt.xlabel('$\omega_r$', fontsize = 22)
plt.ylabel('The wave function', fontsize = 22)
pylab.xticks(fontsize=16)
pylab.yticks(fontsize=16)
plt.legend(['Different values of $\omega_r'], fontsize = 16) 
plt.show()


