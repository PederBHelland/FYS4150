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
        t = np.zeros(n)
        kineticEnergy = np.zeros(n)
        potentialEnergy = np.zeros(n)
        totalEnergy = np.zeros(n)
        angularMomentum = np.zeros(n)
    else:
        t[i-1] = line_[0]
        kineticEnergy[i-1] = line_[1]
        potentialEnergy[i-1] = line_[2]
        totalEnergy[i-1] = kineticEnergy[i-1] + potentialEnergy[i-1]
        angularMomentum[i-1] = line_[3]
    i += 1

print max(kineticEnergy) - min(kineticEnergy)

plt.plot(t[:-1], kineticEnergy[:-1], 'b', t[:-1], potentialEnergy[:-1], 'r', t[:-1], totalEnergy[:-1], 'g', t[:-1], angularMomentum[:-1], 'y')
plt.title('Energy of the Earth: n=%s' % n, fontsize = 22)
plt.xlabel('Time', fontsize = 22)
plt.ylabel('Energy', fontsize = 22)
pylab.xticks(fontsize=16)
pylab.yticks(fontsize=16)
plt.legend(['Kinetic energy', 'Potential energy', 'Mechanical energy', 'Angular momentum'], fontsize = 16) 
plt.show()


