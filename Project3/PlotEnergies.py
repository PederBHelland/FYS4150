from matplotlib import pyplot as plt
import numpy as np
import sys
import pylab

input = sys.argv[1]
file = open(input, "r")

numberOfIterations  = int(file.readline())
dim = int(file.readline()) 
numberOfPlanets = int(file.readline())
"""
t = np.zeros((numberOfIterations/100+1, numberOfPlanets))
kineticEnergy = np.zeros((numberOfIterations/100+1, numberOfPlanets))
potentialEnergy = np.zeros((numberOfIterations/100+1, numberOfPlanets))
totalEnergy = np.zeros((numberOfIterations/100+1, numberOfPlanets))
angularMomentum = np.zeros((numberOfIterations/100+1, numberOfPlanets))
"""
t = np.zeros((numberOfIterations/100, numberOfPlanets))
kineticEnergy = np.zeros((numberOfIterations/100, numberOfPlanets))
potentialEnergy = np.zeros((numberOfIterations/100, numberOfPlanets))
totalEnergy = np.zeros((numberOfIterations/100, numberOfPlanets))
angularMomentum = np.zeros((numberOfIterations/100, numberOfPlanets))



p = 0
i = 0
for line in file: 
    line_ = line.split()

    t[i, p] = line_[0]
    kineticEnergy[i, p] = line_[1]
    potentialEnergy[i, p] = line_[2]
    totalEnergy[i, p] = kineticEnergy[i, p] + potentialEnergy[i, p]
    angularMomentum[i,p] = line_[3]

    p += 1
    if(p == numberOfPlanets):
        i += 1
        p = 0

#amplitude_kinetic = max(kineticEnergy[:,1]) - min(kineticEnergy[:,1])
#print "jasfgk = %e" %amplitude_kinetic
#print max(kineticEnergy[:,1]),  min(kineticEnergy[:,1])

#plt.plot(t[:,1], kineticEnergy[:,1], t[:,1], potentialEnergy[:,1], t[:,1], totalEnergy[:,1])
plt.plot(t[:, 1], kineticEnergy[:,1], 'b', t[:,1], potentialEnergy[:,1], 'r', t[:,1], totalEnergy[:,1], 'g', t[:,1], angularMomentum[:,1], 'y')
#plt.plot(t[:, 0], kineticEnergy[:,0], 'b', t[:,0], potentialEnergy[:,0], 'r', t[:,0], totalEnergy[:,0], 'g', t[:,0], angularMomentum[:,0], 'y')
#plt.title('Energy of the Earth: n=%s' % n, fontsize = 22)
plt.xlabel('Time', fontsize = 22)
plt.ylabel('Energy', fontsize = 22)
pylab.xticks(fontsize=16)
pylab.yticks(fontsize=16)
plt.legend(['Kinetic energy', 'Potential energy', 'Mechanical energy', 'Angular momentum'], fontsize = 16) 
plt.show()


