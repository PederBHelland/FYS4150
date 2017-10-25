from matplotlib import pyplot as plt
import numpy as np
import sys
import pylab

input = sys.argv[1]
file = open(input, "r")

numberOfIterations  = int(file.readline())
dim = int(file.readline()) 
numberOfPlanets = int(file.readline())

print numberOfPlanets

X = np.zeros(((numberOfIterations+1), numberOfPlanets))
Y = np.zeros(((numberOfIterations+1), numberOfPlanets))
Z = np.zeros(((numberOfIterations+1), numberOfPlanets))
"""
X = np.zeros((numberOfIterations+1, numberOfPlanets))
Y = np.zeros((numberOfIterations+1, numberOfPlanets))
Z = np.zeros((numberOfIterations+1, numberOfPlanets))
"""
t = 0
p = 0

for line in file: 
    line_ = line.split()
    x_i = float(line_[0])                
    y_i = float(line_[1])
    z_i = float(line_[2])
    
    X[t,p] = x_i
    Y[t,p] = y_i
    Z[t,p] = z_i 

    p += 1
    
    if(p == numberOfPlanets):
        t += 1
        p = 0

#for i in range(0,numberOfPlanets):
#	X[:,i] -= X[:,0]
#	Y[:,i] -= Y[:,0] 

print X[:,0], Y[:,0]

fig, ax = plt.subplots()


#plt.plot(X[:,0], Y[:,0], 'o', color='gold')
"""
plt.plot(X[:,6], Y[:,6], 'grey')
plt.plot(X[:,4], Y[:,4], 'brown')
"""
plt.plot(X[:-1,1], Y[:-1,1], 'g')
"""
plt.plot(X[:,3], Y[:,3], 'firebrick')
"""
#plt.plot(X[:,2], Y[:,2], 'gold')
"""
plt.plot(X[:,5], Y[:,5])
plt.plot(X[:,7], Y[:,7], 'coral')
plt.plot(X[:,8], Y[:,8], 'dodgerblue')
plt.plot(X[:,9], Y[:,9], 'chocolate')
"""


#plt.plot(X[:,0],t)


sun = plt.Circle((X[0,0],Y[0,0]), 0.005, color='black')
#orbit = plt.Circle((0,0), 1, color='k', fill=False)
plt.title('Position of Earth: n=%s' % numberOfIterations, fontsize = 22)
plt.xlabel('x [AU]', fontsize = 22)
#plt.axis([-1.5, 1.5, -1.5, 1.5])
plt.axis('equal')
plt.ylabel('y [AU]', fontsize = 22)# change to v_tilde(x)
pylab.xticks(fontsize=16)
pylab.yticks(fontsize=16)
plt.legend(['Earth'])
#plt.legend(['Sun', 'Mercury', 'Venus', 'Earth', 'Mars', 'Jupiter', 'Saturn', 'Uranus', 'Neptun', 'Pluto'], fontsize = 16) # Change to Solve solution v_tilde(x) in 1e
ax.add_artist(sun)
#ax.add_artist(orbit)
plt.show()


