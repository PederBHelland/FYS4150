from matplotlib import pyplot as plt
import numpy as np
import sys
import pylab

input = sys.argv[1]
file = open(input, "r")

#numberOfIterations  = 50000
#dim = 2

#r = np.zeros(numberOfIterations)
#t = np.zeros(numberOfIterations)
t = []
r = []
x_sun = []
y_sun = []
z_sun = []
x_M = []
y_M = []
z_M = []

i = 0

t_r_min = [];   t_r_min.append(0.)
r_min = [];     r_min.append(0.39)
x_min_sun = [];     x_min_sun.append(0.)
y_min_sun = [];     y_min_sun.append(0.)
z_min_sun = [];     z_min_sun.append(0.)
x_min_M = [];     x_min_M.append(0.39)
y_min_M = [];     y_min_M.append(0.)
z_min_M = [];     z_min_M.append(0.)

k = []
for line in file: 
    line_ = line.split()
    t_i = float(line_[0])                
    r_i = float(line_[1])
    x_i_sun = float(line_[2])
    y_i_sun = float(line_[3])
    z_i_sun = float(line_[4])
    x_i_M = float(line_[5])
    y_i_M = float(line_[6])
    z_i_M = float(line_[7])

#    t[i] = t_i
#    r[i] = r_i
#    x[i] = x_i
#    y[i] = y_i
#    z[i] = z_i

    t.append(t_i)
    r.append(r_i)
    x_sun.append(x_i_sun)
    y_sun.append(y_i_sun)
    z_sun.append(z_i_sun)
    x_M.append(x_i_M)
    y_M.append(y_i_M)
    z_M.append(z_i_M)
    i += 1
for j in range(len(r)-1):
    if j >= 1:
        if r[j+1] > r[j] and r[j-1] > r[j]:
            t_r_min.append(t[j])
            r_min.append(r[j])
            x_min_sun.append(x_sun[j])
            y_min_sun.append(y_sun[j])
            z_min_sun.append(z_sun[j])
            x_min_M.append(x_M[j])
            y_min_M.append(y_M[j])
            z_min_M.append(z_M[j])
#print r
#print len(t_r_min)
#print 'r_min: ', r_min
#print 't_r_min: ', t_r_min

#plt.plot(t[:-1], r[:-1])
#plt.plot(x_sun[:-1], y_sun[:-1])
plt.plot(x_M[:-1], y_M[:-1], 'r', x_sun[:-1], y_sun[:-1], 'y')
plt.title('Position')
plt.show()

#Finding theta

theta = np.zeros(len(t_r_min))
for i in range(len(t_r_min)): 
    t_orbit = 0.240846
    theta[i] = -(t_r_min[i]-t_orbit*(i))/t_orbit
#print theta
"""

theta = np.zeros(len(r_min))
for i in range(len(r_min)): 
    theta[i] = np.arctan(y_min_M[i]/x_min_M[i])
    #print x_min[i], y_min[i]
#print theta
"""

plt.plot(t_r_min, theta)
plt.title('Theta')
plt.show()

plt.plot(t, r, 'r')
plt.title('r')
plt.show()


#x = np.linspace(-10,10,100)
#print x
#y = np.arctan(x)
#plt.plot(x, y)
#plt.show()

