from matplotlib import pyplot as plt
import numpy as np
import sys
import pylab

input = sys.argv[1]
file = open(input, "r")

t = []
r = []
x_M = []
y_M = []
z_M = []

i = 0

t_r_min = [];   
r_min = [];
x_min_M = [];    
y_min_M = [];  
z_min_M = [];   

k = []

for line in file: 
    line_ = line.split(" ")
    t_i = float(line_[0])                
    r_i = float(line_[1])
    x_i_M = float(line_[2])
    y_i_M = float(line_[3])
    z_i_M = float(line_[4])

    t.append(t_i)
    r.append(r_i)
    x_M.append(x_i_M)
    y_M.append(y_i_M)
    z_M.append(z_i_M)
    i += 1
for j in xrange(len(r)-1):
    if j >= 1:
        if r[j+1] > r[j] and r[j-1] > r[j]:
            t_r_min.append(t[j])
            r_min.append(r[j])
            x_min_M.append(x_M[j])
            y_min_M.append(y_M[j])
            z_min_M.append(z_M[j])

theta = np.zeros(len(r_min))
for i in range(len(r_min)): 
    theta[i] = np.arctan(y_min_M[i]/x_min_M[i])

#plt.plot(t_r_min, theta)
#plt.title('Perihelion angle theta$_p$', fontsize = 22)
#plt.xlabel('Time [s]', fontsize = 22)
#plt.axis([-1.5, 1.5, -1.5, 1.5])
#plt.axis('equal')
#plt.ylabel('Theta$_p$ [rad]', fontsize = 22)# change to v_tilde(x)
#pylab.xticks(fontsize=16)
#pylab.yticks(fontsize=16)
#plt.show()

#plt.plot(t_r_min, theta, 'b')
plt.title('Perihelion angle theta$_p$ with linear regression', fontsize = 22)

plt.xlabel('Time [s]', fontsize = 22)
#plt.axis([-1.5, 1.5, -1.5, 1.5])
#plt.axis('equal')

plt.ylabel('Theta$_p$ [rad]', fontsize = 22)# change to v_tilde(x)
pylab.xticks(fontsize=16)
pylab.yticks(fontsize=16)

fit = np.polyfit(t_r_min, theta,1)
fit_fn = np.poly1d(fit) 
print (fit_fn(t_r_min[-1])- fit_fn(t_r_min[0]))*360*3600/(2*np.pi)
plt.plot(t_r_min, fit_fn(t_r_min)-fit_fn(t_r_min[0]), 'r', t_r_min, theta- fit_fn(t_r_min[0]), 'b')
plt.legend(['Linear regression', 'Perihelion angle'], fontsize = 16, loc=2)
plt.show()

