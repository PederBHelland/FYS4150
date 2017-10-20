from matplotlib import pyplot as plt
import numpy as np
import sys
import pylab

input = sys.argv[1]
file = open(input, "r")

numberOfIterations  = 5000
dim = 2


r = np.zeros(numberOfIterations)
t = np.zeros(numberOfIterations)


i = 0

r_min = []
t_r_min = []

for line in file: 
    line_ = line.split()
    t_i = float(line_[0])                
    r_i = float(line_[1])
    
    t[i] = t_i
    r[i] = r_i
    if r[i+1] > r[i] and r[i-1] > r[i]:
        r_min.append(r[i])
        t_r_min.append(t[i])
    i += 1
'''
for i in range(len(r)):
    a = min(r)
    print r[i]
    if r[i] == a:
        t_r_min.append(t[i])
        print t[i]
'''
print 'r_min: ', r_min
print 't_r_min: ', t_r_min

plt.plot(t[:-1], r[:-1])
plt.axis([0,1,0.3895,0.3905])
plt.show()
#r_min = min(r)
#print r_min
'''

index = []
for i in range(len(r)):
    if r[i] <= 0.39:
        index.append(i)
#        print r[i]
#print index
#print min(r)

for j in range(len(index)):
    print t[j]
'''
