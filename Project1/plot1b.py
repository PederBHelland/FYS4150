from matplotlib import pyplot as plt
import numpy as np
import sys

input = sys.argv[1]
file = open(input, "r")

i = 0
for line in file: 
    line_ = line.split()
    if i==0:
        n=int(line_[0])
        x = np.zeros(n)
        u = np.zeros(n)
        v = np.zeros(n)
    else:
        x[i-1] = line_[0]
        u[i-1] = line_[1]
        v[i-1] = line_[2]
    i += 1

plt.plot(x, u, 'r', x, v, 'b')
plt.title('n=%s' % n)
plt.xlabel('x')
plt.ylabel('u(x),v(x)')
plt.legend(['Exact solution u(x)', 'Our solution v(x)'])
plt.show()


