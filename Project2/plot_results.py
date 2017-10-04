from matplotlib import pyplot as plt
import numpy as np
import sys
import pylab
import sys

def input_file(input):
	file = open(input, "r")
	return file

#input = sys.argv[1]
#file = open(input, "r")

def probability_distribution(file):
	i = 0
	for line in file: 
		line_ = line.split()
		if i==0:
			n = int(line_[0])
			if line_[1] == n:
				w_r = int(line_[1])
			else:
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
	return u1, u2, u3, rho, w_r, n

def plot_figures(file):
	u1, u2, u3, rho, w_r, n = probability_distribution(file)
	if w_r == n:
		plt.plot(rho, u1, 'r')
		plt.legend(['Ground state'], fontsize = 16)
		plt.title('Energy levels without repulsion', fontsize = 20)
	else:
		plt.plot(rho, u1, 'r', rho, u2, 'b', rho, u3, 'g')
		plt.legend(['Ground state', 'First excited state', 'Second excited state'], fontsize = 16)
		plt.title('Energy levels, $w_r$ = %.2f' %w_r, fontsize = 20)
	plt.xlabel('rho', fontsize = 18)
	plt.ylabel('The probability distribution', fontsize = 18)
	pylab.xticks(fontsize=16)
	pylab.yticks(fontsize=16) 
	plt.show()

list_of_wr = []
def plot_comparing(file):
	u1, u2, u3, rho, w_r, n = probability_distribution(file)
	plt.plot(rho, u1)
	if w_r == 200:	
		list_of_wr.append("Without repulsion")
	else:
		list_of_wr.append(w_r)
	plt.hold('on')

input1 = 'results_w0.010000.txt'
input2 = 'results_w0.500000.txt'
input3 = 'results_w1.000000.txt'
input4 = 'results_w5.000000.txt'
input5 = 'results_without_repulsion.txt'

# Plotting energy levels for w_r = 0.01, 0.5, 1.0, 5.0 and without repulsion 
plot_figures(input_file(input1))
plot_figures(input_file(input2))
plot_figures(input_file(input3))
plot_figures(input_file(input4))
plot_figures(input_file(input5))

# plotting the ground state for w_r = 0.01, 0.5, 1.0, 5.0 and without repulsion in the same plot
plot_comparing(input_file(input1))
plot_comparing(input_file(input2))
plot_comparing(input_file(input3))
plot_comparing(input_file(input4))
plot_comparing(input_file(input5))
plt.legend(['$w_r$ = %.2f'%list_of_wr[0], '$w_r$ = %.2f'%list_of_wr[1], '$w_r$ = %.2f'%list_of_wr[2], '$w_r$ = %.2f'%list_of_wr[3], list_of_wr[4]], fontsize = 16)
plt.title('Ground state for different $w_r$', fontsize = 20)
plt.xlabel('rho', fontsize = 18)
plt.ylabel('The probability distribution', fontsize = 18)
pylab.xticks(fontsize=16)
pylab.yticks(fontsize=16)
plt.axis([0, 45, 0, 0.025])
plt.show()

# plotting the ground state for w_r = 0.5, 1.0, 5.0 and without repulsion in the same plot
plot_comparing(input_file(input2))
plot_comparing(input_file(input3))
plot_comparing(input_file(input4))
plot_comparing(input_file(input5))
plt.legend(['$w_r$ = %.2f'%list_of_wr[1], '$w_r$ = %.2f'%list_of_wr[2], '$w_r$ = %.2f'%list_of_wr[3], list_of_wr[4]], fontsize = 16)
plt.title('Ground state for different $w_r$', fontsize = 20)
plt.xlabel('rho', fontsize = 18)
plt.ylabel('The probability distribution', fontsize = 18)
pylab.xticks(fontsize=16)
pylab.yticks(fontsize=16)
plt.axis([0, 5, 0, 0.025])
plt.show()

