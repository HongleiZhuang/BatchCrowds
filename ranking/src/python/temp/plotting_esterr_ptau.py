import numpy as np
import matplotlib.pyplot as plt
import math

def readmodel(fstr):
	fmodel = open(fstr, 'r').readlines()
	l = float(fmodel[0])
	ptau = [float(fmodel[i + 2]) for i in xrange(int(fmodel[1]))]
	return l, ptau

x = []
y = []
for rho in [1, 1.5, 2, 2.5, 3]:
	true_lambda, true_ptau       = readmodel("model_0.5_" + str(rho) + ".txt")
	learned_lambda, learned_ptau = readmodel("learnedmodel_0.5_" + str(rho) + ".txt")
	x.append(rho)
	y.append(sum([(learned_ptau[i] - true_ptau[i]) ** 2 for i in xrange(len(true_ptau))]))


fig = plt.figure()

o = ['yo-', 'bx-', 'r*-', 'gD-', 'ms-']

plt.plot(x, y, 'bx-', markersize=15, linewidth=2)
# plt.plot([0, 1], [0, 0], 'k--', linewidth = 2)
#plt.plot([-100, 100], [0, 0], 'k--', label='_nolegend_')
#plt.tight_layout()

plt.xticks(fontsize=22)
plt.yticks(fontsize=22)
plt.ylim([0, 0.01])
plt.xlim([0.75, 3.25])
#fig.suptitle('Parameter analysis of beta', fontsize=20)
plt.ylabel(r'$\|\|\hat{p}_{\tau} - p_{\tau}\|\|^2$', fontsize=38)
#plt.xticks([1,2,3,4], ['>=1', '>=2', '>=3', '>=4'])
plt.xlabel(r'$\rho$', fontsize=38)
#plt.tight_layout()
#plt.legend(loc=2)
#plt.plot(x, y, 'bo-')
plt.subplots_adjust(top=0.97, bottom=0.14, left=0.22, right=0.96)
plt.show()
