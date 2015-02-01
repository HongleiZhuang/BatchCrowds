import numpy as np
import matplotlib.pyplot as plt
import math

f = open("loglikelihood.txt", "r")
slines = f.readlines()
y = [float(sline)/1000 for sline in slines]
# y = y[1:1]
x = xrange(len(y))


fig = plt.figure()

o = ['yo-', 'bx-', 'r*-', 'gD-', 'ms-']

plt.plot(x, y, 'b-', markersize=10, linewidth=5)

#plt.plot([-100, 100], [0, 0], 'k--', label='_nolegend_')
#plt.tight_layout()
# plt.xlim([-0.5, 5.5])
plt.xticks(fontsize=22)
plt.yticks(fontsize=22)
#plt.ylim([0.8, 2.2])
#fig.suptitle('Parameter analysis of beta', fontsize=20)
plt.ylabel(r'Log-likelihood (10^3)', fontsize=22)
#plt.xticks([1,2,3,4], ['>=1', '>=2', '>=3', '>=4'])
plt.xlabel('Iteration number', fontsize=22)
#plt.tight_layout()
#plt.legend(loc=2)
#plt.plot(x, y, 'bo-')
plt.subplots_adjust(top=0.95, bottom=0.12, left=0.21, right=0.96)
plt.show()
