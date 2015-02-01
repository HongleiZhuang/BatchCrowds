import numpy as np
import matplotlib.pyplot as plt
import math

# f = open("model.txt", "r")
# slines = f.readlines()

# alphasize = int(slines[0].split("\t")[0])
# betasize = int(slines[0].split("\t")[1])

ys = [[0.6705159247532129, 0.16762898118830322, 0.07450176941702365, 0.041907245297075804, 0.026820636990128516, 0.018625442354255912],
	  [0.6581534005952402, 0.1803998458686523, 0.06467879172907545, 0.050469175621349886, 0.026096455180812047, 0.020202331004870143]]
labels = ['Generative parameters', 'Learned parameters']

fig = plt.figure()

o = ['kx--', 'r+-', 'gD-', 'ms-']

for i in xrange(len(ys)):
	# if i < 4:
	# 	continue
	y = ys[i]
	x = xrange(len(y))
	# y = [float(t) for t in slines[i+2+alphasize].strip().split("\t")]
	plt.plot(x, y, o[i], markersize=20, label=labels[i], linewidth=3)

#plt.plot([-100, 100], [0, 0], 'k--', label='_nolegend_')
#plt.tight_layout()
plt.xlim([-0.5, 5.5])
plt.xticks(fontsize=22)
plt.yticks(fontsize=22)
#plt.ylim([0.8, 2.2])
#fig.suptitle('Parameter analysis of beta', fontsize=20)
plt.ylabel(r'$p_\tau$', fontsize=38)
#plt.xticks([1,2,3,4], ['>=1', '>=2', '>=3', '>=4'])
plt.xlabel(r'$\tau$', fontsize=38)
#plt.tight_layout()
plt.legend(loc=1, prop={'size':22})
#plt.plot(x, y, 'bo-')
plt.subplots_adjust(top=0.95, bottom=0.12, left=0.16, right=0.98)
plt.show()
