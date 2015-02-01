import numpy as np
import matplotlib.pyplot as plt
import math

# f = open("model.txt", "r")
# slines = f.readlines()

# alphasize = int(slines[0].split("\t")[0])
# betasize = int(slines[0].split("\t")[1])

ys = [[0.83202,0.83736,0.83096,0.82718],[0.68668,0.8184,0.88362,0.91746],[0.71172,0.78252,0.82394,0.83534],[0.7762,0.85498,0.89774,0.92968]]
ysf1 = [[0.323124787962,0.249819789735,0.173589321553,0.131102479131],[0.491333199938,0.620824544374,0.713606601703,0.78846834819],[0.508444931138,0.598354260472,0.662347127105,0.689109178634],[0.54159292897,0.657618308644,0.739020784969,0.806651864473]]
labels = ['MV', 'MVT', 'PL', 'BAM']
o = ['kx-', 'r+-', 'gD-', 'ms-']
x = [2,5,10,20]

fig = plt.figure()
for i in xrange(len(ys)):
	# if i < 4:
	# 	continue
	y = ys[i]
	# y = [float(t) for t in slines[i+2+alphasize].strip().split("\t")]
	plt.plot(x, y, o[i], markersize=20, label=labels[i], linewidth=3)

#plt.plot([-100, 100], [0, 0], 'k--', label='_nolegend_')
#plt.tight_layout()
plt.xlim([0, 22])
plt.xticks(fontsize=22)
plt.yticks(fontsize=22)
#plt.ylim([0.8, 2.2])
#fig.suptitle('Parameter analysis of beta', fontsize=20)
plt.ylabel('Accuracy', fontsize=31)
#plt.xticks(x, ['2', '5', '10', '20'])
plt.xlabel(r'$m/n$', fontsize=38)
#plt.tight_layout()
plt.legend(loc=4, prop={'size':22})
#plt.plot(x, y, 'bo-')
plt.subplots_adjust(top=0.98, bottom=0.15, left=0.16, right=0.98)

fig = plt.figure()
for i in xrange(len(ys)):
	# if i < 4:
	# 	continue
	y = ysf1[i]
	# y = [float(t) for t in slines[i+2+alphasize].strip().split("\t")]
	plt.plot(x, y, o[i], markersize=20, label=labels[i], linewidth=3)

#plt.plot([-100, 100], [0, 0], 'k--', label='_nolegend_')
#plt.tight_layout()
plt.xlim([0, 22])
plt.xticks(fontsize=22)
plt.yticks(fontsize=22)
#plt.ylim([0.8, 2.2])
#fig.suptitle('Parameter analysis of beta', fontsize=20)
plt.ylabel(r'$F_1$'+"-score", fontsize=31)
#plt.xticks(x, ['2', '5', '10', '20'])
plt.xlabel(r'$m/n$', fontsize=38)
#plt.tight_layout()
plt.legend(loc=4, prop={'size':22})
#plt.plot(x, y, 'bo-')
plt.subplots_adjust(top=0.98, bottom=0.15, left=0.16, right=0.98)
plt.show()


plt.show()
