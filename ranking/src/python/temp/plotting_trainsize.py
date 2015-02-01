import numpy as np
import matplotlib.pyplot as plt
import math

# f = open("model.txt", "r")
# slines = f.readlines()

# alphasize = int(slines[0].split("\t")[0])
# betasize = int(slines[0].split("\t")[1])

ys = [[0.876651305684,0.955760368664,0.959907834101,0.958371735791,0.960061443932,0.962365591398,0.962519201229,0.963748079877,0.964055299539,0.963133640553],[0.962826420891,0.96420890937,0.965437788018,0.964516129032,0.965437788018,0.966052227343,0.965898617512,0.966359447005,0.96712749616,0.966820276498]]
ysf1 = [[0.753356684902,0.871688038735,0.886176171256,0.883663245473,0.883232083984,0.889143806103,0.89003494868,0.892813338666,0.893356454008,0.890064164381],[0.889253987535,0.892976672534,0.897460570728,0.893711144788,0.896946112804,0.898639048815,0.898228378932,0.899750159926,0.902185890868,0.901079412956]]
labels = ['MV', 'MVT', 'PL', 'BAM']
o = ['kx-', 'r+-', 'gD-', 'ms-']
x = [100*(i+1) for i in xrange(10)]

labels=labels[1:4:2]
o=o[1:4:2]

fig = plt.figure()
for i in xrange(len(ys)):
	# if i < 4:
	# 	continue
	y = ys[i]
	# y = [float(t) for t in slines[i+2+alphasize].strip().split("\t")]
	plt.plot(x, y, o[i], markersize=20, label=labels[i], linewidth=3)

#plt.plot([-100, 100], [0, 0], 'k--', label='_nolegend_')
#plt.tight_layout()
#plt.xlim([0, 22])
# plt.xscale('log')
plt.xticks(fontsize=22)
plt.yticks(fontsize=22)
#plt.ylim([0.8, 2.2])
#fig.suptitle('Parameter analysis of beta', fontsize=20)
plt.ylabel('Accuracy', fontsize=31)
#plt.xticks(x, ['2', '5', '10', '20'])
plt.xlabel(r'$m_L$', fontsize=38)
#plt.tight_layout()
plt.legend(loc=4, prop={'size':22})
#plt.plot(x, y, 'bo-')
plt.subplots_adjust(top=0.96, bottom=0.16, left=0.16, right=0.96)

fig = plt.figure()
for i in xrange(len(ys)):
	# if i < 4:
	# 	continue
	y = ysf1[i]
	# y = [float(t) for t in slines[i+2+alphasize].strip().split("\t")]
	plt.plot(x, y, o[i], markersize=20, label=labels[i], linewidth=3)

#plt.plot([-100, 100], [0, 0], 'k--', label='_nolegend_')
#plt.tight_layout()
#plt.xlim([0, 22])
# plt.xscale('log')
plt.xticks(fontsize=22)
plt.yticks(fontsize=22)
#plt.ylim([0.8, 2.2])
#fig.suptitle('Parameter analysis of beta', fontsize=20)
plt.ylabel(r'$F_1$'+"-score", fontsize=31)
#plt.xticks(x, ['2', '5', '10', '20'])
plt.xlabel(r'$m_L$', fontsize=38)
#plt.tight_layout()
plt.legend(loc=4, prop={'size':22})
#plt.plot(x, y, 'bo-')
plt.subplots_adjust(top=0.96, bottom=0.16, left=0.16, right=0.96)
plt.show()


plt.show()
