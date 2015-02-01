#-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.
#
# File Name : avgeval.py
#
# Purpose :
#
# Creation Date : 18-01-2015
#
# Last Modified : Sun Jan 18 17:05:16 2015
#
# Created By : Honglei Zhuang (hzhuang3@illinois.edu) 
#
#_._._._._._._._._._._._._._._._._._._._._.

import sys

if len(sys.argv) < 2:
	print "avgeval.py [filename]"
	sys.exit()

metriclist=["Accuracy", "Precision", "Recall", "F1-score", "Area under ROC"]

slines = open(sys.argv[1], 'r').readlines()
avgmetric = [0 for i in xrange(len(metriclist))]
cnt = 0
for sline in slines:
	sline = sline.strip()
	slist = sline.split("\t")
	if len(slist) < len(metriclist):
		continue
	for i in xrange(len(avgmetric)):
		avgmetric[i] += float(slist[i])
	cnt += 1

for i in xrange(len(avgmetric)):
	print metriclist[i], "=", float(avgmetric[i]) / cnt


