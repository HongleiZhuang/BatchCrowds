import sys
import random

if len(sys.argv) < 4:
	print "sample training data set."
	print "sample_train.py [train.srk] [train.gt] [num_batches]"
	sys.exit()

srklines = open(sys.argv[1], 'r').readlines()
gtlines  = open(sys.argv[2], 'r').readlines()
m = int(sys.argv[3])

gtdict = {}
for gtline in gtlines:
	slist = gtline.strip().split("\t")
	gtdict[long(slist[0])] = int(slist[1])

srklist = []
for srkline in srklines:
	srkline = srkline.strip()
	if srkline == "":
		continue
	srklist.append(srkline)
random.shuffle(srklist)

srklist = srklist[:m]
idset = set()
for srk in srklist:
	idlist = srk.replace(';',' ').replace(',',' ').strip().split(" ")
	for idstr in idlist:
		idset.add(long(idstr))

srkoutput = open(sys.argv[1] + ".sample", 'w')
gtoutput  = open(sys.argv[2] + ".sample", 'w')
for srk in srklist:
	srkoutput.write(srk + "\n")
for gtid in idset:
	gtoutput.write(str(gtid) + "\t" + str(gtdict[gtid]) + "\n")
srkoutput.close()
gtoutput.close()
