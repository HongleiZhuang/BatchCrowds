
# Worker model
import sys
import random

if len(sys.argv) < 9:
	print "Synthetic data generation."
	print "python datagen_beta.py [lambda] [rho] [model_output_file] [n] [m] [rkdata_output_file] [gtlabel_output_file] [gtscore_output_file]"
	sys.exit()

ind_prob = float(sys.argv[1])
rho = float(sys.argv[2])
modelf_name = sys.argv[3]
tot = int(sys.argv[4])
N   = int(sys.argv[5])
rkdata_name = sys.argv[6]
gtlabel_name = sys.argv[7]
gtscore_name = sys.argv[8]

bsize = 5
offset = 1

nullp = [(k + offset) ** (-rho) for k in xrange(bsize + 1)]
null_prob = [x / sum(nullp) for x in nullp]
print null_prob

modelf = open(modelf_name , 'w')
modelf.write(str(ind_prob) + "\n")
modelf.write(str(bsize + 1) + "\n")
for ptau in null_prob:
	modelf.write(str(ptau) + "\n")
modelf.close()

alpha = 2
beta  = 4
scores = [random.betavariate(alpha, beta) for i in xrange(tot)]
label = [1 if x > 0.5 else 0 for x in scores]

gtlabel = open(gtlabel_name, 'w')
gtscore = open(gtscore_name, "w")
itemlist = [(i, scores[i]) for i in xrange(len(scores))]
itemlist = sorted(itemlist, key = lambda x: - x[1])
for item in itemlist:
	gtscore.write(str(item[0]) + "\t" + str(item[1]) + "\n")
	gtlabel.write(str(item[0]) + "\t" + str(1 if item[1] > 0.5 else 0) + "\n")
gtscore.close()
gtlabel.close()




#print scores
print "Positive: ", sum(label)

f = open(rkdata_name, 'w')

for i in xrange(N):
	l = range(len(scores))
	random.shuffle(l)
	selected = l[:bsize]
	w = [scores[i] for i in selected]
	totw = sum(w)
	rl  = []
	
	for j in xrange(bsize):
		p = random.random() * totw
		p0 = 0.0
		for k in xrange(bsize):
			p0 += 0 if selected[k] in rl else w[k]
			if p0 > p:
				rl.append(selected[k])
				totw -= w[k]
				break

	independent = True if random.random() <= ind_prob else False
	if not independent:
		#print "guess"
		p = random.random()
		p0 = 0.0
		for k in xrange(bsize + 1):
			p0 += null_prob[k]
			if p0 > p:
				#print "k=", k
				rlt = rl[:k]
				rlf = rl[k:]
				break
	else:
		rlt = []
		rlf = []
		for i in rl:
			if random.random() < scores[i]:
				rlt.append(i)
			else:
				rlf.append(i)
	f.write(",".join([str(x) for x in rlt]) + ";" + ",".join([str(x) for x in rlf]) + "\n")
f.close()