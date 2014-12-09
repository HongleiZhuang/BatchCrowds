
# Worker model

import random

tot = 100
bsize = 5

offset = 1
rho = 1
nullp = [(k + offset) ** (-rho) for k in xrange(bsize + 1)]
null_prob = [x / sum(nullp) for x in nullp]
print null_prob

alpha = 2
beta  = 6
scores = [random.betavariate(alpha, beta) for i in xrange(tot)]
label = [1 if x > 0.5 else 0 for x in scores]

gt = open("ground_truth_i.txt", 'w')
for i in xrange(len(label)):
	if label[i] == 1:
		gt.write(str(i) + "\n")
gt.close()

gtscore = open("ground_truth_score_i.txt", "w")
itemlist = [(i, scores[i]) for i in xrange(len(scores))]
itemlist = sorted(itemlist, key = lambda x: - x[1])
for item in itemlist:
	gtscore.write(str(item[0]) + "\t" + str(item[1]) + "\n")
gtscore.close()


N = 100000
guess_prob = 1.0

#print scores
print "Positive: ", sum(label)

f = open("rankedlists_beta_i.txt", 'w')

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

	guess = True if random.random() <= guess_prob else False
	if guess:
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
			if label[i] == 1:
				rlt.append(i)
			else:
				rlf.append(i)
	f.write(",".join([str(x) for x in rlt]) + ";" + ",".join([str(x) for x in rlf]) + "\n")
f.close()