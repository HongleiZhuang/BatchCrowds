import random

# scores = [10, 10, 10, 1, 1, 1, 1, 1, 1, 1]

pos = 10
tot = 1000
posv = 10
negv = 1

scores = [posv] * pos + [negv] * (tot - pos)


N = 10000
b = 5
f = open("rankedlists.txt", 'w')
for i in xrange(N):
	l = range(len(scores))
	random.shuffle(l)
	selected = l[:b]
	w = [scores[i] for i in selected]
	totw = sum(w)
	rl  = []
	
	for j in xrange(b):
		p = random.random() * totw
		p0 = 0.0
		for k in xrange(b):
			p0 += 0 if selected[k] in rl else w[k]
			if p0 > p:
				rl.append(selected[k])
				totw -= w[k]
				break
	for num in rl:
		f.write(str(num) + "\t")
	f.write("\n")

f.close()