from collections import defaultdict
import collections

inp = open('C:\data\SOP\gtggttag - promoter\PromoterSeq.txt','r').read()
out = open('C:\data\SOP\gt_result.txt','w')

seq1 = 'gtggttag'
seq2 = 'ggtaat'
length1 = len(seq1)
length2 = len(seq2)
count = 0
positions = []
distance = []
distance_consider = []
spacers = []
a = defaultdict(list)

for i in range(len(inp)):
	
	k1 = inp[i:i+length1]
	k2 = inp[i:i+length2]
	
	if k1 == seq1:
		count = count + 1
		positions.append(i)
		
	if k2 == seq2:
		count = count + 1
		positions.append(i)

for j in range(len(positions)-1):
	if positions[j+1] - positions[j] < 26:
		distance.append(positions[j+1] - positions[j])
		k = inp[positions[j]:positions[j+1]]
		if k.startswith('gtggttag'):
			spacer = k[8:]
			if len(spacer) > 0:
				spacers.append(spacer)
		if k.startswith('ggtaat'):
			spacer = k[6:]
			if len(spacer) > 0:
				spacers.append(spacer)

for i in range(len(distance)):
	if distance[i] < 26:
		distance_consider.append(distance[i])

for i in range(1,26):
	a1 = [x for x in spacers if len(x) == i]
	a[i].append(a1)

spacers.sort(key = lambda s: len(s))

for i in a.keys():
	a_consider = a[i]
	a_consider1 = a_consider[0]
	a_consider1 = ''.join(str(elem) for elem in a_consider1)
	k = collections.Counter(a_consider1)
	out.write(str(i) + "\t" + str(k) + "\n")

print(spacers)
print(a)
out.close()