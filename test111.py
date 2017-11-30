import pandas as pd
import collections
from collections import defaultdict

whole_promoter = open('C:\data\SOP\gtggttag - promoter\/newwhole.txt','r').read()
promoter_dict = dict()
seq1 = 'CCGAC'
#seq2 = 'GGTAAT'
length1 = len(seq1)
#length2 = len(seq2)

count = 0

positions = []

inp = whole_promoter

for i in range(len(inp)):
	
	k1 = inp[i:i+length1]
	#k2 = inp[i:i+length2]

	if k1 == seq1:
		count = count + 1
		positions.append([seq1,i])

	"""
	if k2 == seq2:
		count = count + 1
		positions.append([seq2,i])
	"""
#print(positions)
print(len(positions))

spacers = dict()

for i in range(26):
	spacers[i] = []

for i in range(len(positions) - 1):
	spacer = whole_promoter[positions[i][1] + 5:positions[i+1][1]]
	k = positions[i+1][1] - positions[i][1] - 5
	if(k > 0 and k < 26):
		#spacers.append([k,spacer])
		spacers[k].append(spacer)

#spacers1 = dict((x,spacers.count(x)) for x in set(spacers))

spacer_count_results = []

for i in range(1,26):
	string = ''
	for j in range(len(spacers[i])):
		string = string + spacers[i][j]
	spacer_count_results.append([i,collections.Counter(string)])

print(spacer_count_results)