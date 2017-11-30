"""
Analyzing the distribution of GTGGTTAG and GGTAAT in the promoters of various upregulated and downregulated genes obtained from microarray analysis in case of
heat,osmotic and salinity stress
"""

import pandas as pd
import collections
from collections import defaultdict

def union(a, b):
    """ return the union of two lists """
    return list(set(a) | set(b))

#promoter file
promoter = open('C:\data\SOP\gtggttag - promoter\PromoterSeq.txt','r').read()
promoter_dict = dict()
seq1 = 'gtggttag'
seq2 = 'ggtaat'
length1 = len(seq1)
length2 = len(seq2)
spacers_heat_down = []
spacers_heat_up = []
spacers_osmotic_down = []
spacers_osmotic_up = []
spacers_salinity_down = []
spacers_salinity_up = []

#defining the output file
out = open('C:\Users\Anish\Documents\LOP\element_positions_microarray.txt','w')

#extracting heat gene ids
df_heat_down1 = pd.ExcelFile("C:\Users\Anish\Documents\LOP\heat_result\heat-downreg_acgt_tgac.xls").parse("Sheet1")
genes_heat_down1 = df_heat_down1['Unnamed: 1'].values.tolist()
genes_heat_down1 = genes_heat_down1[2:]
genes_heat_down1 = [x for x in genes_heat_down1 if str(x) != 'nan']

df_heat_down2 = pd.ExcelFile("C:\Users\Anish\Documents\LOP\heat_result\heat-downreg_tgac-acgt.xls").parse("Sheet1")
genes_heat_down2 = df_heat_down2['Unnamed: 1'].values.tolist()
genes_heat_down2 = genes_heat_down2[2:]
genes_heat_down2 = [x for x in genes_heat_down2 if str(x) != 'nan']

df_heat_down3 = pd.ExcelFile("C:\Users\Anish\Documents\LOP\heat_result\heat-downreg_tgac-tgac.xls").parse("Sheet1")
genes_heat_down3 = df_heat_down3['Unnamed: 1'].values.tolist()
genes_heat_down3 = genes_heat_down3[2:]
genes_heat_down3 = [x for x in genes_heat_down3 if str(x) != 'nan']

genes_heat_down = union(union(genes_heat_down1,genes_heat_down2),genes_heat_down3)

df_heat_up1 = pd.ExcelFile("C:\Users\Anish\Documents\LOP\heat_result\heat-upreg_acgt_tgac.xls").parse("Sheet1")
genes_heat_up1 = df_heat_up1['Unnamed: 1'].values.tolist()
genes_heat_up1 = genes_heat_up1[2:]
genes_heat_up1 = [x for x in genes_heat_up1 if str(x) != 'nan']

df_heat_up2 = pd.ExcelFile("C:\Users\Anish\Documents\LOP\heat_result\heat-upreg_tgac-acgt.xls").parse("Sheet1")
genes_heat_up2 = df_heat_up2['Unnamed: 1'].values.tolist()
genes_heat_up2 = genes_heat_up2[2:]
genes_heat_up2 = [x for x in genes_heat_up2 if str(x) != 'nan']

df_heat_up3 = pd.ExcelFile("C:\Users\Anish\Documents\LOP\heat_result\heat-upreg_tgac-tgac.xls").parse("Sheet1")
genes_heat_up3 = df_heat_up3['Unnamed: 1'].values.tolist()
genes_heat_up3 = genes_heat_up3[2:]
genes_heat_up3 = [x for x in genes_heat_up3 if str(x) != 'nan']

genes_heat_up = union(union(genes_heat_up1,genes_heat_up2),genes_heat_up3)

#extracting osmotic gene ids
df_osmotic_down1 = pd.ExcelFile("C:\Users\Anish\Documents\LOP\osmotic_result\osmotic-downreg_acgt-tgac.xls").parse("Sheet1")
genes_osmotic_down1 = df_osmotic_down1['Unnamed: 1'].values.tolist()
genes_osmotic_down1 = genes_osmotic_down1[2:]
genes_osmotic_down1 = [x for x in genes_osmotic_down1 if str(x) != 'nan']

df_osmotic_down2 = pd.ExcelFile("C:\Users\Anish\Documents\LOP\osmotic_result\osmotic-downreg_tgac-acgt.xls").parse("Sheet1")
genes_osmotic_down2 = df_osmotic_down2['Unnamed: 1'].values.tolist()
genes_osmotic_down2 = genes_osmotic_down2[2:]
genes_osmotic_down2 = [x for x in genes_osmotic_down2 if str(x) != 'nan']

df_osmotic_down3 = pd.ExcelFile("C:\Users\Anish\Documents\LOP\osmotic_result\osmotic-downreg_tgac-tgac.xls").parse("Sheet1")
genes_osmotic_down3 = df_osmotic_down3['Unnamed: 1'].values.tolist()
genes_osmotic_down3 = genes_osmotic_down3[2:]
genes_osmotic_down3 = [x for x in genes_osmotic_down3 if str(x) != 'nan']

genes_osmotic_down = union(union(genes_osmotic_down1,genes_osmotic_down2),genes_osmotic_down3)

df_osmotic_up1 = pd.ExcelFile("C:\Users\Anish\Documents\LOP\osmotic_result\osmotic-upreg_acgt-tgac.xls").parse("Sheet1")
genes_osmotic_up1 = df_osmotic_up1['Unnamed: 1'].values.tolist()
genes_osmotic_up1 = genes_osmotic_up1[2:]
genes_osmotic_up1 = [x for x in genes_osmotic_up1 if str(x) != 'nan']

df_osmotic_up2 = pd.ExcelFile("C:\Users\Anish\Documents\LOP\osmotic_result\osmotic-upreg_tgac-acgt.xls").parse("Sheet1")
genes_osmotic_up2 = df_osmotic_up2['Unnamed: 1'].values.tolist()
genes_osmotic_up2 = genes_osmotic_up2[2:]
genes_osmotic_up2 = [x for x in genes_osmotic_up2 if str(x) != 'nan']

df_osmotic_up3 = pd.ExcelFile("C:\Users\Anish\Documents\LOP\osmotic_result\osmotic-upreg_tgac-tgac.xls").parse("Sheet1")
genes_osmotic_up3 = df_osmotic_up3['Unnamed: 1'].values.tolist()
genes_osmotic_up3 = genes_osmotic_up3[2:]
genes_osmotic_up3 = [x for x in genes_osmotic_up3 if str(x) != 'nan']

genes_osmotic_up = union(union(genes_osmotic_up1,genes_osmotic_up2),genes_osmotic_up3)

#extracting salinity gene ids
df_salinity_down1 = pd.ExcelFile("C:\Users\Anish\Documents\LOP\salinity_result\salinity-downreg_acgt-tgac.xls").parse("Sheet1")
genes_salinity_down1 = df_salinity_down1['Unnamed: 1'].values.tolist()
genes_salinity_down1 = genes_salinity_down1[2:]
genes_salinity_down1 = [x for x in genes_salinity_down1 if str(x) != 'nan']

df_salinity_down2 = pd.ExcelFile("C:\Users\Anish\Documents\LOP\salinity_result\salinity-downreg_tgac-acgt.xls").parse("Sheet1")
genes_salinity_down2 = df_salinity_down2['Unnamed: 1'].values.tolist()
genes_salinity_down2 = genes_salinity_down2[2:]
genes_salinity_down2 = [x for x in genes_salinity_down2 if str(x) != 'nan']

df_salinity_down3 = pd.ExcelFile("C:\Users\Anish\Documents\LOP\salinity_result\salinity-downreg_tgac-tgac.xls").parse("Sheet1")
genes_salinity_down3 = df_salinity_down3['Unnamed: 1'].values.tolist()
genes_salinity_down3 = genes_salinity_down3[2:]
genes_salinity_down3 = [x for x in genes_salinity_down3 if str(x) != 'nan']

genes_salinity_down = union(union(genes_salinity_down1,genes_salinity_down2),genes_salinity_down3)

df_salinity_up1 = pd.ExcelFile("C:\Users\Anish\Documents\LOP\salinity_result\salinity-upreg_acgt-tgac.xls").parse("Sheet1")
genes_salinity_up1 = df_salinity_up1['Unnamed: 1'].values.tolist()
genes_salinity_up1 = genes_salinity_up1[2:]
genes_salinity_up1 = [x for x in genes_salinity_up1 if str(x) != 'nan']

df_salinity_up2 = pd.ExcelFile("C:\Users\Anish\Documents\LOP\salinity_result\salinity-upreg_tgac-acgt.xls").parse("Sheet1")
genes_salinity_up2 = df_salinity_up2['Unnamed: 1'].values.tolist()
genes_salinity_up2 = genes_salinity_up2[2:]
genes_salinity_up2 = [x for x in genes_salinity_up2 if str(x) != 'nan']

df_salinity_up3 = pd.ExcelFile("C:\Users\Anish\Documents\LOP\salinity_result\salinity-upreg_tgac-tgac.xls").parse("Sheet1")
genes_salinity_up3 = df_salinity_up3['Unnamed: 1'].values.tolist()
genes_salinity_up3 = genes_salinity_up3[2:]
genes_salinity_up3 = [x for x in genes_salinity_up3 if str(x) != 'nan']

genes_salinity_up = union(union(genes_salinity_up1,genes_salinity_up2),genes_salinity_up3)

#tokenizing promoter for analysis
promoter_tokens = promoter.strip().split('\n')

for i in range(len(promoter_tokens)):
	token = promoter_tokens[i].split('\t')
	token[0] = token[0].upper()
	token[0] = token[0][:-2]
	if token[0] not in promoter_dict.keys():
		promoter_dict[token[0]] = token[1]

#Actual analysis for each type of microarray gene ids - heat, osmotic,salinity - up,down regulated
print("Heat downregulated")
out.write("Heat downregulated\n")
for gene in genes_heat_down:
	if str(gene) in promoter_dict.keys():
		inp = promoter_dict[str(gene)]
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
				positions.append([seq1,i])
		
			if k2 == seq2:
				count = count + 1
				positions.append([seq2,i])
		
		if(len(positions) != 0):
			print(gene.encode('ascii','ignore'),positions)
			out.write(gene.encode('ascii','ignore') + '\t')
			for element in positions:
				out.write(str(element))
				out.write('\t')
			out.write('\n')

		if(len(positions) > 1):
			for j in range(len(positions) - 1):
				spacers_heat_down.append(positions[j + 1][1] - positions[j][1])

count1 = 0
count2 = 0
count3 = 0
count4 = 0
for j in range(len(spacers_heat_down)):
	if spacers_heat_down[j] < 100:
		count1 = count1 + 1
	elif spacers_heat_down[j] > 100 and spacers_heat_down[j] < 500:
		count2 = count2 + 1
	elif spacers_heat_down[j] > 500 and spacers_heat_down[j] < 1000:
		count3 = count3 + 1
	else:
		count4 = count4 + 1

spacers_heat_down_result = dict()
spacers_heat_down_result['spacer less than 100'] = count1
spacers_heat_down_result['spacer between 100 to 500'] = count2
spacers_heat_down_result['spacer between 500 to 1000'] = count3
spacers_heat_down_result['spacer greater than 1000'] = count4

out.write('spacer info\t')
print(spacers_heat_down)
for item in spacers_heat_down:
	out.write("%s\t" %item)
print(spacers_heat_down_result)
out.write('\n')
for item in spacers_heat_down_result:
	out.write("%s\t" %item)
	out.write("%d\n" %spacers_heat_down_result[item])

out.write('\n')
print("Heat upregulated")
out.write("Heat upregulated\n")
for gene in genes_heat_up:
	if str(gene) in promoter_dict.keys():
		inp = promoter_dict[str(gene)]
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
				positions.append([seq1,i])
				
			if k2 == seq2:
				count = count + 1
				positions.append([seq2,i])
				
		if(len(positions) != 0):
			print(gene.encode('ascii','ignore'),positions)
			out.write(gene.encode('ascii','ignore') + '\t')
			for element in positions:
				out.write(str(element))
				out.write('\t')
			out.write('\n')

		if(len(positions) > 1):
			for j in range(len(positions) - 1):
				spacers_heat_up.append(positions[j + 1][1] - positions[j][1])

count1 = 0
count2 = 0
count3 = 0
count4 = 0
for j in range(len(spacers_heat_up)):
	if spacers_heat_up[j] < 100:
		count1 = count1 + 1
	elif spacers_heat_up[j] > 100 and spacers_heat_up[j] < 500:
		count2 = count2 + 1
	elif spacers_heat_up[j] > 500 and spacers_heat_up[j] < 1000:
		count3 = count3 + 1
	else:
		count4 = count4 + 1

spacers_heat_up_result = dict()
spacers_heat_up_result['spacer less than 100'] = count1
spacers_heat_up_result['spacer between 100 to 500'] = count2
spacers_heat_up_result['spacer between 500 to 1000'] = count3
spacers_heat_up_result['spacer greater than 1000'] = count4

out.write('spacer info\t')
print(spacers_heat_up)
for item in spacers_heat_up:
	out.write("%s\t" %item)
print(spacers_heat_up_result)
out.write('\n')
for item in spacers_heat_up_result:
	out.write("%s\t" %item)
	out.write("%d\n" %spacers_heat_up_result[item])

out.write('\n')		
print("Osmotic downregulated")
out.write("Osmotic downregulated\n")
for gene in genes_osmotic_down:
	if str(gene) in promoter_dict.keys():
		inp = promoter_dict[str(gene)]
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
				positions.append([seq1,i])
				
			if k2 == seq2:
				count = count + 1
				positions.append([seq2,i])
				
		if(len(positions) != 0):
			print(gene.encode('ascii','ignore'),positions)
			out.write(gene.encode('ascii','ignore') + '\t')
			for element in positions:
				out.write(str(element))
				out.write('\t')
			out.write('\n')

		if(len(positions) > 1):
			for j in range(len(positions) - 1):
				spacers_osmotic_down.append(positions[j + 1][1] - positions[j][1])

count1 = 0
count2 = 0
count3 = 0
count4 = 0
for j in range(len(spacers_osmotic_down)):
	if spacers_osmotic_down[j] < 100:
		count1 = count1 + 1
	elif spacers_osmotic_down[j] > 100 and spacers_osmotic_down[j] < 500:
		count2 = count2 + 1
	elif spacers_osmotic_down[j] > 500 and spacers_osmotic_down[j] < 1000:
		count3 = count3 + 1
	else:
		count4 = count4 + 1

spacers_osmotic_down_result = dict()
spacers_osmotic_down_result['spacer less than 100'] = count1
spacers_osmotic_down_result['spacer between 100 to 500'] = count2
spacers_osmotic_down_result['spacer between 500 to 1000'] = count3
spacers_osmotic_down_result['spacer greater than 1000'] = count4

out.write('spacer info\t')
print(spacers_osmotic_down)
for item in spacers_osmotic_down:
	out.write("%s\t" %item)
print(spacers_osmotic_down_result)
out.write('\n')
for item in spacers_osmotic_down_result:
	out.write("%s\t" %item)
	out.write("%d\n" %spacers_osmotic_down_result[item])

out.write('\n')
print("Osmotic upregulated")
out.write("Osmotic upregulated\n")
for gene in genes_osmotic_up:
	if str(gene) in promoter_dict.keys():
		inp = promoter_dict[str(gene)]
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
				positions.append([seq1,i])
				
			if k2 == seq2:
				count = count + 1
				positions.append([seq2,i])
				
		if(len(positions) != 0):
			print(gene.encode('ascii','ignore'),positions)
			out.write(gene.encode('ascii','ignore') + '\t')
			for element in positions:
				out.write(str(element))
				out.write('\t')
			out.write('\n')

		if(len(positions) > 1):
			for j in range(len(positions) - 1):
				spacers_osmotic_up.append(positions[j + 1][1] - positions[j][1])

count1 = 0
count2 = 0
count3 = 0
count4 = 0
for j in range(len(spacers_osmotic_up)):
	if spacers_osmotic_up[j] < 100:
		count1 = count1 + 1
	elif spacers_osmotic_up[j] > 100 and spacers_osmotic_up[j] < 500:
		count2 = count2 + 1
	elif spacers_osmotic_up[j] > 500 and spacers_osmotic_up[j] < 1000:
		count3 = count3 + 1
	else:
		count4 = count4 + 1

spacers_osmotic_up_result = dict()
spacers_osmotic_up_result['spacer less than 100'] = count1
spacers_osmotic_up_result['spacer between 100 to 500'] = count2
spacers_osmotic_up_result['spacer between 500 to 1000'] = count3
spacers_osmotic_up_result['spacer greater than 1000'] = count4

out.write('spacer info\t')
print(spacers_osmotic_up)
for item in spacers_osmotic_up:
	out.write("%s\t" %item)
print(spacers_osmotic_up_result)
out.write('\n')
for item in spacers_osmotic_up_result:
	out.write("%s\t" %item)
	out.write("%d\n" %spacers_osmotic_up_result[item])

out.write('\n')
print("Salinity downregulated")
out.write("Salinity downregulated\n")
for gene in genes_salinity_down:
	if str(gene) in promoter_dict.keys():
		inp = promoter_dict[str(gene)]
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
				positions.append([seq1,i])
				
			if k2 == seq2:
				count = count + 1
				positions.append([seq2,i])
				
		if(len(positions) != 0):
			print(gene.encode('ascii','ignore'),positions)
			out.write(gene.encode('ascii','ignore') + '\t')
			for element in positions:
				out.write(str(element))
				out.write('\t')
			out.write('\n')

		if(len(positions) > 1):
			for j in range(len(positions) - 1):
				spacers_salinity_down.append(positions[j + 1][1] - positions[j][1])

count1 = 0
count2 = 0
count3 = 0
count4 = 0
for j in range(len(spacers_salinity_down)):
	if spacers_salinity_down[j] < 100:
		count1 = count1 + 1
	elif spacers_salinity_down[j] > 100 and spacers_salinity_down[j] < 500:
		count2 = count2 + 1
	elif spacers_salinity_down[j] > 500 and spacers_salinity_down[j] < 1000:
		count3 = count3 + 1
	else:
		count4 = count4 + 1

spacers_salinity_down_result = dict()
spacers_salinity_down_result['spacer less than 100'] = count1
spacers_salinity_down_result['spacer between 100 to 500'] = count2
spacers_salinity_down_result['spacer between 500 to 1000'] = count3
spacers_salinity_down_result['spacer greater than 1000'] = count4

out.write('spacer info\t')
print(spacers_salinity_down)
for item in spacers_salinity_down:
	out.write("%s\t" %item)
print(spacers_salinity_down_result)
out.write('\n')
for item in spacers_salinity_down_result:
	out.write("%s\t" %item)
	out.write("%d\n" %spacers_salinity_down_result[item])

out.write('\n')
print("Salinity upregulated")
out.write("Salinity upregulated\n")
for gene in genes_salinity_up:
	if str(gene) in promoter_dict.keys():
		inp = promoter_dict[str(gene)]
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
				positions.append([seq1,i])
				
			if k2 == seq2:
				count = count + 1
				positions.append([seq2,i])
				
		if(len(positions) != 0):
			print(gene.encode('ascii','ignore'),positions)
			out.write(gene.encode('ascii','ignore') + '\t')
			for element in positions:
				out.write(str(element))
				out.write('\t')
			out.write('\n')

		if(len(positions) > 1):
			for j in range(len(positions) - 1):
				spacers_salinity_up.append(positions[j + 1][1] - positions[j][1])

count1 = 0
count2 = 0
count3 = 0
count4 = 0
for j in range(len(spacers_salinity_up)):
	if spacers_salinity_up[j] < 100:
		count1 = count1 + 1
	elif spacers_salinity_up[j] > 100 and spacers_salinity_up[j] < 500:
		count2 = count2 + 1
	elif spacers_salinity_up[j] > 500 and spacers_salinity_up[j] < 1000:
		count3 = count3 + 1
	else:
		count4 = count4 + 1

spacers_salinity_up_result = dict()
spacers_salinity_up_result['spacer less than 100'] = count1
spacers_salinity_up_result['spacer between 100 to 500'] = count2
spacers_salinity_up_result['spacer between 500 to 1000'] = count3
spacers_salinity_up_result['spacer greater than 1000'] = count4

out.write('spacer info\t')
print(spacers_salinity_up)
for item in spacers_salinity_up:
	out.write("%s\t" %item)
print(spacers_salinity_up_result)
out.write('\n')
for item in spacers_salinity_up_result:
	out.write("%s\t" %item)
	out.write("%d\n" %spacers_salinity_up_result[item])