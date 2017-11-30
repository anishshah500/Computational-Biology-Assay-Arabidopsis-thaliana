import os
import numpy as np
import pandas as pd
import collections
from collections import defaultdict
import matplotlib.pyplot as plt
from scipy.stats import ttest_ind

pd.set_option('display.max_rows', 50)
pd.set_option('display.max_columns', 10)
pd.set_option('display.width', 1000)

compileIntersection = False

def union(a, b):
    """ return the union of two lists """
    return list(set(a) | set(b))

spacers_heat_down = []
spacers_heat_up = []
spacers_osmotic_down = []
spacers_osmotic_up = []
spacers_salinity_down = []
spacers_salinity_up = []

debug = False

class graphMaker():

	def __init__(self,config):
		self.typeOfCond = config['condition']
		self.upOrDown   = config['upOrDown']
		self.outFile	= config['self.outFile'] + "_" + config['condition'] + '_' + config['upOrDown'] + ".txt"
		self.degeneracy = config['degeneracy']
		self.genome		= config['genome']

		self.spacers    = []
		self.genes		= []
		
		self.promoter 	= ""
		
		self.promoterTokens = []
		self.promoterDict   = {}
		self.allPos			= []

		self.CSVFile1	= 'C:/Users/Anish/Documents/LOP/' + self.typeOfCond + '_result/' + self.typeOfCond + '-' + self.upOrDown + 'reg_acgt-tgac.xls'
		self.CSVFile2	= 'C:/Users/Anish/Documents/LOP/' + self.typeOfCond + '_result/' + self.typeOfCond + '-' + self.upOrDown + 'reg_tgac-acgt.xls'
		self.CSVFile3	= 'C:/Users/Anish/Documents/LOP/' + self.typeOfCond + '_result/' + self.typeOfCond + '-' + self.upOrDown + 'reg_tgac-tgac.xls'
		
		if(debug):
			print(self.CSVFile1)
			print(self.CSVFile2)
			print(self.CSVFile3)
		
		self.promoterFile = 'C:/data/SOP/gtggttag - promoter/PromoterSeq.txt'
		
		self.promoter_dict = dict()
		self.seq1 		   = config['element1']

		if(self.degeneracy):
			self.seq2		   = config['element2']
		
		self.length1 = len(self.seq1)
		
		if(self.degeneracy):
			self.length2 = len(self.seq2)

		self.spacerDict = {}
		for i in range(30):
			self.spacerDict[i+1] = set()

		self.spacerSeqDict = {}

		self.spacerDF = pd.DataFrame()
		self.tTestDF  = pd.DataFrame()

		self.spacerLengthDict = dict()
		for i in range(30):
			self.spacerLengthDict[i + 1] = []

		self.spacerSequencePredictDict = dict()
		for spacerLength in range(30):
			self.spacerSequencePredictDict[spacerLength + 1] = dict()
			for eachSpacerLength in range(spacerLength + 1):
				self.spacerSequencePredictDict[spacerLength + 1][eachSpacerLength + 1] 		= dict()
				self.spacerSequencePredictDict[spacerLength + 1][eachSpacerLength + 1]['A'] = 0
				self.spacerSequencePredictDict[spacerLength + 1][eachSpacerLength + 1]['T'] = 0
				self.spacerSequencePredictDict[spacerLength + 1][eachSpacerLength + 1]['C'] = 0
				self.spacerSequencePredictDict[spacerLength + 1][eachSpacerLength + 1]['G'] = 0

	def readFiles(self):
		self.promoter = open(self.promoterFile,'r').read()

		self.promoterTokens = self.promoter.strip().split('\n')
		
		self.genes	  += pd.ExcelFile(self.CSVFile1).parse("Sheet1")['Unnamed: 1'].dropna().values.tolist()[2:]
		self.genes    += pd.ExcelFile(self.CSVFile2).parse("Sheet1")['Unnamed: 1'].dropna().values.tolist()[2:]
		self.genes    += pd.ExcelFile(self.CSVFile3).parse("Sheet1")['Unnamed: 1'].dropna().values.tolist()[2:]

		self.genes = list(set(self.genes))

		if(debug):
			print("Genes")
			print(self.genes)
			print("Length: " + str(len(self.genes)))

	def preProcessPromoter(self):
		for i in range(len(self.promoterTokens)):
			
			token 	 = self.promoterTokens[i].split('\t')
			token[0] = token[0].upper()
			token[0] = token[0][:-2]
			
			if token[0] not in self.promoterDict.keys():
				self.promoterDict[token[0]] = token[1]

		if(debug):
			print(self.promoterDict)

	def countNoOfNucleotides(self):

		count = 0

		for gene in self.genes:
			if str(gene) in self.promoterDict.keys():
				count += len(self.promoterDict[str(gene)])

		print(self.typeOfCond + " " + self.upOrDown + ": " + str(count))

	def analysis(self, plotGraph = False, verbose = False, printBins = False, window = 100, saveGraph = False, showGraph = True):
		
		out = open(self.outFile,'w')

		if(verbose):
			print(self.typeOfCond + " " + self.upOrDown + "regulated")

		out.write(self.typeOfCond + " " + self.upOrDown + "regulated")

		for gene in self.genes:
			if str(gene) in self.promoterDict.keys():
				inp = self.promoterDict[str(gene)]
				count = 0
				positions = []
				distance = []
				distance_consider = []
				spacers = []
				a = defaultdict(list)

				for i in range(len(inp)):
			
					k1 = inp[i:i+self.length1]
					
					if(self.degeneracy):
						k2 = inp[i:i+self.length2]
			
					if k1 == self.seq1:
						count = count + 1
						positions.append([self.seq1,i])
						self.allPos.append([self.seq1,gene,i])
				
					if(self.degeneracy):
						if k2 == self.seq2:
							count = count + 1
							positions.append([self.seq2,i])
							self.allPos.append([self.seq1,gene,i])
				
				if(len(positions) != 0):
					if(verbose):
						print(gene.encode('ascii','ignore'),positions)
					
					out.write(str(gene.encode('ascii','ignore')) + '\t')
					
					for element in positions:
						out.write(str(element))
						out.write('\t')
					out.write('\n')

				if(len(positions) > 1):
					for j in range(len(positions) - 1):
						spacers_heat_down.append(positions[j + 1][1] - positions[j][1])

		count1,count2,count3,count4 = 0,0,0,0
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

		if(printBins):
			spacers_heat_down_result['spacer less than 100'] = count1
			spacers_heat_down_result['spacer between 100 to 500'] = count2
			spacers_heat_down_result['spacer between 500 to 1000'] = count3
			spacers_heat_down_result['spacer greater than 1000'] = count4

		out.write('spacer info\t')
		# print(spacers_heat_down)
		for item in spacers_heat_down:
			out.write("%s\t" %item)
		# print(spacers_heat_down_result)
		out.write('\n')
		for item in spacers_heat_down_result:
			out.write("%s\t" %item)
			out.write("%d\n" %spacers_heat_down_result[item])

		if(plotGraph):
			plotData = pd.DataFrame(self.allPos)
			bins 	 = list(range(0,plotData.iloc[:,2].max()+25,window))
			labels   = [i for i in range(len(bins)-1)]

			plotData['bin'] = pd.cut(plotData.iloc[:,2], bins, labels = labels)
			
			plt.figure(figsize=(12, 9))
			ax = plt.subplot(111)
			plt.title(self.typeOfCond + " " + self.upOrDown + "regulated, " + "Window Size = " + str(window))
			plt.xlabel('Position of gene in genome')
			plt.ylabel('Frequency')
			ax.spines["top"].set_visible(False)  
			ax.spines["right"].set_visible(False)
			ax.get_xaxis().tick_bottom()  
			ax.get_yaxis().tick_left()
			plt.yticks(fontsize = 14)
			plt.xticks(range(0,plotData.iloc[:,2].max()+25,window), fontsize = 6, rotation = 60)

			plt.hist(plotData.iloc[:,2],color="#3F5D7D",bins = bins)

			if(saveGraph):
				path = "results/graphResults/" + self.seq1 + "_" + self.typeOfCond + "_" + self.upOrDown + ".png"

				if(not(os.path.exists(path))):
					plt.savefig(path)

			if(showGraph):
				plt.show()

	def spacerAnalysis(self, printNoOfSpacers = True, printZeroSpacer = True):
		if(debug):
			print(self.spacerDict)

		if(printZeroSpacer):
			zeroSpacerCount = 0

		for gene in self.genes:
			if gene.upper() in self.promoterDict:
				promoterSeq 	   = self.promoterDict[gene.upper()]
				currMotifPositions = [n for n in range(len(promoterSeq)) if promoterSeq.find(self.seq1, n) == n]
				currSpacers		   = [currMotifPositions[i+1] - currMotifPositions[i] - len(self.seq1) for i in range(len(currMotifPositions) - 1)]

				for i in currSpacers:
					if i not in self.spacerDict.keys():
						self.spacerDict[i] = set()
					self.spacerDict[i].add(gene)

				for i in range(len(currMotifPositions) - 1):
					currSpacerLength = currMotifPositions[i+1] - currMotifPositions[i] - len(self.seq1)

					if(printZeroSpacer):
						if(currSpacerLength == 30):
							zeroSpacerCount += 1

					if currSpacerLength <= 30:
						if currSpacerLength not in self.spacerSeqDict:
							self.spacerSeqDict[currSpacerLength] = ''
						self.spacerSeqDict[currSpacerLength] += promoterSeq[currMotifPositions[i] + len(self.seq1) : currMotifPositions[i+1]]

		self.spacerDict = sorted(self.spacerDict.items())

		if(debug):
			print(self.spacerDict)
			print(self.spacerSeqDict)

		if(printNoOfSpacers):

			noOfSpacersList = []

			file = "results/no_Of_Spacers_Of_Each_Length/" + self.seq1 + "_" + self.typeOfCond + "_" + self.upOrDown + ".csv"

			noOfSpacersList.append([0,zeroSpacerCount])
			for length in self.spacerSeqDict:
				if(length > 0):
					noOfSpacersList.append([length,int(len(self.spacerSeqDict[length]) / length)])

			if(debug):
				print(pd.DataFrame(noOfSpacersList, columns = ["Length","No of Spacers"]).set_index("Length"))

			pd.DataFrame(noOfSpacersList, columns = ["Length","No of Spacers"]).set_index("Length").to_csv(file)

		return self.spacerDict

	def analyzeFlankingRegions(self, printFlag = False, writeFlag = True):
		
		flankingRegionDict = dict()
		flankingRegionList = []

		for length in range(4,31):
			flankingRegionDict[length] = []

		for gene in self.genes:
			if gene.upper() in self.promoterDict:
				promoterSeq 	   = self.promoterDict[gene.upper()]
				currMotifPositions = [n for n in range(len(promoterSeq)) if promoterSeq.find(self.seq1, n) == n]

				for i in range(len(currMotifPositions) - 1):

					currSpacerLength = currMotifPositions[i+1] - currMotifPositions[i] - len(self.seq1)

					if(currSpacerLength >= 4 and currSpacerLength <= 30):
						
						if(debug):
							print(promoterSeq[currMotifPositions[i] - 3 : currMotifPositions[i+1] + 3 + len(self.seq1)], currSpacerLength)
						
						flankingRegionDict[currSpacerLength].append([promoterSeq[currMotifPositions[i] - 3 : currMotifPositions[i]],\
							promoterSeq[currMotifPositions[i] + len(self.seq1) : currMotifPositions[i] + 3 + len(self.seq1)],\
							promoterSeq[currMotifPositions[i + 1] - 3 : currMotifPositions[i + 1]],\
							promoterSeq[currMotifPositions[i + 1] + len(self.seq1) : currMotifPositions[i + 1] + 3 + len(self.seq1)]])

						if(debug):
							print(flankingRegionDict[currSpacerLength][-1])

		for length in range(4,31):

			leftFirst, rightFirst, leftSecond, rightSecond = '','','',''

			for i in range(len(flankingRegionDict[length])):
				leftFirst   += flankingRegionDict[length][i][0]
				rightFirst  += flankingRegionDict[length][i][1]
				leftSecond  += flankingRegionDict[length][i][2]
				rightSecond += flankingRegionDict[length][i][3]

			if(debug):
				print(leftFirst, collections.Counter(leftFirst))
				print(rightFirst, collections.Counter(rightFirst))
				print(leftSecond, collections.Counter(leftSecond))
				print(rightSecond, collections.Counter(rightSecond))

			leftFirstCollection   = collections.Counter(leftFirst)
			rightFirstCollection  = collections.Counter(rightFirst)
			leftSecondCollection  = collections.Counter(leftSecond)
			rightSecondCollection = collections.Counter(rightSecond)

			flankingRegionList.append([length,"leftFirst",leftFirstCollection['a'],leftFirstCollection['t'],\
				leftFirstCollection['c'],leftFirstCollection['g']])
			flankingRegionList.append([length,"rightFirst",rightFirstCollection['a'],rightFirstCollection['t'],\
				rightFirstCollection['c'],rightFirstCollection['g']])
			flankingRegionList.append([length,"leftSecond",leftSecondCollection['a'],leftSecondCollection['t'],\
				leftSecondCollection['c'],leftSecondCollection['g']])
			flankingRegionList.append([length,"rightSecond",rightSecondCollection['a'],rightSecondCollection['t'],\
				rightSecondCollection['c'],rightSecondCollection['g']])

		flankingRegionDF = pd.DataFrame(flankingRegionList, columns = ['Length','FlankingRegion','A','T','C','G'])

		grouped_df = flankingRegionDF.groupby('Length')

		if(printFlag):
			for key, item in grouped_df:
			    print(grouped_df.get_group(key), "\n\n")

		if(writeFlag):
			file = "results/flanking_Region_Analysis/" + self.seq1 + "_" + self.typeOfCond + "_" + self.upOrDown + ".csv"
			flankingRegionDF.to_csv(file, index = False)

	def printSpacerGenes(self, lowerLimit, UpperLimit, printFlag = True, writeFlag = True):

		if(printFlag):
			for i in range(lowerLimit-1, UpperLimit):
				print("Spacer " + str(i+1) + ": ")
				print(self.spacerDict[i][1])

		if(writeFlag):

			path = "results/geneResults/" + self.genome + "/"

			if(not(os.path.exists(path))):
				os.mkdir(path)

			out = open(path + self.typeOfCond + "_" + self.upOrDown + ".txt", "w")

			for i in range(lowerLimit-1, UpperLimit):
				out.write("\nSpacer " + str(i+1) + ": ")
				out.write(str(self.spacerDict[i][1]))
				out.write('\n')

	def reverseAndPrintDictionary(self, printFlag = False, writeFlag = True):
		reverseDict = {}

		for length in self.spacerDict:

			if length[0] <= 30:
				
				currGenes = length[1]

				for gene in currGenes:
					if gene not in reverseDict:
						reverseDict[gene] = set()
					reverseDict[gene].add(length[0])

		for gene in reverseDict:
			reverseDict[gene] = sorted(reverseDict[gene])

		if(printFlag):
			print(reverseDict)

		if(writeFlag):
			path = "results/geneResults/" + self.genome + "/"

			if(not(os.path.exists(path))):
				os.mkdir(path)

			out = open(path + self.typeOfCond + "_" + self.upOrDown + "_geneToSpacerMapping.txt", "w")

			for i in reverseDict:
				out.write("Gene " + str(i) + ": ")
				out.write(str(reverseDict[i]))
				out.write('\n')

	def performTTest(self, spacerSeries1, spacerSeries2):

		tStat, pValue =  ttest_ind(spacerSeries1, spacerSeries2)

		return tStat, pValue

	def analyzeSpacerDist(self, tTest = False, verbose = False, writeFlag = True, clipBoardFlag = False, printNucleotideCount = True):

		dfList = []
		dfListNormalized = []

		for length in self.spacerSeqDict.keys():
			
			letters 	 = collections.Counter(self.spacerSeqDict[length])
			total		 = len(self.spacerSeqDict[length])

			#To avoid zero division
			if(total < 1):
				total = 1

			if(printNucleotideCount):
				currFreqList = [length,letters['a'],letters['t'], letters['c'], letters['g']]
				dfList.append(currFreqList)

			currFreqListNormalized = [length,letters['a'] / total,letters['t'] / total,letters['c'] / total,letters['g'] / total]
			dfListNormalized.append(currFreqListNormalized)

			if(debug):
				print(length, letters)
				print(currFreqListNormalized)

		if(printNucleotideCount):
			self.spacerDF = pd.DataFrame(dfList, columns = ['length','A','T','C','G'])
			self.spacerDF = self.spacerDF.set_index('length')
			self.spacerDF = self.spacerDF.iloc[1:]			
			
			if(writeFlag):
				self.spacerDF.to_csv('results/nucleotide_Count_For_Spacer_Length/'  + self.typeOfCond + '_' + self.upOrDown + '_' +\
					self.seq1 + '_notNormalized.csv')

		self.spacerDF = pd.DataFrame(dfListNormalized, columns = ['length','A','T','C','G'])
		self.spacerDF = self.spacerDF.set_index('length')
		self.spacerDF = self.spacerDF.iloc[1:]
		if(writeFlag):
				self.spacerDF.to_csv('results/nucleotide_Count_For_Spacer_Length/'  + self.typeOfCond + '_' + self.upOrDown + '_' +\
					self.seq1 + '_normalized.csv')

		if(debug):
			print(self.spacerDF)
			self.spacerDF.to_clipboard()
			input(self.typeOfCond + ', ' + self.upOrDown + ' Dataframe copied')

		if(tTest):
			
			tTestList = []
			
			for seperation in range(2,29,1):
				
				currTTestList = [1, seperation, seperation+1, 30]
				
				for nucleotide in ['A','T','C','G']:
					
					series1 = self.spacerDF[nucleotide].iloc[:seperation]
					series2 = self.spacerDF[nucleotide].iloc[seperation:]

					tStat, pVal = self.performTTest(series1, series2)

					if(verbose):
						print(nucleotide)
						print(series1,series2)
						print(tStat, pVal)

					currTTestList.append(tStat)
					currTTestList.append(pVal)

				tTestList.append(currTTestList)

			if(debug):
				print(tTestList)

			columns = ['Series1 - Start', 'Series1 - End', 'Series2 - Start', 'Series2 - End']
			columns.extend(['t-Stat: A', 'p-value: A'])
			columns.extend(['t-Stat: T', 'p-value: T'])
			columns.extend(['t-Stat: C', 'p-value: C'])
			columns.extend(['t-Stat: G', 'p-value: G'])

			self.tTestDF = pd.DataFrame(tTestList, columns = columns)

			print(self.tTestDF.round(7))

			# if(writeFlag):
			# 	self.tTestDF.round(7).to_csv('results/comparativeSpacerComposition/' + self.typeOfCond + '_' + self.upOrDown + '_tTest_'\
			# 		+ self.seq1 + '_' + self.seq1 + '.csv')

			if(clipBoardFlag):
				self.tTestDF.round(7).to_clipboard()
				raw_input(self.typeOfCond + ', ' + self.upOrDown + 'Dataframe copied')

	def findMax(self, currDict):

		maxCount 	  = 0
		maxNucleotide = 'A'

		for nucleotide in ['A', 'T', 'C', 'G']:
			if currDict[nucleotide] > maxCount:
				maxCount 	  = currDict[nucleotide]
				maxNucleotide = nucleotide

		return maxNucleotide

	def predictSpacer(self, writeFlag = False, writeDistributionOfNucleotides = True):

		if(debug):
			print(spacerDict)
			print(self.spacerSequencePredictDict)

		if(writeFlag):
			path = "results/predictedSpacers/" + self.genome + "/"

			if(not(os.path.exists(path))):
				os.mkdir(path)

			out = open(path + self.typeOfCond + "_" + self.upOrDown + "_predictedSpacers.txt", "w")

		for gene in self.genes:
			if gene.upper() in self.promoterDict:
				promoterSeq 	   = self.promoterDict[gene.upper()]
				currMotifPositions = [n for n in range(len(promoterSeq)) if promoterSeq.find(self.seq1, n) == n]
				currSpacers		   = [currMotifPositions[i+1] - currMotifPositions[i] - len(self.seq1) for i in range(len(currMotifPositions) - 1)]

			for i in range(len(currMotifPositions) - 1):
					currSpacerLength = currMotifPositions[i+1] - currMotifPositions[i] - len(self.seq1)
					if currSpacerLength <= 30 and currSpacerLength > 0:
						currSpacerSeq = promoterSeq[currMotifPositions[i] + len(self.seq1) : currMotifPositions[i+1]]
						self.spacerLengthDict[len(currSpacerSeq)].append(currSpacerSeq)

		for position in range(1,31):
			currAnalyzeSeq = ''
			for length in range(position,31):
				for eachSeq in self.spacerLengthDict[length]:
					
					if(debug):
						print(position, length, eachSeq[position-1])
					
					self.spacerSequencePredictDict[length][position][eachSeq[position-1].upper()] += 1

		if(debug):
			print(self.spacerSequencePredictDict)


		for spacerLength in range(1,31):
			
			sequenceForCurrLength = ''
			
			for eachSpacerPosition in range(1,spacerLength + 1):
				
				predictedNucleotide = self.findMax(self.spacerSequencePredictDict[spacerLength][eachSpacerPosition])
				
				if(debug):
					print(predictedNucleotide)
				
				sequenceForCurrLength += predictedNucleotide
			
			print("Predicted Sequence for length " + str(spacerLength) + ": " + sequenceForCurrLength)

			if(writeFlag):
				out.write("Predicted Sequence for length " + str(spacerLength) + ": " + sequenceForCurrLength + "\n")

		if(writeDistributionOfNucleotides):

			file = "results/count_Of_Nucleotides_For_Each_Position_In_Each_Spacer/"  + self.seq1 + "_" + \
			self.typeOfCond + "_" + self.upOrDown + ".csv"

			outFile = open(file,"w")

			for length in range(1,31):
				outFile.write("\n\nDistribution of spacer length " + str(length) + ":\n")
				currDF = pd.DataFrame(self.spacerSequencePredictDict[length]).T
				currDF.index.names = ['Nucleotide Number']
				currDF.to_csv(outFile)

if __name__ == "__main__":

	if(compileIntersection):
		upSpacer = []
		downSpacer = []

	for condition in ["heat","osmotic","salinity"]:
		for upOrDown in ["down","up"]:
	# for condition in ["heat"]:	
	# 	for upOrDown in ["down"]:
					
			config   				 = {}
			config['condition'] 	 = condition
			config['upOrDown']  	 = upOrDown
			config['self.outFile']   = "C:/Users/Anish/Documents/LOP/element_positions_microarray"
			config['degeneracy']	 = False
			config['element1']		 = 'aaag'
			config['limitLength']	 = False
			config['genome']		 = 'Arabidopsis'

			if(compileIntersection):
				currDict = {}
				for i in range(30):
						currDict[i+1] = set()

			graphObj = graphMaker(config)
			graphObj.readFiles()
			graphObj.preProcessPromoter()
			# graphObj.analysis(plotGraph = False, saveGraph = False, showGraph = False)
			
			graphObj.spacerAnalysis()

			# graphObj.printSpacerGenes(0,30,False, False)

			# graphObj.reverseAndPrintDictionary()

			graphObj.analyzeSpacerDist()

			# graphObj.predictSpacer()

			# graphObj.analyzeFlankingRegions()

			# graphObj.countNoOfNucleotides()

			if(compileIntersection):
				currDict = graphObj.spacerAnalysis()
				if(upOrDown == "up"):
					upSpacer.append(currDict)
				else:
					downSpacer.append(currDict)

	if(compileIntersection):
		finalDictUp, finalDictDown = {}, {}
		
		for i in range(30):
			
			finalDictUp[i+1], finalDictDown[i+1] = set(), set()
			finalDictDown[i+1] = set(downSpacer[0][i+1][1]) & set(downSpacer[1][i+1][1]) & set(downSpacer[2][i+1][1])
			finalDictUp[i+1]   = set(upSpacer[0][i+1][1]) & set(upSpacer[1][i+1][1]) & set(upSpacer[2][i+1][1])
			
		out = open('results/commonGenes/upRegulatedGenes.txt','w')
		for i in range(30):
			out.write("\nSpacer " + str(i+1) + ": ")
			out.write(str(finalDictUp[i+1]))
			out.write('\n')

		out = open('results/commonGenes/downRegulatedGenes.txt','w')
		for i in range(30):
			out.write("\nSpacer " + str(i+1) + ": ")
			out.write(str(finalDictDown[i+1]))
			out.write('\n')