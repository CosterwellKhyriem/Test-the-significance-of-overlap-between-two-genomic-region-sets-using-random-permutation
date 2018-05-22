#! /usr/bin/python:wq

import pybedtools as py
import sys 
import numpy as np
from scipy import stats
import argparse
#import matplotlib #for windows system   
#matplotlib.use('Agg') #for windows system
import matplotlib.pyplot as plt
import random

def getUniqueChromosomeList(chrList):
	inList = set()
	return [x for x in chrList if x not in inList and not inList.add(x)]
	
def suffleRegions(chromosome,start,end,chromSizes):
	chromosome = np.asarray(chromosome)
	start = np.asarray(start)
	end = np.asarray(end)
	for chrom in uniqueSet:
		indexOfChromosomes = np.where(chromosome == chrom)
		chromStart = start[indexOfChromosomes]
		chromEnd = end[indexOfChromosomes]
		randomNumbers = [] 
		randomBinaries = []
		listofGenomicLocations = []
		for index in range(len(chromosome[indexOfChromosomes])):
			randomNumber=random.randint(0,(chromSizes[chrom])/2)
			randomBinary = random.randint(0,1)
			if randomBinary  == 0:
				if chromStart[index]-randomNumber < 0:
					listofGenomicLocations.append([chrom,randomNumber-chromStart[index],randomNumber-chromStart[index]])
					#print chrom+'\t'+str(randomNumber-chromStart[index])+'\t'+str(randomNumber-chromStart[index])
				else:
					listofGenomicLocations.append([chrom,chromStart[index]-randomNumber,chromEnd[index]-randomNumber])
					#print chrom+'\t'+str(chromStart[index]-randomNumber)+'\t'+str(chromEnd[index]-randomNumber)
			
			else:
				listofGenomicLocations.append([chrom,chromStart[index]+randomNumber,chromEnd[index]+randomNumber])
				#print "chr1"+'\t'+str(chromStart[index]+randomNumber)+'\t'+str(chromEnd[index]+randomNumber)
			randomBinaries.append(randomBinary)
			randomNumbers.append(randomNumber)
	return py.BedTool(listofGenomicLocations)

parser = argparse.ArgumentParser()
parser.add_argument('--a',dest='bedA',required=True)
parser.add_argument('--b',dest='bedB',required=True)
parser.add_argument('--n',dest='num_of_iterations',required=True)
parser.add_argument('--g',dest='genome',required=True)
parser.add_argument('--o',dest='output',required=True)


if len(sys.argv)==1:
   parser.print_help(sys.stderr)
   sys.exit(1)

information = parser.parse_args()

bedFileA=py.BedTool(information.bedA)
bedFileB=py.BedTool(information.bedB)

		
chromosome = []
start = []
end = []
 
with open(information.bedA) as bedFile:
	for lines in bedFile:
		row = lines.split('\t')
		chromosome.append(str(row[0]))
		start.append(int(row[1]))
		end.append(int(row[2]))
uniqueSet = getUniqueChromosomeList(chromosome)


fileName = information.genome+".chromInfo.txt"


effectiveGenomeSizes = {}
effectiveGenomeSizes['mm9'] = 2150570000.00
effectiveGenomeSizes['hg19'] = 2451960000.00
effectiveGenomeSizes['dm3'] = 121400000.00
effectiveGenomeSizes['ce1'] = 93260000.00

if information.genome == 'mm9' or information.genome == 'hg19' or information.genome == 'dm3' or information.genome == 'ce1':
	effectiveGenomeSize = effectiveGenomeSizes[information.genome]
	chromSizes = {}
	with open(fileName) as chromFile:
		for lines in chromFile:
			row = lines.split('\t')
			chromSizes[str(row[0])] = int(row[1])
else:
	print "Genome not in database....!!! please enter the effective genome size for your genome!!!\n\n"
	sys.exit()



bedAStart = []
bedAEnd = []
bedBStart = []
bedBEnd = []

for coord in bedFileA:
	bedAStart.append(float(coord[1]))
	bedAEnd.append(float(coord[2]))

for coord in bedFileB:
	bedBStart.append(float(coord[1]))
	bedBEnd.append(float(coord[2]))

bedAStart=np.asarray(bedAStart)
bedAEnd=np.asarray(bedAEnd)
bedBStart=np.asarray(bedBStart)
bedBEnd=np.asarray(bedBEnd)

totallengthOfRegionsBedA=sum(bedAEnd-bedAStart)
totallengthOfRegionsBedB=sum(bedBEnd-bedBStart)




listOflocation = []
listForZscore = []

observedOVerlap=bedFileA.intersect(bedFileB, u=True).count()
normalizedOverlap=observedOVerlap*((totallengthOfRegionsBedA-effectiveGenomeSize)/(totallengthOfRegionsBedB-effectiveGenomeSize))

listForZscore.append(normalizedOverlap)
#listForZscore.append(observedOVerlap)

print "Make "+str(information.num_of_iterations)+" iterations"
	
for shuffle in range(0,int(information.num_of_iterations)):
#	print "Random iterations of "+str(shuffle)+" is going on"
	#randomShuffledRegions=bedFileA.shuffle(genome=information.genome,chrom=True)
	randomShuffledRegions=suffleRegions(chromosome,start,end,chromSizes)
	permOverlap=randomShuffledRegions.intersect(bedFileB, u=True).count()
	normalizedPermOverlap=permOverlap*((totallengthOfRegionsBedA-effectiveGenomeSize)/(totallengthOfRegionsBedB-effectiveGenomeSize))
	listForZscore.append(normalizedPermOverlap)
	#listForZscore.append(permOverlap)

print "Iterations complete"

zScore=stats.zscore(np.asarray(listForZscore))	

print stats.norm.sf(zScore[0])
print zScore[0]
plt.hist(zScore,bins=50,color='b',edgecolor='k', alpha=0.65)
plt.axvline(zScore[0],color='r', linestyle='dashed', linewidth=1)

if zScore[0] >= 1.65:
	plt.title("Distribution of overlapping Z-score\nOverlap enrichment over "+str(information.num_of_iterations)+" random iterations is significant\nP value is "+str(stats.norm.sf(zScore[0]))+"\nRed vertical line shows the Z-score of observed overlap")
elif zScore[0] <= -1.65:
	plt.title("Distribution of overlapping Z-score\nOverlap under-enrichment over "+str(information.num_of_iterations)+" random iterations \nP value is "+str(stats.norm.sf(zScore[0]))+"\nRed vertical line shows the Z-score of observed overlap")
else:
	plt.title("Distribution of overlapping Z-score\nOverlap enrichment over "+str(information.num_of_iterations)+" random iteration is not significant\nP value is "+str(stats.norm.sf(zScore[0]))+"\nRed vertical line shows the Z-score of observed overlap")

plt.xlabel("Z-score of Overlaps")
plt.ylabel("Frequency")
fileName = str(information.output)+'.png'
plt.savefig(fileName)
