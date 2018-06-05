#! /usr/bin/python:wq

import pybedtools as py
import sys 
import numpy as np
from scipy import stats
import argparse
import matplotlib.pyplot as plt



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

effectiveGenomeSizes = {}
effectiveGenomeSizes['mm9'] = 2150570000.00
effectiveGenomeSizes['hg19'] = 2451960000.00
effectiveGenomeSizes['dm3'] = 121400000.00
effectiveGenomeSizes['ce1'] = 93260000.00

if information.genome == 'mm9' or information.genome == 'hg19' or information.genome == 'dm3' or information.genome == 'ce1':
	effectiveGenomeSize = effectiveGenomeSizes[information.genome]
else:
	print "Genome not in database....!!!\n\n"
	sys.exit()
	

bedFileA=py.BedTool(information.bedA)
bedFileB=py.BedTool(information.bedB)

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
normalizedOverlap=observedOVerlap*((totallengthOfRegionsBedA/effectiveGenomeSize)/(totallengthOfRegionsBedB/effectiveGenomeSize))

listForZscore.append(normalizedOverlap)

print "Make "+str(information.num_of_iterations)+" iterations"
	
for shuffle in range(0,int(information.num_of_iterations)):
	randomShuffledRegions=bedFileA.shuffle(genome=information.genome,chrom=True)
	permOverlap=randomShuffledRegions.intersect(bedFileB, u=True).count()
	normalizedPermOverlap=permOverlap*((totallengthOfRegionsBedA/effectiveGenomeSize)/(totallengthOfRegionsBedB/effectiveGenomeSize))
	listForZscore.append(normalizedPermOverlap)

print "Iterations complete"

zScore=stats.zscore(np.asarray(listForZscore))	

#print stats.norm.sf(zScore[0])
#print zScore[0]
plt.hist(zScore,color='b',edgecolor='k', alpha=0.65,bins=20)
plt.axvline(zScore[0],color='r', linestyle='dashed', linewidth=1)

if zScore[0] >= 1.65:
	plt.title("Overlap enrichment over "+str(information.num_of_iterations)+" random iterations is significant\nP value is "+str(stats.norm.sf(zScore[0]))+"\nRed vertical line shows the Z-score of observed overlap",fontsize=10)
elif zScore[0] <= -1.65:
	plt.title("Overlap under-enrichment over "+str(information.num_of_iterations)+" random iterations \nP value is "+str(stats.norm.sf(zScore[0]))+"\nRed vertical line shows the Z-score of observed overlap",fontsize=10)
else:
	plt.title("Overlap enrichment over "+str(information.num_of_iterations)+" random iteration is not significant\nP value is "+str(stats.norm.sf(zScore[0]))+"\nRed vertical line shows the Z-score of observed overlap",fontsize=10)

plt.xlabel("Z-score of Overlaps")
plt.ylabel("Frequency")


fig_size = plt.rcParams["figure.figsize"]
 

fileName=information.output+".png"
plt.savefig(fileName)
