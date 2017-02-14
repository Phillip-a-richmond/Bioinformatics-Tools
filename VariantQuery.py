#!/opt/tools/miniconda/bin/python
import argparse
import os
import math
import pysam
parser = argparse.ArgumentParser()
parser.add_argument("-bam",help="input your bam file",type=str)
parser.add_argument("-bed",help="input the bed file of snps",type=str)
parser.add_argument("-o",help="name of the ouput file",type=str)
parser.add_argument("-ID",help="name of the sample ID",type=str)
args = parser.parse_args()


query = open(args.bed,'r')
outFile = open(args.o,'w')
sampleID = args.ID

class Variant:
	def __init__(self,chromosome,position,samfile):
		self.chromosome = chromosome
		self.position = position
		self.samfile = samfile
	def getInfo(self):
		### Here we are utilizing the pysam pileupcolumn object which is filled with pileupread objects which are filled with aligment objects
		self.nucDict = {'A':0,'C':0,'T':0,'G':0,'N':0}
		self.coverage = 0
		self.indelReads = 0
		#print '#################newVar############################'
		for pileupcolumn in self.samfile.pileup(self.chromosome,self.position+1):
			if pileupcolumn.pos == (self.position-1):
				for pileUpRead in pileupcolumn.pileups:
					if pileUpRead.is_del == 1:
						self.indelReads = self.indelReads + 1
					else:
						self.nucDict[pileUpRead.alignment.seq[pileUpRead.query_position]] = self.nucDict[pileUpRead.alignment.seq[pileUpRead.query_position]] + 1
						self.coverage = self.coverage + 1
				self.ann = self.chromosome+'\t'+str(self.position)+'\t'+'Coverage:'+str(self.coverage)+'\t'
				for i in self.nucDict:
					self.ann = self.ann+i+":"+str(self.nucDict[i])+'\t'
				self.ann = self.ann+"INDELReads:"+str(self.indelReads)
				break
			elif pileupcolumn.pos > (self.position + 10):
				break
		if self.coverage == 0:
			self.ann = self.chromosome+'\t'+str(self.position)+'\t'+'Coverage:'+str(self.coverage)+'\t'+'A:0'+'\t'+'C:0'+'\t'+'T:0'+'\t'+'G:0'+'\t'+'N:0'+'\t'+'INDELReads:'+str(self.indelReads)


samfile = pysam.Samfile(args.bam,"rb")
#ignore header
query.readline()
outFile.write('%s Coverage (wt|variant)\n'%sampleID)
for line in query:
	line = line.strip('\r')
	line = line.strip('\n')
	if line[0] != '#':
		columns = line.split('\t')
		chromosome = columns[7]
		varPos = int(columns[8])
		ref = columns[10]
		alt = columns [11]
		varMan = Variant(chromosome,varPos,samfile)
		varMan.getInfo()
		if (len(ref) > 1) or (len(alt) > 1):
			wt='.'
			variant='.'
			outFile.write("%s|%s\n"%(wt,variant))	
		else:
			wt = varMan.nucDict[ref]
			variant = varMan.nucDict[alt]
			outFile.write("%d|%d\n"%(wt,variant))	
