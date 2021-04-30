'''
Code to perform variant calling. 
'''

import numpy as np
import math
import sys
import getopt




# helper dictionaries mapping bases on number 0-3 and vice versa
index = {'A':0, 'C':1, 'G':2, 'T':3, 'a':0, 'c':1, 'g':2, 't':3}
inverseIndex = {0:'A', 1:'C', 2:'G', 3:'T'}

'''
Function to find the most likely genotype, if we observe bases counts nBases (array of length 4: number
of As, Cs, etc.) , the prior on heterozygous genotype is heteroPrior and the sequencing error rate theta.
'''
def mostLikelyGenotype(nBases, heteroPrior, theta):
	logTheta = math.log(theta/3)
	logOneMinusTheta = math.log(1-theta)
	logHalfMinusTheta = math.log(0.5-theta/3)

	data = np.sort(nBases)
	# coverage
	cov = np.sum(data)
	# we require coverage at least 8
	if cov>=9:
		# probability of the most likely homozygous genotype 
		logProb_homo = data[3]*logOneMinusTheta + (cov-data[3])*logTheta 
		# probability of the most likely heterozygous genotype
		logProb_hetero = (data[2]+data[3])*logHalfMinusTheta + (data[0]+data[1])*logTheta + math.log(heteroPrior)
			
		# if logProb_homo == logProb_hetero, we are not able to decide
		if logProb_homo == logProb_hetero:
			return "nan"
		# homozygous genotype
		elif logProb_homo > logProb_hetero:
			# check that the genotype is unique
			if data[2]==data[3]:
				return "nan"
			# return the most likely base
			else:
				base = inverseIndex[np.argmax(nBases)]
				return base+base
		# heterozygous genotype
		else:
			# check uniqueness
			if data[2]==data[1]:
				return "nan"
			else:
				genotype = ""
				for i in range(4):
					if nBases[i]==data[2] or nBases[i]==data[3]:
						genotype = genotype+inverseIndex[i]
				return genotype
	else:
		return "nan"


'''
# weights from the iterative refinement of cluster assignments and centers
with open("iterative_assignments_cluster0_final", 'r') as f:
	weights_iterative_0 = [float(x) for x in f.read().split(',')]
weights_iterative_1 = [1-x for x in weights_iterative_0]
discardedCells = sum([1 for x in weights_iterative_0 if (x<0.95 and x>0.05)])
print("Discarded cells: " + str(discardedCells))
for i in range(n_cells):
	if weights_iterative_0[i]<0.95 and weights_iterative_0[i]>0.05:
		print(weights_iterative_0[i])
'''



def varCalling(mpileupFile, labels, theta, heteroPrior, ind):
	# file with results
	rf = open("varCalls_chr"+str(ind), 'w')

	# read the mpileupFile line by line and call genotype for each position in the file for both clusters
	mf = open(mpileupFile, 'r') 
	# position of the last mutation
	while True: 
		# Get next line from file
		line = mf.readline()
		# if line is empty, end of file is reached
		if not line:
			break

		# split by tabs
		splitLine = line.split("\t")
		# first column: chromosome
		chrom = splitLine[0]
		# second column: position (1-based!!!)
		position = int(splitLine[1])
		# 3rd column: coverage
		cov = int(splitLine[2])
		# 4th column: bases
		bases = splitLine[3]
		#print(bases)
		# 5th column: cell ids
		cell_ids = splitLine[4].split(",")

		# compute the weighted number of bases
		nBases_0 = [0,0,0,0]
		nBases_1 = [0,0,0,0]
		
		
		for i in range(cov):
			cell_id = int(int(cell_ids[i])/4)
			# iterative
			#w0 = weights_iterative_0[cell_id]
			w0 = 1 if labels[cell_id]<=0.05 else 0
			nBases_0[index[bases[i]]] += w0
			#w1 = weights_iterative_1[cell_id]
			w1 = 1 if labels[cell_id]>=0.95 else 0
			nBases_1[index[bases[i]]] += w1
		
			
			
		### iterative
		genotype_0 = mostLikelyGenotype(nBases_0, heteroPrior, theta)
		genotype_1 = mostLikelyGenotype(nBases_1, heteroPrior, theta)


		################### output
		# output
		if genotype_0!=genotype_1 and genotype_0!="nan" and genotype_1!="nan":
			toWrite = [chrom, str(position), str(nBases_0), genotype_0, str(nBases_1), genotype_1]
			rf.write('\t'.join(toWrite)+"\n")
		
	mf.close() 
	rf.close()


def main(argv):
	heteroPrior = 0.001
	theta=0.01
	# labels file
	labelsFile = ""
	
	# parse the arguments
	try:
		opts, args = getopt.getopt(argv,"l:f:t:h:i:")
	except getopt.GetoptError:
		sys.exit(2)
	for opt, arg in opts:
		if opt == '-l':
			labelsFile = arg
		if opt == '-f':
			mpileupFile = arg
		if opt == '-t':
			theta = arg
		if opt == '-h':
			heteroPrior = arg
		if opt == '-i':
			ind = arg
		
	with open(labelsFile, 'r') as f:
		labels = f.read().split()
	labels = [float(x) for x in labels]
	
	varCalling(mpileupFile, labels, theta, heteroPrior, ind)



if __name__ == "__main__":
    main(sys.argv[1:])