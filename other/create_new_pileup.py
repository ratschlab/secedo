import re
import math
import numpy as np
import sys
import getopt

# log factorial(n) using Stirling formula
# !!! DOES NOT CHECK IF N IS INTEGER !!!
def log_fact(n):
	if n>170:
		return 0.5*math.log(2*math.pi*n)+n*math.log(n/math.e)
	else:
		return math.log(math.factorial(n))


# function to do preprocessing; Bayesian approach, homoOnly
# returns True -- the position is kept, or False -- the position is discarded
### data is a vector of lenght 4 with counts of A,C,G, and T in the pooled data at a fixed position
def preprocessing(data, theta, K):
	logTheta = math.log(theta/3)
	logOneMinusTheta = math.log(1-theta)
	
	# sort data; sorts in ascending order
	data = np.sort(data)
	# total coverage
	cov = np.sum(data)
	# if there are no reads or if we have just one base, do not keep
	if cov>0 and np.count_nonzero(data)>1:
		# the multinomial coefficient
		logMultinom_coef = log_fact(cov) - log_fact(data[0]) - log_fact(data[1]) - log_fact(data[2]) - log_fact(data[3])
		# log P(data|wt) for wt=homo 
		logProb_homo = logMultinom_coef + data[3]*logOneMinusTheta + (cov-data[3])*logTheta 
		# add prior on the most probable genotype (1/4, because all four homozygous genotypes are equally likely)
		logProb_homo += math.log(1/4) 
		# add prior on null hypothesis (0.998)
		logProb_homo += math.log(0.998)
		# the normalizing coefficient (normalizing for read depth)
		logNorm_coef = log_fact(cov+3) - math.log(6) - log_fact(cov)
		# the result
		res = logNorm_coef+logProb_homo
		if res < K:
			return True
		else:
			return False
	else:
		return False


def createNewMpileup(oldF, newF, cells_ids, theta=0.001):
	Ks = [0.872, 0.872, 0.872, -1.238, -1.358, -5.749, -10.139, -13.912, -25.998, -33.936, -34.469, -54.084, -44.453, -56.842, -63.352, -88.297, -90.711, -96.841, -115.730, -112.601]


	of = open(oldF, 'r')
	nf = open(newF, 'w')
	while True:
		line = of.readline()

		if not line: 
			break

		line = line.rstrip("\n")
		# split by tabs
		splitLine = line.split("\t")
		# first field: chromosome
		chrom = splitLine[0]
		# second field: position
		line_pos = splitLine[1]
		# 4th field: bases
		line_bases = splitLine[3]
		# 5th field: cell identifiers
		line_cells = [int(x) for x in splitLine[4].split(",")]
		# 6th field: read identifiers
		line_reads = splitLine[5].split(",")

		takeThisPos = [True if x in cells_ids else False for x in line_cells]
		bases_new = ''.join([line_bases[x] for x in range(len(line_bases)) if takeThisPos[x]==True])

		### pre-processing
		# count number of As, Cs, Gs, Ts
		nAs = bases_new.count('a') + bases_new.count('A')
		nCs = bases_new.count('c') + bases_new.count('C')
		nGs = bases_new.count('g') + bases_new.count('G')
		nTs = bases_new.count('t') + bases_new.count('T')
		#print([nAs, nCs, nGs, nTs])
		# coverage
		cov_new = nAs + nCs + nGs + nTs
		#print(cov)
		# choose K for the closest coverage
		i = round(cov_new/10)-1
		if i==-1:
			i=0
		elif i>19:
			i=19
		K=Ks[i]
		# do pre-processing
		if preprocessing([nAs, nCs, nGs, nTs], theta, K):
			cells_new = ','.join([str(line_cells[x]) for x in range(len(line_cells)) if takeThisPos[x]==True])
			reads_new = ','.join([line_reads[x] for x in range(len(line_reads)) if takeThisPos[x]==True])
			# write the line to the new file
			line = '\t'.join([chrom, line_pos, str(cov_new), bases_new, cells_new, reads_new])
			nf.write(line+"\n")

	of.close()
	nf.close()


###########################
def main(argv):
	theta=0.001
	oldF=""
	newF=""
	labF=""
	label=""

	# parse the arguments
	try:
		opts, args = getopt.getopt(argv,"n:o:l:L:t:c:")
	except getopt.GetoptError:
		sys.exit(2)
	for opt, arg in opts:
		if opt == '-n':
			newF = arg
		elif opt == '-o':
			oldF = arg
		elif opt == '-l':
			labF = arg
		elif opt == '-L':
			label = float(arg)
		elif opt == '-t':
			theta = float(arg)
		elif opt == '-c':
			cells = arg


	#labels = np.genfromtxt(labF, delimiter=",")
	#cells_ids = [x for x in range(len(labels)) if labels[x]==label]
	cells_ids = np.genfromtxt(cells, delimiter=",")
	#with open("cells_cluster1", 'w') as f:
	#	f.write(",".join([str(x) for x in cells_ids]))	

	createNewMpileup(oldF, newF, cells_ids, theta)

if __name__ == "__main__":
	main(sys.argv[1:])