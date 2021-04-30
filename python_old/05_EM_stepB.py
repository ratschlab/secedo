import numpy as np
import sys
import getopt
import math

def EM_stepB(labels):
	n_cells=len(labels)

	# probability of being in cluster 0 and in cluster 1
	prob_cluster_1 = [float(x) for x in labels]
	prob_cluster_0 = [1-x for x in prob_cluster_1]
	# compute the priors
	prior_0 = sum(prob_cluster_0)/n_cells
	prior_1 = sum(prob_cluster_1)/n_cells

	loglikelihood_0 = np.zeros(shape=n_cells)
	loglikelihood_1 = np.zeros(shape=n_cells)

	for chrom in [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,"X"]:
		l0 = np.genfromtxt("loglikelihood_0_"+str(chrom), delimiter=",")
		loglikelihood_0 += l0
		l1 = np.genfromtxt("loglikelihood_1_"+str(chrom), delimiter=",")
		loglikelihood_1 += l1

	np.savetxt("loglikelihood_0", loglikelihood_0, delimiter=",")
	np.savetxt("loglikelihood_1", loglikelihood_1, delimiter=",")

	odds_10 = [math.exp(min(max(loglikelihood_1[x]-loglikelihood_0[x],-100),100)) for x in range(n_cells)]
	new_prob_cluster_0 = [1/(1+odds_10[x]*prior_1/prior_0) for x in range(n_cells)]
	new_prob_cluster_1 = [1-x for x in new_prob_cluster_0]

	with open("new_labels", 'w') as f:
		f.write(','.join([str(x) for x in new_prob_cluster_1]))



######################

def main(argv):
	# labels file
	labelsFile = ""
	
	# parse the arguments
	try:
		opts, args = getopt.getopt(argv,"l:")
	except getopt.GetoptError:
		sys.exit(2)
	for opt, arg in opts:
		if opt == '-l':
			labelsFile = arg
		
	with open(labelsFile, 'r') as f:
		labels = f.read().split()
	
	EM_stepB(labels)



if __name__ == "__main__":
    main(sys.argv[1:])