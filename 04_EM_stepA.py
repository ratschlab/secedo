import numpy as np
import getopt
import sys
import math




# helper dictionaries mapping bases on number 0-3 and vice versa
index = {'A':0, 'C':1, 'G':2, 'T':3, 'a':0, 'c':1, 'g':2, 't':3}
inverseIndex = {0:'A', 1:'C', 2:'G', 3:'T'}

'''
# function to count number of lines in file
def file_len(fname):
    with open(fname) as f:
        for i, l in enumerate(f):
            pass
    return i + 1
'''

'''
# OLD
def read_cell_data(fname, cells):
	n_pos = file_len(fname)
	n_cells = len(cells)
	data=np.zeros(shape=(n_cells, n_pos, 4))
	# process each line
	g = open(fname, 'r') 
	i=-1
	while True: 
		# Get next line from file
		line = g.readline()
		# if line is empty, end of file is reached
		if not line:
			break
		else:
			i += 1

		# split by tabs
		splitLine = line.split("\t")
		# 4th column: bases
		bases = splitLine[3]
		#print(bases)
		# 5th column: cell ids
		cell_ids = splitLine[4].split(",")
		for j in range(len(bases)):
			cell_id = int(cell_ids[j])
			if cell_id in cells:
				cell_index = cells.index(cell_id)
				base = bases[j]
				data[cell_index][i][index[base]] += 1
	g.close()
	return data
'''


'''
# OLD
def compute_cluster_center(cell_data, prob_cluster_0, theta):
	center = np.zeros(shape=(cell_data.shape[1], cell_data.shape[2]))
	# sum the bases over all cells belonging to the cluster (weighted by prob_cluster_0)
	for i in range(cell_data.shape[0]):
		center += prob_cluster_0[i]*cell_data[i]
	# normalize the composition at each position
	for i in range(center.shape[0]):
		s = sum(center[i])
		if s>0:
			center[i] = [x/s if x/s>theta else theta for x in center[i]]
		else:
			center[i] = [theta]*4
		#print(center[i])
		s = sum(center[i])
		center[i] = [math.log(x) for x in center[i]/s]
	return center
'''

def compute_cluster_center2(cell_data, prob_cluster, theta):
	center = np.zeros(shape=(4))
	# sum the bases over all cells belonging to the cluster (weighted by prob_cluster_0)
	for i in range(cell_data.shape[0]):
		center += prob_cluster[i]*cell_data[i]
	# normalize the composition at each position
	s = sum(center)
	if s>0:
		center = [x/s if x/s>theta else theta for x in center]
	else:
		center = [theta]*4
	s = sum(center)
	center = [math.log(x/s) for x in center]
	return center



def prob_cell_in_cluster(cell_data, cluster_center):
	prob = 0
	for i in range(cell_data.shape[0]):
		for j in range(4):
			prob += cell_data[i][j]*cluster_center[i][j]
	return prob

def log_prob_cell_in_cluster(cell_data, cluster_center):
	prob = 0
	for j in range(4):
		prob += cell_data[j]*cluster_center[j]
	return prob





############ main
def EM_stepA(labels, mpileupFile, ind, theta):
	n_cells=len(labels)

	# probability of being in cluster 0 and in cluster 1
	prob_cluster_1 = [float(x) for x in labels]
	prob_cluster_0 = [1-x for x in prob_cluster_1]

	# recompute the priors
	prior_0 = sum(prob_cluster_0)/n_cells
	prior_1 = sum(prob_cluster_1)/n_cells

	log_likelihood_0 = np.zeros(shape=n_cells)
	log_likelihood_1 = np.zeros(shape=n_cells)


	# process position by position
	g = open(mpileupFile, 'r') 
	while True: 
		currentData = np.zeros(shape=(n_cells, 4))
		# Get next line from file
		line = g.readline()
		# if line is empty, end of file is reached
		if not line:
			break

		### parse the data
		# split by tabs
		splitLine = line.split("\t")
		# 4th column: bases
		bases = splitLine[3]
		#print(bases)
		# 5th column: cell ids
		cell_ids = splitLine[4].split(",")
		for j in range(len(bases)):
			cell_id = int(cell_ids[j])
			cell_index = int(cell_id / 4)
			#print(cell_index)
			base = bases[j]
			#print(base)
			currentData[cell_index][index[base]] += 1

		### compute the cluster centers
		center_0 = compute_cluster_center2(currentData, prob_cluster_0, theta)
		center_1 = compute_cluster_center2(currentData, prob_cluster_1, theta)

		### compute the likelihoods
		log_likelihood_0 += np.array([log_prob_cell_in_cluster(currentData[j], center_0) for j in range(n_cells)]) 
		log_likelihood_1 += np.array([log_prob_cell_in_cluster(currentData[j], center_1) for j in range(n_cells)]) 
		#print(log_likelihood_0)
		#print(log_likelihood_1)

	g.close()	

		
	
	# print the resulting soft assignments to file
	with open("loglikelihood_0_"+str(ind), 'w') as f:
		f.write(','.join([str(x) for x in log_likelihood_0]))
	with open("loglikelihood_1_"+str(ind), 'w') as f:
		f.write(','.join([str(x) for x in log_likelihood_1]))







######################

def main(argv):
	theta=0.001
	# mpileupFile
	mpileupFile = ""
	# labels file
	labelsFile = ""
	# index
	ind = 0

	# parse the arguments
	try:
		opts, args = getopt.getopt(argv,"l:m:i:t:")
	except getopt.GetoptError:
		sys.exit(2)
	for opt, arg in opts:
		if opt == '-l':
			labelsFile = arg
		if opt == '-m':
			mpileupFile = arg
		if opt == '-i':
			ind = arg
		if opt == '-t':
			theta = float(arg)

	with open(labelsFile, 'r') as f:
		labels = f.read().split()


	n_cells = len(labels)
	#cells = [x for x in range(n_cells)]

	
	EM_stepA(labels, mpileupFile, ind, theta)



if __name__ == "__main__":
    main(sys.argv[1:])
