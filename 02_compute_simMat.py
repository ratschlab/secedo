'''
Code generating the distance matrix from the pre-processed mpileup file.
'''

import getopt
import sys
import numpy as np
import scipy.special
import math





'''
Functions to retrieve or compute the probability of result (x_s, x_d).
'''
def log_prob_of_result(x_s, x_d, epsilon, h, log_probs, combs_xs_xd, p_same_same, p_same_diff, p_diff_same, p_diff_diff):
	if np.isnan(log_probs[x_s, x_d]):
		p=0
		for i in range(x_s+1):
			for j in range(x_d+1):
				for k in range(x_s-i+1):
					for l in range(x_d-j+1):
						if i+j>0:
							p += scipy.special.binom(x_s, i)*scipy.special.binom(x_d, j)*\
							scipy.special.binom(x_s-i,k)*scipy.special.binom(x_d-j,l)*((1-epsilon-h)**(i+j))*\
							1/2*((p_same_same**i)*(p_same_diff**j) + (p_diff_same**i)*(p_diff_diff**j))*\
							(h**(k+l))*(p_same_same**k)*(p_same_diff**l)*(epsilon**(x_s+x_d-i-j-k-l))*\
							(0.5**(x_s+x_d-i-j-k-l))*((p_same_same+p_diff_same)**(x_s-i-k))*\
							((p_same_diff+p_diff_diff)**(x_d-j-l))
						else:
							p += scipy.special.binom(x_s, k)*scipy.special.binom(x_d, l)*(h**(k+l))*\
							(p_same_same**k)*(p_same_diff**l)*(epsilon**(x_s+x_d-k-l))*(0.5**(x_s+x_d-k-l))*\
							((p_same_same+p_diff_same)**(x_s-k))*((p_same_diff+p_diff_diff)**(x_d-l))
		p = scipy.special.binom(x_s+x_d, x_s)*p
		log_probs[x_s, x_d] = math.log(p)
	combs_xs_xd[x_s, x_d] += 1
	return log_probs[x_s, x_d]



# function to compute mat_same (the matrix giving probabilities of cells i and j given they are in the same cluster)
# and mat_diff (prob. of cells i and j given they are in different clusters)
# propHomo ... estimate of how many positions are actually homozygous germline, were only included because of errors
def computeSimilarityMatrix(mpileupFile, n_cells, cells, cell_groups, epsilon, propHomo, index, theta=0.001, readLen=100, maxInsert=1000):
	# arrays with values of already computed probabilities
	# NaN if not yet computed
	log_probs_same = np.empty((2*readLen,2*readLen))
	log_probs_same[:] = np.nan
	log_probs_diff = np.empty((2*readLen,2*readLen))
	log_probs_diff[:] = np.nan
	# count how many times we have seen a given combination of x_s, x_d
	combs_xs_xd = np.zeros((2*readLen,2*readLen))


	# probability that two same letters will be read as different
	p_same_diff = 2*theta*(1-theta)+2*(theta**2)/3
	# probability that two same letters will be read as same
	p_same_same = 1-p_same_diff
	# probability that two different letters will be read as same
	p_diff_same = 2*(1-theta)*theta/3 + theta*theta/3
	# probability that two different letters will be read as different
	p_diff_diff = 1-p_diff_same

	
	# a dictionary holding the active reads
	# key: read id 
	# value: list with four values: 
	#	- read sequence (string)
	#	- starting position (line number) (int)
	#	- cell id (string)
	#	- starting position in real coordinates (int)
	activeReads = {}
	# distance matrices
	mat_same = np.zeros(shape=(n_cells, n_cells))
	mat_diff = np.zeros(shape=(n_cells, n_cells))


	# line counter
	n_lines = 0
	# process the mpileup file line by line
	f = open(mpileupFile, 'r')
	while True:
		line = f.readline()

		if not line: 
			break

		n_lines += 1
		#print(n_lines)
		line = line.rstrip("\n")
		# split by tabs
		splitLine = line.split("\t")
		# second field: position
		line_pos = int(splitLine[1])
		# 4th field: bases
		line_bases = splitLine[3]
		# 5th field: cell identifiers
		line_cells = splitLine[4].split(",")
		# 6th field: read identifiers
		line_reads = splitLine[5].split(",")

		# for each active read, check if it spans the current position
		for r_id in list(activeReads.keys()):
			r_value = activeReads[r_id]
			# the read continues
			if r_id in line_reads:
				# position of the read in line_reads
				tmp_pos = line_reads.index(r_id)
				# the sequenced base
				tmp_base = line_bases[tmp_pos]
				activeReads[r_id] = [r_value[0]+tmp_base, r_value[1], r_value[2], r_value[3]]
			# if we are not more than maxInsert from start of the mate,
			# it is possible the read still continues, only this particular base is unknown 
			# (deletion, or low quality, or we are inbetween the two mates)
			elif r_value[3]>=line_pos-maxInsert:
				activeReads[r_id] = [r_value[0]+'*', r_value[1], r_value[2], r_value[3]]
			# the read does not continue; compute its overlaps with all other active reads
			else:
				cell_id = r_value[2]
				pos = r_value[1]
				seq = r_value[0]
				# delete the read from activeReads
				del activeReads[r_id]
				#print("Delete " + r_id)

				# compare with all other reads in activeReads
				for r_id_2, r_value_2 in activeReads.items():
					if r_id_2==r_id:
						print("The same read id was present twice in activeReads; something is wrong.")
						quit()
					# read is not from the same cell -> count the number of matches and mismatches in the overlap
					cell_id_2 = r_value_2[2]
					if cell_id_2 != cell_id:
						#print(str(pos)+", "+seq)
						#print("Compare " + r_id + " and " + r_id_2)
						pos_2 = r_value_2[1]
						seq_2 = r_value_2[0]
						#print(str(pos_2)+", "+seq_2)
						x_s = 0
						x_d = 0							
						if pos < pos_2:
							diff = pos_2-pos
							for i in range(len(seq_2)):
								# both read and the same
								if seq_2[i]!='*' and seq[i+diff]!='*' and seq_2[i].capitalize()==seq[i+diff].capitalize():
									x_s +=1
								# both read and different
								elif seq_2[i]!='*' and seq[i+diff]!='*':
									x_d +=1
						else:
							diff = pos - pos_2
							for i in range(len(seq)):
								if seq_2[i+diff]!='*' and seq[i]!='*' and seq[i].capitalize()==seq_2[i+diff].capitalize():
									x_s +=1
								elif seq_2[i+diff]!='*' and seq[i]!='*':
									x_d +=1
						#print("x_s=" + str(x_s) + ", x_d=" + str(x_d))
						#print()

						# update the distance matrices
						if x_s+x_d>0:
							index_1 = int(cell_groups[cells.index(cell_id)] )
							index_2 = int(cell_groups[cells.index(cell_id_2)])
							if index_1 != index_2:
								mat_same[index_1, index_2] += log_prob_of_result(x_s, x_d, 0, propHomo+0.5*epsilon, log_probs_same, combs_xs_xd, p_same_same, p_same_diff, p_diff_same, p_diff_diff)
								mat_same[index_2, index_1] = mat_same[index_1, index_2]
								mat_diff[index_1, index_2] += log_prob_of_result(x_s, x_d, epsilon, propHomo, log_probs_diff, combs_xs_xd, p_same_same, p_same_diff, p_diff_same, p_diff_diff)
								mat_diff[index_2, index_1] = mat_diff[index_1, index_2]


		# for each read in line_reads, check, if it is in activeReads; if not, add it
		for i in range(len(line_reads)):
			if line_reads[i] not in activeReads.keys():
				b = line_bases[i]
				#print("Add " + line_reads[i])
				cell_id = int(line_cells[i])
				activeReads[line_reads[i]] = [b, n_lines, cell_id, line_pos]

	f.close()

	# in the end, process all the remaining active reads
	# TODO: THIS IS UGLY :-(
	for r_id in list(activeReads.keys()):
		r_value = activeReads[r_id]
		cell_id = r_value[2]
		pos = r_value[1]
		seq = r_value[0]
		# delete the read from activeReads
		del activeReads[r_id]
		#print("Delete " + r_id)

		# compare with all other reads in activeReads
		for r_id_2, r_value_2 in activeReads.items():
			if r_id_2==r_id:
				print("The same read id was present twice in activeReads; something is wrong.")
				quit()
			# read is not from the same cell -> count the number of matches and mismatches in the overlap
			cell_id_2 = r_value_2[2]
			if cell_id_2 != cell_id:
				#print("Compare " + r_id + " and " + r_id_2)
				pos_2 = r_value_2[1]
				seq_2 = r_value_2[0]
				x_s = 0
				x_d = 0							
				if pos < pos_2:
					diff = pos_2-pos
					for i in range(len(seq_2)):
						# both read and the same
						if seq_2[i]!='*' and seq[i+diff]!='*' and seq_2[i].capitalize()==seq[i+diff].capitalize():
							x_s +=1
						# both read and different
						elif seq_2[i]!='*' and seq[i+diff]!='*':
							x_d +=1
				else:
					diff = pos - pos_2
					for i in range(len(seq)):
						if seq_2[i+diff]!='*' and seq[i]!='*' and seq[i].capitalize()==seq_2[i+diff].capitalize():
							x_s +=1
						elif seq_2[i+diff]!='*' and seq[i]!='*':
							x_d +=1
				#print("x_s=" + str(x_s) + ", x_d=" + str(x_d))

				# update the distance matrices
				if x_s+x_d>0:
					index_1 = int(cell_groups[cells.index(cell_id)] )
					index_2 = int(cell_groups[cells.index(cell_id_2)])
					if index_1 != index_2:
						mat_same[index_1, index_2] += log_prob_of_result(x_s, x_d, 0, propHomo+0.5*epsilon, log_probs_same, combs_xs_xd, p_same_same, p_same_diff, p_diff_same, p_diff_diff)
						mat_same[index_2, index_1] = mat_same[index_1, index_2]
						mat_diff[index_1, index_2] += log_prob_of_result(x_s, x_d, epsilon, propHomo, log_probs_diff, combs_xs_xd, p_same_same, p_same_diff, p_diff_same, p_diff_diff)
						mat_diff[index_2, index_1] = mat_diff[index_1, index_2]



	# save mat_same and mat_diff into file
	np.savetxt('mat_same_' + str(index) + '.csv', mat_same, delimiter=',')
	np.savetxt('mat_diff_' + str(index) + '.csv', mat_diff, delimiter=',')
	np.savetxt('combs_xs_xd_' + str(index) + '.csv', combs_xs_xd, delimiter=',')




######################

def main(argv):
	theta=0.001
	readLen=100
	# number of cells
	n_cells = 0
	# epsilon
	epsilon = 0
	# chromosome
	chromosome = 0
	# mpileupFile
	mpileupFile = ""
	# file with the relevant cell ids
	cellsFile = ""
	cellGroupFile = ""
	# maximal considered inset size (default 1000)
	maxInset = 1000

	# parse the arguments
	try:
		opts, args = getopt.getopt(argv,"n:e:c:f:t:h:i:g:m:")
	except getopt.GetoptError:
		sys.exit(2)
	for opt, arg in opts:
		if opt == '-n':
			n_cells = int(arg)
		elif opt == '-e':
			epsilon = float(arg)
		elif opt == '-h':
			propHomo = float(arg)
		elif opt == '-i':
			chromosome = arg
		elif opt == '-f':
			mpileupFile = arg
		elif opt == '-t':
			theta = float(arg)
		elif opt == '-c':
			cellsFile = arg
		elif opt== '-g':
			cellGroupFile = arg
		elif opt=='-m':
			maxInsert = int(arg)


	if cellsFile != "":
		with open(cellsFile, 'r') as f:
			cells_ids = f.read().split()
		cells_ids = [int(x) for x in cells_ids]
		#n_cells = len(cells_ids)
	else:
		cells_ids = range(n_cells)

	#print(cells_ids)

	if cellGroupFile != "":
		cell_groups = np.genfromtxt(cellGroupFile, delimiter=",")
	else:
		cell_groups = range(n_cells)

	#print(cell_groups)
	


	computeSimilarityMatrix(mpileupFile,  n_cells, cells_ids, cell_groups, epsilon, propHomo, chromosome, theta, readLen, maxInsert)

if __name__ == "__main__":
	main(sys.argv[1:])
	
