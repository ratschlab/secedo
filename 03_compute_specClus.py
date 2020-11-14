import numpy as np
import math
import scipy.special
from scipy.sparse.csgraph import laplacian as csgraph_laplacian
import random
from sklearn.mixture import GaussianMixture
from sklearn.cluster import SpectralClustering
from sklearn.cluster import KMeans



def compute_AIC_BIC(simMat, index):
	# compute graph Laplacian,the eigenvalues and eigenvectors
	L = csgraph_laplacian(simMat, normed=True)
	eigenvalues, eigenvectors = np.linalg.eig(L)
	idx = eigenvalues.argsort()
	eigenvalues = eigenvalues[idx]
	eigenvectors = eigenvectors[:,idx]
	# save the first 20 eigenvalues
	with open("eigenvalues_"+str(index), 'w') as f:
		f.write(','.join([str(x) for x in eigenvalues[:20]]))
	# print the second and third smallest eigenVectors
	with open("eigenvectors_"+str(index), 'w') as rf:
		for i in range(simMat.shape[0]):
			rf.write(str(eigenvectors[i,1]) + ","+str(eigenvectors[i,2])+"\n")

	# decide if it seems there are two different clusters
	# coordinates of cells ... value of the first eigenvector
	# we assume the coordinates for cells within one cluster are normally distributed
	# -> fit a Gaussian Mixture Model with 1 component and 2 components, compare using AIC, BIC
	cell_coord = eigenvectors[:,[1,2,3,4,5]]
	#cell_coord = eigenvectors[:,1].reshape(-1,1)
	# GMM with 1 component
	gmm_1 = GaussianMixture(n_components=1, n_init=1, random_state=29)
	gmm_1.fit(cell_coord)
	aic_1 = gmm_1.aic(cell_coord)
	bic_1 = gmm_1.bic(cell_coord)
	
	# GMM with 2 components
	gmm_2 = GaussianMixture(n_components=2, n_init=3, random_state=29)
	gmm_2.fit(cell_coord)
	aic_2 = gmm_2.aic(cell_coord)
	bic_2 = gmm_2.bic(cell_coord)

	# GMM with 2 components and tied variance
	gmm_2_tied = GaussianMixture(n_components=2, n_init=3, random_state=29, covariance_type='tied')
	gmm_2_tied.fit(cell_coord)
	aic_2_tied = gmm_2_tied.aic(cell_coord)
	bic_2_tied = gmm_2_tied.bic(cell_coord)


	# k-means on the first five eigenvectors (excluding the trivial one)
	#X = eigenvectors[:,[1,2,3,4,5]]
	#print(X.shape)
	#kmeans = KMeans(n_clusters=2, random_state=0).fit(X)
	#np.savetxt("labels_specClus_5eigen.csv", kmeans.labels_, delimiter=',')
	# return all the AIC's and BIC's
	return aic_1, aic_2, aic_2_tied, bic_1, bic_2, bic_2_tied


def spectralClustering(simMat):	
	sc = SpectralClustering(2, affinity='precomputed', n_init=100, random_state=15)
	sc.fit(simMat)
	#with open("labels", 'w') as f:
	#	f.write(','.join([str(x) for x in sc.labels_]))
	return sc.labels_


def main(argv):
	# parse the arguments
	try:
		opts, args = getopt.getopt(argv,"n:i:")
	except getopt.GetoptError:
		sys.exit(2)
	for opt, arg in opts:
		if opt == '-n':
			n_cells = int(arg)
		elif opt == '-t':
			identifiersFile = arg

	with open(identifiersFile, 'r') as iF:
		identifiers = iF.read().split(",")

	mat_same = np.zeros(shape=(n_cells, n_cells))
	mat_diff = np.zeros(shape=(n_cells, n_cells))


	# read the partial matrices and sum them all together
	for i in identifiers:
		print(i)
		mat_same_i = np.genfromtxt("mat_same_"+str(i)+".csv", delimiter=',')
		mat_diff_i = np.genfromtxt("mat_diff_"+str(i)+".csv", delimiter=',')

		mat_same = mat_same + mat_same_i
		mat_diff = mat_diff + mat_diff_i


	# compute log(P(diff)/P(same))
	# i.e., simMat_diff[i,j] = -w(i,j) for w(i,j) as defined in the draft: https://www.overleaf.com/3934821935gjttfmhfzkyf
	simMat = mat_diff-mat_same
	np.savetxt('simMat_diff.csv', simMat, delimiter=',')


	### option 1: add the abs. value of the min. element
	simMat = -1*simMat_diff
	minEl = np.min(simMat)
	simMat = simMat + abs(minEl)
	np.savetxt('simMat_final_diff.csv', simMat, delimiter=',')

	'''
	### option 2: exponentiate (equivalent to Equation 5.6 in the thesis)
	simMat = 1/(1+np.exp(simMat_diff))
	np.fill_diagonal(simMat, 1)
	np.savetxt('simMat_final_exp.csv', simMat, delimiter=',')
	'''

	'''
	### option 3: scale the exponentiated version so that max. element (excluding diagonal) is 1
	simMat = np.genfromtxt("simMat_final_exp.csv", delimiter=',')
	np.fill_diagonal(simMat, 0)
	maxEl = np.max(simMat)
	simMat = simMat*1/maxEl
	np.fill_diagonal(simMat, 1)
	np.savetxt("simMat_final_exp_scaled.csv", simMat, delimiter=",")
	'''

	'''
	### option 4: take simMat_diff, scale is so that mean of elements (excluding diagonal) = 0,
	### then exponentiate
	mask = np.ones(simMat_diff.shape, dtype=bool)
	np.fill_diagonal(mask, 0)
	mean = np.mean(simMat_diff[mask])
	simMat_diff -= mean
	simMat_diff = 1/(1+np.exp(simMat))
	np.fill_diagonal(simMat, 1)
	np.savetxt('simMat_mean0_exp.csv', simMat, delimiter=',')
	'''

	aicbic_all = compute_AIC_BIC(simMat, "final")
	np.savetxt("aic_bic_final", aicbic_all, delimiter=",")


	labels = spectralClustering(simMat)
	np.savetxt("labels_specClus.csv", labels, delimiter=',')


if __name__ == "__main__":
	main(sys.argv[1:])