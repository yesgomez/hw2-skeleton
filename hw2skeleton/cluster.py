from .utils import Atom, Residue, ActiveSite
import sys #, time
import numpy as np
from rdkit import DataStructs, Chem
from rdkit.Chem.Fingerprints import FingerprintMols
from operator import itemgetter
from scipy.cluster.hierarchy import dendrogram, linkage, cophenet, fcluster
from scipy.spatial.distance import pdist
from matplotlib import pyplot as plt
from Bio.PDB import *
from difflib import SequenceMatcher
from .kmeans import Point, Centroid, Kmeans, makeRandomPoint
from sklearn import metrics
# from sklearn.metrics import pairwise_distances

def make_fp(pdb):
	# converts given pdb to mol object for rdkit use
	site = Chem.rdmolfiles.MolFromPDBFile(pdb, sanitize=False, removeHs=False)
	fp = FingerprintMols.FingerprintMol(site, minPath=1, maxPath=7, fpSize=1024, bitsPerHash=2, useHs=True, tgtDensity=0.0, minSize=256)
	print (fp)
	return fp


def compute_similarity(site_a, site_b):
	"""
	Compute the similarity between two given ActiveSite instances.

	Input: two ActiveSite instances
	Output: the similarity between them (a floating point number)
	"""
	tanSimilarity = rdkit.DataStructs.FingerprintSimilarity(site_a, site_b, metric=DataStructs.TanimotoSimilarity)
	dicSimilarity = rdkit.DataStructs.FingerprintSimilarity(site_a, site_b, metric=DataStructs.AsymmetricSimilarity)
	# can try [ TanimotoSimilarity, DiceSimilarity, CosineSimilarity, SokalSimilarity, *RusselSimilarity*, RogotGoldbergSimilarity,\
	# AllBitSimilarity, KulczynskiSimilarity, McConnaugheySimilarity, *AsymmetricSimilarity*, BraunBlanquetSimilarity]
	
	return tanSimilarity, dicSimilarity


def cluster_by_partitioning(active_sites):
	"""
	Cluster a given set of ActiveSite instances using a partitioning method.
	Implementing kmeans.py from git@github.com:siddheshk/Faster-Kmeans.git [https://github.com/siddheshk/Faster-Kmeans/blob/master/Code/kmeans.py]
	Input: a list of ActiveSite instances
	Output: a clustering of ActiveSite instances
			(this is really a list of clusters, each of which is list of
			ActiveSite instances)
	"""
	clusters = []
	for k in range(2,52): # intialize with different values of k
		# start = time.time()
		x = Kmeans(k, active_sites, 10, initialCentroids=None)
		# print ("Time taken:",time.time() - start)
		clusters.append((x.error, x))
	bestclusters = min(clusters,key=itemgetter(0)) # choose the value of k that gives the lowest error
	num = bestclusters[1].centroidList
	print ("Lowest error was",bestclusters[0],"with",len(num),"clusters.")
	sitelist = []
	idlist = []
	centlist = []
	for i in range(int(len(num))):
		sitelist.append([])
		centlist.append(num[i].point.coordinates)
	for sites in active_sites:
		j = Kmeans.getCentroid(bestclusters[1],Point(sites,2))
		clusterid = j[0] # 0 is cluster id, 1 is distance from centroid
		distid = j[1]
		sitelist[clusterid].append(sites)
		idlist.append(clusterid)
	for c, csite in enumerate(sitelist): # check for empty clusters
		if not csite:
			sitelist.remove(csite)
			centlist.remove(centlist[c])
	print (len(sitelist),"clusters of varying sizes.", len(centlist))
	return sitelist, idlist, centlist


def cluster_hierarchically(active_sites):
	"""
	Cluster the given set of ActiveSite instances using a hierarchical algorithm.                                                                  

	Input: a list of ActiveSite instances
	Output: a list of clusterings
			(each clustering is a list of lists of Sequence objects)
	"""
	
	# scipy has several dist metrics/linkage methods which I ranked based on c->1 to choose this one:
	Z = linkage(active_sites, 'average', 'euclidean') 

	# Max Distance method
	# max_d = 1.0 # determined by looking at dendrogram
	# clusterids = fcluster(Z, max_d, criterion='distance'x)
	# Known No. of Clusters method
	k = int(sys.argv[4]) 
	clusterids = fcluster(Z, k, criterion='maxclust')
	sitelist = []
	for i in range(max(clusterids)+1): # we add 1 because the cluster ids are NOT zero indexed
		sitelist.append([])
	for index, sites in enumerate(active_sites): # map cluster ids to active_site pairs
		cid = clusterids[index] 
		sitelist[cid].append(sites)
	for c, csite in enumerate(sitelist): # check for empty clusters
		if not csite:
			sitelist.remove(csite)
	print (len(sitelist),"clusters of varying sizes.", clusterids)
	return sitelist, clusterids, Z


def graph_clusters(sl, type):
	plt.figure(figsize=(10, 8))
	for i in range(len(sl)):
		cluster = sl[i]
		print (len(cluster))
		plt.plot(*zip(*cluster), marker='o')
		# for j in range(len(clusters)):
	plt.title("%s Clustering" %(type))
	plt.xlabel('Tanimoto Similarity')
	plt.ylabel('Asymmetric Similarity')
	plt.legend()
	plt.show()


def similar(a, b):
	return SequenceMatcher(None, a, b).ratio()

def sim_metric(files):
	strings = []
	aa = { 'TRP': 'W', 'PHE': 'F', 'TYR': 'Y', 'LEU': 'L', 'ILE': 'I', 'MET': 'M', 'VAL': 'V', 'CYS': 'C', 'ALA': 'A', 'GLY': 'G', 'PRO': 'P', 'SER': 'S', 'THR': 'T', 'ASN': 'N', 'GLN': 'Q', 'HIS': 'H', 'ARG': 'R', 'LYS': 'K', 'ASP': 'D', 'GLU': 'E'}
	p=PDBParser()
	for file in files:
		struct = ""
		structure=p.get_structure('FILE', file)
		for model in structure:
			for residue in model.get_residues():
				if residue.get_resname() in aa:
					struct += str(aa[residue.get_resname()])
		strings.append(struct)
	qualmatrix = []
	for i in range(len(strings)): # compute two distance metrics between all pairs of sequences (n^2)
		for j in range(len(strings)):
			ratio = similar(strings[i], strings[j])
			qualmatrix.append(ratio)
	print(len(qualmatrix))
	return qualmatrix


def third_graph(data, labels, gtype):
	# Using sequence similarity as distance
	avgscore = metrics.silhouette_score(data, labels, metric='euclidean')
	print ("Average silhouette score for %s was %s." %(gtype, avgscore))
	silscores = metrics.silhouette_samples(data, labels, metric='euclidean')

	# z = []
	# for i in range(len(qualmatrix)):
	# 	z.append(qualmatrix[i]**2)
	# cm = matplotlib.cm.get_cmap('RdYlBu')
	plt.figure(figsize=(10, 8))
	plt.scatter(labels, silscores, s=None, marker='o') #, c=labels, cmap=cm)
	plt.xlabel('cluster ID')
	plt.ylabel('Silhouette score')
	plt.title("%s Clustering vs Sequence" %gtype)
	plt.show()

	# return silscores

	
def matrix_graph(centroidarr, gtype):
	ca = np.array(centroidarr)
	if gtype == "kmeans":
		cm = plt.cm.Greys
		title = "Visual Representation of All %s Centroids (in 1024 dim)" %len(centroidarr)
		xl = 'Position of coordinate for given cluster centroid'
		yl = 'Centroid index'
	elif gtype == "clustercomp":
		cm = plt.cm.PuBu
		title = "Comparing Clusters by membership and sequence similarity"
		xl = 'K-means Cluster ID'
		yl = 'Hierarchical Cluster ID'
	else:
		print ("Sorry, not a proper graph type. Try kmeans or clustercomp instead.\n")
	plt.matshow(ca, fignum=100, aspect='auto', cmap=cm)
	plt.title(title)
	plt.xlabel(xl)
	plt.ylabel(yl)
	plt.show()

	
def dendrogram_graph(Z):
	plt.figure()
	dn = dendrogram(Z, above_threshold_color='y',orientation='top', leaf_font_size=7.5)
	plt.title("Dendrogram showing relations between clusters")
	plt.show()
