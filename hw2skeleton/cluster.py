from .utils import Atom, Residue, ActiveSite
import os
import time
import numpy as np
import rdkit
import matplotlib
from rdkit import DataStructs
from rdkit.Chem.Fingerprints import FingerprintMols
from operator import itemgetter
from scipy.cluster.hierarchy import dendrogram, linkage, cophenet, fcluster
from scipy.spatial.distance import pdist
from matplotlib import pyplot as plt
from Bio.PDB import *
from difflib import SequenceMatcher
from .kmeans import Point, Centroid, Kmeans, makeRandomPoint

def make_fp(pdb):
	# converts given pdb to mol object for rdkit use
	site = rdkit.Chem.rdmolfiles.MolFromPDBFile(pdb, sanitize=False, removeHs=False)
	fp = FingerprintMols.FingerprintMol(site)
	print (fp)
	return fp


def compute_similarity(site_a, site_b):
	"""
	Compute the similarity between two given ActiveSite instances.

	Input: two ActiveSite instances
	Output: the similarity between them (a floating point number)
	"""
	tanSimilarity = rdkit.DataStructs.FingerprintSimilarity(site_a, site_b, metric=DataStructs.TanimotoSimilarity)
	dicSimilarity = rdkit.DataStructs.FingerprintSimilarity(site_a, site_b, metric=DataStructs.RusselSimilarity)
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
		start = time.time()
		x = Kmeans(k, active_sites, 10, initialCentroids=None)
		print ("Time taken:",time.time() - start)
		clusters.append((x.error, x))
	bestclusters = min(clusters,key=itemgetter(0)) # choose the value of k that gives the lowest error
	num = len(bestclusters[1].centroidList)
	print ("Lowest error was",bestclusters[0],"with",num,"clusters.")
	sitelist = []
	for i in range(int(num)):
		sitelist.append([])
	for sites in active_sites:
		j = Kmeans.getCentroid(bestclusters[1],Point(sites,2))
		clusterid = j[0] # 0 is cluster id, 1 is distance from centroid
		sitelist[clusterid].append(sites)
	print (len(sitelist),"clusters of varying sizes:", len(sitelist[0]),"~", len(sitelist[-1]))
	return sitelist


def cluster_hierarchically(active_sites):
	"""
	Cluster the given set of ActiveSite instances using a hierarchical algorithm.                                                                  

	Input: a list of ActiveSite instances
	Output: a list of clusterings
			(each clustering is a list of lists of Sequence objects)
	"""
	Z = linkage(active_sites, 'average', 'euclidean') # scipy has several dist metrics/linkage methods which I compared based on c->1
 
	# metrics = ['braycurtis', 'canberra', 'chebyshev', 'cityblock', 'cosine', 'dice', 'euclidean', 'hamming', 'jaccard', 'kulsinski', 'mahalanobis', 'matching', 'minkowski', 'rogerstanimoto', 'russellrao', 'seuclidean', 'sokalmichener', 'sokalsneath', 'sqeuclidean', 'yule']
	# methods = ['single', 'complete', 'average', 'weighted', 'ward']
	# for method in methods:
	#     Z = linkage(active_sites, method)
	#     c, coph_dists = cophenet(Z, pdist(active_sites))
	#     print (method, c)
	# for i in range(len(Z)):
	#     j = 0
	#     if Z[i][2] > j:
	#         j = Z[i][2]
	#         print (j)
	#     else:
	#         j = j

	max_d = 0.1 # determined by manually comparing delta(distances)
	clusterids = fcluster(Z, max_d, criterion='distance')
	sitelist = []
	for i in range(max(clusterids)+1): # we add 1 because the cluster ids are NOT zero indexed
		sitelist.append([])
	# for sites in active_sites:
	for index, sites in enumerate(active_sites): # map cluster ids to active_site pairs
		cid = clusterids[index] 
		sitelist[cid].append(sites)
	print (len(sitelist),"clusters of varying sizes:", len(sitelist[0]),"~", len(sitelist[-1]))
	return sitelist


def graph_clusters(sl, type):
	plt.figure(figsize=(10, 8))
	for i in range(len(sl)):
		cluster = sl[i]
		plt.plot(*zip(*cluster), marker='o')
		# for j in range(len(clusters)):
	plt.title("%s Clustering" %(type))
	plt.xlabel('Tanimoto Similarity')
	plt.ylabel('Russel Similarity')
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
				# print (residue)
				if residue.get_resname() in aa:
					struct += str(aa[residue.get_resname()])
		strings.append(struct)
	qualmatrix = []
	for i in range(len(strings)): # compute two distance metrics between all pairs of sequences (n^2)
		for j in range(i, len(strings)):
			ratio = similar(strings[i], strings[j])
			qualmatrix.append(ratio)
	print(len(qualmatrix))


def third_graph(sl, type, qualmatrix):
	cm = matplotlib.cm.get_cmap('RdYlBu')
	colors= qualmatrix
	# xs, ys = zip(*sl)
	plt.figure(figsize=(10, 8))
	for i in range(len(sl)):
		cluster = sl[i]
		xs, ys = zip(*cluster)
		plt.Axes3D.scatter(xs, ys, zs=qualmatrix, zdir='z', s=20, c=colors)
		# plt.plot(*zip(*cluster))
	plt.colorbar(sc)
	plt.title("%s Clustering" %(type))
	plt.xlabel('Tanimoto Similarity')
	plt.ylabel('Russel Similarity')
	plt.show()
