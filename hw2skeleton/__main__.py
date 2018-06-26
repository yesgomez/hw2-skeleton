import sys
import glob
import numpy as np
from rdkit import DataStructs
from operator import itemgetter
from matplotlib import pyplot as plt
from scipy.cluster.hierarchy import dendrogram
from .io import read_active_sites, write_clustering, write_mult_clusterings
from .cluster import cluster_by_partitioning, cluster_hierarchically, compute_similarity, make_fp, graph_clusters, sim_metric, third_graph, matrix_graph, dendrogram_graph
from .kmeans import Point, Centroid, Kmeans, makeRandomPoint


# Some quick stuff to make sure the program is called correctly
if len(sys.argv) < 4:
	print("Usage: python -m hw2skeleton [-P| -H] <pdb directory> <output file> [# of P clusters for H]")
	sys.exit(0)

# make list of pdbs
files = glob.glob('./data/*.pdb') 
fps = []
np_fps = []
# make fingerprints from pdb list
for file in files: 
	fp = make_fp(file)
	fps.append(fp)
# convert the RDKit explicit vectors into numpy arrays
for fp in fps:
  arr = np.zeros((1,))
  DataStructs.ConvertToNumpyArray(fp, arr)
  newa = np.ndarray.tolist(arr)
  np_fps.append(newa)
print (len(fps), len(np_fps), np_fps[0])

# simmatrix = []
# # compute two distance metrics between all pairs of fingerprints (n^2)
# for i in range(len(fps)): 
# 	# j = i + 1
# 	for j in range(i, len(fps)):
# 		simmatrix.append(list(compute_similarity(fps[i], fps[j])))
# print (len(simmatrix)) # 2x9316 matrix of unique (dist1, dist2) for 136 active site pdbs

# Quality metric generation / testing
seqSimMatrix = sim_metric(files)

# Choose clustering algorithm
clusters = []
if sys.argv[1][0:2] == '-P':
	print("Clustering using Partitioning method")  # employing k-means as shown in kmeans.py
	clusters, ids, centroids = cluster_by_partitioning(np_fps)
	matrix_graph(centroids, kmeans)
	# graph_clusters(clusters, 'Partitioning')
	# third_graph(ids, 'Partitioning', seqSimMatrix)
	write_clustering(sys.argv[3], clusters, ids)

if sys.argv[1][0:2] == '-H':
	print("Clustering using Hierarchical method")
	clusters, ids, Z = cluster_hierarchically(np_fps)
	dendrogram_graph(Z)
	# graph_clusters(clusters, 'Hierarchical')
	# third_graph(ids, 'Hierarchical', seqSimMatrix)
	write_clustering(sys.argv[3], clusters, list(ids))

if sys.argv[1][0:2] == '-H':
	print("Generating cluster comparison graphs")
	matrix_graph(centroids, kmeans)
	