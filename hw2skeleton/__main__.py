import sys
import glob
import numpy as np
from rdkit import DataStructs
from operator import itemgetter
from matplotlib import pyplot as plt
from scipy.cluster.hierarchy import dendrogram
from .io import write_clustering, read_active_sites
from .cluster import cluster_by_partitioning, cluster_hierarchically, make_fp, sim_metric, third_graph, matrix_graph, dendrogram_graph
from .kmeans import Point, Centroid, Kmeans, makeRandomPoint


# Some quick stuff to make sure the program is called correctly
if len(sys.argv) < 4:
	print("Usage: python -m hw2skeleton [-P| -H] <pdb directory> <output file> [# of P clusters for H]")
	sys.exit(0)


# Read in active site pdbs, generate fingerprints, and return those fp as numpy vectors	
activesites, np_fps = read_active_sites(sys.argv[2])

# Quality metric generation / testing
seqSimList = sim_metric(activesites)
seqSimMatrix = np.reshape(seqSimList, (136, 136))


# Choose clustering algorithm
clusters = []

if sys.argv[1][0:2] == '-P':
	print("Clustering using Partitioning method")  # employing k-means as shown in kmeans.py
	clusters, ids, centroids = cluster_by_partitioning(np_fps)
	matrix_graph(centroids, 'kmeans')
	third_graph(seqSimMatrix, ids, 'partitioning')
	write_clustering(sys.argv[3], clusters, ids)

if sys.argv[1][0:2] == '-H':
	print("Clustering using Hierarchical method") # employing the scipy average linkage
	clusters, ids, Z = cluster_hierarchically(np_fps)
	dendrogram_graph(Z)
	third_graph(seqSimMatrix, ids, 'hierarchical')
	write_clustering(sys.argv[3], clusters, list(ids))

if sys.argv[1][0:2] == '-C':
	print("Generating cluster comparison graphs")
	hids = read_idarr("smol_h.txt")
	pids = read_idarr("smol_p.txt")
	zip(pids, hids)
	
	matrix_graph(clustids, 'clustercomp')
	