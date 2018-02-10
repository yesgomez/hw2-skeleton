import sys
import glob
from operator import itemgetter
from matplotlib import pyplot as plt
from scipy.cluster.hierarchy import dendrogram
from .io import read_active_sites, write_clustering, write_mult_clusterings
from .cluster import cluster_by_partitioning, cluster_hierarchically, compute_similarity, make_fp, graph_clusters
from .kmeans import Point, Centroid, Kmeans, makeRandomPoint


# Some quick stuff to make sure the program is called correctly
if len(sys.argv) < 4:
	print("Usage: python -m hw2skeleton [-P| -H] <pdb directory> <output file>")
	sys.exit(0)

# active_sites = read_active_sites(sys.argv[2])

files = glob.glob('./data/*.pdb') # make list of pdbs
fps = []
for file in files: # make fingerprints from pdb list
	fp = make_fp(file)
	fps.append(fp)
print (len(fps))

simmatrix = []
for i in range(len(fps)): # compute two distance metrics between all pairs of fingerprints (n^2)
	# j = i + 1
	for j in range(i, len(fps)):
		simmatrix.append(list(compute_similarity(fps[i], fps[j])))
print (len(simmatrix)) # 2x9316 matrix of unique (dist1, dist2) for 136 active site pdbs


# Choose clustering algorithm
clusters = []
if sys.argv[1][0:2] == '-P':
	print("Clustering using Partitioning method")  # employing k-means as shown in kmeans.py
	clusters = cluster_by_partitioning(simmatrix)
	write_clustering(sys.argv[3], clusters)
	graph_clusters(clusters, 'Partitioning')

if sys.argv[1][0:2] == '-H':
	print("Clustering using Hierarchical method")
	clusters = cluster_hierarchically(simmatrix)
	write_clustering(sys.argv[3], clusters)
	graph_clusters(clusters, 'Hierarchical')