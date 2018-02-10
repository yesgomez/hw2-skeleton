import sys
import glob
from operator import itemgetter
from .io import read_active_sites, write_clustering, write_mult_clusterings
from .cluster import cluster_by_partitioning, cluster_hierarchically, compute_similarity, make_fp
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
	# simtemp = []
	for j in range(len(fps)):
		simmatrix.append(list(compute_similarity(fps[i], fps[j])))
	# simmatrix.append(simtemp)
	# print (simtemp[0])
print (len(simmatrix)) # 2x18496 matrix of (dist1, dist2) for 136 active site pdbs


# Choose clustering algorithm
clusters = []
if sys.argv[1][0:2] == '-P':
	print("Clustering using Partitioning method")  # employing k-means as shown in kmeans.py
	clusters = cluster_by_partitioning(simmatrix)
	# write_clustering(sys.argv[3], clusters)

if sys.argv[1][0:2] == '-H':
	print("Clustering using Hierarchical method")
	clusterings = cluster_hierarchically(active_sites)
	write_mult_clusterings(sys.argv[3], clusterings)
