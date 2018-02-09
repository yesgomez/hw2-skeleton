import sys
import glob
from .io import read_active_sites, write_clustering, write_mult_clusterings
from .cluster import cluster_by_partitioning, cluster_hierarchically, compute_similarity, make_fp
from .kmeans import Point, Centroid, Kmeans, makeRandomPoint


# Some quick stuff to make sure the program is called correctly
if len(sys.argv) < 4:
    print("Usage: python -m hw2skeleton [-P| -H] <pdb directory> <output file>")
    sys.exit(0)


# active_sites = read_active_sites(sys.argv[2])

files = glob.glob('./data/*.pdb')
fps = []
for file in files:
	fp = make_fp(file)
	fps.append(fp)
print (len(fps))

simmatrix = []
for i in range(len(fps)):
	simtemp = []
	for j in range(len(fps)):
		simtemp.append(compute_similarity(fps[i], fps[j]))
	simmatrix.append(simtemp)
print (len(simmatrix))

# Choose clustering algorithm
if sys.argv[1][0:2] == '-P':
    print("Clustering using Partitioning method")
    # clustering = cluster_by_partitioning(active_sites)
    clustering = cluster_by_partitioning(5, simmatrix)
    # write_clustering(sys.argv[3], clustering)

if sys.argv[1][0:2] == '-H':
    print("Clustering using Hierarchical method")
    clusterings = cluster_hierarchically(active_sites)
    write_mult_clusterings(sys.argv[3], clusterings)
