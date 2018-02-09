from .utils import Atom, Residue, ActiveSite
from prody import *
import os
import rdkit
from rdkit import DataStructs
from rdkit.Chem.Fingerprints import FingerprintMols
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
    similarity = rdkit.DataStructs.FingerprintSimilarity(site_a, site_b, metric=DataStructs.DiceSimilarity)
    # print (similarity)
    return similarity


def cluster_by_partitioning(k, active_sites):
    """
    Cluster a given set of ActiveSite instances using a partitioning method.
    Implementing kmeans.py from git@github.com:siddheshk/Faster-Kmeans.git [https://github.com/siddheshk/Faster-Kmeans/blob/master/Code/kmeans.py]
    Input: a list of ActiveSite instances
    Output: a clustering of ActiveSite instances
            (this is really a list of clusters, each of which is list of
            ActiveSite instances)
    """
    x = Kmeans(k, active_sites, 10, initialCentroids=None)
    print (x)
    return x

# def cluster_points(X, mu):
#     clusters  = {}
#     for x in X:
#         bestmukey = min([(i[0], np.linalg.norm(x-mu[i[0]])) \
#                     for i in enumerate(mu)], key=lambda t:t[1])[0]
#         try:
#             clusters[bestmukey].append(x)
#         except KeyError:
#             clusters[bestmukey] = [x]
#     return clusters
 
# def reevaluate_centers(mu, clusters):
#     newmu = []
#     keys = sorted(clusters.keys())
#     for k in keys:
#         newmu.append(np.mean(clusters[k], axis = 0))
#     return newmu
 
# def has_converged(mu, oldmu):
#     return (set([tuple(a) for a in mu]) == set([tuple(a) for a in oldmu]))

# def find_centers(X, K):
#     # Initialize to K random centers
#     oldmu = random.sample(X, K)
#     mu = random.sample(X, K)
#     while not has_converged(mu, oldmu):
#         oldmu = mu
#         # Assign all points in X to clusters
#         clusters = cluster_points(X, mu)
#         # Reevaluate centers
#         mu = reevaluate_centers(oldmu, clusters)
#     return(mu, clusters)
    

def cluster_hierarchically(active_sites):
    """
    Cluster the given set of ActiveSite instances using a hierarchical algorithm.                                                                  #

    Input: a list of ActiveSite instances
    Output: a list of clusterings
            (each clustering is a list of lists of Sequence objects)
    """

    return []

