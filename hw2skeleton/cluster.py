from .utils import Atom, Residue, ActiveSite
from prody import *
import os
import time
import numpy as np
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
    tanSimilarity = rdkit.DataStructs.FingerprintSimilarity(site_a, site_b, metric=DataStructs.TanimotoSimilarity)
    dicSimilarity = rdkit.DataStructs.FingerprintSimilarity(site_a, site_b, metric=DataStructs.DiceSimilarity)
    # print (similarity)
    return tanSimilarity, dicSimilarity


def cluster_by_partitioning(k, active_sites):
    """
    Cluster a given set of ActiveSite instances using a partitioning method.
    Implementing kmeans.py from git@github.com:siddheshk/Faster-Kmeans.git [https://github.com/siddheshk/Faster-Kmeans/blob/master/Code/kmeans.py]
    Input: a list of ActiveSite instances
    Output: a clustering of ActiveSite instances
            (this is really a list of clusters, each of which is list of
            ActiveSite instances)
    """
    start = time.time()
    x = Kmeans(k, active_sites, 10, initialCentroids=None)
    print ("Time taken:",time.time() - start)
    print (len(x.centroidList), x.error)
    print (Point(active_sites[0],2), Kmeans.getCentroid(x,Point(active_sites[0],2)))  
    return x.error, x

def return_clusters(kmeansobj, a_s):
    no = len(kmeansobj.centroidList) # WE WERE ASSSUMING KMEANSOBJ = X BUT IT DOES NOT. WE NEED TO IMPORT X OR JUST COMBINE BOTH FUNCTIONS
    clusters = []
    for i in range(no):
        clusters.append([])
    for sites in a_s:
        j = Kmeans.getCentroid(kmeansobj,Point(sites,2))
        clusterid = j[0]
        clusters[clusterid].append(sites)
    print (len(clusters), len(clusters[-1]), len(clusterid))
    return clusters

def cluster_hierarchically(active_sites):
    """
    Cluster the given set of ActiveSite instances using a hierarchical algorithm.                                                                  #

    Input: a list of ActiveSite instances
    Output: a list of clusterings
            (each clustering is a list of lists of Sequence objects)
    """

    return []

