from .utils import Atom, Residue, ActiveSite
from prody import *
import os
import rdkit
from rdkit import DataStructs
from rdkit.Chem.Fingerprints import FingerprintMols

def make_fp(pdb):
    # converts given pdb to mol object for rdkit use
    site = rdkit.Chem.rdmolfiles.MolFromPDBFile(pdb, sanitize=False, removeHs=False)
    print (site)
    # fp = rdkit.Chem.Fingerprints.FingerprintMols.FingerprintsFromMols(site)
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
    print (similarity)
    # similarity = 0.0
    return similarity


def cluster_by_partitioning(active_sites):
    """
    Cluster a given set of ActiveSite instances using a partitioning method.

    Input: a list of ActiveSite instances
    Output: a clustering of ActiveSite instances
            (this is really a list of clusters, each of which is list of
            ActiveSite instances)
    """
    # Fill in your code here!

    return []


def cluster_hierarchically(active_sites):
    """
    Cluster the given set of ActiveSite instances using a hierarchical algorithm.                                                                  #

    Input: a list of ActiveSite instances
    Output: a list of clusterings
            (each clustering is a list of lists of Sequence objects)
    """

    # Fill in your code here!

    return []

