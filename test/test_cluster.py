from hw2skeleton import cluster
from hw2skeleton import io
from prody import *
from rdkit import Chem
import os

def test_fp_generation():
    filename_a = os.path.join("data", "276.pdb")
    site_a = cluster.make_fp(filename_a)
    assert site_a

def test_similarity():
    filename_a = os.path.join("data", "276.pdb")
    filename_b = os.path.join("data", "4629.pdb")
    # activesite_a = io.read_active_site(filename_a)
    # activesite_b = io.read_active_site(filename_b)
    activesite_a = cluster.make_fp(filename_a)
    activesite_b = cluster.make_fp(filename_b)
    # update this assertion
    assert cluster.compute_similarity(activesite_a, activesite_b)
    # assert cluster.compute_similarity(filename_a, filename_b) == 0.0

def test_partition_clustering():
    # tractable subset
    # pdb_ids = [276, 4629, 10701]

    # active_sites = []
    # for id in pdb_ids:
    #     filepath = os.path.join("data", "%i.pdb"%id)
    #     active_sites.append(io.read_active_site(filepath))

    # # update this assertion
    # assert cluster.cluster_by_partitioning(active_sites) == []

def test_hierarchical_clustering():
    # tractable subset
    # pdb_ids = [276, 4629, 10701]

    # active_sites = []
    # for id in pdb_ids:
    #     filepath = os.path.join("data", "%i.pdb"%id)
    #     active_sites.append(io.read_active_site(filepath))

    # # update this assertion
    # assert cluster.cluster_hierarchically(active_sites) == []
