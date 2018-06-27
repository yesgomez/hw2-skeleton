import numpy as np
import glob
import os
from .utils import Atom, Residue, ActiveSite
from .cluster import make_fp
from rdkit import DataStructs


def read_active_sites(dir):
	"""
	Read in all of the active sites from the given directory.

	Input: directory
	Output: list of ActiveSite instances
	"""
	# make list of pdbs
	files = glob.glob("./%s*.pdb" %dir) 

	# make fingerprints from pdb list
	fingerprints = []
	
	# iterate over each .pdb file in the given directory
	for file in files: 
		fingerprint = make_fp(file)
		fingerprints.append(fingerprint)

	print("Read in %d active sites"%len(fingerprints))

	
	# convert the RDKit explicit vectors into numpy arrays
	numpy_fps = []
	for fp in fingerprints:
	  arr = np.zeros((1,))
	  DataStructs.ConvertToNumpyArray(fp, arr)
	  newa = np.ndarray.tolist(arr)
	  numpy_fps.append(newa)
		
	print (len(numpy_fps), numpy_fps[0])
	return files, numpy_fps


def read_active_site(filepath):
	"""
	Read in a single active site given a PDB file

	Input: PDB file path
	Output: ActiveSite instance
	"""
	basename = os.path.basename(filepath)
	name = os.path.splitext(basename)

	if name[1] != ".pdb":
		raise IOError("%s is not a PDB file"%filepath)

	active_site = ActiveSite(name[0])

	r_num = 0

	# open pdb file
	with open(filepath, "r") as f:
		# iterate over each line in the file
		for line in f:
			if line[0:3] != 'TER':
				# read in an atom
				atom_type = line[13:17].strip()
				x_coord = float(line[30:38])
				y_coord = float(line[38:46])
				z_coord = float(line[46:54])
				atom = Atom(atom_type)
				atom.coords = (x_coord, y_coord, z_coord)

				residue_type = line[17:20]
				residue_number = int(line[23:26])

				# make a new residue if needed
				if residue_number != r_num:
					residue = Residue(residue_type, residue_number)
					r_num = residue_number

				# add the atom to the residue
				residue.atoms.append(atom)

			else:  # I've reached a TER card
				active_site.residues.append(residue)
	return active_site


def write_clustering(filename, clusters, ids):
	"""
	Write the clustered ActiveSite instances out to a file.

	Input: a filename and a clustering of ActiveSite instances
	Output: none
	"""

	out = open(filename, 'w')
	np.set_printoptions(threshold=np.inf)
	for i in range(len(clusters)):
		out.write("\nCluster %d\n--------------\n" % i)
		for j in range(len(clusters[i])):
			out.write("%s\n" % clusters[i][j])
	out.write("\nCluster ID Matrix\n")
	out.write("%s\n" % ids)
	out.close()


def write_mult_clusterings(filename, clusterings):
	"""
	Write a series of clusterings of ActiveSite instances out to a file.

	Input: a filename and a list of clusterings of ActiveSite instances
	Output: none
	"""

	out = open(filename, 'w')

	for i in range(len(clusterings)):
		clusters = clusterings[i]

		for j in range(len(clusters)):
			out.write("\nCluster %d\n------------\n" % j)
			for k in range(len(clusters[j])):
				out.write("%s\n" % clusters[j][k])

	out.close()
