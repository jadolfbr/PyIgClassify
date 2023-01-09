import gzip
import os
import sys

from Bio.PDB.PDBParser import PDBParser

from src.modules.restype_definitions import definitions as RestypeDefinitions
from src.tools.numeric import *


### NOTE: All Utility function have been replaced by a Bio Structure wrapper: BioPose.
### Please see this new class for future developments!

########  NEW Biopython utility functions ##########

def peptide_bond_distance(res1, res2):
    """
    Return the bond distance between two residues using Numpy array math.
    :param res1: Bio.PDB.Residue.Residue
    :param res2: Bio.PDB.Residue.Residue
    :rtype: float
    """
    return atomic_distance(res1, res2, 'C', 'N')

def atomic_distance(res1, res2, res1_atom_name, res2_atom_name):
    """
    Return the atomic distance between two arbitrary Bio residues and two arbitrary atom names.
    :param res1: Bio.PDB.Residue.Residue
    :param res2: Bio.PDB.Residue.Residue
    :param res1_atom_name: str
    :param res2_atom_name: str
    :rtype: float
    """
    try:
        return distance_numpy(res1[res1_atom_name].get_vector().get_array(), res2[res2_atom_name].get_vector().get_array())
    except Exception:
        print "Residue does not have the atom name or there is a problem in the vector.  Returning 0"
        raise IndexError



########  OLD Biopython Utility Functions replaced by BIOPose ########

def has_id(model, id):
    """
    Returns true or false if the model has the chain.  Because biopython is not updating it's index that has_id is using.  WTF.
    """
    for i in model:
        if i.id == id:
            return True
    return False

def get_biopython_structure(path, model_id = None):
    structure = None
    path = path.strip()
    parser = PDBParser()
    if not model_id:
        model_id = os.path.basename(path)
    if os.path.basename(path).split('.')[-1] == "pdb":
        structure = parser.get_structure(model_id, path)
    elif os.path.basename(path).split('.')[-1] == 'gz':
        GZ = gzip.open(path, 'rb')
        structure = parser.get_structure(model_id, GZ)
        GZ.close()
    else:
        sys.exit("Unknown extension to read PDB: "+path)

    return structure



def get_seq_from_biostructure(structure, chain_id):
    for biochain in structure[0]:
        if get_chain_length(biochain) == 0:
            continue
        if biochain.id == chain_id:
            return get_seq_from_biochain(biochain)

    print "Chain not found!"
    raise LookupError

def get_seq_from_biochain(bio_chain):

    if get_chain_length(bio_chain) == 0:
            return ""

    seq = ""
    d = RestypeDefinitions()

    for res in bio_chain:
        if res.id[0]==' ':
            aa = d.get_one_letter_from_three(res.resname)
            if not aa:
                print "Skipping non-canonical resname: "+res.resname
                print "This could pose a problem!"
                continue
            seq = seq+aa
    return seq

def get_chain_length(bio_chain):

    l = 0
    for res in bio_chain:
        if res.id[0]==' ':
            l+=1
    return l

def get_num_biochains(model):
    return len(model[0])