#!/usr/bin/env python

#Jared Adolf-Bryfogle (jadolfbr@gmail.com)

import gzip
import os
import re
from collections import defaultdict

from Bio.PDB import MMCIFParser
from Bio.PDB import PDBParser
from Bio.PDB import calc_dihedral
from src.modules.restype_definitions import definitions as RestypeDefinitions
from src.tools import vector1
from src.tools.biopython_util import peptide_bond_distance
from src.tools.path import *

from src.modules.Structure import PDBResInfo as PDBInfo
from src.modules.Structure import Residue as ResidueRecord
from src.tools.numeric import *


class BioPose(object):
    """
    This is my biopython meta class.  Because biopython's interface kinda sucks.
    This is a little cleaner.

    The other way is to sublclass each Biopython class structure, which I'm not ready to do.

    Right now, you need a path as I don't know how we would use this from sequence, etc as you do in Rosetta.
    :path: Is a path to an RCSB file.  PDB (.pdb), mmCIF(.cif), and gzipped (.gz) versions.
    """
    def __init__(self, path, model_num=0):

        self.res_definitions = RestypeDefinitions()
        self.struct, self.header = self.load_from_file(path) #Bio struct, Header dictionary
        self._setup_self(model_num)

    ############ IO ###################
    def load_from_file(self, path):
        """
        Load a file from PDB or mmCIF.  .gz is supported.

        :param path: Path to PDB or mmCIF file
        :rtype: tuple(bio.PDB.Structure.Structure, dict)
        """
        structure = None
        if re.search(".pdb", path):
            parser = PDBParser()
        else:
            parser = MMCIFParser()

        path = path.strip()
        model_id = os.path.basename(path)
        if os.path.basename(path).split('.')[-1] == 'gz':
            GZ = gzip.open(path, 'rb')
            GZ.close()
        else :
            structure = parser.get_structure(model_id, path)
        header = parser.get_header()

        return structure, header

    def reload_from_file(self, path, model_num=0):
        """
        Reload a BioPose from a file path.
        :param path: str
        :param model_num: int
        :return:
        """
        self.struct, self.header = self.load_from_file(path)
        self._setup_self(model_num)

    def _setup_self(self, model_num =0):
        """
        Setup all info for bio pose after loading from file.
        :param model_num: int
        :return:
        """
        self.all_residues = self._setup_all_residues(model_num) #vector1 of Bio Residues
        self.pdb_info = self._setup_pdb_info() #PDBInfo to map the vector1
        self.peptide_bond_distances = self._setup_peptide_bond_distances() #map of bond distances to next residue in pose.

    ############ Getting Components #############
    def structure(self):
        """
        Get the Bio Structure stored in this class.
        :rtype: bio.PDB.Structure.Structure
        """
        return self.struct

    def model(self, model_num = 0):
        """
        Get a Bio Model of the stored structure
        :param id: int
        :rtype: bio.PDB.Model.Model
        """
        return self.struct[model_num]

    def chain(self, chain_id, model_num = 0):
        """
        Get a Bio Chain of the stored structure
        :param chain_id: str
        :param model_num: int
        :rtype: bio.PDB.Chain.Chain
        """
        return self.struct[model_num][chain_id]

    def residue(self, resnum, chain_id, icode=' ', alt=' ', model_num = 0):
        """
        Get a Bio Residue of the stored structure.
        Adds a chain_id attribute.

        :param resnum: int
        :param chain_id: str
        :param icode: str
        :param alt: str
        :param model_num: int
        :rtype: bio.PDB.Residue.Residue
        """
        res = self.struct[model_num][chain_id][(alt, resnum, icode)]
        res.chain_id = chain_id
        return res


    def atom(self, atom_name, resnum, chain_id, icode=' ', alt=' ', model_num=0):
        """
        Get a Bio Atom of the stored structure

        :param atom_name: str
        :param resnum: int
        :param chain_id: str
        :param icode: str
        :param alt: str
        :param model_num: int
        :rtype: bio.PDB.Atom.Atom
        """
        return self.struct[model_num][chain_id][(alt, resnum, icode)][atom_name]


    ### Lists of Structures ###

    def chains(self, model_num = 0):
        """
        Get a list of Bio Chains
        :param model_num: int
        :rtype: list[bio.PDB.Chain.Chain
        """
        return [c for c in self.struct]

    def residues(self, chain_id, model_num = 0, include_alt = False):
        """
        Get residues, including or not including residues with alternate location codes - which can be a PITA
        Adds chain_id attribute to residue.

        :param chain_id: str
        :param model_num: int
        :param include_alt: bool
        :rtype: list[bio.PDB.Residue.Residue]
        """
        resi = []
        for res in self.chain(chain_id, model_num):
            res.chain_id = chain_id
            if res.id[0] ==' ':
                resi.append(res)
            elif include_alt:
                resi.append(res)
            else:
                continue
        return resi

    def atoms(self, resnum, chain_id, icode=' ', alt=' ', model_num = 0):
        """
        Get a list of Bio Atoms
        :param resnum: int
        :param chain_id: str
        :param icode: str
        :param alt: str
        :param model_num: int
        :rtype: list[bio.PDB.Atom.Atom]
        """
        return [atm for atm in self.residue(resnum, chain_id, icode, alt, model_num)]


    ############ Helper Funtions ################
    def total_residue(self):
        return len(self.all_residues)

    def phi(self, resi):
        """
        Get the Phi Angle of i in radians

        :param i: int
        :rtype: float
        """
        if resi == 1 or not self.connected_to_previous(resi):
            return 0.0

        res = self.all_residues[resi]

        try:
            n = res['N'].get_vector()
            ca = res['CA'].get_vector()
            c = res['C'].get_vector()

            res_minus_one = self.all_residues[resi -1]
            cp = res_minus_one['C'].get_vector()
            phi = calc_dihedral(cp, n, ca, c)
            return phi
        except Exception:
            print "Could not get phi for "+repr(resi)
            raise LookupError



    def psi(self, resi):
        """
        Get the Psi Angle of i in radians

        :param i: int
        :rtype: float
        """
        res = self.all_residues[resi]

        if resi == len(self.all_residues) or not self.connected_to_next(resi):
            return 0.0

        try:
            n = res['N'].get_vector()
            ca = res['CA'].get_vector()
            c = res['C'].get_vector()
            res_plus_one = self.all_residues[resi + 1]

            nn = res_plus_one['N'].get_vector()
            psi = calc_dihedral(n, ca, c, nn)
            return psi
        except Exception:
            print "Could not get psi for "+repr(resi)
            raise LookupError


    def omega(self, resi, rosetta_definitions = True):
        """
        Get the Omega Angle of i in radians
        Omega is defined as the dihedral angle between the peptide bond of i and i + 1, as in Rosetta.
        If rosetta_definitions are False, omega is then treated as being between i and i -1

        :param i: int
        :param reverse_rosetta_definitions: bool
        :rtype: float
        """
        res = self.all_residues[resi]

        try:
            n = res['N'].get_vector()
            ca = res['CA'].get_vector()
            c = res['C'].get_vector()

            if rosetta_definitions and resi < len(self.all_residues) -1 and self.connected_to_next(resi):
                res_plus_one = self.all_residues[resi + 1]
                next_n = res_plus_one['N'].get_vector()
                next_ca = res_plus_one['CA'].get_vector()
                omega = calc_dihedral(ca, c, next_n, next_ca)
                return omega

            elif not rosetta_definitions and resi > 1 and self.connected_to_previous(resi):
                res_minus_one = self.all_residues[resi - 1]
                pre_c = res_minus_one['C'].get_vector()
                pre_ca = res_minus_one['CA'].get_vector()
                omega = calc_dihedral(pre_ca, pre_c, n, ca)
                return omega
            else:
                return 0.0

        except BaseException:
            print "Could not get omega for "+repr(resi)
            raise LookupError


    def res_bond_distance(self, resi):
        """
        Get the stored bond distances between residue and residue+1
        :param res: int
        :rtype: float
        """
        return self.peptide_bond_distances[resi]

    def connected_to_next(self, resi, peptide_bond_distance_cutoff=1.8):
        if self.res_bond_distance(resi) <= peptide_bond_distance_cutoff:
            return True
        else:
            return False

    def connected_to_previous(self, resi, peptide_bond_distance_cutoff=1.8):
        if resi == 0:
            return False
        elif self.res_bond_distance(resi-1) <= peptide_bond_distance_cutoff:
            return True
        else:
            return False

    def get_sequence(self, chain_id, model_num = 0):
        """
        Get a sequence of a chain - Not including alternate res locations

        :param chain_id: str
        :param model_num: int
        :rtype: str
        """
        if self.get_chain_length(chain_id, model_num) == 0:
            return ""

        seq = ""
        for res in self.residues(chain_id, model_num):
            aa = self.res_definitions.get_one_letter_from_three(res.resname)
            if aa:
                seq = seq+aa
            else:
                seq = seq+'X'

        return seq

    def get_regional_sequence(self, start, end, model_num =0):
        seq = ""
        for i in range(start, end+1):
            aa = self.res_definitions.get_one_letter_from_three(self.all_residues[i].get_resname())
            if aa:
                seq = seq+aa
            else:
                seq = seq+'X'
        return seq

    def get_chain_length(self, chain_id, model_num = 0):
        """
        Get the number of AA in a chain - Not including alternate res locations
        :param chain_id: str
        :rtype: int
        """
        return len(self.residues(chain_id, model_num))

    def get_chain_ids(self, model_num):
        """
        Get all chain IDS for a model.
        :param model_num: int
        :rtype: list[str]
        """
        ids = []
        for chain in self.model(model_num):
            ids.append(chain.id)
        return ids

    ###### Private Members #######

    def _setup_all_residues(self, model_num=0):
        """
        Setup all residues in the model, for ease of use like a pose.
        :param model_num: int
        :rtype: list[bio.PDB.Residue.Residue]
        """
        all_residues = vector1()

        for chain_id in self.get_chain_ids(model_num):
            print "ChainID: "+chain_id
            residues = self.residues(chain_id, model_num)
            all_residues.extend(residues)

        return all_residues

    def _setup_pdb_info(self):
        """
        Setup the PDBInfo mapping from residues to pdb iformation.
        :param model_num: int
        :rtype: PDBInfo
        """
        pdb_info = PDBInfo()
        i = 1
        for res in self.all_residues:
            s = ResidueRecord(self.res_definitions.get_one_letter_from_three(res.resname), res.id[1], res.chain_id, res.id[2])
            pdb_info.add_residue(s)

        return pdb_info

    def _setup_peptide_bond_distances(self):
        bond_distances = defaultdict()
        for i in range(1, len(self.all_residues)+1):
            res = self.all_residues[i]
            if i == len(self.all_residues):
                bond_distances[i] = 0
            else:
                try:
                    bond_distances[i] = peptide_bond_distance(self.all_residues[i], self.all_residues[i+1])
                except Exception:
                    bond_distances[i] = None

        return bond_distances

########## Testing Functions #############
def test_dihedrals(pose):
    """
    Simple Test for Dihedral output

    :param pose: BioPose
    :rtype: bool
    """
    for i in range(1, pose.total_residue()+1):

        print "\n"+str(pose.pdb_info.pose2pdb(i))
        try:
            print "Phi: "+repr(math.degrees(pose.phi(i)))
            print "Psi: "+repr(math.degrees(pose.psi(i)))
            print "Omega:"+repr(math.degrees(pose.omega(i)))
        except Exception:
            "Print could not get dihedral for resnum "+repr(i)

    return True

def test_pdbinfo(pose):
    """
    Simple Test for pdbinfo output.
    :param pose: BioPose
    :rtype: bool
    """
    for i in range(1, pose.total_residue() +1):
        print repr(i)
        print pose.all_residues[i].id
        print pose.pdb_info.pose2pdb(i)


######## Testing Main ############
if __name__ == "__main__":
    v = vector1([1, 2, 3])
    for i in v:
        print repr(i)

    print len(v)

    print repr(v[3])

    test_pdb = os.path.join(get_testing_inputs_path(),"2j88.pdb")
    pose = BioPose(test_pdb)
    test_dihedrals(pose)






