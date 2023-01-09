#Jared Adolf-Bryfogle
#Quick script for separating out the Fv from the FAB after all Antibodies have been renumbered.  Should really be a class.


import copy
import glob
import os
import sys
from collections import defaultdict

from src import tools as ab_db
from src.modules.PythonPDB import PythonPDB


def run_main(ab_dir, output_dir, only_dimer=True):

    if not os.path.exists(ab_dir):
        sys.exit(ab_dir+" does not exist!  Please check your path and try again.")

    if not os.path.exists(output_dir):
        os.mkdir(output_dir)

    dirs = ["Fv", "FAB", "Fc", "linker", "linker_L", "linker_H"]
    for dir in dirs:
        if not os.path.exists(output_dir+"/"+dir):
            os.mkdir(output_dir+"/"+dir)

    antibodies = glob.glob(ab_dir+"/*.pdb")
    print "Separating Fv, Fc, and linker from FABs. Only_dimer = "+repr(only_dimer)

    for ab in antibodies:
        separate_pdb(ab, output_dir, only_dimer)

def run_split_proto_CDR4(ab_dir, output_dir, overhang = 0, skip_present = False):
    if not os.path.exists(output_dir):
        os.mkdir(output_dir)

    if not os.path.exists(ab_dir):
        sys.exit(ab_dir+" does not exist!  Please check your path and try again.")

    antibodies = glob.glob(ab_dir+"/*.pdb")

    i  = 1
    for ab in antibodies:
        print "Splitting "+repr(i)
        separate_proto_CDR4(ab, output_dir, skip_present)
        i+=1

    print "complete"

def run_split_proto_CDR4_by_gene(db, ab_dir, output_dir, overhang = 0, skip_present = False, res_cutoff = 2.8, rfac_cutoff = .3):

    all_antibodies = glob.glob(ab_dir+"/*.pdb")
    ab_path_dict = defaultdict()

    for ab in all_antibodies:
        fname = os.path.basename(ab)
        fname = fname.split(".")[0]
        fnameSP = fname.split("-")
        for f in fnameSP:
            ab_path_dict[f] = ab

    for gene in ["lambda", "kappa", "heavy", "lambda6"]:
        print "working on: "+gene
        if gene in ["lambda", "kappa", "lambda6"]:
            chain = "L"
        else:
            chain = "H"

        out_dir = output_dir+"/"+gene
        if not os.path.exists(out_dir):
            os.mkdir(out_dir)

        pdb_chains = ab_db.get_pdb_chain_subset(db, gene, True, res_cutoff, rfac_cutoff)

        i = 1
        for pdb_chain in pdb_chains:
            fname = pdb_chain[0].lower() + pdb_chain[1]
            if ab_path_dict.has_key(fname):
                separate_proto_CDR4(ab_path_dict[ fname ], out_dir, chain, overhang, skip_present)
                #if i == 20: break
                i+=1


def separate_proto_CDR4(pdb_path, output_dir, chain, overhang = 0, skip_present = True):
    pdb_name = pdb_path.split("/")[-1]
    parent_PDB = PythonPDB(pdb_path)

    if skip_present and os.path.exists(output_dir+"/"+pdb_name):
        return

    if not os.path.exists(output_dir):
        os.mkdir(output_dir)

    if (has_chain(parent_PDB.get_pdb_map(), chain)):
        cdr4 = split_CDR4(parent_PDB, chain, overhang)
        parent_PDB.set_pdb_map(cdr4)
        parent_PDB.save_PDB(output_dir+"/"+pdb_name)

def separate_pdb(pdb_path, output_dir, only_dimer = True):
    """
    Determine if we have Fv or FAB.  If FAB, split into parts: Fc, Fv, linker.
    """
    pdb_name = pdb_path.split("/")[-1]
    parent_PDB = PythonPDB(pdb_path)
    if not (has_chain(parent_PDB.get_pdb_map(), "L") and has_chain(parent_PDB.get_pdb_map(), "H")):


        if only_dimer:
            print pdb_name +" Missing L or H"
            print "Only splitting dimers - skipping..."
            return

    #split_linker_save(parent_PDB, output_dir, pdb_name)

    if not has_Fc(parent_PDB):
        parent_PDB.save_PDB(output_dir+"/Fv/"+pdb_name)
        #message(pdb_name, "Fv")
    else:
        parent_PDB.save_PDB(output_dir+"/FAB/"+pdb_name)
        #message(pdb_name, "FAB")

        Fv_pdb_map = split_Fv(parent_PDB)
        Fc_pdb_map = split_Fc(parent_PDB)
        li_pdb_map = split_linker(parent_PDB)

        parent_PDB.set_pdb_map(Fc_pdb_map)
        parent_PDB.save_PDB(output_dir+"/Fc/"+pdb_name)


        parent_PDB.set_pdb_map(Fv_pdb_map)
        parent_PDB.save_PDB(output_dir+"/Fv/"+pdb_name)


        parent_PDB.set_pdb_map(li_pdb_map)
        parent_PDB.save_PDB(output_dir+"/linker/"+pdb_name)

        liL_pdb_map = split_linker_L(parent_PDB)
        liH_pdb_map = split_linker_H(parent_PDB)

        if has_chain(liL_pdb_map, "L"):
            parent_PDB.set_pdb_map(liL_pdb_map)
            parent_PDB.save_PDB(output_dir+"/linker_L/"+pdb_name)

        if has_chain(liH_pdb_map, "H"):
            parent_PDB.set_pdb_map(liH_pdb_map)
            parent_PDB.save_PDB(output_dir+"/linker_H/"+pdb_name)

        #message(pdb_name, "Fc, Fv, and linker")

def has_Fc(parent_PDB):
    pdb_map = parent_PDB.get_pdb_map()
    for i in pdb_map:
        if int(pdb_map[i]["residue_number"]) > 200:
            return True
        else:
            continue
    return False


def has_chain(pdb_map, chain):
    for i in pdb_map:
        if pdb_map[i]["chain"] == chain:
            return True
    return False

def message(pdb_name, type):
    print "Saved "+pdb_name+" as "+type


def split_Fc(parent_PDB):
    """
    Split Fc from FAB.  Return new pdb_map to save
    """
    pdb_map = parent_PDB.get_pdb_map()
    pdb_map_cp = copy.deepcopy(parent_PDB.get_pdb_map())
    for i in pdb_map:
        if int(pdb_map[i]['residue_number']) < 153:
            del pdb_map_cp[i]
    return pdb_map_cp

def split_Fv(parent_PDB):
    """
    Split Fv from FAB.  Return new pdb_map to save
    """
    pdb_map = parent_PDB.get_pdb_map()
    pdb_map_cp = copy.deepcopy(parent_PDB.get_pdb_map())
    for i in pdb_map:
        if int(pdb_map[i]['residue_number']) > 149:
            del pdb_map_cp[i]
    return pdb_map_cp

def split_CDR4(parent_PDB, chain, overhang = 0):
    pdb_map = parent_PDB.get_pdb_map()
    pdb_map_cp = copy.deepcopy(parent_PDB.get_pdb_map())
    for i in pdb_map:
        resnum = int(pdb_map[i]['residue_number'])
        if pdb_map[ i ]['chain'] != chain:
            del pdb_map_cp[ i ]
            continue
        if resnum < (82 - overhang) or resnum > (89 + overhang):

            del pdb_map_cp[ i ]
            continue
    return pdb_map_cp

def split_linker(parent_PDB):
    """
    Split 4 Residue linker from FAB.  Linker may be longer than this, but I think this is about it.
    """
    pdb_map = parent_PDB.get_pdb_map()
    pdb_map_cp = copy.deepcopy(parent_PDB.get_pdb_map())
    for i in pdb_map:
        resnum = int(pdb_map[i]['residue_number'])
        if resnum < 149 or resnum > 152:
            del pdb_map_cp[i]
    return pdb_map_cp

def split_linker_L(parent_PDB):
    """
    Split 4 Residue linker from FAB.  Linker may be longer than this, but I think this is about it.
    """
    pdb_map = parent_PDB.get_pdb_map()
    pdb_map_cp = copy.deepcopy(parent_PDB.get_pdb_map())
    for i in pdb_map:
        if pdb_map[i]['chain'] == "H":
            del pdb_map_cp[i]
            continue
        resnum = int(pdb_map[i]['residue_number'])
        if resnum < 149 or resnum > 152:
            del pdb_map_cp[i]
    return pdb_map_cp

def split_linker_H(parent_PDB):
    """
    Split 4 Residue linker from FAB.  Linker may be longer than this, but I think this is about it.
    """
    pdb_map = parent_PDB.get_pdb_map()
    pdb_map_cp = copy.deepcopy(parent_PDB.get_pdb_map())
    for i in pdb_map:
        if pdb_map[i]['chain'] == "L":
            del pdb_map_cp[i]
            continue
        resnum = int(pdb_map[i]['residue_number'])
        if resnum < 149 or resnum > 152:
            del pdb_map_cp[i]
    return pdb_map_cp

if __name__ == '__main__':

    only_dimer = True
    ab_dir = sys.argv[1]
    output_dir = sys.argv[2]

    run_main(ab_dir, output_dir, only_dimer)







