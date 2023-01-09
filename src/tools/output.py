import glob
import sqlite3
import sys
from collections import defaultdict

from src.modules.PythonPDB2 import *
from src.modules.Structure import AntibodyResidue
from src.modules.Structure import PDBResInfo

import AbDbFunctions
import general
from src.modules.chains.AbChain import AbChain
from src.tools.path import *



def output_abchain_for_ben(ab_chain, OUTFILE):
    """
    Outputs data to a text file for input to bens C# renumbering code.
    """
    if not isinstance(ab_chain, AbChain): sys.exit()

    cdrs = defaultdict()
    cdrs["H"] = ["H1", "H2", "H3"]
    cdrs["L"] = ["L1", "L2", "L3"]


    full_line = ab_chain.get_full_id()
    total_resnum = full_line.split()[1]
    exp_method = full_line.split()[2]
    resolution = full_line.split()[3]
    pdbid_chain = ab_chain.get_pdbaa_tag_info().get_pdbid()+ab_chain.get_pdbaa_tag_info().get_original_chain()

    base_line =  pdbid_chain+" "+total_resnum+" "+exp_method+" "+resolution+"\t"

    #Skips bad entry 7fab and 1OCW
    if pdbid_chain.upper()[0:4] in ["7FAB", "1OCW"]:
        return

    #Skip very strange entry that Ben's renumbering code cannot handle.
    if pdbid_chain.upper()[0:4] == "4HJJ":
        return

    cdrs = ab_chain.get_cdrs()
    for cdr in cdrs:
        #print base_line
        base_line = base_line + cdr.get_type()+"@"+str(cdr.get_start())+" "+cdr.get_sequence(ab_chain.get_sequence())+" "
    OUTFILE.write(base_line+"\n")
    output_numbering_db_output(ab_chain, OUTFILE)

def output_abchain_general( ab_chain, OUTFILE):
    """
    Outputs data to a text file which does not use PDBAA header parsing to organize output.  Used for KABAT or general large FASTA analysis.
    """
    cdrs = defaultdict()
    cdrs["H"] = ["H1", "H2", "H3"]
    cdrs["L"] = ["L1", "L2", "L3"]


    #base_line = "ID> "+ab_chain.get_id()+" "+ab_chain.get_description()
    base_line = ab_chain.get_id()+" "+ab_chain.get_description()
    OUTFILE.write(base_line+"\n")

    #cdr_line = "CDRS> "
    cdr_line = ""
    cdrs = ab_chain.get_cdrs()
    for cdr in cdrs:
        #print base_line
        cdr_line = cdr_line + cdr.get_type()+"@"+str(cdr.get_start())+" "+cdr.get_sequence(ab_chain.get_sequence())+" "
    OUTFILE.write(cdr_line+"\n")
    output_numbering_db_output(ab_chain, OUTFILE)

def output_numbering_db_output(ab_chain, OUTFILE):
    """
    Outputs numbering as originally used for database output
    OUTFILE is an opened file
    DOES NOT CURRENTLY WORK WITH ICODES!!!
    """
    assert isinstance(ab_chain, AbChain)

    pdb_res_info = ab_chain.get_pdb_res_info()
    assert isinstance(pdb_res_info, PDBResInfo)

    for i in range(1, pdb_res_info.total_residue()+1):
        res = pdb_res_info.get_residue(i)
        assert isinstance(res, AntibodyResidue)
        new_chain = res.get_chain()
        new_resnum = res.get_pdb_num()
        resname = res.get_aa()
        region_type = res.get_region()
        line = new_chain+"_RENUMBER> "+resname+repr(i)+" "+resname+str(new_resnum)+" "+region_type+"\n"
        #print line
        OUTFILE.write(line)

def write_chain_subset(db_fname, outpath, gene):

    #print db_fname
    OUT = open_file(outpath+"/"+gene+"_chains.txt", 'w')
    con = sqlite3.connect(db_fname)
    con.row_factory = sqlite3.Row
    c = con.cursor()
    c.execute("SELECT DISTINCT PDB, original_chain FROM cdr_data WHERE gene =?", (gene,))
    rows = c.fetchall()
    for row in rows:

        row.keys()
        outname = row['PDB']+row['original_chain']
        #print outname
        OUT.write(outname+"\n")

    OUT.close()

def make_pdblist_from_chains(pices_out, renum_dir, out_name):

    OUTFILE_rel = open_file(renum_dir+"/"+"rel_"+out_name, 'w')
    OUTFILE_full = open_file(renum_dir+"/"+"full_"+out_name, 'w')
    pdbs = glob.glob(renum_dir+"/*.pdb")
    PICES = open_file(pices_out, 'r')
    PICES.readline()

    pices_pdb_chains = defaultdict(str)
    for line in PICES:
        pdb = line.strip().split()[0][0:4]
        pices_pdb_chains[pdb] == line.strip().split()[0][4]


    for pdb_file in pdbs:
        #print pdb_file
        fname = pdb_file.split("/")[-1].split(".")[0]

        fnameSP = fname.split("_")
        pdb = fnameSP[0].upper()[0:4]
        if pices_pdb_chains.has_key(pdb):
            OUTFILE_rel.write(pdb_file.split("/")[-1]+"\n")
            OUTFILE_full.write(pdb_file+"\n")

    OUTFILE_full.close()
    OUTFILE_rel.close()
    PICES.close()

def reorder_and_save_chains(in_path, out_path, remove_het = False):
    blank_pdb = PythonPDB2()
    full_pdb = PythonPDB2(in_path)

    blank_pdb.copy_chain_into_pdb_map(full_pdb, "A")
    blank_pdb.copy_chain_into_pdb_map(full_pdb, "L")
    blank_pdb.copy_all_but_chains_into_pdb_map(full_pdb, ["A", "L"])

    if remove_het:
        blank_pdb.remove_hetatm_atoms()
        blank_pdb.remove_waters()

    blank_pdb.save_PDB(out_path)

def create_center_cluster_dihedral_file(db_path, out_path):
    if not os.path.exists(db_path):
        print "DB path does not exist: "+db_path
        return
    db = sqlite3.connect(db_path)
    data = AbDbFunctions.get_dihedral_string_for_centers(db)

    OUTFILE = open_file(out_path, 'w')
    for fullcluster in sorted(data):
        cdr = fullcluster.split('-')[0]
        length = fullcluster.split('-')[1]
        cluster = "-".join(fullcluster.split('-')[2:])
        length_type = data[fullcluster][0]
        ss = data[fullcluster][1]
        dihedral_string = data[fullcluster][2]

        new_dihedral_string = general.reformat_dihedrals(dihedral_string, ':')
        line = cdr+" "+length+" "+cluster+" "+length_type+" "+fullcluster+" "+ss+" "+new_dihedral_string
        print line
        OUTFILE.write(line+"\n")
    OUTFILE.close()