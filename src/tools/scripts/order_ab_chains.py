#Jared Adolf-Bryfogle
#Reorders PDBFiles in a dirctory according to LH_A in order for Rosetta Antibody Design benchmarking. Removes HetAtm!!!

import sys

from tools.path import *

from src.tools.output import reorder_and_save_chains


if __name__ == "__main__":

    #I don't have time to make this fancy right now.

    #Arguments:
    ### Input Directory
    ### Input PDBList
    ### Output Directory

    #Assumes PDBList has no paths and requires an input directory as if we run Rosetta.

    in_dir = sys.argv[1]
    in_pdblist = sys.argv[2]
    out_dir = sys.argv[3]

    if not os.path.exists(out_dir):
        os.mkdir(out_dir)

    PDBLIST = open_file(in_pdblist, 'r')
    for line in PDBLIST:
        line = line.strip()
        line = in_dir+"/"+line

        print "Reordering "+line
        outpath = out_dir+"/"+os.path.basename(line)
        reorder_and_save_chains(line, outpath)

    print "Done!!"
    PDBLIST.close()