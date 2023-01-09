import os
from optparse import *

from src.tools.output import make_pdblist_from_chains

if __name__ == '__main__':

    p = os.path.split(os.path.abspath(__file__))[0]
    parser = OptionParser()

    parser.add_option("--pices_chains", "-l")

    parser.add_option("--out_name", "-o",
        default="culled_PDBLIST.txt",)

    parser.add_option("--renum_dir", "-d",
        default = p+"/../DBOUT/renumbered_pdbs")

    (options, args) = parser.parse_args()

    make_pdblist_from_chains(options.pices_chains, os.path.abspath(options.renum_dir), options.out_name)