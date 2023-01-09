import os
from optparse import *

from src.tools.output import write_chain_subset
#Script to get list of chains for use in Pices by gene.  Convert to functional program if needed.


if __name__=='__main__':
    p = os.path.split(os.path.abspath(__file__))[0]
    parser = OptionParser()

    parser.add_option("--out_path", "-o", default = ".")
    parser.add_option("--db_path", "-d", default = p+"/../DBOUT/antibody_database_nr_by_cdrs_per_chain.db")
    parser.add_option("--gene", "-g", default = "heavy")

    (options, args) = parser.parse_args()

    write_chain_subset(options.db_path, options.out_path, options.gene)