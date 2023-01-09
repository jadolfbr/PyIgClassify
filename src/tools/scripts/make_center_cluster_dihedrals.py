#!/usr/env python

import os
import sys

from src.tools.output import create_center_cluster_dihedral_file

if __name__ == "__main__":

    if len(sys.argv) < 2:
        print "Arguments are:"
        print " 1) Input DB path"
        print " 2) Output path"
        sys.exit()

    db_path = sys.argv[1]
    out_path = sys.argv[2]

    if not os.path.exists(db_path):
        sys.exit("DB Path does not exist")

    create_center_cluster_dihedral_file(db_path, out_path)
    print "Complete."
