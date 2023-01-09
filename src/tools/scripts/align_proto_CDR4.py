import glob
import os
import sqlite3
import sys
from collections import defaultdict

import src.tools.rosetta_alignment as aln

from rosetta import *
rosetta.init("-ignore_unrecognized_res -ignore_zero_occupancy false -use_input_sc -ex1 -ex2")


def align_and_save(outpath):

    inpath = outpath+"/redun"
    if not os.path.exists(outpath+"/stem_align"):
        os.mkdir(outpath+"/stem_align")

    types = ["redun", "nr"]
    for gene in ["lambda", "kappa", "heavy"]:
        files = glob.glob(inpath+"/"+gene+"/*.pdb")
        print inpath+"/"+gene+"/*.pdb"

        for t in types:
            poses = load_poses(files, t)

            out = outpath+"/stem_align/"+t
            if not os.path.exists(out):
                os.mkdir(out)
            out = out+"/"+gene
            if not os.path.exists(out):
                os.mkdir(out)

            print out
            poses[ 0 ].dump_pdb(out+"/"+os.path.basename(poses[ 0 ].pdb_info().name()))
            for i in range(1, len(poses)):
                print out
                pose_name = os.path.basename(poses[ i ].pdb_info().name())
                print pose_name
                if pose_name == "1oauI.pdb" or pose_name == "1oauO-1oauK.pdb":
                    print "skipping problematic PDB "+pose_name
                    continue
                aln.align_to_second_pose_save_pdb(pose_name, poses[ i ], poses[ 0 ], out, 3, True)

def load_poses(files, t):


    poses = []
    for file in files:
        print "loading "+file
        p = pose_from_pdb(file)
        poses.append(p)

    if t == "redun":
        return poses
    else:
        seq_dic = defaultdict()
        for p in poses:
            seq_dic[p.sequence()] = p
        return seq_dic.values()


if __name__ == "__main__":
    inpath = sys.argv[1]
    outpath = sys.argv[2]
    in_db = sys.argv[3]
    overhang = 3
    skip_current = True
    res_cutoff = 2.8
    rfac_cutoff = .3

    if not os.path.exists(outpath): os.mkdir(outpath)
    if not os.path.exists(inpath): sys.exit(inpath+ " does not exist!")

    db = sqlite3.connect(in_db)
    #ab_splitter.run_split_proto_CDR4_by_gene(db, inpath, outpath+"/redun", overhang, skip_current, res_cutoff, rfac_cutoff)
    align_and_save(outpath)




