#!/usr/bin/env python

#Author Jared Adolf-Bryfogle



#Python Imports
import os
import sys
from argparse import ArgumentParser
import re
from collections import defaultdict

from Bio.PDB import PDBIO
from Bio.PDB.PDBParser import PDBParser

#Append Python Path
p = os.path.split(os.path.abspath(__file__))[0]+"/src"
sys.path.append(p); #Allows all modules to use all other modules, without needing to update pythonpath

#Project Imports
from src.IgClassifyFASTA import IgClassifyFASTA
from src.IgClassifyPDB import IgClassifyPDB
from src.tools.fasta import *




if __name__ == '__main__':

    parser = ArgumentParser("PyIgClassify can identify clusters of a PDB, renumber a PDB, or analyze CDR positions and lengths of a fasta."
                            "Please pass either -s or -l for a PDB file or PDBLIST.  mmCIF and gzipped files are supported.")
    args = sys.argv
    parser.add_argument("--pdb_mode", "-p",
        default = False,
        action="store_true",
        help = "Renumber a PDB")
    
    parser.add_argument("--fasta_mode", "-t",
        default = False,
        action="store_true",
        help = "Identify CDRs in FASTA split by chain")

    parser.add_argument("--single", "-s",
        help = "Path to single file.")

    parser.add_argument("--list", "-l",
        help = "List of PDBs or FASTAs to renumber.")

    parser.add_argument("--outname", "-n",
        default = "PyIgClassifyOutput",
        help = "Prefix to use for renumbered PDB or FASTA identification.")

    parser.add_argument("--outdir", "-o",
        default = "USEROUT",
        help = "Directory to output the results to (Renumbered PDB / FASTA identification")

    parser.add_argument("--input_dir", "-i",
        help = "Directory of PDBs or FASTA file for which to use PDBLIST without full paths.")

    parser.add_argument("--concat", "-c",
        default = False,
        action = "store_true",
        help = "Concat output data to output file instead of overriding each time.  Will create if not exists")

    parser.add_argument("--skip_renumber", "-k",
        default = False,
        action = "store_true",
        help = "Skip renumbering of the input PDB file(s).  Only identify clusters.  Must be already renumbered in the AHO scheme.")

    parser.add_argument("--numbering_scheme", "-a",
        default = "modified_aho",
        help = "Numbering scheme to renumber or identify PDB or FASTA.  Only supporting modified_aho for now.")


    options = parser.parse_args()
    
    #Option Error Checking:
    
    if options.pdb_mode and options.fasta_mode:
        sys.exit("Please pick only one mode.")
    
    if options.pdb_mode and not options.single and not options.list:
        sys.exit("PDB mode requires a PDB or a list!")
    
    if options.fasta_mode and not options.single and not options.list:
        sys.exit("FASTA mode requires a FASTA or a list!")
        
    
    if not os.path.exists(options.outdir):
        os.mkdir(options.outdir)

    if options.list and not os.path.exists(options.list):
            sys.exit("Could not find list!")

    paths = []
    if options.list:
        FILE = open_file(options.list)
        for p in FILE:
            p = p.strip()
            paths.append(p)
        FILE.close()
    else:
        paths.append(options.single)

    #Set a mode based on input, as I am sick of setting --pdb_mode or --fasta_mode.
    if not options.pdb_mode and not options.fasta_mode:
        #Probably a more elegant way to write this by counting what re.search outputs.
        pdb_modes = [".pdb", ".cif", ".mmcif"]
        fasta_modes = [".fasta", ".txt"]

        pdb_matches = [ re.search(x, paths[0]) for x in pdb_modes]
        fasta_matches = [ re.search(x, paths[0]) for x in fasta_modes]
        if len(pdb_matches) - pdb_matches.count(None):
            options.pdb_mode = True
        elif len(fasta_matches) - fasta_matches.count(None):
            options.fasta_mode = True
        else:
            sys.exit("Could not parse mode from file name.  Pleas pass either --fasta_mode or --pdb_mode")

        if options.pdb_mode and options.fasta_mode:
            sys.exit("Cannot have strange extensions (such as .pdb.txt or .pdb.fasta).  "
                     "Please pass either --fasta_mode or --pdb_mode")

    if options.pdb_mode:

        print "Running PDB mode"
        for path in paths:

            if options.input_dir:
                path = os.path.join(options.input_dir, path)

            print "Running "+path
            if not os.path.exists(path):
                sys.exit("Could not find PDB file: "+path)

            pdb_renumber = IgClassifyPDB(path)
            pdb_renumber.set_output_path(options.outdir)
            pdb_renumber.set_output_name(options.outname)
            pdb_renumber.run(options.numbering_scheme, options.concat, options.skip_renumber)
        
    elif options.fasta_mode:
        print "Running FASTA mode"
        for path in paths:
            print "Running "+path
            if not os.path.exists(path):
                sys.exit("Could not find FASTA file:"+path)

            fasta_identifier = IgClassifyFASTA(path)
            fasta_identifier.set_output_name(options.outname)
            fasta_identifier.set_output_path(options.outdir)
            fasta_identifier.run(options.numbering_scheme, options.concat)
        

    print "PyIgClassify complete.  Please see outdir"

    if not options.pdb_mode and not options.fasta_mode:
        sys.exit("Please pick a mode:  Options are --pdb_mode and --fasta_mode")