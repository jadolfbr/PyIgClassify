#!/bin/sh

#PDB Mode
./PyIgClassify.py --pdb_mode -s testing/2j88.pdb

#FASTA Mode
./PyIgClassify.py --fasta_mode --s testing/combined_test.fasta

#PDB List Mode
./PyIgClassify.py --pdb_mode -l testing/PDBLIST.txt

#FASTA List Mode
./PyIgClassify.py --fasta_mode -l testing/FASTALIST.txt
