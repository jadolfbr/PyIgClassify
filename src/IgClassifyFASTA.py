#!/usr/bin/env python

#Author Jared Adolf-Bryfogle



#Python Imports
import os
import sys

#Append Python Path
p = os.path.split(os.path.abspath(__file__))[0]+"/src"
sys.path.append(p); #Allows all modules to use all other modules, without needing to update pythonpath

#Project Imports
from src.modules.chains.ProposedIgChain import ProposedIgChain
from src.modules.chains.IgChainSet import IgChainSet
from src.modules.chains.IgChainFactory import IgChainFactory
from src.modules.chains.AbChainFactory import AbChainFactory

from src.tools.fasta import *


class IgClassifyFASTA:
    """
    Identify CDRs from a Fasta File
    """
    def __init__(self, fasta_path):
        self.fasta_path = fasta_path
        self.outdir = os.getcwd()
        self.fasta_paths = split_fasta_from_fasta(os.path.abspath(self.fasta_path), "user_fasta_split_"+self.outname, self.outdir)
        self.outname = "classified"

    def __exit__(self):
        for path in self.fasta_paths:
            if os.path.exists(path):
                os.remove(path)

    def set_output_path(self, outdir):
        self.outdir = outdir

    def set_output_name(self, outname):
        self.outname = outname

    def run(self, numbering_scheme, concat_out_data = False):
        chain_set = IgChainSet()
        ig_creator = IgChainFactory()
        ab_creator = AbChainFactory()

        if not concat_out_data:
            outpath = self.outdir+"/"+self.outname+"_"+os.path.basename(self.fasta_path).split(".")[0]+".txt"
        else:
            outpath = self.outdir+"/"+self.outname+".txt"

        base_file_name = ".".join(os.path.basename(self.fasta_path).split(".")[0:-1])

        if not os.path.exists(outpath) or not concat_out_data:
            OUTFILE = open(outpath, 'w')
            OUTFILE.write("#file CHAIN original_chain identification is_scfv ig_domains\n")
            OUTFILE.write("#file DOMAIN original_chain identification gene start end score evalue")
        else:
            OUTFILE = open(outpath, 'a')

        for path in self.fasta_paths:
            FASTA = open_file(path)
            for line in FASTA:
                OUTFILE.write("# "+line)
            FASTA.close()
            proposed_chain = ProposedIgChain(path)
            ig_chain = ig_creator.create_ig_vset_chain(proposed_chain)
            if not ig_chain:
                OUTFILE.write(base_file_name+" CHAIN "+proposed_chain.get_id() + " non_ig NA NA\n")
                continue

            OUTFILE.write(ig_chain.print_out(base_file_name))
            OUTFILE.write("\n")
            ab_chain = ab_creator.create_ab_vset_chain(ig_chain, numbering_scheme)
            if ab_chain:
                OUTFILE.write(str(ab_chain))
                OUTFILE.write("\n")

        OUTFILE.close()


        for p in self.fasta_paths:
            if os.path.exists(p):
                os.remove(p)