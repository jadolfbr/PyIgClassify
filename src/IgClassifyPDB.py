#!/usr/bin/env python

#Author Jared Adolf-Bryfogle

import random

#Python Imports
import os
import sys

from Bio.PDB import PDBIO

#Append Python Path
p = os.path.split(os.path.abspath(__file__))[0]+"/src"
sys.path.append(p); #Allows all modules to use all other modules, without needing to update pythonpath

#Project Imports
from src.modules.chains.ProposedIgChain import ProposedIgChain
from src.modules.chains.IgChainSet import IgChainSet
from src.modules.chains.IgChainFactory import IgChainFactory
from src.modules.chains.AbChainFactory import AbChainFactory
from src.modules.CDRs.CDRClusterer import CDRClusterer
from src.modules.CDRs.util import get_norm_distance
from src.modules.CDRs.util import get_norm_distance_deg
from src.modules.BioPose import BioPose
from src.modules.Structure import Antibody_Structure

from src.tools import output
from src.tools.fasta import *



class IgClassifyPDB:
    """
    Identify CDRs from a PDB and renumber
    """

    def __init__(self, pdb_path):
        self.pdb_path = pdb_path


        self.bio_pose = BioPose(pdb_path)
        self.bio_structure = self.bio_pose.structure()
        self.ab_structure = Antibody_Structure()

        #Defaults
        self.outdir = os.getcwd()
        self.outname = "classified"

        self.pwd_dir = os.getcwd()

    def set_output_path(self, outdir):
        self.outdir = outdir

    def set_output_name(self, outname):
        self.outname = outname

    def identify_cdrs(self):
        """
        Identifies clusters of CDRs
        """
        clusterer = CDRClusterer(self.bio_pose)

        data = defaultdict()
        for cdr in self.ab_structure.cdr_names:
            data[cdr] = clusterer.get_fullcluster(cdr)
        return data

    def run(self, numbering_scheme, concat_out_data = False, skip_renumber = False):
        """
        Identifies clusters of CDRs and renumbers the PDB
        """

        citation = "#Please Reference North, B., A. Lehmann, and R.L. Dunbrack, Jr., A new clustering of antibody CDR loop conformations. JMB, 2011. 406(2): p. 228-56.\n"
        cluster_header = "#file CDR type old_resnum old_chain new_resnum new_chain sequence cluster distance normDis normDisDeg\n\n"

        base_file_name = ".".join(os.path.basename(self.pdb_path).split(".")[0:-1])

        outpath = ""

        if skip_renumber:
            print "Skipping renumbering step"

        if not concat_out_data:
            if not skip_renumber:
                outpath = self.outdir+"/"+self.outname+"_"+os.path.basename(self.pdb_path).split(".")[0]+".txt"
            clusters_out = self.outdir+"/"+self.outname+"_"+os.path.basename(self.pdb_path).split(".")[0]+"_clusters.txt"
        else:
            if not skip_renumber:
                outpath = self.outdir+"/"+self.outname+".txt"
            clusters_out = self.outdir+"/"+self.outname+"_clusters.txt"

        if not os.path.exists(clusters_out) or not concat_out_data:
            CLUSTERS = open(clusters_out, 'w')
            CLUSTERS.write(citation)
            CLUSTERS.write("\n")
            CLUSTERS.write(cluster_header)
        else:
            CLUSTERS = open(clusters_out, 'a')

        if not skip_renumber and ( not os.path.exists(outpath) or not concat_out_data):
            OUTFILE = open(outpath, 'w')
            OUTFILE.write(citation)
            OUTFILE.write("\n")
            OUTFILE.write("#file CHAIN original_chain identification is_scfv ig_domains\n")
            OUTFILE.write("#file DOMAIN original_chain identification new_chain gene score evalue\n")
            OUTFILE.write(cluster_header)
        elif not skip_renumber:
            OUTFILE = open(outpath, 'a')


        clusterer = CDRClusterer(self.bio_pose)
        if skip_renumber:
            for cdr_name in self.ab_structure.cdr_names:
                chain = cdr_name[0]
                start = self.ab_structure.get_CDR(cdr_name).get_pdb_start()
                end = self.ab_structure.get_CDR(cdr_name).get_pdb_end()
                seq = self.ab_structure.get_cdr_seq(self.bio_pose, cdr_name)
                clusterer.set_dihedrals_from_cdr(cdr_name, cdr_name[0])
                pair = clusterer.get_fullcluster(cdr_name)
                cluster = pair[0]
                distance = pair[1]
                norm_dis = get_norm_distance(self.ab_structure.get_cdr_length(self.bio_pose, cdr_name), distance)
                norm_dis_deg = get_norm_distance_deg(norm_dis)
                line = "CDR "+cdr_name+"\t"+repr(start)+"\t"+chain+"\t"+repr(start)+"\t"+chain+"\t"+ \
                        seq+"\t"+cluster+"\t%.4f"%distance+"\t%.4f"%norm_dis+"\t%.2f"%norm_dis_deg

                CLUSTERS.write(line+"\n")
            CLUSTERS.close()
            return



        self.fasta_paths = chain_fasta_files_from_biostructure(self.bio_structure, "user_fasta_from_pdb_"+self.outname+"_"+str(random.getrandbits(128)), self.outdir)

        chain_set = IgChainSet()
        ig_creator = IgChainFactory()
        ab_creator = AbChainFactory()



        chains = defaultdict()
        for path in self.fasta_paths:
            """
            FASTA = open(path, 'r')
            for line in FASTA:
                OUTFILE.write("# "+line)
            OUTFILE.write("\n")
            FASTA.close()
            """
            proposed_chain = ProposedIgChain(path)
            #print proposed_chain; #Uncomment this for Debugging of Evalues and Cutoffs!  Will print all of them for each HMM!
            ig_chain = ig_creator.create_ig_vset_chain(proposed_chain)
            if not ig_chain:
                OUTFILE.write(base_file_name+" CHAIN "+proposed_chain.get_id() + " non_ig NA NA\n")
                continue
            else:
                chains[ig_chain] = None

            ab_chain = ab_creator.create_ab_vset_chain(ig_chain, numbering_scheme)
            if ab_chain:


                renumbering.renumber_biopython_chain(self.bio_structure, ab_chain)
                chains[ig_chain] = ab_chain

        temp_pdb_file = self.outdir+"/"+self.outname+"_temp"+os.path.basename(self.pdb_path).split("-")[-1]
        new_pdb_file = self.outdir+"/"+self.outname+"_"+os.path.basename(self.pdb_path)

        dumper = PDBIO()
        dumper.set_structure(self.bio_structure)
        dumper.save( open_file(temp_pdb_file, 'w') )

        #Reorder to be LH_A:
        output.reorder_and_save_chains(temp_pdb_file, new_pdb_file, False)
        os.remove(temp_pdb_file)

        if not os.path.exists(self.outdir+"/"+self.outname+"_"+os.path.basename(self.pdb_path)):
            print "NEW file not found!"
        else:
            print self.outdir+"/"+self.outname+"_"+os.path.basename(self.pdb_path)

        new_pose = BioPose(self.outdir+"/"+self.outname+"_"+os.path.basename(self.pdb_path))
        print "\n"+cluster_header.strip("\n").replace("#", "HEADER ")
        for ig_chain in chains:
            OUTFILE.write(ig_chain.print_out(base_file_name))

            if chains[ig_chain]:
                ab_chain = chains[ig_chain]
                ab_chain.identify_clusters(new_pose)
                for cdr in ab_chain.get_cdrs():
                    outline = "CLUSTER "+base_file_name+" "+cdr.get_pdb_output(ab_chain.sequence)
                    print outline
                    OUTFILE.write(outline+"\n")
                    CLUSTERS.write(outline+"\n")
                OUTFILE.write("\n")
                OUTFILE.write(ab_chain.get_pdb_print(base_file_name))
                OUTFILE.write("\n")

        print ""
        OUTFILE.close()
        CLUSTERS.close()

        for p in self.fasta_paths:
            if os.path.exists(p):
                os.remove(p)

        print "UserPyIgClassify complete"