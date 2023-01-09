import re
#from CreateDesignList import CreateDesignList
#from CreateRelaxList import CreateRelaxList
from src.tools.path import *
from src.tools import general
from collections import defaultdict


class CSVPrune:
    """
    Prunes raw full pdbaa csv file for nr_by_cdr and nr_by_cdr_per_chain data.  Writes two csv files for these, which are the input to CreateCDRDatabase class.
    Started as a script, and never fully got rewritten better.
    """
    
    def __init__(self, outdir, fasta_data, res_cutoff = 2.8, rfac_cutoff=.3):
        """
        fasta_data is dictionary loaded from tools/fasta: pdb_chain: [method, residues, resolution, R factor]
        """
        self.res_cutoff = res_cutoff
        self.rfac_cutoff = rfac_cutoff

        self.fasta_data = fasta_data
        self.resolution_dict = self._convert_fasta_data_to_resolution_dict()
        self.outdir = outdir
        self.logdir = outdir+"/logs"
        if not os.path.exists(self.outdir):
            os.mkdir(self.outdir)
        if not os.path.exists(self.logdir):
            os.mkdir(self.logdir)

        self.ben_input_csv = get_db_path()+"/csv/AbMatchToPaperClusters.csv"
        self.load_centers()

    def _create_new_csv_file(self, to_keep_dict, filename):
        """
        to_keep_dict is [pdb_chain_cdr] = [''
        """

        #If it is a center member, add it always. Add to to_keep_dict if not already present
        for pdb_chain in self.center_data:
            for cdr in self.center_data[pdb_chain]:
                pdb_chain_cdr = pdb_chain+"_"+cdr
                if not to_keep_dict.has_key(pdb_chain_cdr):
                    to_keep_dict[pdb_chain_cdr] = ' '


        #Write data according to pdb_chain_cdr in to_keep_dict
        FILE = open_file(self.ben_input_csv, 'r')
        FILEOUT = open_file(get_db_path()+'/csv/'+filename, 'w')
        for line in FILE:
            if line[0]=='#':
                FILEOUT.write(line)
                continue
            lineSP = line.split(',')
            
            if to_keep_dict.has_key(lineSP[0].upper()):
                FILEOUT.write(line)
        
        FILE.close()
        FILEOUT.close()
        
    def _convert_fasta_data_to_resolution_dict(self):
        """
        Converts fasta_data to resolution dict, which has NA values replaced by numbers.
        fasta_data is dictionary loaded from tools/fasta: pdb_chain: [residues, resolution, R factor]
        """
        
        #print self.fasta_data
        resolution_dict = defaultdict(dict)
        for pdb_chain in self.fasta_data:
            #Uses the best structure, but it seems it could be NMR data/etc.  Need to look into this.
            #If the FASTA is missing residue/resolution/R factor. This is bad.  We don't want this.
            method = self.fasta_data[pdb_chain][0]
            residues = int(self.fasta_data[pdb_chain][1]);
            resolution = self.fasta_data[pdb_chain][2];
            if resolution == 'NA': resolution = 100
            rfactor = self.fasta_data[pdb_chain][3];
            if rfactor == 'NA': rfactor = 100
            
            
            #resolution_dict[pdb_chain] = [resolution, rfactor, residues, method]
            resolution_dict[pdb_chain]["resolution"] = resolution
            resolution_dict[pdb_chain]["rfactor"] = rfactor
            resolution_dict[pdb_chain]["pdb_length"] = residues
            resolution_dict[pdb_chain]["method"] = method

        return resolution_dict

    def load_centers(self, filepath=False):
        """
        Loads center data from Ben's median output file.
        """
        if not filepath:
            FILE = open_file(get_db_path()+"/medians.txt")
        else:
            FILE = open_file(filepath)

        FILE.readline()
        self.center_data = defaultdict(dict)

        #New Style:
        for line in FILE:
            line = line.strip()
            lineSP = line.split(',')
            pdb_chain = lineSP[2][0:5].upper()
            pdb = pdb_chain[0:4].upper()
            cdr = lineSP[2].split("_")[1].upper()
            #Cluster Identification (Cis/Trans):
            SS = lineSP[0].split("-")[2]

            length_name = general.get_loop_key_from_ss(SS)
            cluster = cdr+"-"+length_name+"-"+lineSP[1]
            self.center_data[pdb_chain][cdr] = cluster
        FILE.close()

    def _passes_cutoffs(self, pdb_chain):
        """
        Checks the cutoffs of res and rfactor according to the res dict and returns a boolean
        """
        resolution = float(self.resolution_dict[pdb_chain]["resolution"])
        rfactor = float(self.resolution_dict[pdb_chain]["rfactor"])

        if (resolution <= self.res_cutoff) and (rfactor <= self.rfac_cutoff):
            return True
        else:
            return False

    def prune_for_nr_by_cdr_database(self, paper_only=False):

        outname = ""
        if paper_only:
            outname = "data_for_nr_by_cdr_with_outliers_paper_only.csv"
        else:
            outname = "data_for_nr_by_cdr_with_outliers.csv"

        to_keep = defaultdict(); #Dict so it's easy to search.
        double_check = self.logdir+"/structures_to_double_check_for_nr_by_cdr.txt"
        DATA = open_file(self.ben_input_csv, 'r')
        DATA.readline()
        DOUBLECHECK = open_file(double_check, 'w')
        sequence_dict = defaultdict()
        
        for line in DATA:

            line = line.strip()

            if not line: continue
            if line.startswith('#'): continue

            lineSP = line.split(',')

            paper = lineSP[2]
            pdb_chain = lineSP[0].split('_')[0].upper()

            if paper=="missingSeq":
                continue
            if paper_only and paper != "inPaper":
                continue

            pdb = pdb_chain[0:4].upper(); chain=pdb_chain[5:5]
            cdr = lineSP[0].split('_')[1]
            sequence = lineSP[1]
            sequence = sequence+cdr
            length = lineSP[3].split('-')[1]
            SS = lineSP[3].split('-')[2]
            cluster = lineSP[4]
            distance = lineSP[5]
            
            if lineSP[6]=="?":
                distanceNorm = -1
            else:
                distanceNorm = float(lineSP[6])
            
            #Cluster Fix (Cis/Trans):
            if re.search("C", SS):
                fix = general.get_loop_key_from_ss(SS)
                cluster = fix+'-'+cluster
                
            full_cluster = cdr+'-'+length+'-'+cluster

            if not sequence_dict.has_key(sequence):
                
                if self.resolution_dict.has_key(pdb_chain.upper()):

                    if self._passes_cutoffs(pdb_chain):
                        sequence_dict[sequence] = (self.resolution_dict[pdb_chain]["resolution"], pdb_chain, cdr, distanceNorm)
                    else:
                        continue


                else:
                    #print "pdb_chain not found in fasta List.  Removing."
                    print pdb_chain+ " not found in resolution dict from fasta!."
            #Here we compare, take the lowest one.
            else:
                if self.resolution_dict.has_key(pdb_chain):
                    new_res = self.resolution_dict[pdb_chain]["resolution"]
                    old_res = sequence_dict[sequence][0]
                    #print pdb_chain
                    if new_res<=old_res:
                        #print "Better or Equel to Xray res found.  Checking equivalency."
                        if new_res != old_res:
                            sequence_dict[sequence] = (self.resolution_dict[pdb_chain]["resolution"], pdb_chain, cdr, distanceNorm)
                        else:
                            #print "Resolutions are identical.  Checking R factor."
                            new_r = self.resolution_dict[pdb_chain]["rfactor"]
                            old_r = self.resolution_dict[sequence_dict[sequence][1]]["rfactor"]
                            if new_r < old_r:
                                sequence_dict[sequence] = (self.resolution_dict[pdb_chain]["resolution"], pdb_chain, cdr, distanceNorm)
                            elif new_r ==old_r:
                                #print "R factors are identical!  What to do now?"
                                #print "Checking distanceNorm, taking the structure closer to known structure..."
                                new_d = distanceNorm
                                old_d = sequence_dict[sequence][3]
                                if new_d < old_d:
                                    sequence_dict[sequence] = (self.resolution_dict[pdb_chain]["resolution"], pdb_chain, cdr, distanceNorm)
                                elif new_d==old_d:
                                    #print "distances are Identical.  Need to check if the PDBS have been replaced.  Printing this to a list. KEEPING First Found"
                                    new_pdb = pdb_chain[0:4]
                                    old_pdb = sequence_dict[sequence][1][0:4]
                                    if new_pdb==old_pdb:
                                        new_length = self.resolution_dict[pdb_chain]["pdb_length"]
                                        old_length = self.resolution_dict[sequence_dict[sequence][1]]["pdb_length"]
                                        if new_length> old_length:
                                            sequence_dict[sequence] = (self.resolution_dict[pdb_chain]["resolution"], pdb_chain, cdr, distanceNorm)
                                        elif new_length==old_length:
                                            DOUBLECHECK.write("SAMEID - Same Length. CHECK FOR CRYSTAL CONTACTS AND ENERGY. Last PDB_CDR used. : "+pdb_chain+'_'+cdr+' '+sequence_dict[sequence][1]+'_'+sequence_dict[sequence][2]+"\n")
                                    else:
                                        DOUBLECHECK.write("SAME - CHECK FOR PDB REPLACEMENT : "+pdb_chain+'_'+cdr+' '+sequence_dict[sequence][1]+'_'+sequence_dict[sequence][2]+"\n")
                    else:
                        #print "new_resolution not lower then old."
                        pass
        DATA.close()
        DOUBLECHECK.close()
        
        
        for sequence in sequence_dict:
            pdb_chain = sequence_dict[sequence][1]; cdr = sequence_dict[sequence][2]; distanceNorm = sequence_dict[sequence][3]; res = sequence_dict[sequence][0]
            pdb_chain_cdr = pdb_chain+ '_'+cdr
            to_keep[pdb_chain_cdr] = []
        self._create_new_csv_file(to_keep, outname)
        
    def prune_for_nr_by_cdr_seq_database(self, paper_only = False):

        outname = ""
        if paper_only:
            outname = "data_for_nr_by_cdr_seq_with_outliers_paper_only.csv"
        else:
            outname = "data_for_nr_by_cdr_seq_with_outliers.csv"

        to_keep = defaultdict(); #Dict so it's easy to search.
        double_check = self.logdir+"/structures_to_double_check_for_nr_by_cdr_seq.txt"
        DATA = open_file(self.ben_input_csv, 'r')
        DATA.readline()
        DOUBLECHECK = open_file(double_check, 'w')
        sequence_dict = defaultdict()
        sequences = defaultdict()
        
        for line in DATA:

            line = line.strip()
            if not line: continue
            if line.startswith('#'): continue

            lineSP = line.split(',')
            paper = lineSP[2]
            pdb_chain = lineSP[0].split('_')[0].upper()

            if paper=="missingSeq":
                continue
            if paper_only and paper != "inPaper":
                continue

            pdb = pdb_chain[0:4].upper(); chain=pdb_chain[4]
            cdr = lineSP[0].split('_')[1]
            sequence = lineSP[1]
            #sequence = sequence+cdr
            length = lineSP[3].split('-')[1]
            SS = lineSP[3].split('-')[2]
            cluster = lineSP[4]
            distance = lineSP[5]
            
            if lineSP[6]=="?":
                distanceNorm = -1
            else:
                distanceNorm = float(lineSP[6])
            
            #Cluster Fix (Cis/Trans):
            if re.search("C", SS):
                fix = general.get_loop_key_from_ss(SS)
                cluster = fix+'-'+cluster
                
            full_cluster = cdr+'-'+length+'-'+cluster
            chain_type = cdr[0]
            if not self.resolution_dict.has_key(pdb_chain):
                continue
            if not sequences.has_key(pdb_chain):
                if self._passes_cutoffs(pdb_chain):
                    sequences[pdb_chain]=defaultdict(dict)
                else:
                    continue

                #sequences[pdb_chain][chain_type]=defaultdict()
            
            sequences[pdb_chain][chain_type][cdr]=sequence
            
        light_cdrs = ["L1", "L2", "L3"]
        heavy_cdrs = ["H1", "H2", "H3"]
        
        final_sequences = dict()
        #print sequences
        for pdb_chain in sequences:
            #print pdb_chain
            if sequences[pdb_chain].has_key("L"):
                if all (cdr in sequences[pdb_chain]["L"] for cdr in light_cdrs):
                    sequence = sequences[pdb_chain]["L"]["L1"]+sequences[pdb_chain]["L"]["L2"]+sequences[pdb_chain]["L"]["L3"]
                    if not final_sequences.has_key(sequence):
                        final_sequences[sequence]=("L", pdb_chain, self.resolution_dict[pdb_chain]["resolution"], self.resolution_dict[pdb_chain]["rfactor"])
                    else:
                        new_res = self.resolution_dict[pdb_chain]["resolution"]
                        old_res = final_sequences[sequence][2]
                        if new_res < old_res:
                            final_sequences[sequence]=("L", pdb_chain, self.resolution_dict[pdb_chain]["resolution"], self.resolution_dict[pdb_chain]["rfactor"])
                        elif new_res == old_res:
                            #print "Resolutions same.  Checking  R factors..."
                            new_r = self.resolution_dict[pdb_chain]["rfactor"]
                            old_r = final_sequences[sequence][2]
                            if new_r < old_r:
                                final_sequences[sequence]=("L", pdb_chain, self.resolution_dict[pdb_chain]["resolution"], self.resolution_dict[pdb_chain]["rfactor"])
                            elif new_r==old_r:
                                #print "R_Factors Same.  Need more data. Checking Lengths"
                                new_length = self.resolution_dict[pdb_chain]["pdb_length"]
                                old_length = self.resolution_dict[final_sequences[sequence][1]]["pdb_length"]
                                if new_length > old_length:
                                    final_sequences[sequence]=("L", pdb_chain, self.resolution_dict[pdb_chain]["resolution"], self.resolution_dict[pdb_chain]["rfactor"])
                                elif new_length==old_length:
                                    DOUBLECHECK.write("SAMEID - Same Length. CHECK FOR CRYSTAL CONTACTS AND ENERGY : "+pdb_chain+'_'+sequence+"\n")
                            
                else:
                    for cdr in light_cdrs:
                        if not cdr in sequences[pdb_chain]["L"]:
                            DOUBLECHECK.write("MISSING_CDR : "+pdb_chain+'_'+cdr+'\n')
                    
                    
            if sequences[pdb_chain].has_key("H"):
                if all (cdr in sequences[pdb_chain]["H"] for cdr in heavy_cdrs):
                    sequence = sequences[pdb_chain]["H"]["H1"]+sequences[pdb_chain]["H"]["H2"]+sequences[pdb_chain]["H"]["H3"]
                    if not final_sequences.has_key(sequence):
                        final_sequences[sequence]=("H", pdb_chain, self.resolution_dict[pdb_chain]["resolution"], self.resolution_dict[pdb_chain]["resolution"])
                    else:
                        new_res = self.resolution_dict[pdb_chain]["resolution"]
                        old_res = final_sequences[sequence][2]
                        if new_res < old_res:
                            final_sequences[sequence]=("H", pdb_chain, self.resolution_dict[pdb_chain]["resolution"], self.resolution_dict[pdb_chain]["rfactor"])
                        elif new_res == old_res:
                            #print "Resolutions same.  Checking  R factors..."
                            new_r = self.resolution_dict[pdb_chain]["resolution"]
                            old_r = final_sequences[sequence][2]
                            if new_r < old_r:
                                final_sequences[sequence]=("H", pdb_chain, self.resolution_dict[pdb_chain]["resolution"], self.resolution_dict[pdb_chain]["rfactor"])
                            else:
                                #print "R_Factors Same.  Need more data. Checking Lengths"
                                new_length = self.resolution_dict[pdb_chain]["pdb_length"]
                                old_length = self.resolution_dict[final_sequences[sequence][1]]["pdb_length"]
                                if new_length > old_length:
                                    final_sequences[sequence]=("H", pdb_chain, self.resolution_dict[pdb_chain]["resolution"], self.resolution_dict[pdb_chain]["rfactor"])
                                elif new_length==old_length:
                                    DOUBLECHECK.write("SAMEID - Same Length. CHECK FOR CRYSTAL CONTACTS AND ENERGY : "+pdb_chain+'_'+' '+"\n")
                            
                else:
                    for cdr in heavy_cdrs:
                        if not cdr in sequences[pdb_chain]["H"]:
                            DOUBLECHECK.write("MISSING_CDR : "+pdb_chain+'_'+cdr+'\n')
                    DOUBLECHECK.write("MISSING_CDR : "+pdb_chain+'_'+cdr+'\n')
        
        for sequence in final_sequences:
            if final_sequences[sequence][0]=="L":
                for cdr in light_cdrs:
                    to_keep[final_sequences[sequence][1]+'_'+cdr] = ''
            if final_sequences[sequence][0]=="H":
                for cdr in heavy_cdrs:
                    to_keep[final_sequences[sequence][1]+'_'+cdr] = ''
        DOUBLECHECK.close()
        self._create_new_csv_file(to_keep, outname)