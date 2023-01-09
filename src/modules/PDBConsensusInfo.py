import sys
from collections import defaultdict

from modules.restype_definitions import definitions

from src.modules.Structure import PDBResInfo
from src.modules.Structure import Residue
from src.modules.Structure import ResidueRegion
from src.tools import fasta


class PDBConsensusInfo():
    """
    Class to compute frequency and probability from an array of PDBResInfo classes.
    The sequences within PDBResInfo do not necessarily need to be the same length.
    A given sequence position is identified and stored in the data maps by its [pdb_num, chain, and icode]
    -> Use get_position_from_residue(residue) to get this position from a Residue instance.
    """
    def __init__(self, resinfo_list):

        self.pdb_info_list = resinfo_list
        
        if len(self.pdb_info_list)==0:
            print "No sequences found."
            return

        self.stats = defaultdict(); #Vector of a map of each amino acid
        self.freq = defaultdict(); #Same as stats

        self.aas = definitions().get_all_one_letter_codes()

        #-> If an entry has an unknown residue (most likely X), skip the full sequence or just skip that position
        self.skip_full_sequence_of_unknown_aa = False

        self.init_data_map()
        self.compute_stats()

    def _get_total_entries(self, position):
        totals = 0
        for pdb_info in self.pdb_info_list:
            if pdb_info.total_residue() == 0:
                print "Skipping PDBResInfo for stats.  Has no residues"
                continue;
            for residue in pdb_info.get_all_residues():
                entry_position = self.get_position_from_residue(residue)
                if entry_position == position:
                    totals+=1
        return totals

    def set_sequences(self, pdb_info_list):
        """
        Set a sequence list
        """

        self.pdb_info_list = pdb_info_list
        self.init_data_map()
        self.compute_stats()

    def output_seqlogo(self, outdir, outname, clustalpath = None):

        sequences = []
        for pdb_res_info in self.pdb_info_list:
            assert isinstance(pdb_res_info, PDBResInfo)
            if pdb_res_info.total_residue() == 0:
                #print "PDBResInfo has no residues to get sequence!"
                continue;
            sequences.append(pdb_res_info.get_sequence())

        fasta.output_weblogo_for_sequences(sequences, outdir, outname)

    def output_seqlogo_bt_residues(self, outdir, outname, res1, res2, chain):
        if not isinstance(res1, Residue): sys.exit()
        if not isinstance(res2, Residue): sys.exit()
        sequences = []
        for pdb_res_info in self.pdb_info_list:
            assert isinstance(pdb_res_info, PDBResInfo)
            if pdb_res_info.total_residue() == 0:
                #print "PDBResInfo has no residues to get sequence!"
                continue;

            seq = pdb_res_info.get_sequence_bt_residues(res1, res2, chain)
            if not seq: continue
            sequences.append(seq)

        fasta.output_weblogo_for_sequences(sequences, outdir, outname)

    def output_seqlogo_for_regions(self, regions, outdir, outname, chain):
        """
        Regions is an array of Regions classes.  Basically start/stop points
        """
        sequences = []
        for pdb_res_info in self.pdb_info_list:
            print pdb_res_info.get_sequence()
            assert isinstance(pdb_res_info, PDBResInfo)
            if pdb_res_info.total_residue() == 0:
                #print "PDBResInfo has no residues to get sequence!"
                continue

            seq = ""
            for region in regions:
                #print repr(region)
                #print repr(region.res1)+" "+repr(region.res2)
                if not isinstance(region, ResidueRegion): sys.exit()
                seq = seq + pdb_res_info.get_sequence_bt_residues(region.res1, region.res2, chain)

            if not seq: continue
            sequences.append(seq)

        fasta.output_weblogo_for_sequences(sequences, outdir, outname)

    def get_all_sorted_positions(self):
        return sorted(self.freq.keys())

    def get_consensus_for_position(self, position):

        if not self.stats.has_key(position):
            return ""

        aa_letter = ""
        for aa in self.stats[position]:
            prob = self.stats[position][aa]
            if prob > .9:
                aa_letter = aa.upper()
                break
            elif prob > .2:
                aa_letter = aa.lower()
                break
            else:
                aa_letter = '-'
                continue

        return aa_letter

    def get_consensus(self, residue):
        position = self.get_position_from_residue(residue)
        return self.get_consensus_for_position(position)

    def get_consensus_sequence(self):
        consensus = ""
        for position in sorted(self.freq.keys()):
            aa_letter = self.get_consensus(position)

            consensus = consensus+aa_letter

        return consensus

    def get_consensus_for_residues(self, residue_list):
        """
        Get the consensus for an ORDERED list of Residues
        """
        consensus = ""
        for residue in sorted(residue_list):
            aa_letter = self.get_consensus(residue)

            consensus = consensus+aa_letter

        return consensus

    def get_frequency_for_position(self, position, aa):
        if not self.freq.has_key(position):
            return 0

        if not self.freq[position].has_key(aa):
            return 0

        return self.freq[position][aa]

    def get_frequency(self, residue, aa):
        position = self.get_position_from_residue(residue)
        return self.get_frequency_for_position(position, aa)

    def get_probability_for_position(self, position, aa):
        if not self.stats.has_key(position):
            return 0
        if not self.stats[position].has_key(aa):
            return 0

        return self.stats[position][aa]

    def get_probability(self, residue, aa):
        '''
        Get probability of the current position (starting from 0) and aa
        '''

        position = self.get_position_from_residue(residue)
        return self.get_probability_for_position(position, aa)


    ##Private Methods##
    def compute_stats(self):
        """
        Compute frequency and probability (0-1) for each position for each amino acid
        """

        if len(self.pdb_info_list)==0:return

        #print "Total Unique CDR Sequences: "+repr(len(self.sequence_list))

        for pdb_info in self.pdb_info_list:
            if not isinstance(pdb_info, PDBResInfo): sys.exit()
            if pdb_info.total_residue() == 0:
                print "PDBResInfo has no residues.  Popping"
                self.pdb_info_list.remove(pdb_info)
                continue

            for residue in pdb_info.get_all_residues():
                if not isinstance(residue, Residue): sys.exit()
                position = self.get_position_from_residue(residue)

                if residue.aa not in self.aas:
                    print "Unknown amino acid: "+residue.aa


                    if self.skip_full_sequence_of_unknown_aa:
                        print "Skipping full sequence. Popping. "
                        self.pdb_info_list.remove(pdb_info)
                        break


        #If we have removed sequences with unknown residues or missing residues:
        for pdb_info in self.pdb_info_list:
            if not isinstance(pdb_info, PDBResInfo): sys.exit()
            for residue in pdb_info.get_all_residues():
                if not isinstance(residue, Residue): sys.exit()
                position = self.get_position_from_residue(residue)

                if residue.aa not in self.aas:
                    print "Unknown amino acid: "+residue.aa
                    print "Skipping"
                    continue

                self.freq[position][residue.aa]+=1

        for position in sorted(self.freq.keys()):
            total_entries = self._get_total_entries(position)
            for aa in self.freq[position].keys():
                self.stats[position][aa] = self.freq[position][aa]/float(total_entries)
                #print repr(position)+" "+aa+" "+ repr(self.freq[position][aa]/float(total_entries))

        #print repr(sorted(self.stats.keys()))

    def init_data_map(self):
        """
        Sets all probabilities 0 and appends each map to the stats vector
        """

        self.stats = defaultdict()
        self.freq = defaultdict()
        for pdb_info in self.pdb_info_list:
            if not isinstance(pdb_info, PDBResInfo): sys.exit()
            for residue in pdb_info.get_all_residues():
                if not isinstance(residue, Residue): sys.exit()

                resnum_triple = self.get_position_from_residue(residue)
                if self.stats.has_key(resnum_triple):
                    continue


                prob = defaultdict()
                freq = defaultdict()
                for aa in self.aas:
                    prob[aa] = 0
                    freq[aa] = 0
                self.stats[resnum_triple] = prob
                self.freq[resnum_triple] = freq

    def return_initialized_total_map(self):
        totals = dict()
        for aa in self.aas:
            totals[aa] = 0

        return totals

    def get_position_from_residue(self, residue):
        if not isinstance(residue, Residue): sys.exit()
        triple = (residue.get_pdb_num(), residue.get_chain(), residue.get_icode())
        return triple


    
    