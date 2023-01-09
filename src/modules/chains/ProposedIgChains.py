from collections import defaultdict

from Bio import SeqIO
from src.tools.general import *

from .ProposedIgChain import ProposedIgChain
from src.tools.path import *


class ProposedIgChains(ProposedIgChain):
    """
    Specific class of ProposedIgChain for creating a database from the PDBAA or set of entries in a single fasta.  Holds all results.  Uses specific function in IgChainFactory.
    """
    def __init__(self, pdbaa_path):
        ProposedIgChain.__init__(self, pdbaa_path)
        self.sequences = defaultdict()
        self._create_sequences()
        
        self.organized_results = defaultdict(defaultdict)
        self._organize_hits()
        
        self.ids = []
        self._create_ids()
        
    def _create_ids(self):
        chains = defaultdict()
        
        for hmm in self.get_hmms():
            hits = self.get_hits(hmm)
            for hit in hits:
                if not hit: continue
                tag = hit.id
                chains[tag] = ""
                
        self.ids = sorted(chains.keys())
                
    def _organize_hits(self):
        for hmm in self.get_hmms():
            hits = self.get_hits(hmm)
            for hit in hits:
                self.organized_results[hmm][hit.id]=hit
                
    def _create_sequences(self):
        FASTA = open_file(self.fasta, 'r')
        records = SeqIO.parse(FASTA, 'fasta')
        for record in records:
            self.sequences[record.id] = record
        
        FASTA.close

    def get_ids(self):
        return self.ids

    def get_PDBID_Chains(self):
        """
        Get the id - in the PDBAA, this id (>xxx) is PDBID_CHAIN
        """
        return self.ids
    
    def get_sequence(self, id):
        return str(self.sequences[id].seq)
    
    def get_full_id(self, id):
        return id+" "+self.sequences[id].description
    
    def get_protein_name(self, id):
        return self.sequences[id].description

    def get_description(self, id):
        return self.sequences[id].description

    def get_hits_for_hmms(self, id):
        """
        Returns a dictionary of hmm:hit
        """
        result = defaultdict()
        for hmm in self.get_hmms():
            if self.organized_results[hmm].has_key(id):
                result[hmm] = self.organized_results[hmm][id]
            else:
                continue
            
        return result
    
    def get_hit_for_hmm(self, id, hmm):
        """
        Returns the hit for a given HMM and id
        """
        
        try:
            return self.organized_results[hmm][id]
        except KeyError:
            return None
        
    def get_hsps(self, id, hmm):
        
        hsps = []
        hit = self.get_hit_for_hmm(id, hmm)
        if not hit:
            return []
        
        for hsp in hit:
            hsps.append(hsp)
        return hsps
    
    def get_families(self, id):
        hits = []
        families = self.hmm_family_set.get_hmm_chains()
        per_domain_results = []
        for family_type in families:
            hmm = family_type.get_hmm_name()
            query = self.family_results[hmm]['query_result']
            for hit in query:
                if not hit: continue
                tag = hit.id
                query_id = tag
                if query_id == id:
                    for hsp in hit:
                        #Check Evalue - If not passed, continue
                        if hit.bitscore < family_type.get_score_cutoff(): continue
                        elif hit.evalue > family_type.get_evalue_cutoff(): continue
                        
                        results = defaultdict(dict)
                        results[hmm]['hsp'] = hsp
                        results[hmm]['evalue'] = hit.evalue
                        results[hmm]['score'] = hit.bitscore
                        per_domain_results.append(results)                  
        
        if len(per_domain_results) == 0:
            return None
    
        else:
            best_families = []
            for results in per_domain_results:
            
                best_scoring_hmm = get_best_score_hmm(results)
        
                best_evalue_hmm = get_best_evalue_hmm(results)
                best_families.append(self.hmm_family_set.get_hmm_family(best_scoring_hmm))
            
            return best_families
            