
from Bio.SearchIO._model import Hit, HSP

class IgDomain(object):
    """
    Represents the domain identification from a ProposedIgChain.  Part of an IgChain/AbChain.
    """
    def __init__(self, hsp, gene, new_chain, passing_domain_hmms):
        assert isinstance(hsp, HSP)
        self.hsp = hsp
        self.gene = gene
        self.new_chain = new_chain
        self.passing_domain_hmms = passing_domain_hmms
        
    def get_hsp(self):
        return self.hsp
    
    def get_passing_hsp(self, hmm):
        if not self.passing_domain_hmms.has_key(hmm):
            return None
        else:
            return self.passing_domain_hmms[hmm]
    
    
    def get_sequence(self):
        pass
    
    def get_gene(self):
        return self.gene

    def get_new_chain(self):
        return self.new_chain