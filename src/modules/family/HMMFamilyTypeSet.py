
#Author: Jared Adolf-Bryfogle

from HMMFamilyType import HMMFamilyType

class HMMFamilyTypeSet(object):
    """
    A set of HMMFamilyTypes
    """
    def __init__(self, hmmchains=[]):
        assert type(hmmchains) is  list
        self.HMMFamilies = hmmchains
    
    def add_family_type(self, hmm_family_type):
        assert isinstance(hmm_family_type, HMMFamilyType)
    
        self.HMMFamilies.append(hmm_family_type)
    
    def set_chain_types(self, hmm_chain_type_array):
        assert type(hmm_chain_type_array) is list
        
        self.HMMFamilies = hmm_chain_type_array
        
    def get_hmm_chains(self):
        return self.HMMFamilies
    
    def get_hmm_family(self, hmm_name):
        for fam in self.HMMFamilies:
            if fam.get_hmm_name() == hmm_name:
                return fam.get_family()
    
    
    