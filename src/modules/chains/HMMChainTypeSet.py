#Author: Jared Adolf-Bryfogle

from .HMMChainType import HMMChainType

class HMMChainTypeSet(object):
    """
    A set of HMMChainTypes
    """
    def __init__(self, hmmchains=[]):
        assert type(hmmchains) is list
        
        self.HMM_chain_types = hmmchains
    
    def add_chain_type(self, hmm_chain_type):
        assert isinstance(hmm_chain_type, HMMChainType)
        
        self.HMM_chain_types.append(hmm_chain_type)
    
    def set_chain_types(self, hmm_chain_type_array):
        assert type(hmm_chain_type_array) is list
        self.HMM_chain_types = hmm_chain_type_array
        
    def get_hmm_chains(self):
        return self.HMM_chain_types
    
    def get_hmm_chain(self, hmm_name):
        """
        Get a specific HMMChainType by hmm_name
        """
        for chain_type in self.HMM_chain_types:
            if chain_type.get_hmm_name() == hmm_name:
                return chain_type
            
    def get_ig_ab_chains(self):
        """
        Get all ig_ab ChainTypes
        """
        
        igab_chain_types = []
        
        for chain_type in self.HMM_chain_types:
            if chain_type.get_family() == "ig_ab":
                igab_chain_types.append(chain_type)
        return igab_chain_types
    
    def get_ig_tcr_chains(self):
        """
        Get all ig_tcr ChainTypes
        """
        
        ig_tcr_chain_types = []
        
        for chain_type in self.HMM_chain_types:
            if chain_type.get_family() == "ig_tcr":
                ig_tcr_chain_types.append(chain_type)
        return ig_tcr_chain_types
    
    def get_len(self):
        return len(self.HMM_chain_types)
    
    