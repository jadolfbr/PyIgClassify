#Author: Jared Adolf-Bryfogle


from .IgChain import IgChain


class IgChainSet(object):
    """
    Container and accessor of IgChains
    """
    
    def __init__(self, igchains = []):
        assert type(igchains) is list
        self.igchains = igchains
        
    def add_chain(self, igchain):
        
        
        if not igchain:return
        assert isinstance(igchain, IgChain)
        
        self.igchains.append(igchain)
    
    def set_chains(self, hmm_chain_type_array):
        self.igchains = hmm_chain_type_array
        
    def get_chains(self):
        return self.igchains
    
    def get_ig_ab_chains(self):
        """
        Return all Antibody IG Chains
        """
        
        IgChains = []
        
        for chain in self.igchains:
            if chain.get_type() == "ig_ab":
                IgChains.append(chain)
        return IgChains
    
    def get_ig_tcr_chains(self):
        """
        Return all TCR IG Chains
        """
        
        TcrChains = []
        
        for chain in self.igchains:
            if chain.get_type() == "ig_tcr":
                TcrChains.append(chain)
        return TcrChains
    
    def get_light_chains(self):
        """
        Return All chains with 'L' designation.
        """
        
        chains = []
        
        for chain in self.igchains:
            if chain.get_chain() == 'L':
                chains.append(chain)
        return chains
    
    def get_ScFv_chains(self):
        chains = []
        for chain in self.igchains:
            if chain.is_ScFv():
                chains.append(chain)
        return chains
    
    def get_heavy_chains(self):
        """
        Return All chains with 'H' designation.
        """
        
        chains = []
        
        for chain in self.igchains:
            if chain.get_chain() == 'H':
                chains.append(chain)
        return chains
    
    
    