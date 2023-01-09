




class HMMFamilyType(object):
    """
    A container and accessor class for a chain-based HMM model. Analogous to the HMMChainType class.  These types are used for the identication
    of specific HMM families of a chain or set of chains- such as the ig_vset or ig_constant regions.
    """
    
    def __init__(self, hmm_name, family, score_cutoff, evalue_cutoff):
        self.hmm_model = hmm_name
        self.family = family
        self.score_cutoff = score_cutoff
        self.evalue_cutoff = evalue_cutoff
        
        assert type(self.score_cutoff) is float;
        assert type(self.evalue_cutoff) is float;
        
    def get_hmm_path(self):
        return self.hmm_model;
    
    def get_family(self):
        return self.family;
    
    def get_score_cutoff(self):
        return self.score_cutoff
    
    def get_evalue_cutoff(self):
        return self.evalue_cutoff
    
    def get_hmm_name(self):
        return self.hmm_model
    