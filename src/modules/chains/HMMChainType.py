#Author Jared Adolf-Bryfogle


class HMMChainType(object):
    """
    A container and accessor class for a chain-based HMM model.  The model has a name corresponding to it's relative path in the database.
    It has a type (ig_tcr/ig_ab), chainID for renumbering, a gene, as well as score and evalue cutoffs.
    These values are read by the HMMChainTypeParser from /database/hmm_chain_types.txt
    
    The HMMChainType is then used to create an IgChain from a ProposedIgChain by comparing the HMM Results for each model with the corresponding data in the HMMChainType.
    The IgChain is then created by using the data (gene, chainID, etc.) of the HMMChainType.
    """
    
    def __init__(self, hmm_name, chain_type, chainID, gene, score_cutoff, evalue_cutoff):
        self.hmm_model = hmm_name
        self.chain_type = chain_type
        self.chainID = chainID
        self.gene = gene
        self.score_cutoff = score_cutoff
        self.evalue_cutoff = evalue_cutoff
        
        assert type(self.score_cutoff) is float;
        assert type(self.evalue_cutoff) is float;
        
    def get_hmm_name(self):
        return self.hmm_model;
    
    def get_type(self):
        return self.chain_type;
    
    def get_chain(self):
        return self.chainID
    
    def get_score_cutoff(self):
        return self.score_cutoff
    
    def get_evalue_cutoff(self):
        return self.evalue_cutoff
    
    def get_gene(self):
        return self.gene