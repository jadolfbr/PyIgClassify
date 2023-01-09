#Author Jared Adolf-Bryfogle


import os

from Bio.SearchIO._model import QueryResult, Hit, HSP, HSPFragment

class IgChain(object):
    """
    A chain Identified using the chain hmms in /database.  The IgChain has an identified type (ig_tcr, ig_ab, etc.), as well as gene identifications and the 
    actual result of the HMMER3 search.  The IgChain is then used for CDR identification and renumbering of ig_ab chain types.
    IgChains are part of the IgChainSet, which can be used to get specific types of IgChains of interest.
    """
    
    def __init__(self, chain_sequence, id, chain_type, chainID, ig_domains, description = None, full_id = None):

        #Description and full id needs to fixed for UserPyIgClassify so that it can be required upon creation.  This is currently not possible with the way
        # the fasta is parsed and stored for UserPyIgClassify.

        #Original Chain is also the id
        self.chainID = chainID
        self.sequence = chain_sequence
        self.ig_domains = ig_domains
        
        self.id = id
        self.description = description

        self.chain_type = chain_type
        self.ScFv = False
        self._identify_as_ScFv()

        self.full_id = full_id

        self.pdbaa_tag_info = None

        #Extra data that we really need to clean up
        self.benchmark_data = None
        self.hit_results = None

    def __str__(self):

        if self.is_data_from_pdbaa():
            lines = ""
            lines = lines+ "CHAIN "+self.pdbaa_tag_info.get_original_chain()+" "+self.get_type()+" "+repr(self.is_ScFv())+" "+repr(self.get_num_domains())+"\n"
            for domain in self.get_domains():
                hsp = domain.get_hsp()
                lines = lines+ "DOMAIN " + self.pdbaa_tag_info.get_original_chain()+" "+self.get_type()+" "+domain.get_gene()+" "+repr(hsp.bitscore)+" "+repr(hsp.evalue)+"\n"
        
            return lines


    def print_out(self, base_pdb_file = ""):
            lines = ""
            if base_pdb_file:
                lines = base_pdb_file+" "
            lines = lines+ "CHAIN "+self.get_id()+" "+self.get_type()+" "+repr(self.is_ScFv())+" "+repr(self.get_num_domains())+"\n"
            for domain in self.get_domains():
                hsp = domain.get_hsp()
                if base_pdb_file:
                    lines = lines+base_pdb_file+" "
                lines = lines+ "DOMAIN " + self.get_id()+" "+self.get_type()+" "+domain.get_gene()+" "+repr(hsp.bitscore)+" "+repr(hsp.evalue)+"\n"

            return lines


    def _identify_as_ScFv(self):
        """
        Set this AbChain as an ScFv.  Needed for renumbering code.
        """
        if len(self.ig_domains) == 1:
            self.ScFv = False
        elif len(self.ig_domains) > 2:
            self.ScFv = False
        else:
            heavy_found = False
            light_found = False
            light_genes = ['kappa', 'lambda', 'lambda6']
            for ig_domain in self.ig_domains:
                if ig_domain.get_gene() == 'heavy':
                    heavy_found = True
                elif ig_domain.get_gene() in light_genes:
                    light_found = True
            
            if light_found and heavy_found:
                self.ScFv = True
            else:
                self.ScFv = False
    
    def get_sequence(self):
        return self.sequence
    
    def is_ScFv(self):
        """
        Is this chain part of a full ScFv chain?
        """
        return self.ScFv
    
    def get_chain(self):
        """
        Get the chain that this IgChain should be.  If it's an ScFv, return None as it doesn't make sense.  It should be both. Use the IgDomain instead.
        """
        
        return self.chainID

    def get_id(self):
        return self.id

    def get_description(self):
        return self.description

    def get_type(self):
        return self.chain_type
    
    def get_domains(self):
        """
        Get a list of IgDomains of this chain.
        """
        return self.ig_domains
    
    def get_num_domains(self):
        """
        Get the number of IgDomains in this chain.
        """
        return len(self.ig_domains)

    
    def set_full_id(self, full_id):
        """
        Used mainly for PDBAA
        """
        self.full_id = full_id

    
    def get_full_id(self):
        return self.full_id

##### Benchmarking - Should also be its own class####

    def set_benchmark_data(self, benchmark_data):
        self.benchmark_data = benchmark_data

    def get_benchmark_data(self):
        return self.benchmark_data


#######Specific to PDBAA ######

    def set_pdbaa_tag_info(self, pdbaa_tag_info):
        assert isinstance(pdbaa_tag_info, PDBAATagInfo)
        self.pdbaa_tag_info = pdbaa_tag_info

    def get_pdbaa_tag_info(self):
        return self.pdbaa_tag_info

    def is_data_from_pdbaa(self):
        if self.pdbaa_tag_info == None:
            return False
        else:
            return True


class PDBAATagInfo(object):
    def __init__(self, pdbid, original_chain, protein_name):
        self.pdbid = pdbid
        self.orginal_chain = original_chain
        self.protein_name = protein_name

    def get_pdbid(self):
        return self.pdbid

    def get_original_chain(self):
        return self.orginal_chain

    def get_protein_name(self):
        return self.protein_name

class PDBAA(object):
    def __init__(self):
        os.path.abspath("/usr/bin")