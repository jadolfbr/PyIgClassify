#Author Jared Adolf-Bryfogle

#Project Imports
from .HMMChainTypeParser import HMMChainTypeParser
from .IgChain import IgChain
from .IgChain import PDBAATagInfo
from .IgChainSet import IgChainSet
from .ProposedIgChain import ProposedIgChain
from .ProposedIgChains import ProposedIgChains
from src.tools import fasta

import src.tools.general as gen_tools
from src.modules.chains.IgDomain import IgDomain

#Python Imports
from collections import defaultdict


class IgChainFactory(object):
    """
    Create a list of IgChains from a ProposedIgChain based on cutoff values/score/etc stored in HMMChainTypes.
    """
    
    def __init__(self):
        self.chain_parser = HMMChainTypeParser()
        self.HMM_chain_set = self.chain_parser.get_HMMChainTypeSet()
        self.domain_overlap_cutoff = 35
    
    def _does_hmm_pass_cutoff(self, hsp, hmm):
        hmms = []
        
        chaintype = self.HMM_chain_set.get_hmm_chain(hmm)
        if (hsp.bitscore >= chaintype.get_score_cutoff()) and (hsp.evalue <= chaintype.get_evalue_cutoff()):
            return True
        else:
            return False
    
    def _are_seperate_domains(self, hsp_pair_list):
        """
        Rudemntarily, we want to know if they are indeed seperate domains.
        """
        
        if len(hsp_pair_list) ==1:
            return True
        
        for pair1 in hsp_pair_list:
            for pair2 in hsp_pair_list:
                if pair1 == pair2:
                    continue
                
                if self._are_seperate_domain(pair1[1], pair2[1]):
                    continue
                else:
                    return False
        return True
    
    def _are_seperate_domain(self, hsp1, hsp2):
        overlap = self._get_overlap(hsp1, hsp2)
        
        if overlap <= self.domain_overlap_cutoff:
            return True
        else:
            return False
        
    def _are_same_domain(self, hsp1, hsp2):
        overlap = self._get_overlap(hsp1, hsp2)
        
        if overlap >= self.domain_overlap_cutoff:
            return True
        else:
            return False
    
    def _are_same_domains(self, hsp1, hsp_pair_list):
        
        for pair in hsp_pair_list:
            if self._are_same_domain(hsp1, pair[1]):
                continue
            else:
                return False
        
        return True
    
    def _are_same_domains_from_avg_overlap(self, hsp1, hsp_pair_list):
        pass
    
    def _get_avg_overlap(self, hsp1, hsp_pair_list):
        pass
    
    def _get_overlap(self, hsp1, hsp2):
        """
        Get the number of overlapping alignment residues from these two hsps.
        """
        start1 = hsp1.hit_start
        start2 = hsp2.hit_start
        end1 = hsp1.hit_end
        end2 = hsp2.hit_end
        
        #There is probably a one line way to do this, but...
        overlap = 0
        for i in range(start1, end1+1):
            if i in range(start2, end2+1):
                overlap+=1
        return overlap
    
    def _sort_hmm_pair_list(self, hmm_pair_list):
        """
        Sorts the list into a list of list.
        Similar domains are in one list.  Length of this list of list is the number of domains found.
        """
        sorted_lists = []
        sorted_lists.append([hmm_pair_list[0]])
        
        for pair in hmm_pair_list:
            pair_identified= False
            for tentative_domain_pair_list in sorted_lists:
                if self._are_same_domains(pair[1], tentative_domain_pair_list):
                    tentative_domain_pair_list.append(pair); #Not sure if this will work - is it a copy or a reference?
                    pair_identified= True
                    break
                else:
                    continue
            if not pair_identified:
                    new_domain = []
                    new_domain.append(pair)
                    sorted_lists.append(new_domain)
        
        #print "Total domains identified with a "+repr(self.domain_overlap_cutoff)+"res domain overlap cutoff: "+repr(len(sorted_lists))
        return sorted_lists
                    
        
    
    def _sort_hmm_pair_list_avg_mode(self, hmm_pair_list):
        pass
    
    def _create_ig_domain(self, hmm, hsp, passing_domain_hmms):
        
        chain_type = self.HMM_chain_set.get_hmm_chain(hmm)
        #Get the subfamily type for what will be the ig_chain
        subfamily_type = chain_type.get_type()
        
        #Create the IgChain.
        ig_domain = IgDomain(hsp, chain_type.get_gene(), chain_type.get_chain(), passing_domain_hmms)
        return ig_domain
    
    def _create_benchmark_data(self, all_hmms, passing_hmms, proposed_ig_chain):
        """
        Used for benchmarking pdbaa in order to find the best cutoffs + wrong data
        all_hmms is [hmm] = [hsp, hsp, hsp]
        passing_hmms is [hmm] = [hsp, hsp, hsp]
        
        Returns a dictionary:
          [domain_num][hmm]['hsp'] = hsp
          [domain_num][hmm]['best']= boolean
          [domain_num][hmm]['passed_cutoffs'] = boolean
        """
        
        benchmark_data = defaultdict(lambda: defaultdict(dict))
        hmm_pairs = []
        for hmm in all_hmms:
            for hsp in all_hmms[hmm]:
                hmm_pairs.append([hmm, hsp])
        sorted_hmms = self._sort_hmm_pair_list(hmm_pairs)
        
        domain_num = 0
        for hmm_pair_list in sorted_hmms:
            domain_num+=1
            best_pair = gen_tools.get_best_score_pair(hmm_pair_list)
            
            
            for hmm_pair in hmm_pair_list:
                hmm = hmm_pair[0]
                hsp = hmm_pair[1]
                benchmark_data[domain_num][hmm]['hsp'] = hsp
                
                #Does the pair pass cutoffs?
                if hsp in passing_hmms[hmm]:
                    benchmark_data[domain_num][hmm]['passed_cutoffs'] = True
                else:
                    benchmark_data[domain_num][hmm]['passed_cutoffs'] = False
                
                #Is it the best pair within this pseudo domain?
                if hmm_pair == best_pair:
                    benchmark_data[domain_num][hmm]['best'] = True
                else:
                    benchmark_data[domain_num][hmm]['best'] = False
        
        return benchmark_data
    
    def _create_ig_chain(self, passing_hmms, proposed_ig_chain, sequence, new_id = None):
        """
        Gets the best scoring HMM of dictionary of a list of [hmm] = [hsp, hsp, hsp]
        Returns an IgChain
        """
        
        #Return results if we have them already:
        ig_domains = []
        chain = ""; #Chain from the chain_type.  If it's an ScFv, we will deal with that in the renumbering
        original_chain = ""; #This comes from the FASTA query.  This is why we want per-chain in the FASTA
        subfamily_type = ""; #This is a subfamily for the chain. Not the domain.  Domain has seperate genes.
        
        if len(passing_hmms) == 0:
            return None
      
        else:
            hmm_pairs = []
            for hmm in passing_hmms:
                for hsp in passing_hmms[hmm]:
                    hmm_pairs.append([hmm, hsp])
            sorted_hmms = self._sort_hmm_pair_list(hmm_pairs)
            for domain_hmms in sorted_hmms:
                best_pair = gen_tools.get_best_score_pair(domain_hmms)
                
                #For benchmarking passed hmms:
                passing_domain_hmms = defaultdict()
                for pair in domain_hmms:
                    passing_domain_hmms[pair[0]] = pair[1]
                
                ig_domain = self._create_ig_domain(best_pair[0], best_pair[1], passing_domain_hmms)
                
                ig_domains.append(ig_domain)
                if not new_id:
                    id = proposed_ig_chain.get_query_result(best_pair[0])[0].id
                else:
                    id = new_id
                
                #What if multiple subfamilies based on domain split due to wierd TCR/Etc.
                subfamily_type = self.HMM_chain_set.get_hmm_chain(best_pair[0]).get_type()
                chain = self.HMM_chain_set.get_hmm_chain(best_pair[0]).get_chain()
        """
        elif len(passing_hmms)==1:
            hmm = passing_hmms.keys()[0]
            if self._are_seperate_domains(passing_hmms[hmm]):
                #Very unlikely, however..
                for hsp in passing_hmms[hmm]:
                    ig_domain = self._create_ig_domain(hmm, hsp)
                    ig_domains.append(ig_domain)
                    original_chain = hsp.id
                    subfamily_type = self.HMM_chain_set.get_hmm_chain(hmm).get_type()
            else:
                sorted_hmms = self._sort_hmm_pair_list(passing_hmms[hmm])
                for domain_hmms in sorted_hmms:
                    best_pair = gen_tools.get_best_score_pair(domain_hmms)
                    ig_domain = self._create_ig_domain(best_pair[0], best_pair[1])
                    ig_domains.append(ig_domain)
                    original_chain = hsp.id
                    subfamily_type = self.HMM_chain_set.get_hmm_chain(hmm).get_type()
        """  
        
        #If we make it this far, we create an IgChain and return it.
        
        ig_chain = IgChain(sequence, id, subfamily_type, chain, ig_domains)
        return ig_chain

    def create_ig_vset_chain_set_from_multi_fasta(self, proposed_ig_chain, pdbaa = True, add_benchmark_data = True):
        """
        Create an IgChainSet from the PDBAA-specific proposed_ig_chain
        """
        
        assert isinstance(proposed_ig_chain, ProposedIgChains)
        
        ig_chains = []
        
        ids = proposed_ig_chain.get_ids()
        print "# of FASTA Entries: "+repr(len(ids))
        for id in ids:
            family_ids = proposed_ig_chain.get_families(id)
            if not family_ids: continue
            if "ig_vset" not in family_ids: continue
            
            
            #print "Checking "+id
            #Check all hmm chain types for passing of cutoffs.
            all_hmms = defaultdict(list)
            passing_hmms = defaultdict(list)
            for hmm in proposed_ig_chain.get_hmms():
                query_result = proposed_ig_chain.get_query_result(hmm)
                if not query_result: continue
                if len(query_result) == 0: continue               
                for hsp in proposed_ig_chain.get_hsps(id, hmm):
                    if self._does_hmm_pass_cutoff(hsp, hmm):
                        passing_hmms[hmm].append(hsp)
                    if add_benchmark_data:
                        all_hmms[hmm].append(hsp)
                        


            ig_chain = self._create_ig_chain(passing_hmms, proposed_ig_chain, proposed_ig_chain.get_sequence(id), id)

            #If nothing passes: Don't keep it!
            if not ig_chain:
                continue

            #For benchmarking:
            if add_benchmark_data:
                benchmark_data = self._create_benchmark_data(all_hmms, passing_hmms, proposed_ig_chain)
                ig_chain.set_benchmark_data(benchmark_data)
                
            #print "Creating: "+proposed_ig_chain.get_full_id(id)

            #These are hacky - but cannot currently be set via UserPyIgClassify.  Fix if possible.

            ig_chain.hit_results = proposed_ig_chain.get_hits_for_hmms(id) #Just a hack to get all data.
            ig_chain.full_id = id+" "+" ".join(proposed_ig_chain.get_description(id).split()[1:])
            ig_chain.description = " ".join(proposed_ig_chain.get_description(id).split()[1:])

            if pdbaa:
                pdbid_chain = id
                pdbaa_tag = PDBAATagInfo(pdbid_chain[0:4], pdbid_chain[4], proposed_ig_chain.get_protein_name(id))
                ig_chain.set_pdbaa_tag_info(pdbaa_tag)
            ig_chains.append(ig_chain)


        ig_set = IgChainSet(ig_chains)
        return ig_set

    def create_ig_vset_chain(self, proposed_ig_chain):
        """
        Create a regular IgChain or return None if chain is NOT an igchain or a tcr igchain.
        """
        
        #If it passes all tests, set the family and the original chain id from the ID.  For sequence only, this may take some refactoring already!
        assert isinstance(proposed_ig_chain, ProposedIgChain)
        
        

        #Determine the Family.  If it does not have an ig_vset domain, do not create and return an IgChain.
        #If it didn't match anything, even routinely, we have no results to create anything.

        #Determine family.  Make sure there is an ig_vset in the chain.
        family_ids = proposed_ig_chain.get_families()
        #print family_ids
        if "ig_vset" not in family_ids: return None
    
        
        #Check all hmm chain types for passing of cutoffs. 
        passing_hmms = defaultdict(list)
        for hmm in proposed_ig_chain.get_hmms():
            #print hmm
            query_result = proposed_ig_chain.get_query_result(hmm)
            if not query_result: continue
            if len(query_result) == 0: continue
            
            for hsp in query_result[0]:
                if self._does_hmm_pass_cutoff(hsp, hmm):
                    passing_hmms[hmm].append(hsp)
                    
        if len(passing_hmms) == 0:
            return None
        
        original_chain = proposed_ig_chain.get_id()
        sequence = fasta.get_sequence_from_fasta(proposed_ig_chain.fasta, original_chain)
        #print "Orginal_chain: "+original_chain+" "+sequence
        ig_chain = self._create_ig_chain(passing_hmms, proposed_ig_chain, sequence)
        return ig_chain
