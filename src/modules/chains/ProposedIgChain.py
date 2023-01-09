#Author: Jared Adolf-Bryfogle
import random

#Python Imports
from collections import defaultdict
#Project Imports
from src.modules.family.HMMFamilyType import HMMFamilyType
from src.modules.family.HMMFamilyTypeParser import HMMFamilyTypeParser
from src.tools.path import *
from src.tools import fasta
from src.tools import general

#BioPython Imports
from Bio import SearchIO

from src.modules.hmmer import HMMERRunner

class ProposedIgChain(object):
    """
    A collection of HMMER3 result data for a specific fasta file.  Could be multiple chains.
    Runs HMMER3, uses BioPython and SimpleHMMER to do the dirty work.
    There should only be ONE FASTA entry in file.
    """
    def __init__(self, fasta_path):
        self.fasta = fasta_path
        self.hmms = [];
        self.hmm_results = defaultdict(dict);
        self.families = []
        self.family_results = defaultdict(dict)
        self._initialize_data()
        
        self.family_parser = HMMFamilyTypeParser()
        self.hmm_family_set = self.family_parser.get_HMMFamilyTypeSet()
        self._determine_family_types()

    def __str__(self):
        
        #For Debugging.
        print self.fasta
        for hmm in self.hmm_results:
            print hmm
            for hit in self.get_query_result(hmm):
                for hsp in hit:
                    print hsp
                    print "\n"
                    
        print "\nFamily: \n" 
        for hmm in self.family_results:
            print hmm
            for hit in self.family_results[hmm]['query_result']:
                for hsp in hit:
                    print hsp
                    print "\n"
        

        return ""
    
    def _initialize_data(self):
        chain_type_path = get_db_path()+"/hmm_chain_types.txt"
        FILE = open_file(chain_type_path, 'r')
        for line in FILE:
            if line[0]=="#":continue
            line = line.strip()
            lineSP = line.split()
            if len(lineSP)!=6:continue
            self.hmms.append(lineSP[0])
        FILE.close()
        
        
        for hmm in self.hmms:
            
            #Run the HMM.
            #Add BioPython Result and output text file name
            result_name = 'hmm_result_'+os.path.basename(self.fasta)+"_"+str(random.getrandbits(128))

            HR = HMMERRunner(prefix=result_name)

            
            HR.search(get_db_path()+"/"+hmm, self.fasta, os.getcwd())
            result_name = result_name+"_out"
            result_path = os.getcwd()+"/"+result_name
            RESULT = open_file(result_path+".hmmer3", 'r')
            
            #Roundabout way until biopython is fixed.
            try:
                query_result = SearchIO.read(RESULT, 'hmmer3-text')
            except AttributeError:
                query_result = None
            self.hmm_results[hmm]['query_result'] = query_result
            RESULT.close()
        
            os.remove(result_path+".hmmer3")
            os.remove(result_path+".txt")
    
    def get_hmms(self):
        """
        Return list of hmms used to propose the chain.
        """
        
        return self.hmms
    
    def get_query_result(self, hmm_name):
        """
        Return the BioPython Qeury Result
        """
        return self.hmm_results[hmm_name]['query_result']
    
    def get_num_hits(self, hmm_name):
        """
        Get the num of hits, aka, number of chains in FASTA
        """
        if not self.hmm_results[hmm_name]['query_result']:
            return 0
        
        return len(self.hmm_results[hmm_name]['query_result'])
    
    def get_hits(self, hmm_name):
        """
        Return a list of BioPython Hits for specific hmm_name
        """
        if not self.hmm_results[hmm_name]['query_result']:
            return None
        
        return self.hmm_results[hmm_name]['query_result'].hits
        
    def get_id(self):
        """
        Return the ID of query.  Since this is chain based, it should, but not nessessarily be the original chain ID used in creation of fasta file.
        This is a problem if no hit was found, with the BioPython bug!
        """
        """
        if not self.hmm_results[self.hmms[0]]['query_result']:
            return ""
        
        return self.hmm_results[self.hmms[0]]['query_result'].id
        """
        
        label = fasta.get_label_from_fasta(self.fasta)
        return label


#########Families##################################

    def _determine_family_types(self):
        
        
        hmm_family_set = self.hmm_family_set
        families = hmm_family_set.get_hmm_chains()
        
        per_domain_results = []
        
        for family_type in families:
            assert isinstance(family_type, HMMFamilyType)
            
            hmm = family_type.get_hmm_path()
            result_name = 'hmm_result_family_'+os.path.basename(self.fasta)+"_"+str(random.getrandbits(128))

            HR = HMMERRunner(prefix=result_name)

            HR.search(get_db_path()+"/"+hmm, self.fasta, os.getcwd())
            
            result_name = result_name+"_out"
            result_path = os.getcwd()+ "/"+result_name
            RESULT = open_file(result_path+".hmmer3", 'r')
            
            
            try:
                query_result = SearchIO.read(result_name+".hmmer3", 'hmmer3-text')
                self.family_results[hmm]['query_result']=query_result
            except AttributeError:
                continue
            
            per_domain_results = []
            if len(query_result)==1:
                for hsp in query_result[0]:
                    #Check Evalue - If not passed, continue
                    if hsp.bitscore < family_type.get_score_cutoff(): continue
                    elif hsp.evalue > family_type.get_evalue_cutoff(): continue
                    results = defaultdict(dict)
            
                    results[hmm]['hsp'] = hsp
                    results[hmm]['evalue'] = hsp.evalue
                    results[hmm]['score'] = hsp.bitscore
                    per_domain_results.append(results)
                    
            RESULT.close()
        
            os.remove(result_path+".hmmer3")
            os.remove(result_path+".txt")
        
        if len(per_domain_results)==0:
            return
        else:
            best_families = []
            for results in per_domain_results:
                
                best_scoring_hmm = general.get_best_score_hmm(results)
            
                best_evalue_hmm = general.get_best_evalue_hmm(results)
                best_families.append(hmm_family_set.get_hmm_family(best_scoring_hmm))
                
            self.families = best_families
        
    def get_families(self):
        return self.families
