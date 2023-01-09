import sqlite3
import sys
from collections import defaultdict

from src.modules.chains.HMMChainTypeParser import HMMChainTypeParser
from src.modules.chains.IgChain import IgChain
from src.modules.chains.IgDomain import IgDomain

from src.modules.chains.AbChain import AbChain
from src.tools.path import *


class InsertTestTables(object):
    """
    Weakly object-oriented class for getting benchmark data for pdbaa DatabasePyIgClassify.
    """
    def __init__(self, outdir, outname):
        self.outname = outname
        self.outdir = outdir
        self.outdb = self.outdir+"/"+self.outname+".db" #Same output DB, different
        print "Testing database will be: "+self.outdb
        
        if os.path.exists(self.outdb):
            os.remove(self.outdb)
        self.parser = HMMChainTypeParser()
        self.chain_set = self.parser.get_HMMChainTypeSet()
        self.chaintypes = sorted(self.chain_set.get_hmm_chains())
        
    def set_outdb(self, outdb):
        self.outdb = outdb
        
    def set_outname(self, outname):
        self.outname = outname
    
    def set_outdir(self, outdir):
        self.outdir = outdir
    
    
    def compare_to_cdr_db(self, comparison_db, ab_chains):
        """
        Output database that compares an old ab db to new data
        """
        
        if not self.compare_db: return
        for ab_chain in ab_chains:
            cdrs = ab_chain.get_cdrs()
    
    def create_hmm_result_outputs(self, ig_chains):
        """
        Creates a txt file and sqlite3 db of hmm results passing cutoff values.
        """
        self._create_main_identification_results(ig_chains, self.outname+"_identifications.txt")
        self._create_main_identification_results_db(ig_chains, "pdbaa_ident")
        
        txt_output_path = self._create_hmm_results_all_domains_txt_output(ig_chains, self.outname+"_all_domains_identified.txt")
        self._create_hmm_results_sqlite3_db(txt_output_path, "all_domains_no_cutoffs")
        
        txt_output_path = self._create_hmm_results_txt_output(ig_chains, self.outname+"_all_domains_identified_by_passing_cutoffs.txt")
        self._create_hmm_results_sqlite3_db(txt_output_path, "cutoff_domains")

    def create_pdb_cdr_results_from_abchains_pdbaa(self, ab_chains, table_name="cdr_data"):
        """
        Acts as the equivalent of Ben's output via a db.
        """
        
        db = sqlite3.connect(self.outdb)

        with db:
            cur = db.cursor()
            cur.execute("Drop table if exists "+table_name)
            cur.execute("CREATE table "+table_name+"(id INT, tag TEXT, original_chain TEXT, cdr TEXT, sequence TEXT, gene TEXT, cdr_start INT)")
            id = 0;
            for ab_chain in ab_chains:
                if isinstance(ab_chain, AbChain):
                    pass
                else:
                    sys.exit("Not an AbChain")

                cdrs = defaultdict()
                cdrs["H"] = ["H1", "H2", "H3"]
                cdrs["L"] = ["L1", "L2", "L3"]
                
                full_line = ab_chain.get_full_id()
                total_resnum = full_line.split()[1]
                exp_method = full_line.split()[2]
                resolution = full_line.split()[3]
                
                cdrs = ab_chain.get_cdrs()
                for cdr in cdrs:
                    id+=1
                    cur.execute("INSERT INTO "+table_name+ " VALUES(?,?,?,?,?,?,?)", \
                        (id, self.get_pdbid_or_tag(ab_chain), self.get_chain_id(ab_chain), cdr.get_type(), cdr.get_sequence(ab_chain.get_sequence()), cdr.get_gene(), cdr.get_start()))

    def create_pdb_cdr_results_from_abchains_fasta(self, ab_chains, table_name="cdr_data"):
        db = sqlite3.connect(self.outdb)

        with db:
            cur = db.cursor()
            cur.execute("Drop table if exists "+table_name)
            cur.execute("CREATE table "+table_name+"(id INT, entry_id TEXT, full_id TEXT, cdr TEXT, sequence TEXT, gene TEXT, cdr_start INT)")
            id = 0;
            for ab_chain in ab_chains:
                if isinstance(ab_chain, AbChain):
                    pass
                else:
                    sys.exit("Not an AbChain")

                cdrs = defaultdict()
                cdrs["H"] = ["H1", "H2", "H3"]
                cdrs["L"] = ["L1", "L2", "L3"]

                full_line = ab_chain.get_full_id()
                total_resnum = full_line.split()[1]
                exp_method = full_line.split()[2]
                resolution = full_line.split()[3]

                cdrs = ab_chain.get_cdrs()
                for cdr in cdrs:
                    id+=1
                    cur.execute("INSERT INTO "+table_name+ " VALUES(?,?,?,?,?,?,?)", \
                        (id, ab_chain.get_id(), ab_chain.get_full_id(), cdr.get_type(), cdr.get_sequence(ab_chain.get_sequence()), cdr.get_gene(), cdr.get_start()))

    def _create_main_identification_results(self, ig_chains, outname):
        outpath = self.outdir+"/"+outname
        OUTFILE = open_file(outpath, 'w')
        print "Creating main identification results: "+outpath
        header = "#tag description original_chain total_domains domain_num new_chain gene_id gene_score gene_evalue"
        
        OUTFILE.write(header+"\n")
        for ig_chain in ig_chains:
            assert isinstance(ig_chain, IgChain)
            i = 0
            
            for domain in ig_chain.get_domains():
                i+=1
                assert isinstance(domain, IgDomain)
                
                hsp = domain.get_hsp()

                base_line = self.get_full_id(ig_chain)
                base_line = base_line + " "+self.get_chain_id(ig_chain)
                base_line = base_line + " "+repr(ig_chain.get_num_domains())
                base_line = base_line+" "+repr(i)+ " "+domain.get_new_chain()+" "+domain.get_gene()+" "+repr(hsp.bitscore)+" "+repr(hsp.evalue)

                OUTFILE.write(base_line+"\n")
        OUTFILE.close()
        
        return outpath
    
    def _create_hmm_results_txt_file_header(self):
        """
        Returns the header to the main txt file result.
        """

        header = "#tag description original_chain total_domains domain_num gene gene_score gene_evalue "
        for chaintype in self.chaintypes:
            header = header+chaintype.get_gene()+"_evalue "
            header = header+chaintype.get_gene()+"_score "
            header = header+chaintype.get_gene()+"_passed "
            #header = header+chaintype.get_gene()+"_best "
        
        return header
    
    def _create_hmm_results_all_domains_txt_output(self, ig_chains, outname):
        """
        Creates a text file listing all database results - Not nessessarily all that pass cutoffs. 
        Returns path to text file for input into sqlite3 database.
        Identified Gene is best without cutoffs
        """
        
        outpath = self.outdir+"/"+outname
        
        print "Creating hmm results all domains: "+outpath
        OUTFILE = open_file(outpath, 'w')
        
        header = self._create_hmm_results_txt_file_header()
        OUTFILE.write(header+ "\n")
        
        
        
        for ig_chain in ig_chains:
            
            benchmark = ig_chain.get_benchmark_data()
            
            for i in range(1, len(benchmark)+1):
                best_hmm = self._get_best_hmm(benchmark, i)
                best_hsp = benchmark[i][best_hmm]['hsp']

                base_line = self.get_full_id(ig_chain)

                base_line = base_line + " "+self.get_chain_id(ig_chain)+" "+repr(len(benchmark))
                base_line = base_line + " "+repr(i)+ " "+self.chain_set.get_hmm_chain(best_hmm).get_gene()+" "+repr(best_hsp.bitscore)+" "+repr(best_hsp.evalue)
                for chaintype in self.chaintypes:
                    hmm = chaintype.get_hmm_name()
                    if benchmark[i].has_key(hmm):
                        hsp = benchmark[i][hmm]['hsp']
                        #best = benchmark[i][hmm]['best']
                        passed = benchmark[i][hmm]['passed_cutoffs']
                        
                        base_line = base_line+" "+str(hsp.evalue)+" "+str(hsp.bitscore)+" "+repr(int(passed))
                    else:
                        base_line = base_line+" 1000 0 0"
                OUTFILE.write(base_line+"\n")
                
        OUTFILE.close()
        return outpath
    
    def _get_best_hmm(self, benchmark, i):
        """
        Gets the best hmm of a domain from the benchmark results.
        """
        for chaintype in self.chaintypes:
            hmm = chaintype.get_hmm_name()
            if benchmark[i].has_key(hmm):
                if benchmark[i][hmm]['best']:
                    return hmm
                else:
                    continue
                
        return ""
    
    def _create_hmm_results_txt_output(self, ig_chains, outname):
        """
        Creates a text file listing all database results (passing cutoffs) for easy input into graphing software
        Returns path to text file for input into sqlite3 database.
        """

        
        outpath = self.outdir+"/"+outname
        print "Creating hmm results text: "+outname
        OUTFILE = open_file(outpath, 'w')
        
        header = self._create_hmm_results_txt_file_header()
        OUTFILE.write(header+ "\n")
        
        
        for ig_chain in ig_chains:
            domain_i = 0
            
            for domain in ig_chain.get_domains():


                base_line = self.get_full_id(ig_chain)
                base_line = base_line + " "+self.get_chain_id(ig_chain)+" "+repr(ig_chain.get_num_domains())
                domain_i+=1
                hsp = domain.get_hsp()
                base_line = base_line+" "+repr(domain_i)+" "+domain.get_gene()+" "+repr(hsp.bitscore)+" "+repr(hsp.evalue)
                for chaintype in self.chaintypes:
                    hmm = chaintype.get_hmm_name()                
                    hsp  = domain.get_passing_hsp(hmm)
                    if hsp:
                        base_line = base_line+" "+str(hsp.evalue)+" "+str(hsp.bitscore)+" "+"1"
                    else:
                        base_line = base_line+" 1000 0 1"
                
                OUTFILE.write(base_line+"\n")
        
        
        
        OUTFILE.close()
        return outpath
    
    def _create_main_identification_results_db(self, ig_chains, table_name):
        db = sqlite3.connect(self.outdb)
              
        with db:
            cur = db.cursor()
            cur.execute("DROP table IF EXISTS "+table_name)
            cur.execute("CREATE TABLE "+table_name+"(id INT, tag TEXT, description TEXT, original_chain TEXT, total_domains INT, domain_num INT, \
            new_chain TEXT, gene TEXT, gene_score FLOAT, gene_evalue FLOAT)")
            
            row_num = 1
            for ig_chain in ig_chains:
                assert isinstance(ig_chain, IgChain)
                i = 0
                for domain in ig_chain.get_domains():
                    i+=1
                    assert isinstance(domain, IgDomain)
                    
                    hsp = domain.get_hsp()
                    cur.execute("INSERT INTO "+table_name+" VALUES(?,?,?,?,?,?,?,?,?,?)", \
                    (row_num, self.get_pdbid_or_tag(ig_chain), self.get_name_or_description(ig_chain), self.get_chain_id(ig_chain), ig_chain.get_num_domains(), \
                    i, domain.get_new_chain(), domain.get_gene(), hsp.bitscore, hsp.evalue))          
                    row_num+=1
    
    def _create_hmm_results_sqlite3_db(self, input_txt_file, table_name):
        """
        Uses the text file generated in create_hmm_results_txt_output to create an sqlite3 db with that information.
        """
        
        db = sqlite3.connect(self.outdb)
        header = "id INT, tag TEXT, description TEXT, original_chain TEXT, total_domains INT, domain_num INT, gene TEXT, gene_score FLOAT, gene_evalue FLOAT"
        
        for chaintype in self.chaintypes:
            header = header+","+chaintype.get_gene()+"_evalue FLOAT"
            header = header+","+chaintype.get_gene()+"_score FLOAT"
            header = header+","+chaintype.get_gene()+"_passed INT"
            #header = header+","+chaintype.get_gene()+"_best INT"
         
        with db:
            cur = db.cursor()
            cur.execute("DROP table IF EXISTS "+table_name)
            cur.execute("CREATE table "+table_name+"("+header+")")
       
            
            FILE = open_file(input_txt_file, 'r')
            row_num= 1
            for line in FILE:
                line = line.strip()
                if line[0]=='#':continue
                if not line: continue
                
                lineSP = line.split()
                
                result = [row_num]
                for i in range(0, len(lineSP)):
                    result.append(lineSP[i])
                    
                cur.execute("INSERT INTO "+table_name+" VALUES("+self._get_question_mark_string(len(lineSP)+1)+")", \
                            (result))
            
            FILE.close()
        
    def _get_question_mark_string(self, length):
        
        if not length:return
        
        s = '?'
        for i in range(1, length):
            s = s+',?'
        return s
    
    def _insert_dot_to_protein_name(self, protein_name):
        """
        Makes the protein name a connected string.
        """
        nameSP = protein_name.split()
        together = '.'
        together = together.join(nameSP[6:])
        return together

    def get_pdbid_or_tag(self, chain):
        """
        Return the PDBID of the AbChain or IgChain if PDBAA or the tag if a different fasta.
        """
        assert isinstance(chain, IgChain)
        if chain.is_data_from_pdbaa():
            return chain.pdbaa_tag_info.get_pdbid()
        else:
            return chain.get_id()

    def get_full_id(self, chain):
        assert isinstance(chain, IgChain)
        if chain.is_data_from_pdbaa():
            return chain.pdbaa_tag_info.get_pdbid()+" "+self._insert_dot_to_protein_name(chain.pdbaa_tag_info.get_protein_name())
        else:
            return chain.get_id()+" "+".".join(chain.get_description().split())

    def get_chain_id(self, chain):
        assert isinstance(chain, IgChain)
        if chain.is_data_from_pdbaa():
            return chain.pdbaa_tag_info.get_original_chain()
        else:
            return "~"

    def get_name_or_description(self, chain):
        assert isinstance(chain, IgChain)
        if chain.is_data_from_pdbaa():
            return self._insert_dot_to_protein_name(chain.pdbaa_tag_info.get_protein_name())
        else:
            return ".".join(chain.get_description().split())