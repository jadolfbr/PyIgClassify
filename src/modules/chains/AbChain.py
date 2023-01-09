#Author: Jared Adolf-Bryfogle

import sys

from src.modules.CDRs.CDRInfo import CDRInfo
from src.modules.Structure import AntibodyResidue
from src.modules.Structure import PDBResInfo

from .IgChain import IgChain
from src.modules.CDRs.CDRClusterer import CDRClusterer


class AbChain(IgChain):
    """
    A Specific type of IgChain.  Requires a NumberingMap dictionary, which is given to the AbChain during creation through the AbChainFactory.
    This class is then used to renumber a PDB, print numbering + CDR info, or determine a CDR Cluster from the CDR sequence.
    """
    def __init__(self, pdb_res_info, cdrs_info_array, ig_chain):

        assert isinstance(ig_chain, IgChain)
        IgChain.__init__(self, ig_chain.get_sequence(), ig_chain.get_id(), "ig_ab", ig_chain.get_chain(), ig_chain.get_domains(), ig_chain.get_description(), ig_chain.get_full_id())
        
        self.cdrs = cdrs_info_array
        self.pdb_res_info = pdb_res_info
        if not isinstance(self.pdb_res_info, PDBResInfo): sys.exit()

    def get_fasta_print(self, base_file_name = ""):
        
        #Used ONLY for FASTA output
        full_string = "#original_chain new_chain resname old_resnum new_resnum region_type cdr_type\n"
        #print self.numbering_map

        for i in range(1, self.pdb_res_info.total_residue()+1):
            
            """
            print self.get_original_chain()
            print self.get_chain()
            print self.get_resname(i)
            print str(self.get_new_resnum(i))
            print self.get_region_type(i)
            print self.get_cdr_type(i)
            """
            res = self.pdb_res_info.get_residue(i)
            if not isinstance(res, AntibodyResidue): sys.exit()

            line = ""
            if base_file_name:
                line = base_file_name+" "

            line = line+self.get_id() + "\t"+res.get_chain()+"\t"+res.get_aa()+"\t"+str(i).ljust(5)+" "+str(res.get_pdb_num()).ljust(5)+"  "+res.get_region()+"\t"+res.get_cdr_type()+"\n"
            full_string = full_string+line
        
        return full_string
    
    def get_pdb_print(self, base_file_name = ""):
        """"
        Print function for PDB UserPyIgClassify
        """
        full_string = "#original_chain chain_type new_chain resname old_resnum old_icode new_resnum new_icode region_type cdr_type cluster\n"


        for i in range(1, self.pdb_res_info.total_residue()+1):
            
            """
            print self.get_original_chain()
            print self.get_chain()
            print self.get_resname(i)
            print str(self.get_new_resnum(i))
            print self.get_region_type(i)
            print self.get_cdr_type(i)
            """
            res = self.pdb_res_info.get_residue(i)
            if not isinstance(res, AntibodyResidue): sys.exit()
            line = ""
            if base_file_name:
                line = base_file_name+" "

            line = line+self.get_id()+"  "+res.get_chain_type()+"  "+res.get_chain()+" \t"+res.get_aa()+"  "+str(res.get_old_resnum()).rjust(5)+ " "+str(res.get_old_icode())+"  "+\
            str(res.get_pdb_num()).rjust(5)+" "+str(res.get_icode())+"\t"+res.get_region()+"\t"+res.get_cdr_type()+"\t"+res.get_cluster()+"\n"
            
            full_string = full_string+line
        
        return full_string

    def num_residues(self):
        """
        Return the number of residues in the map.
        """
        return len(self.sequence)

    ###########################################################
    # Numbering
    #
    #
    def get_pdb_res_info(self):
        """
        Get the newer PDBResInfo, made of AntibodyResidues.
        Much easier to use.
        """
        return self.pdb_res_info

    def get_new_chain(self, resnum):
        return self.pdb_res_info.pose_to_record_map[resnum].get_chain()

    ###########################################################
    # CDRs
    #
    #
    def has_cdr_type(self, cdr_type):
        found = False
        for i in range(1, self.pdb_res_info.total_residue()+1):

            if self.pdb_res_info.pose_to_record_map[ i ].get_cdr_type() == cdr_type:
                found = True
                break
        
        return found

    def get_cdrs(self):
        return self.cdrs
    
    def get_cdrs_of_type(self, cdr_type):
        cdrs = []
        for cdr in self.cdrs:
            if cdr.get_type()==cdr_type:
                cdrs.append(cdr)
        return cdrs
    
    def identify_clusters(self, pose):
        """
        ReCreates CDRInfo objects, with Cluster Info.  Poplulates self.pdb_res_info with cluster info
        """
        
        CDRs = []
        i = 1
        clusterer = CDRClusterer(pose)
        while i < self.pdb_res_info.total_residues()+1:

            res = self.pdb_res_info.get_residue(i)
            if not isinstance(res, AntibodyResidue): sys.exit()

            if res.get_region() == "CDR":
                cdr_type = res.get_cdr_type()
                gene = res.get_gene()

                new_chain = res.get_chain()
                old_chain = res.get_old_chain()
                new_start = res.get_pdb_num()
                cdr_start = i
                cdr_end = i
                x = i
                #Get length of CDR
                for x in range(cdr_start, self.pdb_res_info.total_residues()+1):
                    if not self.pdb_res_info.pose_to_record_map[x].is_cdr():
                        break
                    else:
                        cdr_end+=1
                        i+=1
                i+=1
                
                #print "CDR: "+repr(cdr_start)+" "+repr(cdr_end)
                cdr_info = CDRInfo(cdr_type, cdr_start, cdr_end -1, gene)
                cdr_info.set_new_chain(new_chain)
                clusterer.set_dihedrals_from_cdr(cdr_type, new_chain)
                cluster = clusterer.get_fullcluster(cdr_type, new_chain)

                for x in range(cdr_start, cdr_end):
                    self.pdb_res_info.pose_to_record_map[x].set_extra_info('cluster', cluster[0])
                    self.pdb_res_info.pose_to_record_map[x].set_extra_info('cluster_dis', cluster[1])
                    
                cdr_info.set_cluster(cluster[0])
                cdr_info.set_distance(cluster[1])
                cdr_info.new_start = new_start
                cdr_info.set_original_cdr_info(self.pdb_res_info.pose_to_record_map[cdr_start].get_old_resnum(), self.pdb_res_info.pose_to_record_map[cdr_end].get_old_resnum(), self.get_id())
                CDRs.append(cdr_info)
            else:
                i+=1
                
        self.cdrs = CDRs


    
