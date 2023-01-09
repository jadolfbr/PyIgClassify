import math
import sys
from collections import defaultdict

from src.modules.CDRs.CDRInfo import CDRInfo
from src.modules.RenumberMapTypeParser import RenumberMapTypeParser
from src.modules.Structure import AntibodyResidue
from src.modules.Structure import PDBResInfo

import src.tools.renumbering as renum_tools
from src.modules.chains.AbChain import AbChain
from src.modules.chains.IgChain import IgChain
from src.modules.CDRs.CDRDefinitionParser import CDRDefinitionParser


class AbChainFactory(object):
    """
    Creates an AbChain from an IgChain.  The AbChain has renumbering information and is the final Class that is used for renumbering.
    Mutates the raw seq_to_renumbering dictionary to the final seq_to_renumbering that is housed and accessed by AbChain.
    
    Caveats:
      1) numbering_map attached to AbChain is the full chain.  For now, it will renumber the Vset and if a Cset is attached to the domain, will renumber to n
      2) An ScFv is dealt with properly.  However, a double Light or double Heavy AbChain will give very wrong numbering.
      
    """
    
    def __init__(self):
        self.numbering_parser = RenumberMapTypeParser()
        self.cdr_types = ["L1", "L2", "L3", "H1", "H2", "H3"]
        self.icodes = icodes = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz"

        self.pdb_info = PDBResInfo()
            
    def _debug_print(self):
        """
        Prints the numbering map in an easy to see way for debugging.
        """
        for i in range(1, self.pdb_info.total_residues()+1):
            print repr(self.pdb_info.residue( i ))
            
    def _create_cdr_numbering_in_renumbering_map(self, numbering_scheme, cdr_type):
        
        if not self._has_cdr_type(cdr_type):
            return
        
        numbering_parser = CDRDefinitionParser(numbering_scheme)
        max_length = numbering_parser.get_cdr_end(cdr_type) - numbering_parser.get_cdr_start(cdr_type) + 1

        if numbering_scheme == 'modified_aho' or numbering_scheme == 'aho':
            for i in range(1, self.pdb_info.total_residues()+1):
                res_i = self.pdb_info.residue( i )
                if not isinstance(res_i, AntibodyResidue): sys.exit()

                if res_i.pdb_num == cdr_type:
                    cdr_start = i
                    cdr_length = 0
                    for cdr_position in range(i, self.pdb_info.total_residues()+1):
                        if self.pdb_info.residue( cdr_position ).pdb_num != cdr_type: break
                        cdr_length+=1
                    
                    cdr_end = cdr_start+cdr_length -1
                    #print cdr_type+ " "+repr(cdr_length)+ " "+repr(cdr_start) + " "+ repr(cdr_end)
                    
                    max_start_residues  = int(math.ceil(float(max_length)/2))
                    max_end_residues = int(math.floor(float(max_length)/2))

                    start_residues = int(math.ceil(float(cdr_length)/2))
                    end_residues = int(math.floor(float(cdr_length)/2))

                    floor_start = cdr_start+start_residues
                    floor_new_start = numbering_parser.get_cdr_end(cdr_type) - end_residues + 1

                    #print "start: "+repr(cdr_start)+" end: "+repr(cdr_end)+" start_res: "+repr(start_residues)+" end_res:"+repr(end_residues)

                    
                    #Renumber the starting residues:
                    offset = 0
                    last_resnum = 0
                    for res in range(cdr_start, cdr_start+start_residues):
                        self.pdb_info.residue( res ).set_pdb_num( numbering_parser.get_cdr_start( cdr_type ) + offset )
                        offset+=1

                    #Use insertion codes for longer CDRs that are normally allowed by the numbering scheme.
                    if cdr_length > max_length:

                        #print "CDR Length > Max Length"
                        #print repr(cdr_length)+" > "+repr(max_length)

                        floor_new_start = numbering_parser.get_cdr_end(cdr_type) - max_end_residues + 1

                        floor_start = cdr_end - max_end_residues + 1

                        #print "start: "+repr(cdr_start)+" max_start_res: "+repr(max_start_residues)

                        i_code_start = cdr_start + max_start_residues
                        i_code_end = cdr_end - max_end_residues
                        i_code_resnum = self.pdb_info.residue( i_code_start-1 ).get_pdb_num()
                        i = 0

                        for res in range(i_code_start, i_code_end+1):
                            #print repr(res)

                            if i > len(self.icodes): break

                            self.pdb_info.residue(res).set_pdb_num(i_code_resnum )
                            self.pdb_info.residue(res).set_icode( self.icodes[i] )
                            i+=1

                    #Renumber the ending residues:
                    offset = 0
                    for res in range(floor_start, cdr_end +1):
                        #print repr(res)+" "+repr(floor_new_start+offset)
                        self.pdb_info.residue(res).set_pdb_num( floor_new_start+offset )
                        offset+=1
                else:
                    continue
        
                
        
        else:
            sys.exit("Numbering scheme not yet supported.")
        
        #self._debug_print()
        
    def _identify_cdrs_from_anchors(self, numbering_scheme, cdr_type):
        
        #If this is slow, fix it up.
        
        cdr_n_anchor = cdr_type + "_N"
        cdr_c_anchor = cdr_type + "_C"
        #print "Identifying CDRs from anchors"
        for i in range(1, self.pdb_info.total_residues()+1):
            if self.pdb_info.pose_to_record_map[ i ].get_meta() == cdr_n_anchor:
                for anchor_position in range(i , self.pdb_info.total_residues()+1):
                    #print repr(anchor_position)
                    
                    #Reached the C-term anchor
                    if self.pdb_info.pose_to_record_map[ anchor_position ].get_meta() == cdr_c_anchor: break
                    #Reached regular residues.  No C-term found.  This is a problem.
                    elif self.pdb_info.pose_to_record_map[ anchor_position ].get_meta() == "":
                        pass
                    #Another N_anchor.  Good.  Should be fine.
                    elif self.pdb_info.pose_to_record_map[ anchor_position ].get_meta() == cdr_n_anchor: continue
                    
                    
                    if self.pdb_info.pose_to_record_map[ anchor_position ].pdb_num == "-2":
                        self.pdb_info.pose_to_record_map[ anchor_position ].pdb_num = cdr_type
                    elif self.pdb_info.pose_to_record_map[ anchor_position ].pdb_num == "+1":
                        self.pdb_info.pose_to_record_map[ anchor_position ].pdb_num = cdr_type
                    else:
                        continue
                    

            
    def _check_for_insertions(self, numbering_scheme):
        for i in range(1, self.pdb_info.total_residues()+1):
            if self.pdb_info.pose_to_record_map[ i ].get_pdb_num() == "+1":
                "Insertions found outside of CDR.  Will use insertion codes in PDB renumbering!"
    
    def _check_for_deletions(self, numbering_scheme):
        for i in range(1, self.pdb_info.total_residues()+1):
            if self.pdb_info.pose_to_record_map[ i ].get_pdb_num() == "-1":
                print "Deletions still found.  Please fix."
                
    def _has_cdr_type(self, cdr_type):
        found = False
        for i in range(1, self.pdb_info.total_residues()+1):
            if self.pdb_info.pose_to_record_map[ i ].get_pdb_num() == cdr_type:
                found = True
                break
        
        return found
    
    def _has_metadata(self, meta_string):
        found = False
        for i in range(1, self.pdb_info.total_residues()+1):
            if self.pdb_info.res( i ).get_meta() == meta_string:
                found = True
                break
        return found

    def _make_pdb_nums_as_integers(self):
        for i in range(1, self.pdb_info.total_residue()+1):
            self.pdb_info.res( i ).pdb_num =  int(self.pdb_info.res( i ).pdb_num)

    def _expand_c_term_num(self):
        """
        Extends the numbering of the C-terminal end or patches of -1.  Mainly for constant regions.  Will also number the ScFv.
        """
        
        expand_from = 0
        insert_chain = ' '
        chain_type = ' '
        for i in range(1, self.pdb_info.total_residues()+1):
            if self.pdb_info.res( i ).get_pdb_num() == '-1':
                if expand_from == 0:
                    continue
                else:
                    new_resnum = expand_from
                    for x in range(i, self.pdb_info.total_residues()+1):
                        if self.pdb_info.res( x ).get_pdb_num() != '-1':
                            break
                        new_resnum+=1
                        self.pdb_info.res( x ).set_pdb_num(new_resnum)
                        self.pdb_info.res( x ).set_region("")
                        self.pdb_info.res( x ).set_chain(insert_chain)
                        self.pdb_info.res( x ).set_chain_type(chain_type)
                        
            else:
                expand_from = int(self.pdb_info.res( i ).get_pdb_num())
                insert_chain = self.pdb_info.res( i ).get_chain()
                chain_type = self.pdb_info.res( i ).get_chain_type()
        
        #self._debug_print()
        
    def _expand_n_term_num(self):
        """
        Extends the numbering of the N-terminal end or patches of -1.  Mainly for constant regions.  Will also number the ScFv.
        """

        expand_from = self.pdb_info.total_residues()
        insert_chain = ' '
        chain_type = ' '
        for i in range(self.pdb_info.total_residues(), 0, -1):
            #print repr(i)
            #print self.numbering_map[i]['new_resnum']
            res_i = self.pdb_info.pose_to_record_map[ i ]
            if not isinstance( res_i, AntibodyResidue ): sys.exit()
            if res_i.pdb_num == '-1':
                if expand_from == self.pdb_info.total_residues():
                    continue
                else:
                    #print "Attempting N-terminal insertion expansion."
                    new_resnum = expand_from
                    for x in range(i, 0, -1):
                        res_x = self.pdb_info.pose_to_record_map[ x ]
                        if not isinstance( res_x, AntibodyResidue ): sys.exit()

                        if res_x.pdb_num != '-1':
                            break
                        if new_resnum != 1:
                            new_resnum-=1
                        
                        #Add insertion codes for end of normal Fv chain, so we don't have two res number 1s
                        if new_resnum == 1 == expand_from:
                            res_x.set_icode( "A" )
                            
                        res_x.set_pdb_num( new_resnum )
                        res_x.set_region( "")
                        res_x.set_chain( insert_chain )
                        res_x.set_chain_type( chain_type )
            else:
                expand_from = int(res_i.get_pdb_num())
                insert_chain = res_i.get_chain()
                chain_type = res_i.get_chain_type()
                
    def _expand_insertions_to_insertion_code(self):
        """
        Expands any remaining insertions that did not match the numbering scheme into same resnum as the last, with insertion residues.
        """

        codes = ["-1", "+1", "-2"]
        for i in range(2, self.pdb_info.total_residues()+1):
            res_i = self.pdb_info.pose_to_record_map[ i ]
            if not isinstance( res_i, AntibodyResidue ): sys.exit()

            if res_i.pdb_num == "+1":
                new_resnum = self.pdb_info.pose_to_record_map[i - 1].pdb_num
                if new_resnum in codes:
                    continue


                x = i
                icode_index = 0

                print "Extra insertion code found.  Inserting"
                res_x = self.pdb_info.pose_to_record_map[ x ]
                if not isinstance( res_x, AntibodyResidue ): sys.exit()

                while res_x.pdb_num == "+1":
                    icode = self.icodes[ icode_index ]

                    res_x.set_pdb_num( new_resnum )
                    res_x.set_icode( icode )

                    icode_index+=1
                    x+=1

    
    def _create_cdr_objects(self, numbering_scheme):
        """
        Creates CDRInfo objects. Used mainly for DatabasePyIgClassify.  
        CDRInfo Objects are used because a chain CAN have multiple CDRs of the same type.  Due to current antibody engineering methods(diabody, triabody, etc.).
        """
        
        CDRs = []
        i = 1
        while i < self.pdb_info.total_residues()+1:
            res_i = self.pdb_info.pose_to_record_map[ i ]
            if not isinstance(res_i, AntibodyResidue): sys.exit()

            if res_i.get_region() == "CDR":
                cdr_type = res_i.get_cdr_type()
                gene = res_i.get_gene()
                cdr_start = i
                cdr_end = i
                x = i
                #Get length of CDR
                for x in range(cdr_start, self.pdb_info.total_residues()+1):
                    if not self.pdb_info.pose_to_record_map[ x ].get_region() == "CDR":
                        break
                    else:
                        cdr_end+=1
                        i+=1
                i+=1
                
                #print "CDR: "+repr(cdr_start)+" "+repr(cdr_end)
                cdr_info = CDRInfo(cdr_type, cdr_start, cdr_end - 1, gene)
                CDRs.append(cdr_info)
            else:
                i+=1
                
        return CDRs
                    
                
    def create_ab_vset_chain(self, ig_chain, numbering_scheme):
        
        if not ig_chain: return None
        assert isinstance(ig_chain, IgChain)
        
        if ig_chain.get_type() != 'ig_ab':return None
        
        domain_to_numbering_map = defaultdict()
        for domain in ig_chain.get_domains():
            map_path = self.numbering_parser.get_full_map_path(numbering_scheme, domain.get_gene())
            domain_to_numbering_map[domain]=map_path
            
        self.pdb_info = renum_tools.get_seq_to_pdb_info(ig_chain, domain_to_numbering_map)
        sequence = ig_chain.sequence
        
        #Identify the CDR positions - add 'type' to dictionary.
        #Identify the CDR Types from the positions - add 'cdr_type' to dictionary.
        
        for cdr in self.cdr_types:
            self._identify_cdrs_from_anchors(numbering_scheme, cdr)
            for i in range(1, len(sequence)+1):
                res_i = self.pdb_info.pose_to_record_map[ i ]
                if not isinstance(res_i, AntibodyResidue): sys.exit()

                if res_i.pdb_num == cdr:
                    res_i.set_cdr_type( cdr )
                    res_i.set_region( "CDR" )

        for i in range(1, len(sequence)+1):
            res_i = self.pdb_info.pose_to_record_map[ i ]
            if not isinstance(res_i, AntibodyResidue): sys.exit()

            if not res_i.get_region() == "CDR":
                res_i.set_region("FRAME"); #Will be set to empty string in expand_c_term
        
        #Insert the true numbering to the seq_to_numbering map
        for cdr in self.cdr_types:
            self._create_cdr_numbering_in_renumbering_map(numbering_scheme, cdr)
        
        #Insert numbering of N and C terminal residues outside of the numbering maps.
        self._expand_c_term_num()
        self._expand_n_term_num()
        self._expand_insertions_to_insertion_code()

        """
        #Replace all new_resnum strings as integers
        for i in range(1, len(self.numbering_map)+1):
            if type(self.numbering_map[i]['new_resnum'] == str) and self.numbering_map[i]['new_resnum'] not in self.cdr_types:
                self.numbering_map[i]['new_resnum']  = int(self.numbering_map[i]['new_resnum'])
            else:
                print "Something is very wrong.  The CDR has not been identified and the resnum not converted to the numbering schemes number.  Please figure out what is wrong."
                print "This message comes from modules.chains.AbChainFactory."
                sys.exit()
        """
        cdrs = self._create_cdr_objects(numbering_scheme)
        self._make_pdb_nums_as_integers()

        ab_chain = AbChain(self.pdb_info, cdrs, ig_chain)


        #Values Optionally set (Mainly for PDBAA) - still needs clean up:
        if ig_chain.is_data_from_pdbaa():
            ab_chain.set_pdbaa_tag_info(ig_chain.get_pdbaa_tag_info())

        if ig_chain.hit_results:
            ab_chain.hit_results = ig_chain.hit_results

        if ig_chain.full_id:
            ab_chain.full_id = ig_chain.full_id
            ab_chain.pdb_res_info.set_extra_data("fullname", ab_chain.full_id)
            ab_chain.pdb_res_info.set_extra_data("description", ab_chain.description)

        ab_chain.pdb_res_info.set_extra_data("name", ab_chain.id)
        return ab_chain

            