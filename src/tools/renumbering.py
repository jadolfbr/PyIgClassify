#Author: Jared Adolf-Bryfogle
#A collection of functions for renumbering withing PyIgClassify.  A forewarning: Renumbering is not fun, and the code is not elegant.  If you can
# make it elegant, by all means go for it.

#System
import copy
import sys
from collections import defaultdict

#Biopthon
from Bio.PDB.Chain import Chain as BioChain

#PyIgClassify
from src.modules.chains.AbChain import AbChain
from src.modules.chains.IgChain import IgChain
from src.modules.Structure import PDBResInfo
from src.modules.Structure import AntibodyResidue
from src.modules.restype_definitions import definitions
from src.tools.path import *

#class QueryToHMM(object):
    #def __init__(sequence, hmm_sequence, hit_sequence, query_start):

def get_seq_to_pdb_info(ig_chain, domain_to_numbering_map_file):
    """
    Creates and returns a NumberingMap dictionary for an IgChain.  Reads the numbering_map_file and the query_result in the ig_chain.

    Raw seq to PDBResInfo class.  Other methods find the CDRs(AbChainFactory) or other regions of interest
    Returns PDBResInfo.  Meta data such as L1_N for an L1 Nterminus is is 'meta'
      Residue of is the single letter amino acid code from the query
      Resnum is a string that shows the match to the renumbering_map.
      -1 if that residue does not have an alignment (Like some wierd N-terminal addition on a variable chain or the C-terminal Constant region)
      +1 if that residue is an addition (Such as in a CDR)
      -2 if that residue is in the CDR and is identified by the renumbering map (At least for the vset renumber maps)
      
       
    """
    #This is a bitch.  Andreas, yea, respect man.
    
    sequence = ig_chain.sequence
    pdb_info = PDBResInfo()
    assert isinstance(ig_chain, IgChain)
    
    #print ":"+sequence+":"
    for i in range(1, len(sequence)+1):

        res = AntibodyResidue(sequence[i-1], "-1", '', ' ')
        res.set_extra_info('cluster_dis', ' ')
        res.set_extra_info('cluster', ' ')
        pdb_info.add_residue(res)

    assert len(sequence) == pdb_info.total_residue()

    for domain in domain_to_numbering_map_file:
        #print domain_to_numbering_map_file[domain]
        numbering_map = read_numbering_map(domain_to_numbering_map_file[domain])
        hsp = domain.get_hsp()
        if not hsp:continue
        
        hmm_sequence = hsp.query_all[0].seq
        hit_sequence = hsp.hit_all[0].seq
        
        assert len(hmm_sequence) == len(hit_sequence)
        
        new_seq_position = hsp.hit_start+1
        num_map_position = hsp.query_start
        
        """
        print hmm_sequence
        print hit_sequence
        print repr(num_map_position)
        print repr(new_seq_position)
        """
        
        for i in range(0, len(hmm_sequence)):
            #print repr(new_seq_position)
            hmm_res = hmm_sequence[i]
            hit_res = hit_sequence[i]

            res = pdb_info.residue(new_seq_position)
            if not isinstance(res, AntibodyResidue): sys.exit()

            res.set_extra_info('gene', domain.get_gene())

            if hit_res == '-':
                #print "Deletion detected."
                num_map_position+=1
                continue
            
            elif hmm_res == '.':
                #print "Addition detected."
                res.set_pdb_num("+1")
                res.set_chain(domain.get_new_chain())
                res.set_chain_type(domain.get_new_chain())

                new_seq_position+=1
            else:
                res.set_pdb_num(numbering_map[num_map_position][1])
                res.set_chain(domain.get_new_chain())
                res.set_meta(numbering_map[num_map_position][2])
                res.set_chain_type(domain.get_new_chain())

                num_map_position+=1
                new_seq_position+=1
    
    #print repr(result_numbering_map)
    assert len(sequence) == pdb_info.total_residue()

    return pdb_info

def has_id(model, id):
    """
    Returns true or false if the model has the chain.  Because biopython is not updating it's index that has_id is using.  WTF.
    """
    for i in model:
        if i.id == id:
            return True
    return False

def renumber_biopython_chain(structure, abchain):
    model = structure[0]

    chain = model[abchain.get_id()]
    
    print "Old chain: "+abchain.get_id()

    chain_index = get_new_chain_index_for_old_chain(abchain)
    
    #print abchain.get_original_chain()
    assert(get_chain_length(chain) == abchain.num_residues())
    assert(abchain.sequence == get_biochain_sequence(chain))
    
    alphabet = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"

    if not isinstance(abchain, AbChain): sys.exit()

    if len(chain_index)==1:
        seq_position = 1
        
        #This makes sure you don't have multiple H or Ls in the resulting structure - the structure is still renumbered, with other letters - A,B,C,D, etc.
        if not has_id(model, abchain.get_new_chain(1)) or chain.id == abchain.get_new_chain(1):
            chain.id = abchain.get_new_chain(1)
        else:
            for letter in alphabet:
                if not has_id(model, letter):
                    chain.id = letter
                    change_chain_in_numbering_map(abchain.get_new_chain(1), letter, abchain)
                    break
        
        #print "Renumbering biopython chain "+chain.id

        pdb_res = abchain.get_pdb_res_info().residue(seq_position)
        if not isinstance(pdb_res, AntibodyResidue): sys.exit()
        for res in chain:
            #print seq_position
            if not res.id[0]==' ':
                #Skipping any heteroatoms
                continue

            pdb_res = abchain.get_pdb_res_info().residue(seq_position)
            pdb_res.set_old_resnum(res.id[1], res.id[2])
            
            res.id = tuple([' ', int(pdb_res.get_pdb_num()), pdb_res.get_icode()])
            seq_position+=1
    else:
        #print "Renumbering biopython chain "+chain.id
        new_chains = defaultdict()
        for chain_id in chain_index:
            new_chain = create_new_chain_from_old(chain, chain_id, chain_index[chain_id][0], chain_index[chain_id][1], abchain)
            new_chains[new_chain.id] = new_chain
            
        fix_model_with_old_chain(model, chain, chain_index, new_chains)
        
        #Here, we first check the chain id in the model.  If it already exists, lets say multiple L or multiple H's in the model, we give it a new letter.  We then fix up the numbering map.
        for id in new_chains:
            new_chain = new_chains[id]
            if has_id(model, id):
                for letter in alphabet:
                    if not has_id(model, letter):
                        new_chain.id = letter
                        change_chain_in_numbering_map(id, letter, abchain)
                        break
            model.add(new_chain)
    
    print "Chain Renumbered via biopython"
    
    
###########################################

def change_chain_in_numbering_map(old_chain, new_chain, abchain):
    
    for i in range(1, abchain.get_pdb_res_info().total_residue()+1):
        res = abchain.get_pdb_res_info().get_residue(i)
        if not isinstance(res, AntibodyResidue): sys.exit()
        if res.chain == old_chain:
            res.chain = new_chain


def renumber_pose_chains(pose, abchain_array):
    """
    Renumber the PDBInfo object in a pose using an array of AbChains.
    """
    
    for abchain in abchain_array:
        renumber_pose_chain(pose, abchain)

def get_new_chain_index_for_old_chain(abchain):
    """
    Gets an index for new chains for biopython renumbering.  Mainly for SCFV.
    """
    
    chain_index = defaultdict(); #[new_chain] = [start, end] (Inclusive)
    start = 1
    end = 0
    start_chain = abchain.get_new_chain(1)
    seq_position = 0
    #print "Residues: "+repr(abchain.num_residues())
    if start_chain == abchain.get_new_chain(abchain.num_residues()):
        chain_index[start_chain]=[start, abchain.num_residues()]
        return chain_index
    
    for seq_position in range(1, abchain.num_residues()+1):
        #print repr(seq_position)+" "+abchain.get_new_chain(seq_position)

        if seq_position == abchain.num_residues():
            end = seq_position
            chain_index[start_chain]=[start, end]
            #print "Last residue"
        elif start_chain == abchain.get_new_chain(seq_position):
            end = seq_position
        else:
            end = seq_position
            chain_index[start_chain]=[start, end-1]
            start = seq_position
            
            start_chain = abchain.get_new_chain(seq_position)
            #print "SCFV Switch"
            
    return chain_index

def get_chain_length(bio_chain):
    
    l = 0
    for res in bio_chain:
        if res.id[0]==' ':
            l+=1
    return l

def get_biochain_sequence(bio_chain):
    seq = ""
    d = definitions()
    
    for res in bio_chain:
        if res.id[0]==' ':
            aa = d.get_one_letter_from_three(res.resname)
            if not aa:
                print "Using X for non-canonical resname: "+res.resname
                seq = seq+'X'
            seq = seq+aa
    return seq

def create_new_chain_from_old(old_chain, new_chain_id, start, end, abchain):
    """
    Creates a new chain as a subset of old chain for biopython renumbering and returns the new chain.
    """
    
    new_chain = BioChain(new_chain_id)
    chain = copy.deepcopy(old_chain)
    seq_position =  1
    i = 1
    #print "Domain start: "+repr(start)+" End: "+repr(end)
    for res in chain:
        #print "Old res: "+repr(res)
        if start <= seq_position <= end:
            #print seq_position
            new_res = res
            if new_res.id[0] == ' ':
                new_chain.add(new_res)

                pdb_residue = abchain.get_pdb_res_info()[seq_position]
                if not isinstance(pdb_residue, AntibodyResidue): sys.exit()

                pdb_residue.set_old_resnum(res.id[1], res.id[2])
                new_res.id = tuple([' ', int(pdb_residue.get_pdb_num()), pdb_residue.get_icode()])
            #print "New res: "+repr(res)
            
            
                
            
            i+=1
        if res.id[0]==' ':
            seq_position+=1
    
    return new_chain
            
def fix_model_with_old_chain(model, chain, chain_index, new_chains):
    """
    Fixes the model with hetatms and the correct residues.  Creates a new_chain with those residues and hetatms that are not renumbered.  If the chainid exists already, we add the info to the chain (L -> L, etc.)
    Then, this chain is detached from the model, and if the chainid does not exist in the new_chains, we attach the chain to the model.
    """
    
    def find_residue(i, chain_index):
        """
        Returns true or false if i in found in between start/end of entries in chain index.
        """
        for id in chain_index:
            if chain_index[id][0] <= i <= chain_index[id][1]:
                return True
            else:
                continue
        return False
    
    
    i=0
    new_chain = BioChain(chain.id)
    for res in chain:
        
        #Copy hetatm from old_chain to another new_chain or add it to a chain from new_chains, which contain the renumbered residues.
        if res.id[0] != ' ':
            #new_res = copy.deepcopy(res)
            new_res = res
            if not new_chains.has_key(new_chain.id):
                new_chain.add(new_res)
            else:
                new_chains[new_chain.id].add(new_res)
                
            continue
        else:
            i+=1
        
        #If the residue is not a hetatm and is not between the new_chains, it adds the residue to the new_chain, or adds it to a chain from new_chains, which contain the renumbered residues.
        if not find_residue(i, chain_index):
            #new_res = copy.deepcopy(res)
            #print repr(i)
            new_res = res
            if not new_chains.has_key(new_chain.id):
                new_chain.add(new_res)
            else:
                new_chains[new_chain.id].add(new_res)
    
    
    model.detach_child(chain.id)
    if not new_chains.has_key(new_chain.id):
        model.add(new_chain)
                
def read_numbering_map(numbering_map_file):

    numbering_map = []
    FILE = open_file(numbering_map_file, 'r')
    for line in FILE:
        if line[0] == "#": continue
        lineSP = line.split()
        if len(lineSP) < 2: continue

        query = lineSP[0]
        num = lineSP[1]
        meta = ""
        if len(lineSP) ==3:
            meta = lineSP[2]
        tripple = [query, num, meta]
        numbering_map.append(tripple)

    return numbering_map
    
    