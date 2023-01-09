import re
import sys
from collections import defaultdict
from copy import deepcopy

from modules.Structure import AntibodyResidue
from src.modules.Structure import PDBResInfo
from src.modules.chains.IgDomain import IgDomain

from src.modules.chains.AbChain import AbChain
from src.tools.path import *


def get_conensus_bt_pdb_nums(pdb_consensus_info, start, stop, chain):
    """
    Assumes no icode
    """
    con = ""
    for i in range(start, stop+1):
        position = (i, chain, " ")
        con = con + pdb_consensus_info.get_consensus_for_position(position)
    return con

def get_seq_between_pdb_nums(pdb_res_info, start, stop, chain):
    seq = ""
    for i in range(1, pdb_res_info.total_residues()+1):
        res = pdb_res_info.pose_to_record_map[ i ]
        if res.get_pdb_num() >= start and res.get_pdb_num() <= stop and res.get_chain() == chain:
            seq = seq+res.get_aa()

    return seq

def get_nr_framework_res_info_array(pdb_res_info_array):

    #First, we get all unique sequences.  These may still be the same, but have more residues on either end of the antibody.
    unique_seq = defaultdict()
    for pdb_info in pdb_res_info_array:
        frame_seq = get_raw_framework_sequence(pdb_info)
        unique_seq[frame_seq] = pdb_info

    #Now we have to compare the overlap to account for deletions in the beginning or end of the CDR
    #If they both have residue 1 ->, we take the longest framework.

    #### I cannot figure this out right now.  For now, we skip this step. ###
    final_array = [unique_seq[ seq ] for seq in unique_seq]
    return final_array


    final_dict = defaultdict()
    overlaps = defaultdict()

    for seq_i in unique_seq.keys():
        similar_seq_found = 0
        for seq_x in unique_seq.keys():

            similar, overlap_seq = compare_frame_overlap(unique_seq[ seq_i ], unique_seq[ seq_x ])
            if similar:
                similar_seq_found+=1
                if not overlaps[ overlap_seq ].has_key(overlap_seq):
                    overlaps[ overlap_seq ] = [ unique_seq[ seq_i ]] #This is still wrong

                overlaps[ overlap_seq ].append( unique_seq[ seq_x ])

            #Do something if no similar seq are found
            if similar_seq_found == 0:
                final_dict[ seq_i ] = unique_seq[ seq_x ]

def compare_frame_overlap(pdb_info_i, pdb_info_j):
    pass


    #Not right to compare at index.  What about insertion codes?

    #What are we returning here?

def get_overlap_start_indexes(pdb_info_i, pdb_info_x):
    def get_residue_start(pdb_info):
        for i in range(1, pdb_info.total_residues()+1):
            if pdb_info.pose_to_record_map[ i ].get_region() == "FRAME":
                return i, pdb_info.pose_to_record_map[ i ]

    i_index, i_start_res = get_residue_start(pdb_info_i)
    x_index, x_start_res = get_residue_start(pdb_info_x)

    if i_start_res.get_pdb_num() == x_start_res.get_pdb_num():
        return i_index, x_index

    elif i_start_res.get_pdb_num() > x_start_res.get_pdb_num():
        pass


def get_missing_n_term_residues(pdb_res_info):

    start_num = pdb_res_info.pose_to_record_map[1].get_pdb_num()
    missing_residues = start_num - 1
    return missing_residues

def get_raw_framework_sequence(pdb_res_info, stop_at_res_148 = True):
    """
    Get the concatonated framework sequence (without constant or cdrs)
    """
    if not isinstance(pdb_res_info, PDBResInfo): sys.exit()

    sequence = ""
    for res in pdb_res_info.get_all_residues():
        if not isinstance(res, AntibodyResidue): sys.exit()

        if stop_at_res_148:
            if res.get_region() == "FRAME" and int(res.get_pdb_num()) < 149:
                sequence = sequence + res.get_aa()

        elif res.get_region() == "FRAME":
            sequence = sequence + res.get_aa()

    return sequence

def get_raw_framework_sequence_in_frame(pdb_res_info, index_start, index_end):
    if not isinstance(pdb_res_info, PDBResInfo): sys.exit()

    sequence = ""
    for i in range(index_start, index_end+1):
        res = pdb_res_info.pose_to_record_map[ i ]
        if not isinstance(res, AntibodyResidue): sys.exit()


        if res.get_region() == "FRAME":
            sequence = sequence + res.get_aa()

    return sequence


def get_new_pdb_infos_for_gene(ab_chain_array, gene, skip_scFv = True):
    """
    Get an array of PDBResInfos that have the gene.  This works with SCFvs, so we get chain-specific arrays.
    Makes a deep copy of the residue for the new array
    """
    genes = ['lambda', 'kappa', 'heavy', 'lambda6']

    new_infos = []
    for abchain in ab_chain_array:
        if skip_scFv and abchain.is_ScFv():
            continue

        pdb_info = abchain.get_pdb_res_info()

        #Skip BAD PDB entries if parsing pdbaa.
        if pdb_info.get_extra_data('name')[0:4] in ["4HJJ", "7FAB", "1OCW"]: continue

        if not isinstance(pdb_info, PDBResInfo): sys.exit()
        new_info = PDBResInfo()
        new_info.extra_data = pdb_info.extra_data

        #print repr(pdb_info.extra_data)
        for res in pdb_info.get_all_residues():
            if not isinstance(res, AntibodyResidue): sys.exit()
            if res.get_gene() == gene:
                #print res
                res_copy = deepcopy(res)
                new_info.add_residue(res_copy)


        new_infos.append(new_info)

    return new_infos

def replace_lambda6_with_lambda_gene(ab_chain_array):
    """
    Replace the variant lambda6 with lambda.  Lambda6 is only different in one region, and is still a lambda chain.
    Lambda6 is mainly used for renumbering, since there is an insertion in a non-cdr loop.
    """
    for abchain in ab_chain_array:
        if not isinstance(abchain, AbChain): sys.exit()

        for domain in abchain.get_domains():
            if not isinstance(domain, IgDomain): sys.exit()
            if domain.gene == 'lambda6':
                domain.gene = 'lambda'

        for i in range(1, abchain.get_pdb_res_info().total_residues()+1):
            if abchain.get_pdb_res_info().pose_to_record_map[ i ].get_gene() == 'lambda6':
                abchain.get_pdb_res_info().pose_to_record_map[ i ].set_gene('lambda')


def load_centers(filepath=False):
    """
    Loads center data from Ben's median output file.
    """
    if not filepath:
        FILE = open_file(get_db_path()+"/medians.txt")
    else:
        FILE = open_file(filepath)

    FILE.readline()
    center_data = defaultdict(dict)

    #New Style:
    for line in FILE:
        line = line.strip()
        lineSP = line.split(',')
        pdb_chain = lineSP[2][0:5].upper()
        pdb = pdb_chain[0:4].upper()
        cdr = lineSP[2].split("_")[1].upper()
        #Cluster Identification (Cis/Trans):
        SS = lineSP[0].split("-")[2]

        length_name = get_loop_key_from_ss(SS)
        cluster = cdr+"-"+length_name+"-"+lineSP[1]
        center_data[pdb_chain][cdr] = cluster

    FILE.close()
    return center_data


def get_best_evalue_hmm(result):
    """
    Returns hmm name of best Evalue from a result.
    Result is a dict of dict with hmm or hmm_pair [hmm, hsp] as first key and 'evalue' as second key.
    """
    
    assert isinstance(result, dict)
    
    hmm = min(result.keys(), key=lambda x:result[x]['evalue'])
    return hmm

def get_best_score_hmm(result):
    """
    Returns hmm name of best score from a result.
    Result is a dict of dict with hmm or hmm_pair [hmm, hsp] key and 'score' as second key.
    """
    
    assert isinstance(result, dict)
    
    hmm  = min(result.keys(), key=lambda x:result[x]['score'])
    return hmm
    
def get_best_score_pair(hmm_pair_list):
    """
    Get best pair [hmm,hsp] in list of pairs based on score.
    """
    
    best_pair = hmm_pair_list[0]
    for pair in hmm_pair_list:
        if pair[1].bitscore > best_pair[1].bitscore:
            best_pair = pair
    
    return best_pair

def get_best_evalue_pair(hmm_pair_list):
    """
    Get best pair [hmm, hsp] in a list of pairs based on evalue.
    """
    
    best_pair = hmm_pair_list[0]
    for pair in hmm_pair_list:
        if pair[1].bitscore < best_pair[1].evalue:
            best_pair = pair
    
    return best_pair

def get_loop_key_from_ss(SS):
    """
    Gets the cluster key from secondary structure.
    Example: "TTTCTT", returns cis4
    Example: "TTTCTC", returns cis4,6
    Example: "TTTTTT", returns 6
    """
    name = ''
    if re.search("C", SS):
        result = re.search("C", SS)
        name = "cis"+repr(result.start()+1)
        for i in range(result.start()+1, len(SS)):
            position = i+1
            if SS[i]=="C":
                name = name+','+repr(position)
    else:
        name = str(len(SS))
        
    return name

def get_rama_type(dih_array):
    """
    Gets the rama type as a charactor from an array [phi, psi]
    Does not do upper/lower based on cis/trans residues.
    B is the beta-sheet region, P is polyproline II, A is the alpha-helix, D is the delta region (near alpha-helix but at more negative values of phi),
    L is the left-handed helix, and G is the gamma region (phi > 0, excluding the L and B regions). 
    """
    phi = dih_array[0]
    psi = dih_array[1]
    
    if phi > 180.0:
        phi = -360.0+phi
    
    if psi > 180.0:
        psi = -360.0+psi
    
    if phi == -180:
        phi = 180
        
    if psi == -180:
        psi = 180
        
    regions = defaultdict(list)
    
    #B, P, G, B, L, D, A
    #Regions denotes Regionstring + list of rama array:  -phi, +phi, -psi, +psi
    regions["B"].append([-180, -100, 50, 180])
    regions["B"].append([-180, -100, -180, -100])
    regions["D"].append([-180, -100, -100, 50])
    regions["P"].append([-100, 0, -180, -100])
    regions["P"].append([-100, 0, 50, 180])
    regions["A"].append([-100, 0, -100, 50])
    regions["G"].append([0, 150, -180, -50])
    regions["G"].append([0, 150, 100, 180])
    regions["B"].append([150, 180, 100, 180])
    regions["B"].append([150, 180, -180, -50])
    regions["L"].append([0, 180, -50, 100])
    
    for region in regions:
        for rama in regions[region]:
            
            if (rama[0] < phi <= rama[1]) and (rama[2] < psi <= rama[3]):
                return region

def reformat_dihedrals(dihedral_string, split_char=':'):
    """
    This script simply reformats a string of dihedrals that are phi:psi:omega:phi:psi:omega to be phi,phi psi,psi omega,omega.
    This is used for cluster_center_dihedrals.txt within Rosetta to make reading the file easier.  Eventually this format will be changed.
    """

    dihSP = dihedral_string.split(split_char)

    print "Dihedral length: "+repr(len(dihSP))

    phis=[]; psis=[];omegas=[]
    phi_index = 0; psi_index = 1; ome_index=2;
    delim=3


    phis.append(dihSP[phi_index])
    psis.append(dihSP[psi_index])
    omegas.append(dihSP[ome_index])

    for i in range(1, (len(dihSP)/3)):
        #print i
        phi_index = phi_index+delim
        phis.append(dihSP[phi_index])

        psi_index = psi_index+delim
        psis.append(dihSP[psi_index])

        ome_index = ome_index+delim
        omegas.append(dihSP[ome_index])


    #print "Phis: "+repr(len(phis)) + "Psis: "+repr(len(psis))+"omes: "+repr(len(omegas))

    s = ","
    newString = s.join(phis) +" "+ s.join(psis)+" "+s.join(omegas)

    return newString