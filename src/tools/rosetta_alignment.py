

##### These are all for same length poses #####

from rosetta import *

def align_to_second_pose_save_pdb( pose_name, pose, second_pose, outdir, overhang=0, stem_align = False):


    rosetta.core.pose.remove_lower_terminus_type_from_pose_residue(pose, 1)
    rosetta.core.pose.remove_upper_terminus_type_from_pose_residue(pose, pose.total_residue())

    rosetta.core.pose.remove_lower_terminus_type_from_pose_residue(second_pose, 1)
    rosetta.core.pose.remove_upper_terminus_type_from_pose_residue(second_pose, second_pose.total_residue())

    if stem_align:
        id_map = get_mask_for_stem_alignment(pose, second_pose, overhang)
    else:
        id_map = get_mask_for_alignment(pose, second_pose, overhang)


    superimpose_pose(pose, second_pose, id_map)
    pose.dump_pdb(str(outdir+"/"+pose_name+".pdb"))


def get_rmsd( pose, second_pose, overhang = 0):
    """
    Get RMSD assuming they are both the same length!
    """

    #id_map = get_mask_for_alignment(pose, second_pose, cdr, overhang)

    #rms = rms_at_corresponding_atoms_no_super(pose, second_pose, id_map)

    start = 1 + overhang
    end = pose.total_residue() - overhang
    l = Loop(start, end)
    loops = Loops()
    loops.push_back(l)

    rms = loop_rmsd(pose, second_pose, loops, False, True)

    return rms

def get_mask_for_alignment(pose, second_pose, overhang = 0):
    """
    Get mask assuming they are both the same length!
    """
    id_map = AtomID_Map_AtomID()
    rosetta.core.pose.initialize_atomid_map_AtomID(id_map, pose, AtomID(0,0))

    start = 1 + overhang
    end = pose.total_residue() - overhang

    for i in range(start, end+1):
        #print repr(i)
        for ii in range(1, 4+1):
            #print repr(ii)

            atom_cdr = AtomID(ii, i)
            atom_center = AtomID(ii, i)

            id_map.set(atom_cdr, atom_center)

    return id_map

def get_map_for_rmsd(pose, second_pose, overhang=3):
    pass
    m = rosetta.utility.rosetta.utility.map_string_Real()

def get_mask_for_stem_alignment(pose, second_pose, stem_size):
    id_map = AtomID_Map_AtomID()
    rosetta.core.pose.initialize_atomid_map_AtomID(id_map, pose, AtomID(0,0))

    start = 1

    for i in range(start, start+stem_size):
        #print repr(i)
        for ii in range(1, 4+1):
            #print repr(ii)

            atom_cdr = AtomID(ii, i)
            atom_center = AtomID(ii, i)

            id_map.set(atom_cdr, atom_center)

    for i in range(0, stem_size):
        #print repr(i)
        for ii in range(1, 4+1):
            #print repr(ii)

            #print "pose"+ repr(pose.total_residue() - i)+" second_pose "+repr(second_pose.total_residue() - i)
            atom_pose = AtomID(ii, pose.total_residue() - i)
            atom_second_pose = AtomID(ii, second_pose.total_residue() - i)

            id_map.set(atom_pose, atom_second_pose)

    return id_map
