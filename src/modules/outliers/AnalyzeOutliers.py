


import os
import sys

p = os.path.abspath("../../../")
sys.path.append(p); #Allows all modules to use all other modules, without needing to update pythonpath




import sqlite3
from collections import defaultdict
from optparse import OptionParser
import shutil

from src.modules import Structure
from src.modules import ClusterData
from src.tools.path import *

from rosetta import *
rosetta.init("-ignore_unrecognized_res -ignore_zero_occupancy false")



class AnalyzeOutliers:
    def __init__(self, db_path, cdr_dir, res_cutoff = 2.8, rfac_cutoff = .3):
        self.db_path = db_path
        self.db = sqlite3.connect(db_path)

        self.cluster_data = ClusterData.CDRData(self.db, cdr_dir, True)
        self.cluster_data.set_res_cutoff(res_cutoff)
        self.cluster_data.set_rfac_cutoff(rfac_cutoff)

        self.cluster_data.load_data()

        self.ab = Structure.Antibody_Structure()

    def analyze_rms(self, outdir, load_aligned = True, overhang = 3, add_data_to_current_db_table = False, db_name="outlier_analysis.db3"):

        rms_data = defaultdict()
        if not os.path.exists(outdir):
            os.mkdir(outdir)

        if add_data_to_current_db_table:
            shutil.copyfile(self.db_path, os.path.basename(self.db_path+"_temp.db3"))
            self.db.execute("ALTER TABLE cdr_data ADD COLUMN bb_rmsd_cdr_align")
            self.db.execute("ALTER TABLE cdr_data ADD COLUMN bb_rmsd_stem_align")

            #Set a default at an arbitrarily large number
            self.db.execute("UPDATE cdr_data SET bb_rmsd_cdr_align = -1 ")
            self.db.execute("UPDATE cdr_data SET bb_rmsd_stem_align = -1 ")
        else:
            outdb = sqlite3.connect(outdir+"/"+db_name)
            outdb.execute("DROP TABLE if exists cdr_data")
            outdb.execute("CREATE TABLE cdr_data( id INT, PDB TEXT, CENTER TEXT, original_chain TEXT, CDR TEXT, length INT, fullcluster TEXT, DistDegree FLOAT, bb_rmsd_cdr_align FLOAT, bb_rmsd_stem_align FLOAT)")

        by_cdr_dir = "by_cdr_culled"
        by_stem_dir = "by_stem_culled"
        if not os.path.exists(outdir+"/"+by_cdr_dir):
            os.mkdir(outdir+"/"+by_cdr_dir)
        if not os.path.exists(outdir+"/"+by_stem_dir):
            os.mkdir(outdir+"/"+by_stem_dir)

        if not os.path.exists(outdir+"/aligned_cdrs_stem"):
            os.mkdir(outdir+"/aligned_cdrs_stem")

        id = 1
        LOG = open_file(outdir+"/ALIGN_LOG.txt", 'w')
        print "Calculating RMSDs"
        for cdr in self.ab.get_CDRs():
            for length in self.cluster_data.get_lengths(cdr.name):
                for cluster in self.cluster_data.get_clusters(cdr.name, length):

                    data = self.cluster_data.get_cluster_data(cdr.name, length, cluster)
                    if isinstance(data, ClusterData.CDRClusterData): pass
                    #print "Loading "+cluster
                    ## Load center pose - which is required.  Skip if not found!
                    center_path = data.get_center_path()
                    if not center_path:
                        print "Could not find center for: "+cluster
                        LOG.write(cluster+" missing center data\n")
                        continue
                        #sys.exit("Temp")

                    if not os.path.exists(center_path):
                        print "Could not load center pose for cluster: "+cluster
                        print "CDR Pose does not exist: "+data.get_center_path()
                        LOG.write(cluster+" missing center pdb\n")
                        continue

                    center_pose = Pose()
                    print "Loading Center Pose "+data.get_center_name()
                    pose_from_pdb(center_pose, str(data.get_center_path()))

                    infos = data.get_infos()
                    for info in infos:
                        #print info.name
                        if isinstance(info, ClusterData.Info): pass

                        if info.name == data.get_center_name():
                            continue


                        #PDBs that fail for some fucked up reason of course.
                        ignored = ["1oayM_L2", "1oayO_L2"]
                        if info.name in ignored:
                            print "Ignoring "+info.name
                            LOG.write(cluster+" "+info.name+" ignored\n")
                            continue

                        #Strange antibody that Ben's renumbering code totally screws up on.
                        if info.name[0:4].lower() == "4hjj":
                            continue


                        print cluster +" "+info.name
                        if load_aligned and os.path.exists(outdir+"/"+by_cdr_dir+"/"+info.name+".pdb"):
                            all_path = outdir+"/"+by_cdr_dir+"/"+info.name+".pdb"
                            stem_path = outdir+"/"+by_stem_dir+"/"+info.name+".pdb"
                            if not os.path.exists(all_path):
                                print "Could not locate: "+all_path
                                print "Skipping... "
                                LOG.write(cluster+" "+info.name+" missing pdb\n")
                                continue

                            cdr_pose_all_aligned =  Pose()
                            pose_from_pdb(cdr_pose_all_aligned, str(all_path))
                            cdr_pose_stem_aligned = Pose()
                            pose_from_pdb(cdr_pose_stem_aligned, str(stem_path))
                        else:


                            if not os.path.exists(info.path):
                                print "Could not locate: "+info.path
                                print "Skipping... "
                                LOG.write(cluster+" "+info.name+" missing pdb\n")
                                continue

                            cdr_pose = pose_from_pdb(str(info.path))
                            cdr_pose_all_aligned = Pose()
                            cdr_pose_all_aligned.assign(cdr_pose)

                            cdr_pose_stem_aligned = Pose()
                            cdr_pose_stem_aligned.assign(cdr_pose)

                            if cdr_pose.total_residue() != length+overhang+overhang:
                                print cdr_pose
                                print cluster+" "+info.name+" missing residues in PDB "+repr(cdr_pose.total_residue()) +" vs "+repr(length+overhang+overhang)
                                LOG.write(cluster+" "+info.name+" missing residues in PDB "+repr(cdr_pose.total_residue()) +" vs "+repr(length+overhang+overhang)+"\n")
                                continue

                            if cdr_pose.total_residue() != center_pose.total_residue():
                                print cdr_pose
                                print cluster+" "+info.name+" not equal center pose length "+repr(cdr_pose.total_residue()) +" vs "+repr(center_pose.total_residue())
                                LOG.write(cluster+" "+info.name+" not equal center pose length "+repr(cdr_pose.total_residue()) +" vs "+repr(center_pose.total_residue())+"\n")
                                continue


                            align_to_cluster_center_save_pdb(cdr, info.name, cdr_pose_all_aligned, center_pose, outdir+"/"+by_cdr_dir, overhang)
                            align_to_cluster_center_save_pdb(cdr, info.name, cdr_pose_stem_aligned, center_pose, outdir+"/"+by_stem_dir, overhang, True)
                        if cdr_pose_all_aligned.total_residue() != length+overhang+overhang:
                            print cdr_pose_all_aligned
                            print cluster+" "+info.name+" missing residues in PDB "+repr(cdr_pose_all_aligned.total_residue())+" vs "+ repr(length+overhang+overhang)
                            LOG.write(cluster+" "+info.name+" missing residues in PDB "+repr(cdr_pose_all_aligned.total_residue())+" vs "+ repr(length+overhang+overhang)+"\n")
                            continue

                        if cdr_pose_all_aligned.total_residue() != center_pose.total_residue():
                            print cdr_pose_all_aligned
                            print cluster+" "+info.name+" not equal center pose length "+repr(cdr_pose_all_aligned.total_residue())+" vs "+ repr(center_pose.total_residue())
                            LOG.write(cluster+" "+info.name+" not equal center pose length "+repr(cdr_pose_all_aligned.total_residue())+" vs "+ repr(center_pose.total_residue())+"\n")
                            continue


                        all_rms = get_rmsd(cdr, cdr_pose_all_aligned, center_pose, overhang)
                        stem_rms = get_rmsd(cdr, cdr_pose_stem_aligned, center_pose, overhang)

                        print "All: "+repr(all_rms)
                        print "Stem: "+repr(stem_rms)


                        if add_data_to_current_db_table:
                            columns = [all_rms, cdr.name, length, info.PDB, info.original_chain, info.cluster]

                            cur = self.db.cursor()
                            cur.execute("UPDATE cdr_data SET bb_rmsd_cdr_align = ? WHERE CDR=? AND length=? AND PDB=? AND original_chain=? and fullcluster=?", columns)
                            cur.close()

                            columns = [stem_rms, cdr.name, length, info.PDB, info.original_chain, info.cluster]

                            cur = self.db.cursor()
                            cur.execute("UPDATE cdr_data SET bb_rmsd_stem_align = ? WHERE CDR=? AND length=? AND PDB=? AND original_chain=? and fullcluster=?", columns)
                            cur.close()

                        else:
                            columns = [id, info.PDB, data.get_center_name()[0:5].upper(), info.original_chain, cdr.name, length, info.cluster, info.dihedral_distance, all_rms, stem_rms]

                            with outdb:
                                cur = outdb.cursor()
                                cur.execute("INSERT INTO cdr_data VALUES(?,?,?,?,?,?,?,?,?,?)", columns)
                            id+=1
                    #sys.exit()
        LOG.close()
        if add_data_to_current_db_table:
            os.remove(os.path.basename(self.db_path+"_temp.db3"))


    ######################################################


    def plot_data(self):
        pass

def attempt_align_and_save(cdr, cdr_path, center_path ,out_path, stem_align = False, overhang = 3):
    """
    Aligns and saves both by CDR and Overhnage alignments.  Returns False if not successfull
    """

    if not os.path.exists(cdr_path):
        print "Could not locate: "+cdr_path
        print "Skipping... "
        return False

    center_pose = pose_from_pdb(str(center_path))

    cdr_pose = pose_from_pdb(str(cdr_path))

    if cdr_pose.total_residue() != center_pose.total_residue():
        print "CDR Pose not equal to center pose length!  Ignoring!!"
        return False


    align_to_cluster_center_save_pdb(cdr, os.path.basename(cdr_path).split(".")[0:-1], cdr_pose, center_pose, out_path, overhang, stem_align)

def align_to_cluster_center_save_pdb( cdr, pose_name, pose, center_pose, outdir, overhang=3, stem_align = False):


    rosetta.core.pose.remove_lower_terminus_type_from_pose_residue(pose, 1)
    rosetta.core.pose.remove_upper_terminus_type_from_pose_residue(pose, pose.total_residue())

    if stem_align:
        id_map = get_mask_for_stem_alignment(pose, center_pose, cdr, overhang)
    else:
        id_map = get_mask_for_alignment(pose, center_pose, cdr, overhang)


    superimpose_pose(pose, center_pose, id_map)
    pose.dump_pdb(str(outdir+"/"+pose_name+".pdb"))

def get_rmsd( cdr, pose, center_pose, overhang = 3):


    #id_map = get_mask_for_alignment(pose, center_pose, cdr, overhang)

    #rms = rms_at_corresponding_atoms_no_super(pose, center_pose, id_map)

    start = pose.pdb_info().pdb2pose(cdr.get_pdb_chain(), cdr.get_pdb_start())
    end = pose.pdb_info().pdb2pose(cdr.get_pdb_chain(), cdr.get_pdb_end())
    l = Loop(start, end)
    loops = Loops()
    loops.push_back(l)

    rms = loop_rmsd(pose, center_pose, loops, False, True)

    return rms

def get_mask_for_alignment(cdr_pose, center_pose, cdr, overhang = 3):
    id_map = AtomID_Map_AtomID()
    rosetta.core.pose.initialize_atomid_map_AtomID(id_map, cdr_pose, AtomID(0,0))

    start = cdr_pose.pdb_info().pdb2pose(cdr.get_pdb_chain(), cdr.get_pdb_start())
    end = cdr_pose.pdb_info().pdb2pose(cdr.get_pdb_chain(), cdr.get_pdb_end())

    for i in range(start, end+1):
        #print repr(i)
        for ii in range(1, 4+1):
            #print repr(ii)

            atom_cdr = AtomID(ii, i)
            atom_center = AtomID(ii, i)

            id_map.set(atom_cdr, atom_center)

    return id_map

def get_map_for_rmsd(cdr_pose, center_pose, cdr, overhang=3):
    pass
    m = rosetta.utility.rosetta.utility.map_string_Real()

def get_mask_for_stem_alignment(cdr_pose, center_pose, cdr, overhang = 3):
    id_map = AtomID_Map_AtomID()
    rosetta.core.pose.initialize_atomid_map_AtomID(id_map, cdr_pose, AtomID(0,0))

    start = cdr_pose.pdb_info().pdb2pose(cdr.get_pdb_chain(), cdr.get_pdb_start())
    end = cdr_pose.pdb_info().pdb2pose(cdr.get_pdb_chain(), cdr.get_pdb_end())

    for i in range(1, start):
        #print repr(i)
        for ii in range(1, 4+1):
            #print repr(ii)

            atom_cdr = AtomID(ii, i)
            atom_center = AtomID(ii, i)

            id_map.set(atom_cdr, atom_center)

    for i in range(end+1, cdr_pose.total_residue()+1):
        #print repr(i)
        for ii in range(1, 4+1):
            #print repr(ii)

            atom_cdr = AtomID(ii, i)
            atom_center = AtomID(ii, i)

            id_map.set(atom_cdr, atom_center)

    return id_map


if __name__ == "__main__":

    parser = OptionParser()

    ############################
    ## Required Options
    ############################
    parser.add_option("--cdr_dir",
                      help = "Structures of renumbered split CDRs",
                      default="/home/jadolfbr/Documents/modeling/databases/antibody_databases/PyIgClassify/DBOUT/cdr_pdbs_redun_by_cdr_overhang_3")

    parser.add_option("--ab_db",
                      help = "rel path to log dir for benchmark",
                      default="/home/jadolfbr/Documents/modeling/databases/antibody_databases/PyIgClassify/DBOUT/website/antibody_database_redundant.db")


    ############################
    ## Optional
    ############################
    parser.add_option("--skip_alignment",
                      help = "Skip the alignment and dumping part.  Load pre-aligned structure",
                      default = False,
                      action = "store_true")

    parser.add_option("--outdir",
                      help = "Main output directory",
                      default = os.path.split(os.path.abspath(__file__))[0])

    parser.add_option("--add_data_to_current_db",
                      help = "Add the data as a new column in the current DB",
                      default = False,
                      action = "store_true")

    parser.add_option("--out_db_name",
                      help = "If NOT adding data to current db, use this name as the output database",
                      default = "outlier_analysis.db3")



    ############################
    ## R Factor and Resolution
    ############################
    parser.add_option("--res_cutoff",
                      help = "Cutoff for resolution",
                      default = 2.8)

    parser.add_option("--rfac_cutoff",
                      help = "Cutoff for Rfactor",
                      default = .3)


    (options, args) = parser.parse_args(sys.argv)


    analyzer = AnalyzeOutliers(options.ab_db, options.cdr_dir, float(options.res_cutoff), float(options.rfac_cutoff))
    analyzer.analyze_rms(options.outdir, options.skip_alignment, 3, options.add_data_to_current_db, options.out_db_name)