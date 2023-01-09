#!/usr/bin/python

p = os.path.abspath("../../../../")
sys.path.append(p); #Allows all modules to use all other modules, without needing to update pythonpath

import os

from modules import ClusterData
from modules import Structure
from modules.PyMolScriptWriter import PyMolScriptWriter


from src.tools import AbDbFunctions
from src.modules.outliers import AnalyzeOutliers

import sys
from optparse import OptionParser
from collections import defaultdict
import sqlite3

class VisualizeCDRs:
    """
    Create A PyMol script for visualizing within cluster distances
    This was more of a script than an actual PyIgClassify module, and it hasn't been refactored to reflect its inclusion
    """
    def __init__(self, ab_db, cdr_directory, res_cutoff=2.8, rfac_cutoff = .3):
        self.ab_db = ab_db
        self.db = sqlite3.connect(ab_db)

        self.cdr_directory = cdr_directory
        self.pymol_writer = PyMolScriptWriter(os.getcwd())
        self.ab_structure = Structure.Antibody_Structure()
        self.base_dir = os.path.split(os.path.abspath(__file__))[0]


        self.cluster_data = ClusterData.CDRData(self.db, cdr_directory, True)
        self.cluster_data.set_res_cutoff(res_cutoff)
        self.cluster_data.set_rfac_cutoff(rfac_cutoff)

        self.cluster_data.load_data()

    def _make_directories(self, align_to_overhang, name):
        outdir = "cdr_visualizations"
        if not os.path.exists(outdir):
            os.mkdir(outdir)

        outdir = outdir+"/"+name
        if not os.path.exists(outdir):
            os.mkdir(outdir)

        if align_to_overhang:
            outdir = outdir+"/"+"align_overhang_bb"

        else:
            outdir = outdir+"/"+"align_all_bb"

        if not os.path.exists(outdir):
            os.mkdir(outdir)

        return outdir

    def group_names_to_degree(self, cdr, data, step_size):
        """
        1) Group the Data by DistDegree via steps.
        2) Add the commands to the current PyMolWriter class
        3) Color the groups

        Return data used to create
        """
        degrees = []
        names_to_deg = defaultdict()
        groups = defaultdict()
        groups_array = []

        for row in data:
            #print row
            distDegree = float(row[2])
            name = self.get_cdr_name(cdr, str(row[0]), str(row[1]))
            degrees.append(distDegree)
            names_to_deg[name] = distDegree

        min_deg = min(degrees)
        max_deg = max(degrees)

        ### Create Groups ###
        left = 0
        right = step_size
        while left <= max_deg:
            group = (left, right)
            groups[group] = []
            groups_array.append(group)

            left = left + step_size
            right = right + step_size

        ### Group the data ###
        for name in names_to_deg:
            degree = names_to_deg[name]

            for group in groups:
                left = group[0]
                right = group[1]

                if (degree >= left) and (degree < right):
                    groups[group].append(name)
                    break

        #print repr(groups)
        groups_class = DegreeGroups(groups_array, names_to_deg, groups)
        return groups_class

    def get_cdr_name(self, cdr, pdb, original_chain):
        name = pdb.lower()+original_chain.upper()+"_"+cdr
        return name

    def get_cdr_path(self, cdr, pdb, original_chain):
        """
        Get the CDR path using the cdr directory
        """
        name = self.get_cdr_name(cdr, pdb, original_chain)
        path = self.cdr_directory+"/"+name+".pdb"
        return path

    def get_all_paths(self, cdr, data):
        """
        Get all paths.  Data is tuple of pdb, original_chain
        """

        paths = []
        for row in data:
            path = self.get_cdr_path(cdr, row[0], row[1])
            paths.append(path)
        return paths

    def get_center_path_and_name(self, cdr, length, cluster):
        """
        Get a tuple of the path to the center cluster member and the name
        """

        data_types = ["PDB", "original_chain", "DistDegree"]
        data = AbDbFunctions.get_center_for_cluster_and_length(self.db, cdr, length, cluster, data_types)
        if not data:
            #sys.exit("No center for cluster: "+cluster)
            return False

        path = self.get_cdr_path(cdr, data[0], data[1])

        name = "center_"+data[0].lower()+data[1].upper()+"_"+cdr
        return (path, name)

    def visualize_outliers(self, align_to_overhang, alignment_path, use_pyrosetta_alignment= True, skip_empty_groups = True, step_size = 10.0, ):


        outdir = self._make_directories(align_to_overhang, "dihedral_outliers")

        outdir_sessions = self.base_dir+"/"+outdir+"/sessions"
        outdir_scripts = self.base_dir+"/"+outdir+"/scripts"

        pdb_path = alignment_path
        if not os.path.exists(pdb_path):
            os.mkdir(pdb_path)


        if not os.path.exists(outdir_sessions):
            os.mkdir(outdir_sessions)
        if not os.path.exists(outdir_scripts):
            os.mkdir(outdir_scripts)


        self.pymol_writer.set_outdir(outdir_scripts)
        self.pymol_writer.reset_script()


        for cdr in self.ab_structure.get_CDRs():
            for length in self.cluster_data.get_lengths(cdr.name):
                for cluster in self.cluster_data.get_clusters(cdr.name, length):
                    self.pymol_writer.reset_script()

                    names = [cdr, str(length), cluster]
                    outname = "dist_groups_"+cluster

                    data = self.cluster_data.get_cluster_data(cdr.name, length, cluster)
                    if isinstance(data, ClusterData.CDRClusterData): pass

                    paths = data.get_cdr_paths()

                    paths2 = []
                    for path in paths:
                        if os.path.exists(path):
                            paths2.append(path)
                    paths = paths2

                    ########## Build PyMol Script #############
                    center_name = data.get_center_name()

                    ### Load and Align BB Only###
                    if not use_pyrosetta_alignment:
                        if data.has_center_data():
                            self.pymol_writer.add_load_pdbs([data.get_center_path()], data.get_center_name())
                            self.pymol_writer.add_color(data.get_center_name(), "purple")
                            self.pymol_writer.add_load_pdbs(paths)
                            sele1 = ""
                            sele2 = ""

                            if align_to_overhang:
                                sele1 = self.get_overhang_sele(cdr, 3)
                                sele2 = sele1

                            self.pymol_writer.add_align_all_to(data.get_center_name(), sele1, sele2, True)

                        else:
                            self.pymol_writer.add_load_pdbs(paths)
                            sele1 = ""
                            sele2 = ""


                            if align_to_overhang:
                                sele1 = self.get_overhang_sele(cdr, 3)
                                sele2 = sele1

                            self.pymol_writer.add_align_all(sele1, sele2, True)

                    else:
                        if data.has_center_data():
                            self.pymol_writer.add_load_pdbs([data.get_center_path()], data.get_center_name())
                            self.pymol_writer.add_color(data.get_center_name(), "purple")

                        new_paths = []
                        for info in data.get_infos():
                            if isinstance(info, ClusterData.Info): pass
                            if info.name == data.get_center_name():
                                continue

                            new_path = pdb_path+"/"+info.name+".pdb"
                            if os.path.exists(new_path):
                                #print "Exists in : "+new_path
                                new_paths.append(new_path)
                            else:
                                #print "Could not find: "+new_path
                                if data.has_center_data():
                                    success = AnalyzeOutliers.attempt_align_and_save(cdr, info.path, data.get_center_path(), pdb_path, align_to_overhang)
                                else:
                                    print "Aligning to first cluster member as no center member is found:"
                                    success = AnalyzeOutliers.attempt_align_and_save(cdr, info.path, paths[0], pdb_path, align_to_overhang)
                                if success:
                                    new_paths.append(new_path)

                        self.pymol_writer.add_load_pdbs(new_paths)


                    ### Group ###
                    groups_class = self.group_names_to_degree(cdr.name, data.get_data(), step_size)

                    groups = groups_class.get_groups()
                    for group in groups_class.get_group_list():
                        if (skip_empty_groups and len(groups[group]) == 0):
                            continue

                        self.pymol_writer.add_group_objects(groups[group], groups_class.get_str_group_name(group))

                    ### Color ###

                    i = 0
                    for group in groups_class.get_group_list():

                        name = groups_class.get_str_group_name(group)

                        if len(groups[group]) == 0:
                            continue

                        color = self.pymol_writer.colors[i]
                        if color == "purple":
                            i+=2
                            color = self.pymol_writer.colors[i]


                        self.pymol_writer.add_color(name, color)

                        i+=2

                    ### Finish ###
                    self.pymol_writer.add_show("cartoon")
                    self.pymol_writer.add_line("center")
                    self.pymol_writer.add_save_session(outdir_sessions+"/"+outname+".pse")

                    ### Write out Script ###
                    fname = outname+".pml"
                    self.pymol_writer.write_script(fname)
                    self.pymol_writer.reset_script()

                    ### Run the script to get the session ###
                    run_pymol_script(outdir_scripts+"/"+fname)



    def visualize_clusters_per_length(self, align_to_overhang, skip_empty_groups = True):
        pass
        """
        outdir = self._make_directories(align_to_overhang, "all_clusters")

        outdir_sessions = self.base_dir+"/"+outdir+"/sessions"
        outdir_scripts = self.base_dir+"/"+outdir+"/scripts"




        if not os.path.exists(outdir_sessions):
            os.mkdir(outdir_sessions)
        if not os.path.exists(outdir_scripts):
            os.mkdir(outdir_scripts)


        self.pymol_writer.set_outdir(outdir_scripts)
        self.pymol_writer.reset_script()


        cdrs = ["L1", "L2", "L3", "H1", "H2", "H3"]

        for cdr in cdrs:
            lengths = AbDbFunctions.get_all_lengths(self.db, cdr, True)


            for length in lengths:
                clusters = AbDbFunctions.get_all_clusters_for_length(self.db, cdr, length, True)

                groups = defaultdict()
                outname = "cluster_groups_"+cdr+"_"+repr(length)
                self.pymol_writer.reset_script()
                for cluster in clusters:




                    names = [cdr, str(length), cluster]


                    data_types = ["PDB", "original_chain", "DistDegree"]

                    data = AbDbFunctions.get_data_for_cluster_and_length(self.db, cdr, length, cluster, data_types, True)

                    paths = self.get_all_paths(cdr, data)

                    groups[cluster] = paths


                    ########## Build PyMol Script #############
                    center_path_name = self.get_center_path_and_name(cdr, length, cluster)

                    ### Load and Align BB Only###
                    if center_path_name:
                        self.pymol_writer.add_load_pdbs([center_path_name[0]], center_path_name[1])
                        self.pymol_writer.add_color(center_path_name[1], "purple")
                        self.pymol_writer.add_load_pdbs(groups[cluster])
                        sele1 = ""
                        sele2 = ""

                        if align_to_overhang:
                            sele1 = self.get_overhang_sele(cdr, 3)
                            sele2 = sele1

                        self.pymol_writer.add_align_all_to(center_path_name[1], sele1, sele2, True)
                        self.pymol_writer.add_group_objects([center_path_name[0]], cluster)

                    else:
                        self.pymol_writer.add_load_pdbs(groups[cluster])
                        sele1 = ""
                        sele2 = ""

                        if align_to_overhang:
                            sele1 = self.get_overhang_sele(cdr, 3)
                            sele2 = sele1

                        self.pymol_writer.add_align_all(sele1, sele2, True)


                ### Group ###

                sorted_clusters = sorted(groups.keys())
                for group in sorted_clusters:
                    if (skip_empty_groups and len(groups[group]) == 0):
                        continue

                    self.pymol_writer.add_group_objects(self.get_names_from_paths(groups[group]), group)

                ### Color ###

                i = 0
                for group in sorted_clusters:

                    name = group

                    if len(groups[group]) == 0:
                        continue

                    color = self.pymol_writer.colors[i]
                    if color == "purple":
                        i+=2
                        color = self.pymol_writer.colors[i]


                    self.pymol_writer.add_color(name, color)

                    i+=2

                ### Finish ###
                self.pymol_writer.add_show("cartoon")
                self.pymol_writer.add_line("center")
                self.pymol_writer.add_save_session(outdir_sessions+"/"+outname+".pse")

                ### Write out Script ###
                fname = outname+".pml"
                self.pymol_writer.write_script(fname)
                self.pymol_writer.reset_script()

                ### Run the script to get the session ###
                run_pymol_script(outdir_scripts+"/"+fname)

            """


    def visualize_species(self, align_to_overhang):
        pass

    def get_names_from_paths(self, paths):
        names = []
        for path in paths:
            name = os.path.basename(path)
            name = "".join(name.split(".")[0:-1])
            names.append(name)
        return names

    def get_overhang_sele(self, cdr, overhang):
        """
        Gets the selection for terminal superposition in PyMol
        """
        start = cdr.get_pdb_start() - 1
        end = cdr.get_pdb_end() + 1

        chain = cdr.get_pdb_chain()

        start_sele = start - overhang
        end_sele = end + overhang

        resi = []
        resi.append((start_sele, start))
        resi.append((end, end_sele))

        sele = self.pymol_writer.get_sele(chain, resi)
        print sele
        return sele

class DegreeGroups:
    """
    Very simple class to hold the groups and explain what each member is.
    """
    def __init__(self, group_array, names_to_degrees, groups):
        self.group_array = group_array
        self.names_to_degrees = names_to_degrees
        self.groups = groups

    def get_group_list(self):
        """
        Get the list of all groups.  Groups are tuple: (Left, Right)
        """
        return self.group_array

    def get_names_to_deg(self):
        """
        Get a map of CDR names to their distance
        """
        return self.names_to_degrees

    def get_groups(self):
        """
        Get the main group map.  Each key is a tuple of (left, right) according to step size.
        Each value is a list of CDR names belonging to that group.
        """
        return self.groups

    def get_str_group_name(self, tuple_name):
        return repr(tuple_name[0]) +"-"+repr(tuple_name[1])

    def get_group_names(self):
        return self.group_names

    ####################################################################################################################
    ##                                                  OPTIONS
    ####################################################################################################################

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
                      default="/home/jadolfbr/Documents/modeling/databases/antibody_databases/PyIgClassify/DBOUT/website/antibody_database_rosetta_design.db")

    parser.add_option("--overhang",
                      help = "The expected number of final decoys",
                      default = 3)

    parser.add_option("--step_size",
                      help = "Step size used for outlier spits",
                      default = 5)

    parser.add_option("--show_empty_groups",
                      help = "Show the empty groups for outliers",
                      default = False,
                      action = "store_true")

    parser.add_option("--run_outlier_visualizations",
                      default = False,
                      action = "store_true")

    parser.add_option("--run_cluster_per_length_visualizations",
                      help = "Show extra visualizations - IE - groups of CDRs, clusters, lengths, etc.",
                      default = False,
                      action = "store_true")


    ############################
    ## Cutoffs
    ############################
    parser.add_option("--res_cutoff",
                      default = 2.8)

    parser.add_option("--rfac_cutoff",
                      default = .3)


    ############################
    ## Pre-aligned Structures
    ############################
    parser.add_option("--aligned_all_path",
                      default = os.path.split(os.path.abspath(__file__))[0]+"/aligned_cdrs_all")

    parser.add_option("--aligned_stem_path",
                      default = os.path.split(os.path.abspath(__file__))[0]+"/aligned_cdrs_stem")

    parser.add_option("--disable_pyrosetta_align",
                      default = False,
                      action = "store_true")
    (options, args) = parser.parse_args(sys.argv)


    visualizer = VisualizeCDRs(options.ab_db, options.cdr_dir, options.res_cutoff, options.rfac_cutoff)

    if options.run_outlier_visualizations:
        print "Running outlier vis"
        visualizer.visualize_outliers(True, options.aligned_stem_path, not bool(options.disable_pyrosetta_align), not bool(options.show_empty_groups), int(options.step_size))
        visualizer.visualize_outliers(False, options.aligned_all_path, not bool(options.disable_pyrosetta_align), not bool(options.show_empty_groups), int(options.step_size) )

    if options.run_cluster_per_length_visualizations:
        print "Running per-length vis"
        visualizer.visualize_clusters_per_length(True, not bool(options.show_empty_groups))
        visualizer.visualize_clusters_per_length(False, not bool(options.show_empty_groups))