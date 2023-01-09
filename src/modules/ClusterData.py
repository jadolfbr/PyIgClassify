import os
import sqlite3
import sys
from collections import defaultdict

import Structure
from src.tools import AbDbFunctions


class Data:
    def __init__(self, ab_db, cdr_dir, limit_to_known = True):
        self.ab_db = ab_db
        self.cdr_dir = cdr_dir
        self.limit_to_known = limit_to_known
        self.ab_struct = Structure.Antibody_Structure()
        self.res_cutoff = 2.8
        self.rfac_cutoff = .3


    def set_res_cutoff(self, res_cutoff):
        self.res_cutoff = res_cutoff

    def set_rfac_cutoff(self, rfac_cutoff):
        self.rfac_cutoff = rfac_cutoff

    def load_data(self, extra_sele = ""):
        pass



#An attempt at organizing final cluster data by cdr, length, and clusters.

class CDRData(Data):
    def load_data(self, extra_sele = []):
        self.data = defaultdict(CDRLengths)

        for cdr in self.ab_struct.CDRS:
            lengths = CDRLengths(self.ab_db, self.cdr_dir, self.limit_to_known)
            lengths.load_data(cdr.name, extra_sele)
            self.data[cdr.name] = lengths

        self._populate_all_data()

    def _populate_all_data(self):
        self.all_data = defaultdict()

        for cdr_name in self.data:
            self.all_data[cdr_name] = defaultdict()
            print "ClusterData "+cdr_name+": loaded "+repr(len(self.data[cdr_name].get_lengths()))+" lengths"
            for length in self.data[cdr_name].get_lengths():
                self.all_data[cdr_name][length] = defaultdict()
                #print "ClusterData "+cdr_name+": loaded "+repr(len(self.data[cdr_name].get_clusters_data(length).get_clusters())+" clusters"
                for cluster in self.data[cdr_name].get_clusters_data(length).get_clusters():
                    self.all_data[cdr_name][length][cluster] = self.data[cdr_name].get_clusters_data(length).get_cluster_data(cluster)
                    #print "ClusterData "+cdr_name+": loaded "+repr(len(self.data[cdr_name][length][cluster]))+" cdrs"

    def get_all_data(self):
        return self.all_data

    def get_lengths(self, cdr_name):
        '''
        :param cdr_name: string
        :return: list
        '''
        return sorted(self.all_data[cdr_name].keys())

    def get_lengths_data(self, cdr_name):
        '''
        :param cdr_name: string
        :return: CDRLengths
        '''
        return self.data[cdr_name]

    def get_clusters(self, cdr_name, length):
        return sorted(self.all_data[cdr_name][length].keys())

    def get_clusters_data(self, cdr_name, length):
        '''
        :param cdr_name: string
        :param length: int
        :return: CDRClusters
        '''
        return self.data[cdr_name].get_clusters_data(length)

    def get_cluster_data(self, cdr_name, length, cluster):
        '''
        :param cdr_name: string
        :param length: int
        :param cluster: string
        :return: CDRClusterData
        '''

        return self.data[cdr_name].get_clusters_data(length).get_cluster_data(cluster)

class CDRLengths(Data):

    def load_data(self, cdr, extra_sele = []):
        self.cdr_name = cdr
        self.lengths = AbDbFunctions.get_all_lengths(self.ab_db, self.cdr_name, self.limit_to_known, self.res_cutoff, self.rfac_cutoff)
        self.data = defaultdict(CDRClusters)

        for length in self.lengths:
            #print repr(length)
            clusters = CDRClusters(self.ab_db, self.cdr_dir, self.limit_to_known)
            clusters.load_data(self.cdr_name, length, extra_sele)
            self.data[length] = clusters

    def get_lengths(self):
        return self.lengths

    def get_clusters_data(self, length):
        '''
        :param length: int
        :return: CDRCLusters
        '''
        return self.data[length]

class CDRClusters(Data):
    def load_data(self, cdr, length, extra_sele = []):
        self.cdr_name = cdr
        self.length = length
        self.clusters = AbDbFunctions.get_all_clusters_for_length(self.ab_db, self.cdr_name, self.length, self.limit_to_known, self.res_cutoff, self.rfac_cutoff)
        self.data = defaultdict(CDRClusterData)

        for cluster in self.clusters:
            #print cluster
            cluster_data = CDRClusterData(self.ab_db, self.cdr_dir, self.limit_to_known)
            cluster_data.load_data(self.cdr_name, self.length, cluster, extra_sele)
            self.data[cluster] = cluster_data

    def get_clusters(self):
        return self.clusters

    def get_cluster_data(self, cluster):
        '''
        :param cluster: string
        :return: CDRClusterData
        '''
        return self.data[cluster]



class CDRClusterData(Data):
    def load_data(self, cdr, length, cluster, extra_sele = []):
        self.cdr_name = cdr
        self.length = length
        self.cluster = cluster
        self.sele = ["PDB", "original_chain", "DistDegree"]

        self.extra_data = []

        self.center_data = AbDbFunctions.get_center_for_cluster_and_length(self.ab_db, self.cdr_name, self.length, self.cluster, self.sele)
        self.data = AbDbFunctions.get_data_for_cluster_and_length(self.ab_db, self.cdr_name, self.length, self.cluster, self.sele, self.limit_to_known, self.res_cutoff, self.rfac_cutoff)

        if extra_sele:
            self.extra_data = AbDbFunctions.get_data_for_cluster_and_length(self.ab_db, self.cdr_name, self.length, self.cluster, extra_sele, self.limit_to_known, self.res_cutoff, self.rfac_cutoff)

        self._populate_infos()

    def _populate_infos(self):
        self.infos = []
        for row in self.data:
            #print row
            path = self._get_cdr_path(row[0], row[1])
            name = self._get_cdr_name(row[0], row[1])
            info = Info(path, name, row[0], row[1], self.cluster, row[2])
            self.infos.append(info)

    def get_infos(self):
        return self.infos

    def get_data(self):
        return self.data

    def get_pdb(self):
        return self.data[0]

    def get_original_chain(self):
        return self.data[1]

    def get_extra_data(self):
        return self.extra_data

    def get_center_data(self):
        return self.center_data

    def has_center_data(self):
        '''
        :return: boolean
        '''
        if len(self.center_data) > 0:
            return True
        else:
            return False

    def get_center_name(self):
        if not self.center_data:
            #sys.exit("No center for cluster: "+cluster)
            return False

        name = self._get_cdr_name(self.center_data[0], self.center_data[1])
        return name

    def get_center_path(self):

        if not self.has_center_data():
            return False

        name = self.get_center_name()
        path = self.cdr_dir+"/"+name+".pdb"
        return path

    def get_cdr_names(self):
        names = []
        for row in self.data:
            name = self._get_cdr_name(row[0], row[1])
            names.append(name)
        return names

    def get_cdr_paths(self):
        """
        Get the CDR path using the cdr directory
        """
        paths = []
        for row in self.data:

            path = self._get_cdr_path(row[0], row[1])
            paths.append(path)

        return paths

    def _get_cdr_path(self, pdb, original_chain):
        name = self._get_cdr_name(pdb, original_chain)

        path = self.cdr_dir+"/"+name+".pdb"
        return path
    def _get_cdr_name(self, pdb, original_chain):
        return pdb.lower()+original_chain.upper()+"_"+self.cdr_name


class Info:
    def __init__(self, path, name, PDB, original_chain, cluster, dihedral_distance):
        self.path = path
        self.name = name
        self.PDB = PDB
        self.cluster = cluster
        self.original_chain = original_chain
        self.dihedral_distance = dihedral_distance


if __name__ == "__main__":

    ab_db ="/home/jadolfbr/Documents/modeling/databases/antibody_databases/PyIgClassify/DBOUT/website/antibody_database_redundant.db"
    cdr_dir = "/home/jadolfbr/Documents/modeling/databases/antibody_databases/PyIgClassify/DBOUT/cdr_pdbs_redun_by_cdr_overhang_3"
    if not os.path.exists(ab_db):
        sys.exit("ab_db does not exist!")
    if not os.path.exists(cdr_dir):
        sys.exit("cdr dir does not exist!")

    ClusData = CDRData(sqlite3.connect(ab_db), cdr_dir)
    ClusData.load_data()