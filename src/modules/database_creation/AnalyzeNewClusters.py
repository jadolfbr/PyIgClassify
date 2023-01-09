
#Jared Adolf-Bryfogle
import sqlite3
import sys
from collections import defaultdict

from src.tools.path import *


class AnalyzeNewClusters:
    """
    Class responsible for analyzing new cluster data in databases.
    CAN connect to multiple databases, and will add this data to the output - so you can compare all at once (NR, redundant, NR_by_cdr, NR_by_seq_cdrs)
    Will eventually run the reclustering on a fullcluster IF the number of sequences surpasses some threshold value.
    """
    
    def __init__(self, db_paths):
        self.connections = defaultdict()
        self.db_paths = db_paths
        for db_path in db_paths:
            if not os.path.exists(db_path):
                sys.exit(db_path+" does not exist!")
            self.connections[db_path] = sqlite3.connect(db_path)
            
        #Number of new fullclusters
        #Number of cdrs for each fullcluster
    
    def output_new_cluster_data(self, outpath):
        OUTFILE = open_file(outpath, 'w')
        OUTFILE.write("#DB new_clusters fullcluster num_entries\n")
        for db_path in self.connections:
            db_name = os.path.basename(db_path)
            db_name = "_".join(db_name.split(".")[0].split("_")[2:])
            new_clusters = self.get_new_fullclusters(self.connections[db_path])
            base = db_name+" "+repr(len(new_clusters))
            num_dict = self.get_all_num_for_fullclusters(new_clusters, self.connections[db_path])
            
            key_num_tup = sorted(num_dict.items(), key=lambda x: x[1], reverse=True)
            for pair in key_num_tup:
                #print pair
                out = base+" "+pair[0]+" "+repr(pair[1])
                OUTFILE.write(out+"\n")
        print "Basic new cluster data written"
        OUTFILE.close()
        
    def recluster_new_lengths(self):
        pass
    
    def recluster_all_lengths(self):
        pass

###################################
    
    def get_all_num_for_fullclusters(self, new_clusters, db):
        
        num_dict = defaultdict(int)
        for new_cluster in new_clusters:
            num = self.get_num_cdrs_for_cluster(db, new_cluster)
            num_dict[new_cluster] = num
        
        return num_dict
    
    def get_new_fullclusters(self, db):
        cur = db.cursor()
        fullclusters = []
        datatag1 = 'loopKeyNotInPaper'
        datatag2 = 'loopKeyNotInPaper_clustered'
        for row in cur.execute("SELECT distinct fullcluster FROM cdr_data WHERE datatag=? or datatag=?", [datatag1, datatag2]):
            fullclusters.append(row[0])
        cur.close()
        return fullclusters
    
    def get_num_cdrs_for_cluster(self, db, cluster):
        cur = db.cursor()
        num = 0
        for row in cur.execute("SELECT * FROM cdr_data WHERE fullcluster=?", [cluster]):
            num+=1
        cur.close()
        return num