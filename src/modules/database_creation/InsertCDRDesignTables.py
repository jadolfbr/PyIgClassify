#Author: Jared Adolf-Bryfogle
#Creates a design table from existing cdr_data table.  Either in a new database or the one passed
import os
import sqlite3
import sys

p = os.path.dirname(os.path.abspath("../../"+__file__))
sys.path.append(p);

from src.tools.path import *

from src.modules.SequenceStats import *
from optparse import OptionParser


class InsertCDRDesignTables:
    """
    Class for creating sequence probability data table from CDR data table of a database.
    """
    
    def __init__(self, db_filename):
        self.db_filename = db_filename
        self.db = sqlite3.connect(db_filename)
        self.set_output_db(db_filename)
        self.aas = definitions().get_all_one_letter_codes()
        self.outlier_definitions = ["conservative", "liberal"]
        self.outlier_definition = "conservative"
        self.include_outliers = False

    def set_include_outliers(self, include_outliers = False):
        self.include_outliers = include_outliers

    def set_outlier_definition(self, definition):
        if definition not in self.outlier_definitions:
            sys.exit("Definition does not exist!!")
        self.outlier_definition = definition

    def set_use_distinct_seq(self, use_distinct_seq):
        self.use_distinct_seq = use_distinct_seq

    def set_output_db(self, out_db_name):
        """
        Sets the output database.  Will create the database if it does not exist.
        """
        
        self.outname = out_db_name
        if self.outname == self.db_filename:
            self.outdb = self.db

    def remove_design_table(self, name="cdr_residue_probabilities"):
        c = self.outdb.cursor()
        c.execute("DROP table IF EXISTS "+name)
        c.close()


    def create_design_table(self, name="cdr_residue_probabilities"):
        print "Creating the Design table"
        
        #Get unique clusters.
        c = self.outdb.cursor()
        c.execute("DROP table IF EXISTS "+name)
        c.execute("CREATE table "+name+"(id INT, fullcluster TEXT, length_type TEXT, position INT, aa TEXT, probability FLOAT, frequency INT, total_seq INT)")
        c.close()
        clusters = self.get_unique_clusters()
        
        #Added due to wierdness with database
        cluster_map = dict()
        for cluster in clusters:
            cluster_map[cluster] = self.get_type(cluster)
        
        with self.outdb:
            c = self.outdb.cursor()
            for cluster in clusters:
                #print cluster
                sequences = self.get_sequences(cluster)
                if len(sequences)==0:
                    continue
                stats = SequenceStats(sequences)
                id = 1
                
                for i in range(0, len(sequences[0])):
                    for aa in self.aas:
                        type = cluster_map[cluster]
                        c.execute("INSERT INTO "+name+" VALUES(?,?,?,?,?,?,?, ?)",(id, cluster,type , i+1, aa, stats.get_probability(i, aa), stats.get_frequency(i, aa), len(sequences)))
                        id+=1
            
    def print_data(self, num=20, prob_data_name = "cdr_residue_probabilities", outdir=""):
        """
        Outputs data for dealing with cutoff values.  Should output more pretty stats at some point.  Could also add it to the database to use in Rosetta.
        """
        
        total_out = open_file(outdir+"/seq_cutoff_totals.txt", 'w')
        total_out.write("#fullclusters not including loopKeyNotInPaper\n")
        total_out.write("#total_seq total_clusters total_h3 total-h3 total_1 total_2 total_3 \n")
        
        total_out.write("#fullclusters not including loopKeyNotInPaper\n")
        cum_out = open_file(outdir+"/seq_cutoff_cumulative_totals.txt", 'w')
        cum_out.write("#<=total_seq total_clusters total_h3 total-h3 total_1 total_2 total_3 \n")
        
        cum_total = cum_h3 = cum_1 = cum_2 = cum_3 = 0
        for i in range(1, num+1):
            total = self.get_num_clusters("select DISTINCT fullcluster from "+prob_data_name+" where total_seq = ?", [i])
            cum_total = cum_total + total
            
            total_h3 = self.get_num_clusters("select DISTINCT fullcluster from "+prob_data_name+" WHERE length_type=? and total_seq = ?",["H3-NA", i])
            cum_h3 = cum_h3 + total_h3
            
            total_1 =  self.get_num_clusters("select DISTINCT fullcluster from "+prob_data_name+" WHERE length_type=? and total_seq = ?",["1", i])
            cum_1 = cum_1+total_1
            
            total_2 = self.get_num_clusters("select DISTINCT fullcluster from "+prob_data_name+" WHERE length_type=? and total_seq = ?",["2", i])
            cum_2 = cum_2+total_2
            
            total_3 = self.get_num_clusters("select DISTINCT fullcluster from "+prob_data_name+" WHERE length_type=? and total_seq = ?",["3", i])
            cum_3 = cum_3+total_3
            
            total_out.write("\t".join([repr(i), repr(total), repr(total_h3), repr(total-total_h3), repr(total_1), repr(total_2), repr(total_3)]) +"\n")
            cum_out.write("\t".join([repr(i), repr(cum_total), repr(cum_h3), repr(cum_total-cum_h3), repr(cum_1), repr(cum_2), repr(cum_3)])+"\n")
        
        total_out.close()
        cum_out.close()
                
            
    ##Private Methods##
    def get_outlier_string(self, add_AND = True):
        addline = ""
        if not self.include_outliers:
            addline = src.tools.outliers.get_outlier_definition_string(self.outlier_definition)
            if add_AND:
                addline = " AND "+addline

        #print addline
        return addline

    def get_unique_clusters(self):
        c = self.db.cursor()
        clusters = []
        datatag = 'loopKeyNotInPaper'

        for row in c.execute("SELECT DISTINCT fullcluster from cdr_data WHERE datatag!=? ", [datatag]):
            clusters.append(row[0])

        c.close()
        return clusters
    
    def get_type(self, fullcluster):
        """
        Get the length type of the cluster
        """

        c = self.db.cursor()
        c.execute("SELECT DISTINCT length_type from cdr_data WHERE fullcluster=? ",[fullcluster])
        type = c.fetchone()[0]
        c.close()
        return type
        
    def get_sequences(self, fullcluster):
        """
        Get unique sequences.  If non-redundant by CDR should already be unique!
        """

        c = self.db.cursor()
        sequences = []

        distinct_str = ""
        if self.use_distinct_seq:
            distinct_str = "DISTINCT"

        for row in c.execute('SELECT '+distinct_str+' seq FROM cdr_data where fullcluster=? '+self.get_outlier_string(), [fullcluster]):
            sequences.append(row[0])
        c.close()
        return sequences
    
    def get_num_clusters(self, statement, l):
        c = self.outdb.cursor()
        c.execute(statement, l)
        num = len(c.fetchall())
        c.close()
        return num

if __name__ == '__main__':
    
    if os.path.exists("test.db"):
        os.remove("test.db")
        
    parser = OptionParser()
    args = sys.argv
    parser.add_option("--input_db", "-i",
        default = None,
        help = "Path to input database with cdr_data table"
    )
    
    parser.add_option("--output_db", "-o",
        default = "test.db",
        help = "Path to output database.  New one if it doesn't exist.  New table if it does.")
    
    parser.add_option("--table_name", "-t",
        default = "CDR_residue_probabilities",
        help = "Table name of sequence data")
    
    parser.add_option("--data_only", "-d",
        action="store_true",
        default = False,
        help = "Use this option to only print data from existing tables in the database. Will print more information is both tables are available.")
    
    parser.add_option("--data_dir", "-p",
        default = "RESULTS",
        help = "Directory to output the sequence stats.  This is used to keep track of cutoff and adjust in Rosetta as more data is added to the PDB.")
    
    (options, args) = parser.parse_args(args=args[1:])
    if not options.input_db:
        sys.exit("Input database needed to calculate sequence probabilities.")
    
    creator = CreateDesignData(options.input_db)
    creator.set_output_db(options.output_db)
    if not options.data_only:
        creator.create_design_table(options.table_name)
    creator.print_data(20, options.table_name, options.data_dir)
    
    print "Complete"
    
    
