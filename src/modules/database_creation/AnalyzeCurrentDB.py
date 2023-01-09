
#Jared Adolf-Bryfogle

import glob
import sqlite3
import sys
import time
from collections import defaultdict

from src.tools import AbDbFunctions
from src.tools import outliers
from weblogolib import *

from src.modules import SequenceStats
from src.tools.path import *


class AnalyzeCurrentDB:
    """
    Run on a NR set of CDRs with the design_table created - antibody_database_rosetta_design
    Qifang already will have a summary table of each cluster.  So this is done.  Here, we are just getting data for Qifang to add.
    """
    
    def __init__(self, outdir, dbpath, DBOUT_dir, table='cdr_data'):
        self.outdir = outdir
        if not os.path.exists(dbpath):
            sys.exit(dbpath+ " does not exist!")
        if not os.path.exists(self.outdir):
            os.mkdir(self.outdir)

        self.dbout_dir = DBOUT_dir
        self.unknown_loop = 'loopKeyNotInPaper'
        self.cdrs = ["L1", "L2", "L3", "H1", "H2", "H3"]

        self.db = sqlite3.connect(dbpath)
        self.table = table
        self.clusters = self._get_cluster_list(self.db)
        self.data = defaultdict(dict)

        self.outlier_definitions = ["conservative", "liberal"]
        self.outlier_definition = "conservative"
        self.include_outliers = False
        self.use_distinct_seq = True

    def set_include_outliers(self, include_outliers = False):
        """
        Include Outliers for Consensus data computation (seq, rama, seqLogos?

        :param include_outliers: boolean
        :return: None
        """
        self.include_outliers = False

    def set_outlier_definition(self, definition):
        """
        Definition to use if removing outliers for consensus data computation

        :param definition: string
        :return:
        """
        if definition not in self.outlier_definitions:
            sys.exit("Definition does not exist!!")
        self.outlier_definition = definition

    def set_use_distinct_seq(self, use_distinct_seq):
        """
        Use the DISTINCT keyword when analyzing CDR seq from probablity (Rama, logos, etc.)
        """
        self.use_distinct_seq(use_distinct_seq)

    def set_main_db(self, db_path):
        self.db = sqlite3.connect(db_path)

    def set_design_table(self, design_table):
        self.design_table = design_table

    def set_main_table(self, main_table):
        self.table = main_table

    ####################################################
    def _get_rama_data(self):
        for cluster in self.clusters:
            self.data[cluster]['rama'] = ""
            self.data[cluster]['rama'] = self.get_consensus_rama(self.db, cluster)
    
    def _get_consensus_data(self):
        print "Generating consensus sequence data for known clusters"
        for cluster in self.clusters:

            #result = self._get_aa_probs(cluster)
            #consensus = self._calculate_consensus(result)

            consensus = self._calculate_consensus_via_seq(cluster)
            self.data[cluster]['consensus'] = consensus
            
    def _output_data(self, outname):
        OUTFILE = open_file(self.outdir+"/"+outname, 'w')
        OUTFILE.write("#cluster consensus rama\n")
        for cluster in sorted(self.data):
            OUTFILE.write(cluster+" "+self.data[cluster]['consensus']+" "+self.data[cluster]['rama']+'\n')
        
        print "Current cluster data written"
        OUTFILE.close()
    
    def get_and_output_data(self, outname):
        self._get_rama_data()
        self._get_consensus_data()
        self._output_data(outname)

    def write_db_totals(self, db_paths, outname, limit_known = True, dih_cutoff = None):

        connections = self._get_connections(db_paths)

        OUTFILE = open_file(self.outdir+"/"+outname, 'w')


        OUTFILE.write("# limit_to_current_clusters total_entries total_lambda_entries total_kappa_entries total_chains total_unique_pdbs total_clusters db\n")
        for db_path in connections:

            db_name = os.path.basename(db_path)
            db_name = "_".join(db_name.split(".")[0].split("_")[2:])
            total_clusters = len(self._get_cluster_list(connections[db_path], limit_known))
            total_entries = self._get_total_cdrs(connections[db_path], limit_known)
            total_pdbs = self._get_total_unique_PDBs(connections[db_path], limit_known)
            total_chains = self._get_total_chains(connections[db_path], limit_known)
            total_lambda = self._get_total_gene(connections[db_path], "lambda", limit_known)
            total_kappa = self._get_total_gene(connections[db_path], "kappa", limit_known)

            data = [limit_known, total_entries, total_lambda, total_kappa, total_chains, total_pdbs, total_clusters]
            out_str = "\n"
            for d in data:
                out_str = out_str+" "+repr(d)
            out_str = out_str+"\t"+db_name
            OUTFILE.write(out_str+"\n")

            limit_known = False
            total_clusters = len(self._get_cluster_list(connections[db_path], limit_known))
            total_entries = self._get_total_cdrs(connections[db_path], limit_known)
            total_pdbs = self._get_total_unique_PDBs(connections[db_path], limit_known)
            total_chains = self._get_total_chains(connections[db_path], limit_known)
            total_lambda = self._get_total_gene(connections[db_path], "lambda", limit_known)
            total_kappa = self._get_total_gene(connections[db_path], "kappa", limit_known)

            data = [limit_known, total_entries, total_lambda, total_kappa, total_chains, total_pdbs, total_clusters]
            out_str = ""
            for d in data:
                out_str = out_str+" "+repr(d)
            out_str = out_str+"\t"+db_name
            OUTFILE.write(out_str+"\n")

        OUTFILE.close()


    def write_unique_seq_totals(self, db_paths, outname, limit_to_known = True, dih_cutoff=None):

        connections = self._get_connections(db_paths)
        OUTFILE = open_file(self.outdir+"/"+outname, 'w')
        head = "db\ttype\ttotal\tratio_to_cdr\tratio_to_length\tmedian\trama\tseq"
        OUTFILE.write(head+"\n")

        for db_path in db_paths:
            db_name = os.path.basename(db_path)
            db_name = "_".join(db_name.split(".")[0].split("_")[2:])
            base = db_name

            for cdr in self.cdrs:
                total_cdr = self._get_total_unique_entries_per_cdr(connections[db_path], cdr, limit_to_known, dih_cutoff)
                OUTFILE.write(base+"\t"+cdr+"\t"+repr(total_cdr)+"\n")

                lengths = self._get_length_list_for_cdr(connections[db_path], cdr, limit_to_known, dih_cutoff)
                for length in sorted(lengths):
                    total_length = self._get_total_unique_entries_per_length(connections[db_path], cdr, length, limit_to_known, dih_cutoff)
                    OUTFILE.write(base+"\t"+cdr+"_"+repr(length)+"\t"+repr(total_length)+"\t"+"%.2f"%((total_length/float(total_cdr))*100)+"\n")

                    clusters = self._get_cluster_list_for_cdr_and_length(connections[db_path], cdr, length, limit_to_known, dih_cutoff)
                    for cluster in sorted(clusters):
                        total_cluster = self._get_total_unique_entries_per_cluster(connections[db_path], cluster, limit_to_known, dih_cutoff)
                        rama = self.get_consensus_rama(connections[db_path], cluster, dih_cutoff)
                        seq = self._get_consensus_seq(connections[db_path], cluster)
                        median = self._get_median_pdb(connections[db_path], cluster)

                        ratio_to_cdr = total_cluster/float(total_cdr)
                        ratio_to_length = total_cluster/float(total_length)

                        #Temp:
                        OUTFILE.write(base+"\t"+cluster+"\t"+repr(total_cluster)+"\t"+"%.2f"%(ratio_to_cdr*100)+"\t"+"%.2f"%(ratio_to_length*100)+"\t"+median+"\t"+rama+"\t"+seq+"\n")

        OUTFILE.close()

    def write_cluster_totals(self, db_paths, outname):

        connections = self._get_connections(db_paths)
        OUTFILE = open_file(self.outdir+"/"+outname, 'w')

        for db_path in db_paths:
            db_name = os.path.basename(db_path)
            db_name = "_".join(db_name.split(".")[0].split("_")[2:])
            base = db_name
            clusters = self._get_cluster_totals(connections[db_path], True)
            key_num_tup = sorted(clusters.items(), key=lambda x: x[1], reverse=True)
            for pair in key_num_tup:
                #print pair
                out = base+" "+pair[0]+" "+repr(pair[1])
                OUTFILE.write(out+"\n")
        OUTFILE.close()

    def write_date_stamp_data_simple(self, db_paths, outname, dih_cutoff = None):
        ext = 'w'

        #If found, we append the file
        outpath = self.outdir+"/"+outname
        if os.path.exists(outpath):
            ext = 'a'

        OUTFILE = open_file(self.outdir+"/"+outname, ext)
        if ext == 'w':
            OUTFILE.write("#date PDBIDs renum_ab chains lambda kappa heavy cdrs cdrs_known cdrs_known_dih_cutoff db\n")
        else:
            OUTFILE.write("\n")

        connections = self._get_connections_by_name(db_paths)

        limit_known = False

        for db_name in ["nr_by_cdr", "redundant"]:


            total_entries = self._get_total_cdrs(connections[db_name], False)
            total_entries_pruned = self._get_total_cdrs(connections[db_name], True)
            total_entries_pruned_dih = self._get_total_cdrs(connections[db_name], True, 40.0)

            total_pdbs = self._get_total_unique_PDBs(connections[db_name], limit_known)
            total_chains = self._get_total_chains(connections[db_name], limit_known)
            total_lambda = self._get_total_gene_chains(connections[db_name], "lambda", limit_known)
            total_kappa = self._get_total_gene_chains(connections[db_name], "kappa", limit_known)
            total_heavy = self._get_total_gene_chains(connections[db_name], "heavy", limit_known)
            total_renumbered = self._get_renumbered_antibody_totals()
            unique_len_cis = len(self._get_cluster_list(connections[db_name], False))
            data = [total_pdbs, total_renumbered, total_chains, total_lambda, total_kappa, total_heavy, total_entries, total_entries_pruned, total_entries_pruned_dih]
            out_str = time.strftime("%m/%d/%Y")
            for d in data:
                out_str = out_str+"\t"+repr(d)
            out_str = out_str+"\t"+db_name
            OUTFILE.write(out_str+"\n")

        OUTFILE.close()

    def plot_date_stamp_data_simple(self, in_name, out_name, only_redundant=True):
        try:
            import matplotlib as mpl
        except ImportError:
            print "Could not import matplotlib. Skipping."

        if only_redundant:
            types = ["redundant", "nr_by_cdr"]
        else:
            types = ["redundant"]

        new_dir = self.outdir+"/plots"
        if not os.path.exists(new_dir):
            os.mkdir(new_dir)

        if not os.path.exists(in_name):
            sys.exit("Path does not exist:"+in_name)

        INFILE = open_file(in_name, 'r')

        data = defaultdict()

        #data_types = ["PDBIDs", "renum_ab chains", "lambda", "kappa" "heavy" "cdrs" "cdrs_known" "cdrs_known_dih_cutoff"]
        for line in INFILE:
            line = line.strip()
            if not line: continue
            if line[0] == "#":continue

            lineSP = line.split()
            if lineSP[-1] in types:
                data[lineSP[-1]] = defaultdict()

    def write_date_stamp_data_by_cdr(self, db_paths, outname, only_current = True):
        ext = 'w'

        #If found, we append the file
        outpath = self.outdir+"/"+outname
        if os.path.exists(outpath):
            ext = 'a'

        OUTFILE = open_file(self.outdir+"/"+outname, ext)
        if ext == 'w':
            OUTFILE.write("#date L1 L2 L3 H1 H2 H3 db\n")
        else:
            OUTFILE.write("\n")
        cdrs = ["L1", "L2", "L3", "H1", "H2", "H3"]

        connections = self._get_connections(db_paths)

        for db_path in sorted(connections.keys()):
            db_name = os.path.basename(db_path)
            db_name = "_".join(db_name.split(".")[0].split("_")[2:])
            total_entries = self._get_total_cdrs(connections[db_path], only_current)
            line = time.strftime("%m/%d/%Y")
            line = line +"\t"+repr(total_entries)
            for cdr in cdrs:
                t = self._get_total_entries_per_cdr(connections[db_path], cdr, only_current)
                line = line +"\t"+repr(t)
            line = line + "\t"+db_name
            OUTFILE.write(line+"\n")

        OUTFILE.close()

    def get_time_stamp_entries(self, outname, backlog = 15):

        INFILE = open_file(self.outdir+"/"+outname)

        lines = []
        for line in INFILE:
            #line = line.strip()
            #if not line: continue
            lines.append(line)

        s = ""

        back = 0
        if len(lines) <= backlog:
            back = 0
        else:
            back = len(lines)-backlog
        for i in range(len(lines)-1, back, -1):
            s = s+lines[i]
        INFILE.close()

        return s

    def run_weblogo(self):
        outdir = self.outdir+"/weblogos"
        if not os.path.exists(outdir):
            os.mkdir(outdir)
        
        print "Outputting new weblogos"
        
        for cluster in self.clusters:
            #print "Working on cluster "+cluster
            msa = outdir+"/"+cluster+"_msa.txt"
            self._output_msa_fasta(outdir, cluster)
            
            MSA_IN = open_file(msa)
            seqs = read_seq_data(MSA_IN) 
            data = LogoData.from_seqs(seqs)
            options = LogoOptions()
            
            #Weblogo Options
            options.logo_title = cluster
            #options.fineprint = datetime.datetime.now().strftime('%b-%d-%G')
            options.creator_text = "Dunbrack Lab - Antibody Database Team"
            options.show_fineprint = False
            options.resolution = 900
            options.number_interval = 1
            options.scale_width = True
            options.unit_name='probability'
            options.color_scheme=std_color_schemes["charge"]
            #options.color_scheme=std_color_schemes["hydrophobicity"]
            
            format = LogoFormat(data, options)
            LOGO_OUT = open_file(outdir+"/"+cluster+"_weblogo.png", 'w')
            png_formatter( data, format, LOGO_OUT)
            
            MSA_IN.close()
            LOGO_OUT.close()
        
        print "WebLogos Written"

####################################
    def _get_connections(self, db_paths):
        connections = defaultdict()
        db_paths = db_paths
        for db_path in db_paths:
            if not os.path.exists(db_path):
                sys.exit(db_path+" does not exist!")
            connections[db_path] = sqlite3.connect(db_path)
        return connections

    def _get_connections_by_name(self, db_paths):
        connections = defaultdict()
        for db_path in db_paths:
            if not os.path.exists(db_path):
                sys.exit(db_path+" does not exist!")
            db_name = os.path.basename(db_path)
            db_name = "_".join(db_name.split(".")[0].split("_")[2:])
            connections[db_name] = sqlite3.connect(db_path)

        return connections

    def _add_datatag_and_dist(self, query, sele, limit_to_known, dih_cutoff):

        if limit_to_known :
            if (len(sele) == 0):
                query = query+" WHERE datatag != ? "
            else:
                query = query+" AND datatag != ? "

            sele.append(self.unknown_loop)

        if dih_cutoff:
            if (len(sele) == 0):
                query = query +" WHERE DistDegree <= ?"
            else:
                query = query +" AND DistDegree <= ?"


            sele.append(dih_cutoff)
        return query, sele

####################################
    def _get_cluster_list(self, db, limit_to_known = True, dih_cutoff = None):
        """
        Get a list of all fullclusters
        """
        c = db.cursor()
        clusters = []
        sele = []
        query = "SELECT DISTINCT fullcluster from "+self.table

        query, sele = self._add_datatag_and_dist(query, sele, limit_to_known, dih_cutoff)

        for row in c.execute(query, sele):
                clusters.append(row[0])

        c.close()
        return clusters

    def _get_length_list_for_cdr(self, db, cdr, limit_to_known = True, dih_cutoff=None):
        c = db.cursor()
        lengths = []
        sele = [cdr]
        query = "SELECT DISTINCT length from "+self.table+" WHERE CDR = ?"

        query, sele = self._add_datatag_and_dist(query, sele, limit_to_known, dih_cutoff)

        for row in c.execute(query, sele):
                lengths.append(row[0])

        c.close()
        return lengths

    def _get_cluster_list_for_cdr(self, db, cdr, limit_to_known = True, dih_cutoff=None):

        c = db.cursor()
        clusters = []
        sele = [cdr]
        query = "SELECT DISTINCT fullcluster from "+self.table+" WHERE CDR=?"

        query, sele = self._add_datatag_and_dist(query, sele, limit_to_known, dih_cutoff)

        for row in c.execute(query, sele):
                clusters.append(row[0])

        c.close()
        return clusters

    def _get_cluster_list_for_cdr_and_length(self, db, cdr, length, limit_to_known=True, dih_cutoff=None):

        length = int(length)
        c = db.cursor()
        clusters = []
        sele = [cdr, length]
        query = "SELECT DISTINCT fullcluster from "+self.table+" WHERE CDR=? AND length=?"

        query, sele = self._add_datatag_and_dist(query, sele, limit_to_known, dih_cutoff)

        for row in c.execute(query, sele):
            clusters.append(row[0])

        c.close()
        return clusters

    def _get_renumbered_antibody_totals(self):
        files = glob.glob(self.dbout_dir+"/renumbered_pdbs/*.pdb")
        return len(files)

    def _get_cluster_totals(self, db, limit_to_known = True, dih_cutoff = None):
        """
        Get a dictionary of all fullclusters and their totals.
        """
        c = db.cursor()
        clusters = defaultdict()
        sele = []
        query = "SELECT fullcluster from "+self.table

        query, sele = self._add_datatag_and_dist(query, sele, limit_to_known, dih_cutoff)


        for row in c.execute(query, sele):
            if not row[0] in clusters:
                clusters[row[0]] = 1
            else:
                clusters[row[0]]+=1

        c.close()
        return clusters

    def _get_consensus_seq(self, db, cluster):

        c = db.cursor()
        sele = [cluster]
        cons =[]
        query = "SELECT ConsSeq from CDRClusterSum WHERE Loop = ?"

        for row in c.execute(query, sele):
            cons.append(row[0])

        c.close

        if cons[0]:
            return cons[0]
        else:
            return "..."

    def _get_median_pdb(self, db, cluster):

        c = db.cursor()
        sele = [cluster]
        cons =[]
        query = "SELECT MedianPDB from CDRClusterSum WHERE Loop = ?"

        for row in c.execute(query, sele):
            cons.append(row[0])

        c.close

        if cons[0]:
            return cons[0]
        else:
            return "..."

    def _get_cluster_totals_for_cdr(self, db, cdr, limit_to_known=True, dih_cutoff = None):

        c = db.cursor()
        clusters = defaultdict()
        sele = [cdr]
        query = "SELECT fullcluster from "+self.table+" WHERE CDR=?"

        query, sele = self._add_datatag_and_dist(query, sele, limit_to_known, dih_cutoff)

        for row in c.execute(query, sele):
            if not row[0] in clusters:
                clusters[row[0]] = 1
            else:
                clusters[row[0]]+=1

        c.close()
        return clusters

    def _get_total_cdrs(self, db, limit_to_known = True, dih_cutoff = None):

        c = db.cursor()
        clusters = []
        sele = []
        query = "SELECT fullcluster from "+self.table

        query, sele = self._add_datatag_and_dist(query, sele, limit_to_known, dih_cutoff)


        for row in c.execute(query, sele):
            clusters.append(row[0])

        c.close()

        return len(clusters)


    def _get_total_entries_per_cdr(self, db, cdr, limit_to_known = True, dih_cutoff = None):
        c = db.cursor()
        clusters = []
        sele = [cdr]
        query = "SELECT fullcluster from "+self.table+" where CDR=?"

        query, sele = self._add_datatag_and_dist(query, sele, limit_to_known, dih_cutoff)

        for row in c.execute(query, sele):
            clusters.append(row[0])

        return len(clusters)

    def _get_total_unique_entries_per_cdr(self, db, cdr, limit_to_known=True, dih_cutoff = None):
        c = db.cursor()
        sequences = []
        sele = [cdr]
        query = "SELECT DISTINCT seq from "+self.table+" where CDR=?"
        query, sele = self._add_datatag_and_dist(query, sele, limit_to_known, dih_cutoff)


        for row in c.execute(query, sele):
            sequences.append(row[0])

        return len(sequences)

    def _get_total_unique_entries_per_length(self, db, cdr, length, limit_to_known=True, dih_cutoff = None):
        length = int(length)
        c = db.cursor()
        sequences = []
        sele = [cdr, length]
        query = "SELECT DISTINCT seq from "+self.table+" where CDR=? AND length=?"
        query, sele = self._add_datatag_and_dist(query, sele, limit_to_known, dih_cutoff)


        for row in c.execute(query, sele):
            sequences.append(row[0])

        return len(sequences)

    def _get_total_unique_entries_per_cluster(self, db, fullcluster, limit_to_known=True, dih_cutoff = None):
        c = db.cursor()
        sequences = []
        sele = [fullcluster]
        query = "SELECT DISTINCT seq from "+self.table+" where fullcluster=?"
        query, sele = self._add_datatag_and_dist(query, sele, limit_to_known, dih_cutoff)


        for row in c.execute(query, sele):
            sequences.append(row[0])

        return len(sequences)

    def _get_total_unique_PDBs(self, db, limit_to_known = True, dih_cutoff = None):

        c = db.cursor()
        clusters = []
        sele = []
        query = "SELECT DISTINCT PDB from "+self.table
        query, sele = self._add_datatag_and_dist(query, sele, limit_to_known, dih_cutoff)


        for row in c.execute(query, sele):
            clusters.append(row[0])

        c.close()
        return len(clusters)

    def _get_total_gene(self, db, gene, limit_to_known = True, dih_cutoff = None):

        c = db.cursor()
        clusters = []
        sele = [gene]
        query = "SELECT PDB from "+self.table + " where gene=?"
        query, sele = self._add_datatag_and_dist(query, sele, limit_to_known, dih_cutoff)

        for row in c.execute(query, sele):
            clusters.append(row[0])

        c.close()
        return len(clusters)

    def _get_total_chains(self, db, limit_to_known = True, dih_cutoff = None):
        c = db.cursor()
        chains = defaultdict()
        sele = []

        query = "SELECT PDB, original_chain from "+self.table
        query, sele = self._add_datatag_and_dist(query, sele, limit_to_known, dih_cutoff)

        for row in c.execute(query, sele):
            pdb_chain = row[0] + row[1]
            chains[pdb_chain] = 1

        c.close()
        return len(chains)

    def _get_total_gene_chains(self, db, gene, limit_to_known = True, dih_cutoff = None):
        c = db.cursor()
        chains = defaultdict()
        sele = [gene]

        query = "SELECT PDB, original_chain from "+self.table + " where gene=?"
        query, sele = self._add_datatag_and_dist(query, sele, limit_to_known, dih_cutoff)

        for row in c.execute(query, sele):
            pdb_chain = row[0] + row[1]
            chains[pdb_chain] = 1

        c.close()
        return len(chains)

    def get_consensus_rama(self, db, cluster, dih_cutoff = None):
        """
        Gets the most represented rama.  If there is a tie, then there is not preference.  
        """
        
        temp = defaultdict();
        cur = db.cursor()

        sele = [cluster]

        #Get DISTINCT sequences using the Group By command.  Ignore the tyring to use the better resolution or dih distance as only very
        #Small changes to rama will happen (Roland).

        #Determine whether to include outliers or not here

        query = "SELECT rama FROM "+self.table+" WHERE fullcluster=? " + outliers.get_outlier_string(self.include_outliers, self.outlier_definition) + " group by seq"


        for row in cur.execute(query, sele):
            if not temp.has_key(str(row[0])):
                temp[str(row[0])] = 0
            temp[str(row[0])]+=1


        max_rama= ""
        last_num = 0
        for rama in temp:
            if temp[rama] >= last_num:
                max_rama = rama
                last_num = temp[rama]

        return max_rama

    def _calculate_consensus_via_seq(self, cluster):
        sequences = AbDbFunctions.get_unique_sequences_for_cluster(self.db, cluster, self.include_outliers, self.outlier_definition)
        stats = SequenceStats.SequenceStats(sequences)
        consensus = stats.get_consensus_sequence()
        return consensus

    def _output_msa_fasta(self, outdir, cluster):
        sequences = AbDbFunctions.get_unique_sequences_for_cluster(self.db, cluster, self.include_outliers, self.outlier_definition)
        OUTFILE = open_file(outdir+"/"+cluster+"_msa.txt", 'w')

        n = 0
        for seq in sequences:
            n+=1
            OUTFILE.write("> "+cluster+" "+repr(n)+"\n")
            OUTFILE.write(seq+"\n")


        OUTFILE.close()


    """
    def _get_aa_probs(self, cluster):
        data = defaultdict(dict); #[position][aa] = [prob]
        cur = self.db.cursor()
        for row in cur.execute("SELECT position, aa, probability FROM "+self.design_table+" WHERE fullcluster=?", [cluster]):
            data[row[0]][row[1]] = row[2]
        cur.close()
        return data

    def _calculate_consensus(self, data):
        consensus = ""
        for position in range(0, len(data)):
            pos = position+1
            aa_letter = ""
            for aa in data[pos]:
                prob = data[pos][aa]
                if prob > .9:
                    aa_letter = aa.upper()
                    break
                elif prob > .2:
                    aa_letter = aa.lower()
                    break
                else:
                    aa_letter = '-'
                    continue

            consensus = consensus+aa_letter

        return consensus
    """

if __name__ == "__main__":

    #Getting data for paper:

    DBOUT = "/Users/jadolfbr/Documents/modeling/rosetta/projects/antibody_databases/PyIgClassify/DBOUT"
    WEBSITE = DBOUT+"/"+"website"
    log_dir = DBOUT+"/"+"logs"

    databases = [WEBSITE+"/antibody_database_redundant.db", WEBSITE+"/antibody_database_nr_by_cdr.db", WEBSITE+"/antibody_database_nr_by_cdrs_per_chain.db"]

    analyzer = AnalyzeCurrentDB(log_dir, WEBSITE+"/antibody_database_redundant.db")
    analyzer.write_unique_seq_totals(databases, "db_analysis_unique_sequences.txt")
    analyzer.write_unique_seq_totals(databases, "db_analysis_unique_sequences_cutoff_40_deg.txt", True, 40.0)
