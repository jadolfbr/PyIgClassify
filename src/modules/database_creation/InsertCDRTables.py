#from modules.Structure import *
#from modules.PythonPDB import *
import glob
import math
import re
from collections import defaultdict

from src.modules.SQLPDB import *
from src.tools import AbDbFunctions
from src.tools import fasta
from src.tools import general

from src.modules.outliers import AnalyzeOutliers as outlier_control
from src.tools.path import *


class InsertCDRTables:
    def __init__(self, inputDIR, outputDIR, renumDIR, new_dbfilename=False, new_clustertypefilename=False, pdbaa_path = 'database/pdbaa'):
        '''
        Main class responsible for creating the antibody databases including for Rosetta Design.
        Initialize.  Connect to database.  Set directories.  Make RESULTS folder.
        inputDIR should have all the files with specific filenames.  Yes, hardcoded for now.
        renumDIR should have all renumbered antibodies.  Would like both chains together.
        '''

        self.outputDIR = outputDIR
        self.outputLOG = outputDIR+"/logs"
        self.renumDIR = renumDIR
        self.inputDIR = inputDIR
        self.pdbaa_path = pdbaa_path

        if not os.path.exists(self.outputDIR):
            os.mkdir(self.outputDIR)
        if not os.path.exists(self.outputLOG):
            os.mkdir(self.outputLOG)
            
        #Connect to Database:
        if new_dbfilename:
            self.db = sqlite3.connect(self.outputDIR+'/'+new_dbfilename)
            self.dbPath = self.outputDIR+'/'+new_dbfilename
        else:
            self.db = sqlite3.connect(self.outputDIR+'/antibody_database')
            self.dbPath = self.outputDIR+'/antibody_database.db'
        
        #Load Cluster information
        if new_clustertypefilename:
            clusdata = self.inputDIR+'/'+new_clustertypefilename
        else:
            clusdata = get_db_path()+"/TypeDef.txt"
        FILE = open_file(clusdata, 'r')
        self.cluster_data = dict()
        for line in FILE:
            line = line.strip()
            lineSP = line.split()
            #print lineSP
            self.cluster_data[lineSP[0]]=lineSP[1]
        FILE.close()
        
        #Load Center Information
        self.center_data = general.load_centers()
        
        #Load FASTA information that has xray/resolution/Rfactor
        self.load_fasta_data()
        self.rmsd_db_path = self.outputDIR+"/aligned_cdrs/rmsd_data.db"
        if os.path.exists(self.rmsd_db_path):
            self.rmsd_db = sqlite3.connect(self.rmsd_db_path)

    def __exit__(self):
        self.db.close()
        
    def change_DB(self, filename, dir=False):
        '''
        Changes the database file being worked on.
        Can specify directory
        '''
        if not dir:
            self.db = sqlite3.connect(self.outputDIR+'/'+filename)
        else:
            self.db = sqlite3.connect(dir+'/'+filename)
        return
        
        
        """
        Old Style: path.get_db_path()+"/paper_supplemental/supplClusterData.txt"

        for line in FILE:
            line = line.strip()
            lineSP = line.split(',')
            pdb = lineSP[2][0:4]
            pdb_chain = lineSP[2][0:5]
            cdr = lineSP[0].split('-')[0]
            if not self.center_data.has_key(pdb_chain):
                self.center_data[pdb_chain]=dict()
            cluster = lineSP[0]+"-"+lineSP[1]
            #print pdb; print cdr
            self.center_data[pdb_chain][cdr]=cluster
        FILE.close()
        """
        
    def load_csv_data(self, paper_only=True, set_new_filename=False):
        '''
        Parses the newctoold.txt file from Ben.
        if paper_only = True, returns only results from the paper.
        Better written as a dictionary of tuples.  One for each column...
        '''
        if set_new_filename:
            path = self.inputDIR+'/'+set_new_filename
        else:
            path = self.inputDIR+'/'+'AbMatchToPaperClusters.csv'
            
            
        FILE = open_file(path, 'r')
        data = defaultdict()
        self.all_pdb_chain = dict()
        for line in FILE:
            #print line
            if line[0]=='#':continue
            line = line.strip()
            if not line:continue
            lineSP = line.split(',')
            #Array identifiers:
            pdbchain_cdr = lineSP[0]; seq=lineSP[1]; paper=lineSP[2]; loopkey = lineSP[3]; cluster = lineSP[4]
            
            #Fix for different versions of Ben's output Csv:
            clusterSP = cluster.split("_")
            if len(clusterSP)==2:
                cluster = clusterSP[1]
            elif cluster == '?':
                cluster = '-1'
                
            dis = lineSP[5]; normdis = lineSP[6]
            if normdis == "?":
                normdis = "-1"
            if dis == "?":
                dis = "-1"
            missing_seq = lineSP[7]
            if (paper_only==True and paper != "inPaper"):
                continue
            elif paper=="missingSeq":
                continue
            
            cdr = pdbchain_cdr.split('_')[1]
            pdb_chain = pdbchain_cdr[0:5].upper()
            original_chain = pdbchain_cdr[4]
            self.all_pdb_chain[pdb_chain]=''
            pdb = (pdbchain_cdr.split('_')[0])[0:4].upper()
            #key = pdb+'-'+cdr
            key = pdbchain_cdr
            #print key
            loopkeySP = loopkey.split('-')
            length = loopkeySP[1]; SS = loopkeySP[2]
            if not data.has_key(key):
                data[key]=dict()
            
            #Cluster Identification (Cis/Trans):
            if re.search("C", SS):
                fix = general.get_loop_key_from_ss(SS)
                cluster = fix+'-'+cluster
                
            
            data[key]['original_chain']=original_chain
            data[key]['datatag']=paper
            data[key]['cdr']=cdr
            data[key]['length']=length
            data[key]['cluster']=cluster
            
            #Makes full cluster --1 to -* if new.
            if cluster == '-1':
                fullcluster = cdr+'-'+length+'-'+'*'
            elif cluster[-3:len(cluster)] == '--1':
                fullcluster = cdr+'-'+length+'-'+cluster[:-3]+'-*'
            else:
                fullcluster = cdr+'-'+length+'-'+cluster
                
            data[key]['fullcluster']=fullcluster
            data[key]['seq']=seq
            #Find out how all the distances compare and what the hell they mean
            
            data[key]['center']=0
            
            #Check for Center cluster
            if self.center_data.has_key(pdb_chain):
                if self.center_data[pdb_chain].has_key(cdr):
                    data[key]['center']=1
            data[key]['normDis']=float(normdis)        
            data[key]['dis']=float(dis)
            data[key]['ss']=SS
            #print dis
            #print repr(data[key]['dis'])
        self.data = data
        FILE.close()
        
        
        #Makes sure multiple PDB's are not loaded into database.
        all_pdb_temp = dict()
        all_pdbs = []
        for pdb_cdr in sorted(self.data):
            pdb = pdb_cdr.split('-')[0]
            if not all_pdb_temp.has_key(pdb):
                all_pdb_temp[pdb]=1
                all_pdbs.append(pdb)
        self.all_pdbs = all_pdbs
        
        self.load_dihedrals()
        self.load_extra_db_data()
        
        return
    
    def load_fasta_data(self):
        """
        Loads data from fasta.  Replaces values NA with numbers.
        fasta_data is dictionary loaded from tools/fasta: pdb_chain: [method, residues, resolution, R factor]
        """
        self.fasta_data = fasta.read_header_data_from_fasta(self.pdbaa_path)
        
        resolution_dict = defaultdict()
        for pdb_chain in self.fasta_data:
            #Uses the best structure, but it seems it could be NMR data/etc.  Need to look into this.
            #If the FASTA is missing residue/resolution/R factor. This is bad.  We don't want this.
            method = self.fasta_data[pdb_chain][0]
            residues = int(self.fasta_data[pdb_chain][1]);
            resolution = self.fasta_data[pdb_chain][2];
            if resolution == 'NA': resolution = 100
            rfactor = self.fasta_data[pdb_chain][3];
            if rfactor == 'NA': rfactor = 100
            
            
            resolution_dict[pdb_chain] = [method, residues, resolution, rfactor]
        
        self.fasta_data = resolution_dict
        
    def load_extra_db_data(self, dbname='testing_pdbaa.db'):
        """
        Loads gene data from the testing database.  This database is output by default by PyIgClassify upon pdbaa updates.
        It has a lot of information, but also acts as a place for any extra data to go that we need for the resultant database.  In this case gene
        """
        con = sqlite3.connect(self.outputDIR+"/"+dbname)
        

        for pdb_cdr in self.data:
            
            pdb = pdb_cdr[0:4].upper()
            cdr = pdb_cdr.split("_")[1]
            
            cur = con.cursor()
            
            #This try block is only to account for if we are using a mismatched CSV with the test DB (while we wait for bens code to complete).  It will be removed.
            try:
                cur.execute("SELECT gene from cdr_data WHERE tag = ? and original_chain = ? and cdr = ?", [pdb, self.data[pdb_cdr]['original_chain'], cdr])
                gene = cur.fetchone()[0]
                self.data[pdb_cdr]['gene'] = gene
            except TypeError:
                self.data[pdb_cdr]['gene'] = "NA"
            cur.close()
        con.close()

    def load_dihedrals(self, new_file=False):
        """
        Loads dihedral data from Ben's dihedral output file.
        """
        if new_file:
            file = self.inputDIR+'/'+new_file
        else:
            file = get_db_path()+"/csv/dihedralOutputMatchToPaperClusters.csv"
        self.dihedrals = defaultdict()
        DIH = open_file(file, 'r')
        for line in DIH:
            line = line.strip()
            lineSP = line.split(',')
            pdb_cdr = lineSP[0]
            #pdb_cdr = pdb_chain_cdr[0:4]+pdb_chain_cdr[5:]
            dih_string = ''
            for dih in lineSP[1:]:
                dih_string = dih_string+dih+':'
            self.dihedrals[pdb_cdr]=dih_string
            #print "DIHEDRAL :"+ pdb_cdr
        #print self.dihedrals
        DIH.close()
        self.calculate_ramachandran()

    def update_rmsd_data(self):

        data = []
        #self.db = sqlite3.connect(self.dbPath)
        cur = self.db.cursor()
        for row in cur.execute("SELECT PDB, original_chain, CDR, length, fullcluster FROM cdr_data"):
            data.append(row)
        cur.close()

        with self.db:
            for d in data:

                pdb = str(d[0]); original_chain = str(d[1]); cdr = str((d[2])); length = d[3]; fullcluster = str(d[4])

                #print "UPDATE RMSD "+repr(d)
                bb_rmsd_cdr_align = AbDbFunctions.get_cdr_rmsd_for_entry(self.rmsd_db, d[0], d[1], d[2], d[3], d[4])
                bb_rmsd_stem_align = AbDbFunctions.get_stem_rmsd_for_entry(self.rmsd_db, d[0], d[1], d[2], d[3], d[4])

                columns = [bb_rmsd_cdr_align, cdr, length, pdb, original_chain, fullcluster]
                #sanity = [cdr, length, pdb, original_chain, fullcluster]


                #print repr(columns)

                cur = self.db.cursor()
                #for row in cur.execute("SELECT * from cdr_data WHERE CDR=? AND length=? AND PDB=? AND original_chain=? and fullcluster=?", sanity):
                    #print row

                cur.execute("UPDATE cdr_data SET bb_rmsd_cdr_align = ? AND CDR=? AND length=? AND PDB=? AND original_chain=? and fullcluster=?", columns)

                columns = [bb_rmsd_stem_align, cdr, length, pdb, original_chain, fullcluster]

                cur.execute("UPDATE cdr_data SET bb_rmsd_stem_align = ? AND CDR=? AND length=? AND PDB=? AND original_chain=? and fullcluster=?", columns)
                cur.close()

                #print self.dbPath

    def create_rmsd_data(self):
        """
        Align using PyRosetta to cluster center (Which is now always included in the databases)
          and calculate both the RMSDs for stem alignment and full CDR alignment.
        Adds two columns for both types of alignments.  Sets both at 1000 to start.
        Set to always generate alignments due to possible changing center cluster membership.
        """
        res_cutoff = 1000 #Attempt analysis of all databases.
        rfac_cutoff = 1000 #Attempt analysis of all databases.

        self.db.close()
        analyzer = outlier_control.AnalyzeOutliers(self.dbPath, self.outputDIR+"/cdr_pdbs_redun_by_cdr_overhang_3/", res_cutoff, rfac_cutoff)
        print "Analyzing RMSD"

        analyzer.analyze_rms(self.outputDIR+"/aligned_cdrs", False, 3, False, "rmsd_data.db")
        #analyzer.analyze_rms(path.get_DBOUT_path_full()+"/aligned_cdrs", True, 3, True ) #Temporarily do not realign to speed debugging

        self.db = sqlite3.connect(self.dbPath)
        self.rmsd_db = sqlite3.connect(self.rmsd_db_path)
        print "Creation of RMSD database complete"

    def calculate_dist_degree(self, normDis):
        """
        Calculate the normalized distance in degrees from the normalized distance
        """

        if normDis == -1.0:
            return -1.0

        dist_degree = math.acos(1.0 - ( float( normDis )/2.0)) * (180/ math.pi)
        return dist_degree

    def calculate_ramachandran(self):
        """
        Adds a 'rama' part to data from dihedrals.
        """
        print "Calculating ramachandran for each entry."
        if not self.dihedrals:
            print "No dihedral data loaded.  Returning."

        OUT = open_file(self.outputLOG+"/MISSING_DIHEDRALS.txt", 'w')
        for pdb_cdr in self.data:
            #print self.data[pdb_cdr]['seq']
            if not self.dihedrals.has_key(pdb_cdr):
                print pdb_cdr +" missing dihedrals!  Skipping..."
                OUT.write(pdb_cdr+"\n")
                continue

            dih = self.dihedrals[pdb_cdr].split(":")
            dih_len = len(dih)/3
            
            rama_string = ""
            
            ss = self.data[pdb_cdr]['ss']
            #print ss
            index = 0
            for pos in range(0, dih_len):
                
                type_string = general.get_rama_type([float(dih[index]), float(dih[index + 1]), float(dih[index + 2])])
                #print pdb_cdr+" "+repr(pos)+" "+type_string+" "+dih[index]+" "+dih[index+1]+" "+dih[index+2]+" "
                #Here, we use lower case for cis residues
                if ss[pos]=="C":
                    type_string = type_string.lower()
                rama_string = rama_string + type_string
                #print rama_string
                index+=3
            self.data[pdb_cdr]['rama'] = rama_string
            
                
            
    def create_data_table(self, update_rmsd = False):
        '''
        Creates the data table part of the database.  columns are:
        PDB|CDR|Length|Cluster|FullCluster|Type|Center|Sequence|Dis|SS
        Update RMSD is for recursion.  Using the UPDATE keyword for RMSD data is not working, so we are skipping it.
        '''
        #TEST = open(self.outputDIR+'/TESTING.txt', 'w')
        #TEST.write(repr(self.dihedrals))

        if update_rmsd:
            self.rmsd_db = sqlite3.connect(self.rmsd_db_path)

        with self.db:
            cur = self.db.cursor()
            cur.execute("DROP TABLE IF EXISTS cdr_data")
            cur.execute("CREATE TABLE cdr_data(id INT, datatag TEXT, PDB TEXT, original_chain TEXT, CDR TEXT, length INT, cluster TEXT, length_type TEXT, fullcluster TEXT, center INT, seq TEXT, dis FLOAT, normDis FLOAT, DistDegree FLOAT, bb_rmsd_cdr_align FLOAT, bb_rmsd_stem_align FLOAT, ss TEXT, rama TEXT, dihedrals TEXT, gene TEXT, species TEXT, method TEXT, resolution FLOAT, rfactor FLOAT)")
            #10 Fields
            id = 1
            for pdb_cdr in sorted(self.data):
                pdb = pdb_cdr.split('_')[0]
                pdb = pdb[0:4]
                fullcluster = self.data[pdb_cdr]['fullcluster']
                #print "Adding "+pdb_cdr +" to table"
                #if not self.cluster_data.has_key(fullcluster):
                
                #This WILL need refactoring once we have NEW clusters.  Should have a NEW datatag.  loopKeyNotPublished loopKeyNew?
                if self.data[pdb_cdr]['cdr']=="H3":
                    type = "H3-NA"
                elif self.data[pdb_cdr]['datatag']=='loopKeyNotInPaper':
                    #print "LoopKeyNotInPaper: "+fullcluster
                    type = 'NA'
                else:
                    type = self.cluster_data[fullcluster]
                
                pdb_chain = pdb.upper()+self.data[pdb_cdr]['original_chain']

                ## Fuck 1OCW.  It is solved incorrectly.  We don't need it or want it in the database.
                ## Until we use data from PDB Redo or re-build antibody structures from the density, leave this.
                if pdb.upper() == "1OCW": continue;

                #######################################################################################################

                if not self.dihedrals.has_key(pdb_cdr):
                    print pdb_cdr+ " missing dihedrals. Placing nothing for now"
                    self.dihedrals[pdb_cdr]= -1

                #Get RMSD Data
                if self.data[pdb_cdr]['center']:
                    bb_rmsd_cdr_align = 0.0
                    bb_rmsd_stem_align = 0.0
                elif update_rmsd:

                    bb_rmsd_cdr_align = AbDbFunctions.get_cdr_rmsd_for_entry(self.rmsd_db, pdb.upper(),
                                                                             self.data[pdb_cdr]['original_chain'], self.data[pdb_cdr]['cdr'], self.data[pdb_cdr]['length'],
                                                                             self.data[pdb_cdr]['fullcluster'])

                    bb_rmsd_stem_align = AbDbFunctions.get_stem_rmsd_for_entry(self.rmsd_db, pdb.upper(),
                                                                               self.data[pdb_cdr]['original_chain'], self.data[pdb_cdr]['cdr'], self.data[pdb_cdr]['length'],
                                                                               self.data[pdb_cdr]['fullcluster'])
                else:

                    bb_rmsd_cdr_align = -1
                    bb_rmsd_stem_align = -1

                gene = self.data[pdb_cdr]['gene']
                if gene == "lambda6":
                    gene = "lambda"

                try:
                    cur.execute("INSERT INTO cdr_data VALUES(?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)", \
                    (id,self.data[pdb_cdr]['datatag'],pdb.upper(), self.data[pdb_cdr]['original_chain'],self.data[pdb_cdr]['cdr'], \
                    self.data[pdb_cdr]['length'],self.data[pdb_cdr]['cluster'],type,self.data[pdb_cdr]['fullcluster'],\
                    self.data[pdb_cdr]['center'],self.data[pdb_cdr]['seq'],self.data[pdb_cdr]['dis'],self.data[pdb_cdr]['normDis'], self.calculate_dist_degree(self.data[pdb_cdr]['normDis']), bb_rmsd_cdr_align, bb_rmsd_stem_align, self.data[pdb_cdr]['ss'], self.data[pdb_cdr]['rama'], self.dihedrals[pdb_cdr], \
                    gene, "TBD", self.fasta_data[pdb_chain][0], self.fasta_data[pdb_chain][2], self.fasta_data[pdb_chain][3]))
                
                    id+=1
                    
                except KeyError:
                    print pdb_cdr +" skipped for missing a piece of needed data."
                    continue


        #End WithDB
        if update_rmsd:
            return

        if os.path.exists(self.rmsd_db_path):
            #print "Updating RMSD in cdr_data table from: "+self.rmsd_db_path
            self.create_data_table(True)
        else:
            #print "Creating RMSD data table from: "+self.rmsd_db_path
            self.create_rmsd_data()

            #print "Updating RMSD in cdr_data table from "+self.rmsd_db_path
            self.create_data_table(True)

    def create_coord_table(self, newFILE=False):
        '''
        Creates the CDR Coordinate Table.
        PDB|CDR|---pdb information---
        InewFILE is the new file name.  If false, uses the one already opened by closing it first, and reopening it at the end. 
        '''
        
        self.db.close()
        if not newFILE:
            file = self.dbPath
        else:
            file = self.outputDIR+'/'+newFILE
        
        MISSING = open_file(self.outputLOG+'/missing_renum_structures.txt', 'w')

        match = self.renumDIR+"/*.pdb"
        all_pdbs = glob.glob(match)

        PDBLIST = open_file(self.renumDIR+'/PDBLIST.txt', 'w')
        print self.pdbaa_path
        if self.pdbaa_path == "database/pdbaa_for_tests":
            print "Writing structures for tests"
            for i in range(0, 20):

                PDBLIST.write(all_pdbs[i]+"\n")
        else:
            print "Writing structures for all pdbs"
            for path in all_pdbs: PDBLIST.write(path + "\n")

        PDBLIST.close()
        
        PDBLIST = open_file(self.renumDIR+'/PDBLIST.txt', 'r')
        
        #Changing it to work with ScFvs where you don't have uniqe pdb_chains!
        chain_dict = defaultdict()
        chain_dict['single'] = defaultdict()
        chain_dict['light'] = defaultdict()
        chain_dict['heavy'] = defaultdict()
        
        chain_list = []
        
        for line in PDBLIST:
            #print line
            line = line.strip()
            renumbered_pdb_name = os.path.basename(line)
            nameSP = renumbered_pdb_name.split('-')
            if len(nameSP)==1:
                pdb_chain = nameSP[0].split('.')[0]
                chain_dict['single'][pdb_chain.lower()]=[line, False]
                chain_list.append(pdb_chain.lower())
                #print x[0].lower()
            else:
                pdb_chain1 = nameSP[0]
                pdb_chain2 = nameSP[1].split('.')[0]
                chain_dict['light'][pdb_chain1.lower()]=[line, 'L']
                #print "LH : "+x[0].lower()+' : '+x[1].lower()
                chain_dict['heavy'][pdb_chain2.lower()]=[line, 'H']
                chain_list.append(pdb_chain1.lower()); chain_list.append(pdb_chain2.lower())
                
        PDBLIST.close()
        i=1
        x=1
        for pdb_chain in sorted(self.all_pdb_chain):
            '''
            Keep opening and closing the database, this is not quite what I want, but ok.
            '''
            #Taking outstructID of 1 because we need the original chain....
            #print pdb_chain
            pdb = pdb_chain[0:4]
            chain = pdb_chain[4]
            #pdb_chain = pdb+chain
            p = SQLPDB(pdb, "x", chain, False, file)
            #We only want to load the specific chain.
            x+=1
            
            if pdb_chain.lower() in chain_list:
                #print "Loading Chain into coordinate DB: "+pdb_chain

                
                types = ['single', 'light', 'heavy']
                
                for chain_type in types:
                    if chain_dict[chain_type].has_key(pdb_chain.lower()):
                        p.read_pdb_into_database_flat(chain_dict[chain_type][pdb_chain.lower()][0], chain_dict[chain_type][pdb_chain.lower()][1])
                        i+=1
            else:
                print "ERROR: Cannot Find PDB chain - "+pdb_chain
                MISSING.write(pdb_chain+"\n")
                
        print repr(x)+ " original chains in database. "+repr(i)+" Chains Added to coord db.  More due to SCFVs."
        db_util = PDB_database(sqlite3.connect(file))
        table = db_util.scrub('pdb')
        db_util.query_all(table)
        db_util.update_modelID_CDRS(table)
        db_util.db.close()
          
                
        MISSING.close()
        
        
        self.db = sqlite3.connect(self.dbPath)    
            
        return
    
    def parse_cdrs_to_PDB(self, overhang=3, filePath=False, dir_base = "cdr_pdbs_nr_by_cdr_overhang_", copy_from_redun = True):
        """
        Connects to a database with Renumbered PDB structures (coord table) - created using create_coord_table (SQLPDB) and
         saves all CDRs to the output directory according to overhang.
        """

        redunDIR = self.outputDIR+"/cdr_pdbs_redun_by_cdr_overhang_"+repr(overhang)

        outDIR = self.outputDIR+"/"+dir_base+repr(overhang)
        if not os.path.exists(outDIR):
            os.mkdir(outDIR)
            
        if not filePath:
            db_util = PDB_database(self.db)
        else:
            db_util = PDB_database(sqlite3.connect(filePath))
        db_util.set_output_DIR(outDIR)
        db_util.set_output_occupancy_1(True)
        FILE = open_file(outDIR+'/'+'PDBLIST.txt', 'w')
        table = db_util.scrub("pdb")
        
        #Setup Antibody to help in parsing
        from src.modules.Structure import Antibody_Structure
        ab = Antibody_Structure()
        
        #Go through each pdb_cdr found, and output.
        for pdbchain_cdr in sorted(self.data):
            pdbchain = pdbchain_cdr[0:5]
            cdr = pdbchain_cdr.split("_")[1]
            pdb = pdbchain[0:4]
            original_chain = pdbchain[4]
            #if os.path.exists(self.renumDIR+'/'+pdb.lower()+'.pdb'):
            #print "Writing Out CDR: "+pdbchain_cdr

            filename = pdb.lower()+original_chain+"_"+cdr+".pdb"
            outpath = outDIR+"/"+filename
            redunpath = redunDIR+"/"+filename

            if copy_from_redun and os.path.exists(redunpath) and redunpath != outpath:
                os.system("cp "+redunpath+" "+outpath)
                FILE.write(outpath +"\n")
                continue

            db_util._reset_cursor()

                
            #Then, query the piece
            db_util.query_piece_pdbID_and_strucID(table, pdb.upper(), ab.CDR[cdr].Nter-overhang, ab.CDR[cdr].Cter+overhang, ab.CDR[cdr].chain, original_chain)
            
            #If the query of the chain was successful, we write to the PDBLIST file.

            if db_util.save_cur_as_pdb(outpath, False):
                FILE.write(outpath+"\n")

        FILE.close()    
        return
    
    def convert_pdbs_to_rosetta_db(self, jd2path, rosetta_database, overhang=3, database_name="", dir_base = "cdr_pdbs_nr_by_cdr_overhang_"):
        """
        Quick and dirty call to rosetta.
        """
        
        if not database_name:
            database_name = self.dbPath
        else:
            database_name = self.outputDIR+"/"+database_name
        
        pdbs = self.outputDIR+"/"+dir_base+repr(overhang)+"/PDBLIST.txt"
        call = jd2path+" -database "+rosetta_database+" -l "+pdbs+ " -out:use_database "+"-inout:dbms:database_name "+database_name+" -ignore_unrecognized_res -ignore_zero_occupancy false"
        print "RosettaDB command: "+call
        os.system(call)
        print "CDR Coordinates added to Rosetta database."
        
    