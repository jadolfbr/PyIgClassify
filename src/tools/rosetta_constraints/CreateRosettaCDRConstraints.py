import math
import re
import sqlite3

from src.tools import AbDbFunctions

from src.modules.Structure import *
from src.tools.path import *


class CreateRosettaCDRConstraints:
    """
    Runs the associated R script on the SQLITE3 AntibodyDatabase and creates Rosetta Constraint files for each dihedral.
    Constraint means are now at the cluster Center to be more physically-relevant.
    Means and SDs produced use the circular R package for correct statistics.
    -- If NO center member is found in the database, will use means instead -- .
    """
    def __init__(self, outdir, constraint_type="CIRCULARHARMONIC"):
        self.constraint_type = constraint_type
        self.set_outdir(outdir)
        if not os.path.exists(self.outdir+"/PLOTS"):
            os.mkdir(self.outdir+"/PLOTS")
        
        self.type = type
        self.phi_atoms = ['C','N','CA','C']
        self.psi_atoms = ['N','CA','C','N']
        self.omega_atoms=['CA','C','N','CA']
        self.cluster_list = get_db_path()+"/TypeDef.txt"

        self.set_remove_liberal_outliers(True)
        self.center_data = defaultdict()
        self.use_center_as_mean = True

    def set_use_center_as_mean(self, use_center_as_mean):
        self.use_center_as_mean = use_center_as_mean

    def read_centers_set_dihedrals(self, database_path, mean_sd_txt):

        db = sqlite3.connect(database_path)

        RAW = open_file(mean_sd_txt, 'r')
        #cluster|type|length|n|position|mean_phi|sd_phi|mean_psi|sd_psi|mean_omega|sd_omega
        RAW.readline()

        antibody = Antibody_Structure()
        for line in RAW:
            lineSP = line.split()
            cluster = lineSP[0];
            if not re.search("-", cluster):
                print "Skipping line "+ line
                continue

            clusterSP = cluster.split("-")
            cdr = clusterSP[0]
            length = clusterSP[1]

            self.center_data[cluster] = AbDbFunctions.get_center_dih_degrees_for_cluster_and_length(db, cdr, int(length), cluster)
        RAW.close()

    def set_remove_liberal_outliers(self, remove_liberal_outliers):
        self.remove_liberal_outliers = remove_liberal_outliers

    def set_outdir(self, outdir):
        self.outdir = outdir
        if not os.path.exists(self.outdir):
            os.mkdir(self.outdir)

    def run_(self, database_path):

        self._create_mean_sd_output(database_path)
        self.read_centers_set_dihedrals(database_path, self.outdir+"/MEAN_SD.txt")

        print self.outdir+"/MEAN_SD.txt"

        self._create_rosetta_constraints(self.outdir+"/MEAN_SD.txt")


    def _create_mean_sd_output(self, database_path):
        if self.remove_liberal_outliers:
            os.system("Rscript "+os.path.abspath(os.path.dirname(__file__))+"/get_dihedral_constraints_outliers_false_liberal.R "+database_path+" "+self.cluster_list+" "+self.outdir)
        else:
            os.system("Rscript "+os.path.abspath(os.path.dirname(__file__))+"/get_dihedral_constraints.R "+database_path+" "+self.cluster_list+" "+self.outdir)
        print "Mean/SD data for Rosetta dihedral constraints created"


    def _create_rosetta_constraints(self, mean_sd_txt):


        RAW = open_file(mean_sd_txt, 'r')
        #cluster|type|length|n|position|mean_phi|sd_phi|mean_psi|sd_psi|mean_omega|sd_omega
        RAW.readline()
        
        antibody = Antibody_Structure()
        for line in RAW:
            lineSP = line.split()
            cluster = lineSP[0];
            if not re.search("-", cluster):
                print "Skipping line "+ line
                continue
            
            #print line
            type = lineSP[1]; length = int(lineSP[2]); n = int(lineSP[3]); position = int(lineSP[4]) #Position starts at 1!
            phi_mean = float(lineSP[5]); phi_sd = lineSP[6]; psi_mean = float(lineSP[7]); psi_sd = lineSP[8]; omega_mean = float(lineSP[9]); omega_sd = lineSP[10]

            #Skip constraint if N=1
            if n == 1:
                continue

            sds = [phi_sd, psi_sd, omega_sd]
            means = [phi_mean, psi_mean, omega_mean] #Use means if center member is unavailable.
            if not self.center_data[cluster]:
                print "No center data for cluster "+cluster+" using means instead!"
                centers = means
            elif not self.use_center_as_mean:
                centers = means
            else:
                centers = [ self.center_data[cluster]["phis"][position-1], self.center_data[cluster]["psis"][position-1], self.center_data[cluster]["omegas"][position-1] ]

            #Make Angles on -pi <-> pi scale
            for i in range(0, len(means)):
                if means[i] > 180.0:
                    means[i] = -360.0+means[i]
                
            #Fix SD of NA, which indicates N=1.  SD in this case would be 180.
            for i in range(0, len(sds)):
                if sds[i]=="NA":
                    sds[i]=float(180)
                else:
                    sds[i] = float(sds[i])
            
            cdr = cluster.split("-")[0]
            
            chain = cdr[0]
            out = self.outdir+'/'+cluster+'.txt'
            if os.path.exists(out):
                OUTFILE = open_file(out, 'a')
            else:
                OUTFILE = open_file(out, 'w')
            
            resnums = []
            
            #Not elegent, but too late for that.
            start_residues = int(math.ceil(float(length)/2.0))
            end_residues = int(math.floor(float(length)/2.0))
            
            start = int(antibody.get_CDR(cdr).Nter)
            for i in range(0, start_residues):
                resnums.append(start+i)
                
            start = int(antibody.get_CDR(cdr).Cter)-end_residues+1
            
            for i in range(0, end_residues):
                resnums.append(start+i)

            reschain = ' '+str(resnums[position-1])+chain+' '

            #First element will hit -1 index, which in Python is the last element. This indicates that its the first cst in
            #  the CDR, which goes from C-N of the residue number BEFORE the start of the CDR.  This is how we fix that.
            pre_reschain_index = position-2;
            if pre_reschain_index == -1:
                pre_reschain = ' '+str(resnums[0]-1)+chain+' '
            else:
                pre_reschain = ' '+str(resnums[pre_reschain_index])+chain+' '


            try:
                #Mid CDR - might be at break point.  Hence no simple addition, and the resnums list
                post_reschain = ' '+str(resnums[position-1+1])+chain+' '
            except IndexError:
                #End of CDR
                post_reschain = ' '+str(resnums[position-1]+1)+chain+' '
                
                
            line = "Dihedral "+self.phi_atoms[0]+pre_reschain+self.phi_atoms[1]+reschain+self.phi_atoms[2]+reschain+self.phi_atoms[3]+reschain+self.constraint_type+' '+"%.3f"%math.radians(centers[0])+' '+"%.3f"%math.radians(sds[0])+"\n"
            OUTFILE.write(line)
            line = "Dihedral "+self.psi_atoms[0]+reschain+self.psi_atoms[1]+reschain+self.psi_atoms[2]+reschain+self.psi_atoms[3]+post_reschain+self.constraint_type+' '+"%.3f"%math.radians(centers[1])+' '+"%.3f"%math.radians(sds[1])+"\n"
            OUTFILE.write(line)
            #line = "Dihedral "+self.omega_atoms[0]+pre_reschain+self.omega_atoms[1]+pre_reschain+self.omega_atoms[2]+post_reschain+self.omega_atoms[3]+post_reschain+self.constraint_type+' '+"%.3f"%math.radians(means[2])+' '+"%.3f"%math.radians(sds[2])+"\n"
            #OUTFILE.write(line)
            OUTFILE.close()
        
        RAW.close()