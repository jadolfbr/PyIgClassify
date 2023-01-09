import math
import util

class CDRInfo(object):
    """
    A simple object for holding and accessing CDR information.
    """
    def __init__(self, cdr_name, start, end, gene="unk", cluster = "unk", distance=None):
        self.cdr_type = cdr_name

        # I know this is really stupid and needs to be refactored.  It basically existed as a grab bag of data for output at some point.
        self.new_start = start
        self.old_start = start

        self.cdr_start = start; #Actual CDR start
        self.cdr_end = end

        self.old_end = end

        self.old_chain = cdr_name[0]
        self.new_chain = cdr_name[0]

        self.cluster = cluster
        self.gene = gene

        if distance:
            self.set_distance(distance)

    def __str__(self):
        line = "CDR "+self.cdr_type +"\t"+ repr(self.cdr_start)+"\t"+repr(self.cdr_end)+"\t"+repr(self.cluster)
        return line
    
    def get_pdb_output(self, sequence):
        line = "CDR "+self.cdr_type+"\t"+repr(self.old_start)+"\t"+self.old_chain+"\t"+repr(self.new_start)+"\t"+self.new_chain+"\t"+ \
        self.get_sequence(sequence)+"\t"+self.cluster+"\t%.4f"%self.distance+"\t%.4f"%(self.norm_dis)+"\t%.2f"%(self.norm_dis_deg)
        return line
    
    def set_new_chain(self, chain):
        self.new_chain = chain
        
    def get_new_chain(self):
        return self.new_chain
    
    def get_gene(self):
        return self.gene
    
    def get_type(self):
        return self.cdr_type
    
    def get_start(self):
        """
        Starting residue of the CDR in chain.  1 - n
        """
        return self.cdr_start
    
    def get_end(self):
        """
        Ending residue of the CDR in chain. 1 - n
        """
        return self.cdr_end
    
    def get_length(self):

        #print self.cdr_type+" "+repr(self.cdr_end)+" : "+repr(self.cdr_start)
        l = self.cdr_end-self.cdr_start + 1
        #print "New Length: "+repr(l)
        return l
    
    def set_original_cdr_info(self, start, end, chain):
        """
        Set the original info for a CDR with a PDB.  Used after renumbering.
        """
        
        self.old_start = start
        self.old_end = end
        self.old_chain = chain
    
    def get_original_start(self):
        return self.old_start
    
    def get_original_end(self):
        return self.old_end
    
    def get_original_chain(self):
        return self.old_chain
    
    ##################
    #Sequence
    #  
    def get_sequence(self, sequence):
        """
        Get the sequence from a sequence string.
        """
        seq = sequence[(self.cdr_start-1):self.cdr_end]
        return seq
        
    ##################
    #Clusters
    #  
    
    def set_cluster(self, cluster):
        self.cluster = cluster
    
    def set_distance(self, distance):
        if distance != 1000:
            self.distance = distance
            self.norm_dis = util.get_norm_distance(self.get_length(), distance)
            self.norm_dis_deg = util.get_norm_distance_deg(self.norm_dis)
        else:
            self.distance = distance
            self.norm_dis = distance
            self.norm_dis_deg = distance

    def get_cluster(self):
        return self.cluster
    
    def get_distance(self):
        return self.distance
    