#Author Jared Adolf-Bryfogle
from .HMMChainType import HMMChainType
from .HMMChainTypeSet import HMMChainTypeSet
from src.tools.path import *



class HMMChainTypeParser(object):
    """
    Class for parsing HMMChainType data from chain type file.  Creates a HMMChainTypeSet
    """
    
    def __init__(self):
        self.chain_type_path = get_db_path()+"/hmm_chain_types.txt"
        
        self.HMMchainset = HMMChainTypeSet()
        
        #Fix for some really fucked up error - where something in this class is global for no reason..
        if not self.HMMchainset.get_len():
            self._read_data()
        
    
    def get_HMMChainTypeSet(self):
        return self.HMMchainset
    
    def _read_data(self):
        FILE = open_file(self.chain_type_path, 'r')
        for line in FILE:
            if line[0]=="#":continue
            
            line = line.strip()
            lineSP = line.split()
            if len(lineSP)!= 6:continue
            HMMchain = HMMChainType(lineSP[0], lineSP[1], lineSP[2], lineSP[3], float(lineSP[4]), float(lineSP[5]))
            self.HMMchainset.add_chain_type(HMMchain)
            #print repr(self.HMMchainset.get_len())
        FILE.close()
            
