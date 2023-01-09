from HMMFamilyType import HMMFamilyType
from HMMFamilyTypeSet import HMMFamilyTypeSet
from src.tools.path import *


class HMMFamilyTypeParser(object):
    """
    Class for parsing HMMFamilyType data from family type file. Creates a HMMFamilySet.
    """
    
    def __init__(self):
        self.chain_type_path = get_db_path()+"/hmm_family_types.txt"
        self.HMMfamilyset = HMMFamilyTypeSet()
        self._read_data()
        
    
    def get_HMMFamilyTypeSet(self):
        return self.HMMfamilyset
    
    def _read_data(self):
        FILE = open_file(self.chain_type_path, 'r')
        for line in FILE:
            if line[0]=="#":continue
            
            line = line.strip()
            lineSP = line.split()
            if len(lineSP)!= 4:continue
            
            HMMfamily = HMMFamilyType(lineSP[0], lineSP[1], float(lineSP[2]), float(lineSP[3]))
            self.HMMfamilyset.add_family_type(HMMfamily)
        
        FILE.close()
