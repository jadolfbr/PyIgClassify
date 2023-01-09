from collections import defaultdict

from src.tools.path import *


class RenumberMapTypeParser(object):
    """
    Simple class for getting the correct numbering map.  See /database/hmm_numbering_maps.txt
    """
    
    def __init__(self):
        self.map_type_path = get_db_path()+"/hmm_numbering_maps.txt"
        self.numbering_types = defaultdict(dict)
        self._read_map_types()
        
    def _read_map_types(self):
        FILE = open_file(self.map_type_path, 'r')
        for line in FILE:
            if line[0] =='#': continue
            line = line.strip()
            lineSP = line.split()
            if len(lineSP) == 0: continue

            self.numbering_types[lineSP[0]]['scheme'] = lineSP[1]
            self.numbering_types[lineSP[0]]['type'] = lineSP[2]
    
    def get_numbering_types(self):
        """
        Get a dictionary of numbering_types (relative db filenames).
        """
        
        return self.numbering_types.keys()
    
    def get_full_map_path(self, numbering_scheme, numbering_type):
        """
        Get the full path of the numbering map file. type is kappa/lambda/etc.
        """
        for number_map in self.numbering_types:
            if self.numbering_types[number_map]['scheme'] == numbering_scheme and self.numbering_types[number_map]['type'] == numbering_type:
                return get_db_path()+"/"+number_map
    
    def get_all_numbering_schemes(self):
        """
        Get all of the numbering schemes included in hmm_numbering_types
        """
        v = dict()
        for number_map in self.numbering_types:
            v[self.numbering_types[number_map]['scheme']] = 0
        
        return v.keys()
    
    def get_all_numbering_types(self):
        """
        Get all fo the numbering types included in hmm_numbering_types
        """
        
        v = dict()
        for number_map in self.numbering_types:
            v[self.numbering_types[number_map]['type']] = 0
        
        return v.keys()
        