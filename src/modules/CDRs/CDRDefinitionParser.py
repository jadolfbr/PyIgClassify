from collections import defaultdict

from src.tools.path import *


class CDRDefinitionParser(object):
    def __init__(self, cdr_definition_type):
        self.cdr_definition_type = cdr_definition_type
        self.cdr_definition = defaultdict(dict)
        self._read_data()
    def _read_data(self):
        FILE = open_file(get_db_path()+"/cdr_definitions/"+self.cdr_definition_type+".txt")
        for line in FILE:
            if line[0]=="#":continue
            line = line.strip()
            lineSP = line.split()
            if len(lineSP) != 3: continue
            self.cdr_definition[lineSP[0]]['start'] = int(lineSP[1])
            self.cdr_definition[lineSP[0]]['end'] = int(lineSP[2])
        FILE.close()

    def get_cdr_start(self, cdr_name):
        return self.cdr_definition[cdr_name]['start']
    
    def get_cdr_end(self, cdr_name):
        return self.cdr_definition[cdr_name]['end']
    
    def get_cdr_definition(self):
        return self.cdr_definition