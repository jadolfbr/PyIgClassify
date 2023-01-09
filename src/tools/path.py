#Author: Jared Adolf-Bryfogle
#Purpose: A collection of path tools

import os
import gzip

def get_PyIgClassify_root():
    return os.path.abspath(os.path.dirname(__file__)+"/../..")

def get_db_path():
    return get_PyIgClassify_root()+"/database"

def get_db_path_full():
    return get_PyIgClassify_root()+"/database"

def get_DBOUT_path_full():
    return get_PyIgClassify_root()+"/DBOUT"

def get_main_path_full():
    return get_PyIgClassify_root()

def open_file(file, mode='r'):
    if file.split(".")[-1] =="gz":
        #print "opening gzipped file"
        INFILE = gzip.open(file, mode+'b')
    else:
        INFILE = open(file, mode)

    return INFILE