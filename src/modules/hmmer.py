###############################################################################
#                                                                             #
#    simpleHMMER.py                                                           #
#                                                                             #
#    Run HMMER and provide an API for Parsing output                          #
#                                                                             #
#    Copyright (C) Michael Imelfort                                           #
#                                                                             #
###############################################################################
#                                                                             #
#            .d8888b.  d8b                        888                         #
#           d88P  Y88b Y8P                        888                         #
#           Y88b.                                 888                         #
#            "Y888b.   888 88888b.d88b.  88888b.  888  .d88b.                 #
#               "Y88b. 888 888 "888 "88b 888 "88b 888 d8P  Y8b                #
#                 "888 888 888  888  888 888  888 888 88888888                #
#           Y88b  d88P 888 888  888  888 888 d88P 888 Y8b.   .                #
#            "Y8888P"  888 888  888  888 88888P"  888  "Y8888                 #
#                                        888                                  #
#                                        888                                  #
#                                        888                                  #
#                                                                             #
#        888    888 888b     d888 888b     d888 8888888888 8888888b.          #
#        888    888 8888b   d8888 8888b   d8888 888        888   Y88b         #
#        888    888 88888b.d88888 88888b.d88888 888        888    888         #
#        8888888888 888Y88888P888 888Y88888P888 8888888    888   d88P         #
#        888    888 888 Y888P 888 888 Y888P 888 888        8888888P"          #
#        888    888 888  Y8P  888 888  Y8P  888 888        888 T88b           #
#        888    888 888   "   888 888   "   888 888        888  T88b          #
#        888    888 888       888 888       888 8888888888 888   T88b         #
#                                                                             #
###############################################################################
#                                                                             #
#    This program is free software: you can redistribute it and/or modify     #
#    it under the terms of the GNU General Public License as published by     #
#    the Free Software Foundation, either version 3 of the License, or        #
#    (at your option) any later version.                                      #
#                                                                             #
#    This program is distributed in the hope that it will be useful,          #
#    but WITHOUT ANY WARRANTY; without even the implied warranty of           #
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the            #
#    GNU General Public License for more details.                             #
#                                                                             #
#    You should have received a copy of the GNU General Public License        #
#    along with this program. If not, see <http://www.gnu.org/licenses/>.     #
#                                                                             #
###############################################################################

__author__ = "Michael Imelfort"
__copyright__ = "Copyright 2012"
__credits__ = ["Michael Imelfort"]
__license__ = "GPL3"
__version__ = "0.2.3"
__maintainer__ = "Michael Imelfort"
__email__ = "mike@mikeimelfort.com"
__status__ = "Development"

###############################################################################

from re import split as re_split
from os import system, makedirs
from os.path import join as osp_join
import errno
import sys

###############################################################################
###############################################################################
###############################################################################
###############################################################################

class FormatError(BaseException): pass
class HMMERError(BaseException): pass
class HMMMERModeError(BaseException): pass

###############################################################################
###############################################################################
###############################################################################
###############################################################################

class HMMERRunner():
    """Wrapper for running HMMER3"""
    def __init__(self, mode="dom", prefix=''):
        # make sure HMMER is installed!
        # checkForHMMER()

        # set the mode
        if mode == "dom":
            self.mode = 'domtblout'
        elif mode == "tbl":
            self.mode = 'tblout'
        elif mode == 'align':
            self.mode = 'align'
        elif mode == 'fetch':
            self.mode = 'fetch'
        else:
            raise HMMMERModeError("Mode %s not understood" % mode)
        # make the output file names
        (self.txtOut, self.hmmOut) = makeOutputFNs(prefix, mode=self.mode)

    def search(self, db, query, outputDir):
        """Run hmmsearch"""
        # make the output dir and files
        if self.mode != 'domtblout' and self.mode != 'tblout':
            raise HMMMERModeError("Mode %s not compatible with search" % self.mode)
        makeSurePathExists(outputDir)
        txt_file = osp_join(outputDir, self.txtOut)
        hmm_file = osp_join(outputDir, self.hmmOut)

        # run hmmer!
        system('hmmsearch --%s %s %s %s > %s' % (self.mode, txt_file, db, query, hmm_file))

    def align(self, db, query, outputFile, writeMode='>',
            outputFormat='PSIBLAST', allcol=True, trim=True):
        """Run hmmalign"""
        if self.mode != 'align':
            raise HMMMERModeError("Mode %s not compatible with align" % mode)
        #build up the option string
        opts=''
        if allcol:
            opts+=' --allcol '
        if trim:
            opts+=' --trim '
        # run hmmer!
        cmd = 'hmmalign%s--outformat %s %s %s %s %s' % (opts, outputFormat, db, query, writeMode, outputFile)
        system(cmd)

    def fetch(self, db, key, fetchFileName, emitFileName=''):
        """Run hmmfetch and possible hmmemit"""
        if self.mode != 'fetch':
            raise HMMMERModeError("Mode %s not compatible with fetch" % self.mode)
        # run hmmer!
        system('hmmfetch %s %s > %s' % (db, key, fetchFileName))
        # run emit if we;ve been told to
        if emitFileName != '':
            system('hmmemit -c %s > %s' % (fetchFileName, emitFileName))

def makeSurePathExists(path):
    try:
        makedirs(path)
    except OSError as exception:
        if exception.errno != errno.EEXIST:
            raise

def makeOutputFNs(prefix='', mode='domtblout'):
    """Consistent interface for making output filenames"""
    if prefix == '':
        txtOut = 'hmmer_out.txt'
        if mode == 'align':
            hmmOut = 'hmmer_out.align'
        else:
            hmmOut = 'hmmer_out.hmmer3'
    else:
        txtOut = '%s_out.txt' % prefix
        if mode == 'align':
            hmmOut = '%s_out.align' % prefix
        else:
            hmmOut = '%s_out.hmmer3' % prefix
    return (txtOut, hmmOut)

def checkForHMMER():
    """Check to see if HMMER is on the system before we try fancy things

    We assume that a successful hmmsearch -h returns 0 and anything
    else returns something non-zero
    """
    # redirect stdout so we don't get mess!
    try:
        exit_status = system('hmmsearch -h > /dev/null')
    except:
      print "Unexpected error!", sys.exc_info()[0]
      raise

    if exit_status != 0:
        raise HMMERError("Error attempting to run hmmsearch, is it in your path?")

###############################################################################
###############################################################################
###############################################################################
###############################################################################

class HMMERParser():
    """Parses tabular output"""
    def __init__(self, fileHandle, mode='dom'):
        """Give this guy an open file!"""
        self.handle = fileHandle
        if mode == 'dom':
            self.mode = 'domtblout'
        elif mode == 'tbl':
            self.mode = 'tblout'
        else:
            raise HMMERError("Mode %s not understood, please use 'dom' or 'tbl'" % mode)

    def next(self):
        """Get the next result in the file"""
        while 1:
            if self.mode == 'domtblout':
                hit = self.readHitsDOM()
            elif self.mode == 'tblout':
                self.hit = self.readHitsTBL()
            else:
                raise HMMERError("Mode %s not understood" % self.mode)

            if hit == {}:
                return None
            else:
                return hit

    #-------------------------------------------------
    # TBL mode
    def readHitsTBL(self):
        """Look for the next hit in tblout format, package and return"""
        """
We expect line to look like:
NODE_110054_length_1926_cov_24.692627_41_3 -          Ribosomal_S9         PF00380.14   5.9e-48  158.7   0.0   6.7e-48  158.5   0.0   1.0   1   0   0   1   1   1   1 # 1370 # 1756 # 1 # ID=41_3;partial=00;start_type=ATG;rbs_motif=None;rbs_spacer=None
        """
        while (1):
            line = self.handle.readline().rstrip()
            try:
                if line[0] != '#' and len(line) != 0:
                    dMatch = re_split( r'\s+', line.rstrip() )
                    if len(dMatch) < 19:
                        raise FormatError( "Something is wrong with this line:\n%s" % (line) )
                    refined_match = dMatch[0:18] + [" ".join([str(i) for i in dMatch[18:]])]
                    return HmmerHitTBL(refined_match)
            except IndexError:
                return {}

    #-------------------------------------------------
    # DOM mode
    def readHitsDOM(self):
        """Look for the next hit in domtblout format, package and return"""
        """
We expect the line to look like:
NODE_925902_length_6780_cov_18.428171_754_2 -            399 PGK                  PF00162.14   384  2.2e-164  543.7   0.1   1   1  1.3e-167  2.5e-164  543.5   0.1     1   384     9   386     9   386 1.00 # 1767 # 2963 # -1 # ID=754_2;partial=00;start_type=ATG;rbs_motif=AGGA;rbs_spacer=5-10bp
        """
        while (1):
            line = self.handle.readline().rstrip()
            try:
                if line[0] != '#' and len(line) != 0:
                    dMatch = re_split( r'\s+', line.rstrip() )
                    if len(dMatch) < 23:
                        raise FormatError( "Something is wrong with this line:\n%s" % (line) )
                    refined_match = dMatch[0:22] + [" ".join([str(i) for i in dMatch[22:]])]
                    return HmmerHitDOM(refined_match)
            except IndexError:
                return {}

###############################################################################
###############################################################################
###############################################################################
###############################################################################

class HmmerHitTBL():
    """Encapsulate a hmmer hit given in tblout format"""
    def __init__(self, values):
        """Load em' in"""
        if len(values) == 19:
            self.target_name = values[0]
            self.target_accession = values[1]
            self.query_name = values[2]
            self.query_accession = values[3]
            self.full_e_value = float(values[4])
            self.full_score = float(values[5])
            self.full_bias = float(values[6])
            self.best_e_value = float(values[7])
            self.best_score = float(values[8])
            self.best_bias = float(values[9])
            self.exp = float(values[10])
            self.reg = int(values[11])
            self.clu = int(values[12])
            self.ov = int(values[13])
            self.env = int(values[14])
            self.dom = int(values[15])
            self.rep = int(values[16])
            self.inc = int(values[17])
            self.target_description = values[18]

    def __str__(self):
        """when we need to print"""
        return "\t".join([self.target_name,
                          self.target_accession,
                          self.query_name,
                          self.query_accession,
                          str(self.full_e_value),
                          str(self.full_score),
                          str(self.full_bias),
                          str(self.best_e_value),
                          str(self.best_score),
                          str(self.best_bias),
                          str(self.exp),
                          str(self.reg),
                          str(self.clu),
                          str(self.ov),
                          str(self.env),
                          str(self.dom),
                          str(self.rep),
                          str(self.inc),
                          self.target_description
                          ]
                         )

###############################################################################
###############################################################################
###############################################################################
###############################################################################

class HmmerHitDOM():
    """Encapsulate a hmmer hit given in domtblout format"""
    def __init__(self, values):
        """Load em' in"""
        if len(values) == 23:
            self.target_name = values[0]
            self.target_accession = values[1]
            self.target_length = int(values[2])
            self.query_name = values[3]
            self.query_accession = values[4]
            self.query_length = int(values[5])
            self.full_e_value = float(values[6])
            self.full_score = float(values[7])
            self.full_bias = float(values[8])
            self.dom = int(values[9])
            self.ndom = int(values[10])
            self.c_evalue = float(values[11])
            self.i_evalue = float(values[12])
            self.dom_score = float(values[13])
            self.dom_bias = float(values[14])
            self.hmm_from = int(values[15])
            self.hmm_to = int(values[16])
            self.ali_from = int(values[17])
            self.ali_to = int(values[18])
            self.env_from = int(values[19])
            self.env_to = int(values[20])
            self.acc = float(values[21])
            self.target_description = values[22]

    def __str__(self):
        """when we need to print"""
        return "\t".join([self.target_name,
                          self.target_accession,
                          str(self.target_length),
                          self.query_name,
                          self.query_accession,
                          str(self.query_length),
                          str(self.full_e_value),
                          str(self.full_score),
                          str(self.full_bias),
                          str(self.dom),
                          str(self.ndom),
                          str(self.c_evalue),
                          str(self.i_evalue),
                          str(self.dom_score),
                          str(self.dom_bias),
                          str(self.hmm_from),
                          str(self.hmm_to),
                          str(self.ali_from),
                          str(self.ali_to),
                          str(self.env_from),
                          str(self.env_to),
                          str(self.acc),
                          self.target_description]
                         )

###############################################################################
###############################################################################
###############################################################################
###############################################################################


