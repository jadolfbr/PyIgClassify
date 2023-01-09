import os
import sys
import multiprocessing



class ClustalRunner:
    """
    A very simple class wrapper to run clustal omega.
    """
    def __init__(self, fasta_path, clustal_name="clustal_omega", clustal_dir = None):

        self.fasta_path = fasta_path

        self.clustal_name = clustal_name
        self.clustal_dir = clustal_dir

        self.output_format = "clu"
        self.accepted_output_formats = ['fasta', 'clustal', 'msf', 'phylip', 'selex', 'stockholm', 'vienna',
                                        'a2m', 'fa', 'clu', 'phy', 'st', 'vie']

        self.hard_wrap = 100
        self.threads = multiprocessing.cpu_count()
        self.extra_options = ""

    def set_fasta_path(self, fasta_path):
        """
        Set the fasta path for alignment.
        """
        self.fasta_path = fasta_path

    def set_output_format(self, output_format):
        """
        Set the output format
        """
        if not output_format in self.accepted_output_formats:
            print "Unrecognized clustal output format: "+output_format
            return

        self.output_format = output_format

    def set_hard_wrap(self, hard_wrap):
        """
        Set the number of charactors before clustal will wrap.  Usually 60-80.
        """
        self.hard_wrap = hard_wrap

    def set_threads(self, threads):
        """
        Limit the number of threads for Clustal
        """
        self.threads = threads

    def set_extra_options(self, extra_options = ""):
        """
        Set any extra options as a string which will be added to the end of the command line.
        """

        self.extra_options = extra_options

    def output_alignment(self, out_dir, out_name):
        """
        Configure command line and Run Clustal Omega
        """
        if not os.path.exists(self.fasta_path):
            sys.exit("Fasta Path not good!  Cannot continue!")
        if not os.path.exists(out_dir): os.mkdir(out_dir)

        print "Running alignment on: "+self.fasta_path
        print "Outputting to: "+out_dir
        print "As - "+out_name

        if self.clustal_dir:
            clustal = self.clustal_dir+"/"+self.clustal_name
        else:
            clustal = self.clustal_name

        out_path = out_dir+"/"+out_name

        cmd_line = clustal+ \
                  " -i "+self.fasta_path+ \
                  " -o "+out_path+ \
                  " --outfmt "+self.output_format+ \
                  " --wrap "+repr(self.hard_wrap)+ \
                  " "+self.extra_options

        if self.threads != multiprocessing.cpu_count():
            cmd_line = cmd_line+ " --threads "+repr(self.threads)

        print "Running Clustal Omega: "+cmd_line

        os.system(cmd_line)
