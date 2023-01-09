from IgClassifyFASTA import IgClassifyFASTA
from IgClassifyPDB import IgClassifyPDB




class ModelAntibody:
    def __init__(self, fasta_path):
        self.FASTAClass = IgClassifyFASTA(fasta_path)