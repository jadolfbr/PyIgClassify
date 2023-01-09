import sys

from src.modules.Structure import FrameworkRegions
from src.modules.Structure import PDBResInfo
from src.modules.restype_definitions import definitions

import src.tools.general as general_tools
from src.modules.PDBConsensusInfo import PDBConsensusInfo


class InsertFrameworkTables(object):
    def __init__(self, pdb_info_array, gene):
        self.set_data(pdb_info_array, gene)
        self.aas = definitions().get_all_one_letter_codes()

    def set_data(self, pdb_info_array, gene):
        self.pdb_info_array = pdb_info_array
        self.gene = gene

        self.regions = self._setup_regions(gene)

    def _setup_regions(self, gene):

        if gene == 'heavy':
            self.chain = 'H'
        else:
            self.chain = 'L'

        return FrameworkRegions(self.chain)


    def create_raw_table(self, con, dataset_name, table_name = "framework_raw", drop_table = False):

        with con:
            cur = con.cursor()
            if drop_table:
                cur.execute("Drop table if exists "+table_name)

            cur.execute("CREATE TABLE IF NOT EXISTS "+table_name+
                        " (ID INTEGER PRIMARY KEY, name TEXT, description TEXT, dataset TEXT, gene TEXT, species TEXT, " +
                          "seq TEXT, f1_seq TEXT, f1_2_seq TEXT, f2_3_seq TEXT, f3_seq TEXT)")

            species = "NA"

            i = 1
            for pdb_res_info in self.pdb_info_array:
                if not isinstance(pdb_res_info, PDBResInfo): sys.exit()
                seq = general_tools.get_raw_framework_sequence(pdb_res_info)
                f1_seq = general_tools.get_seq_between_pdb_nums(pdb_res_info, self.regions.get_start("F1"), self.regions.get_stop("F1"), self.chain)
                f1_2_seq = general_tools.get_seq_between_pdb_nums(pdb_res_info, self.regions.get_start("F1_2"), self.regions.get_stop("F1_2"), self.chain)
                f2_3_seq = general_tools.get_seq_between_pdb_nums(pdb_res_info, self.regions.get_start("F2_3"), self.regions.get_stop("F2_3"), self.chain)
                f3_seq = general_tools.get_seq_between_pdb_nums(pdb_res_info, self.regions.get_start("F3"), self.regions.get_stop("F3"), self.chain)

                data = [None, pdb_res_info.get_extra_data('name'), pdb_res_info.get_extra_data('description'), dataset_name, self.gene, species, seq, f1_seq, f1_2_seq, f2_3_seq, f3_seq]

                cur.execute("INSERT INTO "+table_name+" VALUES(?,?,?,?,?,?,?,?,?,?,?)", data)

                i+=1

    def create_consensus_table(self, pdb_consensus_info, con, dataset_name, table_name = "framework_consensus", drop_table = False):
        if not isinstance(pdb_consensus_info, PDBConsensusInfo): sys.exit()
        with con:
            cur = con.cursor()
            if drop_table:
                cur.execute("Drop table if exists "+table_name)
            cur.execute("CREATE TABLE IF NOT EXISTS "+table_name+" " +
                        "(ID INTEGER PRIMARY KEY, dataset TEXT, gene TEXT, species TEXT, consensus_seq TEXT, " +
                         "f1_consensus TEXT, f1_2_consensus TEXT, f2_3_consensus TEXT, f3_consensus TEXT)")

            f1_con = general_tools.get_conensus_bt_pdb_nums(pdb_consensus_info, self.regions.F1()[0], self.regions.F1()[1], self.chain)
            f1_2_con = general_tools.get_conensus_bt_pdb_nums(pdb_consensus_info, self.regions.F1_2()[0], self.regions.F1_2()[1], self.chain)
            f2_3_con = general_tools.get_conensus_bt_pdb_nums(pdb_consensus_info, self.regions.F2_3()[0], self.regions.F2_3()[1], self.chain)
            f3_con = general_tools.get_conensus_bt_pdb_nums(pdb_consensus_info, self.regions.F3()[0], self.regions.F3()[1], self.chain)

            consensus = f1_con+f1_2_con+f2_3_con+f3_con

            species = "NA"
            data = [None,dataset_name, self.gene, species, consensus, f1_con, f1_2_con, f2_3_con, f3_con]

            cur.execute("INSERT INTO "+table_name+" VALUES(?,?,?,?,?,?,?,?,?)", data)


    def create_prob_table(self, pdb_consensus_info, con, dataset_name, table_name = "framework_prob", drop_table = False):
        if not isinstance(pdb_consensus_info, PDBConsensusInfo): sys.exit()
        with con:
            cur = con.cursor()
            if drop_table:
                cur.execute("Drop table if exists "+table_name)

            cur.execute("CREATE TABLE IF NOT EXISTS "+table_name+ \
                        " (ID INTEGER PRIMARY KEY, dataset TEXT, gene TEXT, species TEXT," +
                          "pdb_num INT, chain TEXT, icode TEXT, aa TEXT, probability FLOAT, frequency INT, total_seq INT)")
            regions = self.regions.get_regions()
            species = "NA"
            for region in regions:
                for i in range(region[0], region[1]+1):
                    position = (i, self.chain, " ")
                    for aa in self.aas:
                        prob = pdb_consensus_info.get_probability_for_position(position, aa)
                        freq = pdb_consensus_info.get_frequency_for_position(position, aa)
                        total_seq = pdb_consensus_info._get_total_entries(position)

                        data = [None, dataset_name, self.gene, species, i, self.chain, " ",aa, prob, freq, total_seq ]
                        cur.execute("INSERT INTO "+table_name+" VALUES(?,?,?,?,?,?,?,?,?,?,?)", data)