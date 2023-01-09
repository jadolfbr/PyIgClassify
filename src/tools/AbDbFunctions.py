import sqlite3
from collections import defaultdict

from src.tools import outliers


#A Collection of functions to query the AbDb


def get_pdb_chain_subset(db, gene, use_cutoffs = False, res_cutoff = 2.8, rfac_cutoff = .3):
    """
    Return a list of tuples of [pdb, chain] of the particular gene
    """
    #print db_fname

    db.row_factory = sqlite3.Row
    c = db.cursor()
    if not use_cutoffs:
        c.execute("SELECT DISTINCT PDB, original_chain FROM cdr_data WHERE gene =?", (gene,))
    else:
        c.execute("SELECT DISTINCT PDB, original_chain FROM cdr_data WHERE gene = ? and resolution <= ? and rfactor <= ?", [gene, res_cutoff, rfac_cutoff])
    rows = c.fetchall()
    entries = []
    for row in rows:

        row.keys()
        entries.append([row['PDB'],row['original_chain']])
    return entries

def get_all_lengths(db, cdr, limit_to_known = True, res_cutoff = 2.8, rfac_cutoff = .3):
    """
    Get all unique lengths for a CDR
    """


    c = db.cursor()
    lengths = []

    if limit_to_known:
        in_data = ['loopKeyNotInPaper', cdr, res_cutoff, float(rfac_cutoff)]

        for row in c.execute("SELECT DISTINCT length from "+"cdr_data"+" WHERE datatag!=? and CDR=? and resolution<=? and rfactor<=?", in_data):
            lengths.append(row[0])
    else:
        in_data = [cdr, res_cutoff, rfac_cutoff]
        for row in c.execute("SELECT DISTINCT length from "+"cdr_data WHERE CDR=? and resolution<=? and rfactor<=?", in_data):
            lengths.append(row[0])
    c.close()
    return lengths

def get_all_clusters_for_length(db, cdr, length, limit_to_known = True, res_cutoff = 2.8, rfac_cutoff = .3):
    """
    Get all unique clusters for a length and a cdr
    """
    c = db.cursor()
    clusters = []
    if limit_to_known:
        data = ['loopKeyNotInPaper', cdr, length, res_cutoff, rfac_cutoff]
        for row in c.execute("SELECT DISTINCT fullcluster from cdr_data"+" WHERE datatag!=? and CDR=? and length=? and resolution <= ? and rfactor <= ?", data):
            clusters.append(row[0])
    else:
        data = [cdr, length, res_cutoff, rfac_cutoff]
        for row in c.execute("SELECT DISTINCT fullcluster from cdr_data WHERE CDR=? and length=? and resolution<=? and rfac_cutoff<=?", data):
            clusters.append(row[0])
    c.close()

    return clusters

def get_center_for_cluster_and_length(db, cdr, length, cluster, data_names_array):

    sele = ", ".join(data_names_array)

    #print sele
    data = []
    c = db.cursor()

    #print  "Getting data for: "+" ".join([cdr, repr(length), cluster])
    in_data = ['loopKeyNotInPaper', cdr, length, cluster]
    for row in c.execute("SELECT "+sele+" FROM cdr_data WHERE center=1 and datatag!=? and CDR=? and length=? and fullcluster=?",in_data):
        data.append(row)

    #print data
    if len(data) == 0: return data
    else: return data[0]

def get_center_dih_degrees_for_cluster_and_length(db, cdr, length, cluster):
    """
    Returns a dictionary of center dihedral angles in positional order.  Or returns False if not found.
    result["phis'] = [phis as floats]
    result["psis"] = [Psis as floats]
    result["omegas"] = [Omegas as floats]
    """


    result = defaultdict()
    result["phis"] = []
    result["psis"] = []
    result["omegas"] = []

    data_names_array = ["dihedrals"]
    data = get_center_for_cluster_and_length(db, cdr, length, cluster, data_names_array)
    if len(data) == 0: return False
    dihedrals = data[0]
    dihSP = dihedrals.split(":")

    phii = 0; psii = 1; omega_i = 2
    for i in range(0, length):
        phi = dihSP[phii]; psi = dihSP[psii]; omega = dihSP[omega_i]

        result["phis"].append(float(phi)); result["psis"].append(float(psi)); result["omegas"].append(float(omega))

        phii+=3; psii+=3; omega_i+=3

    assert(len(result["phis"]) == len(result["psis"]) == len(result["omegas"]))
    return result

def get_dihedral_string_for_centers(db, limit_to_known = True):
    sele = "SELECT fullcluster, length_type, ss, dihedrals from cdr_data where center=1"

    data = defaultdict()
    c = db.cursor()
    for row in c.execute(sele):
        data[row[0]] = [row[1], row[2], row[3]]
    return data

def get_data_for_cluster_and_length(db, cdr, length, cluster, data_names_array, limit_to_known = True, res_cutoff = 2.8, rfac_cutoff = .3):
    """

    Get a set of data of a particular length, cdr, and cluster.
    data_names_array is a list of the types of data.  Can include DISTINCT keyword
      Example: data_names_array = ["PDB", "original_chain", "new_chain", "sequence"]


    """
    sele = ", ".join(data_names_array)
    #print sele
    data = []

    c = db.cursor()
    if limit_to_known:
        in_data = ['loopKeyNotInPaper', cdr, length, cluster, res_cutoff, rfac_cutoff]
        for row in c.execute("SELECT "+sele+" FROM cdr_data WHERE datatag!=? and CDR=? and length=? and fullcluster=? and resolution<=? and rfactor <=?", in_data):
            data.append(row)
    else:
        in_data = [cdr, length, cluster, res_cutoff, rfac_cutoff]
        for row in c.execute("SELECT "+sele+" FROM cdr_data WHERE CDR=? and length=? and fullcluster=? and resolution<=? and rfactor<=?", in_data):
            data.append(row)

    c.close()

    return data


def get_unique_sequences_for_cluster(db, cluster, include_outliers, outlier_definition = "conservative"):
    sequences = []
    cur = db.cursor()

    for row in cur.execute("select DISTINCT seq FROM cdr_data WHERE fullcluster=? "+ outliers.get_outlier_string(include_outliers, outlier_definition), [cluster]):
        sequences.append(row[0])
    cur.close()
    return sequences


def get_cdr_rmsd_for_entry(db, pdb, original_chain, cdr, length, fullcluster):
    #Make sure string to to unicode bullshit.
    sele = [str(pdb), str(original_chain), str(cdr), length, str(fullcluster)]

    cur = db.cursor()
    rmsds = []
    for row in cur.execute("select bb_rmsd_cdr_align FROM cdr_data WHERE PDB=? AND original_chain = ? AND CDR=? AND length=? and fullcluster=?", sele):
        rmsds.append(row[0])

    cur.close()
    if len(rmsds) < 1:
        return -1
    else:
        return rmsds[0]

def get_stem_rmsd_for_entry(db, pdb, original_chain, cdr, length, fullcluster):
    #Make sure string to to unicode bullshit.
    sele = [str(pdb), str(original_chain), str(cdr), length, str(fullcluster)]

    cur = db.cursor()
    rmsds = []
    for row in cur.execute("select bb_rmsd_stem_align FROM cdr_data WHERE PDB=? AND original_chain = ? AND CDR=? AND length=? and fullcluster=?", sele):
        rmsds.append(row[0])

    cur.close()
    if len(rmsds) < 1:
        return -1
    else:
        return rmsds[0]