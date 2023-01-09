import sys

def get_outlier_string(include_outliers, outlier_definition, add_AND = True):
    addline = ""
    if not include_outliers:
        addline = get_outlier_definition_string(outlier_definition)
        if add_AND:
            addline = " AND "+addline
    return addline

def get_outlier_definition_string(outlier_definition, rmsd_cutoff = 1.5, dihdist_cutoff = 40):
    """
    Returns a string for adding to a database query which removes outliers.  Need to add AND manually to the string.
    """
    if outlier_definition == "conservative":
        s = " (bb_rmsd_cdr_align < "+repr(rmsd_cutoff)+" AND DistDegree < "+repr(dihdist_cutoff)+") AND DistDegree != -1 AND bb_rmsd_cdr_align != -1"
        return s
    elif outlier_definition == "liberal":
        s = " (bb_rmsd_cdr_align < "+repr(rmsd_cutoff)+" OR DistDegree < "+repr(dihdist_cutoff)+")  AND DistDegree != -1 AND bb_rmsd_cdr_align != -1"
        return s
    else:
        sys.exit("Could not understand outlier definition:"+outlier_definition)

