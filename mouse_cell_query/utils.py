# utils for querying...

import numpy as np
import scipy.stats

def gene_overlap_indices(input_gene_names, db_gene_names):
    """
    This returns two arrays of integers: an array of indices in input_genes,
    and an array of indices in db_genes, such that they refer to the same genes.

    Args:
        input_gene_names: list or array of strings
        db_gene_names: list or array of strings
    """
    # match genes between the db and the input data
    db_genes_set = set(db_gene_names)
    # map of gene names to IDs
    data_gene_ids_map = {}
    for i, gene in enumerate(input_gene_names):
        data_gene_ids_map[gene] = i
    db_gene_ids = []
    data_gene_ids = []
    for index, gene in enumerate(input_gene_names):
        if gene in db_genes_set:
            db_gene_ids.append(index)
            data_gene_ids.append(data_gene_ids_map[gene])
    db_gene_ids = np.array(db_gene_ids)
    data_gene_ids = np.array(data_gene_ids)
    return data_gene_ids, db_gene_ids

def spearman_search(input_data, db_data, db_gene_ids):
    """
    Search using Spearman correlation. assumes that the input data has already been aligned
    by gene name.

    Args:
        input_data (array): 1d array of shape genes
        db_data (dict): dict of cell type name : 1d array
        db_gene_ids (array): 1d array of ints

    Returns:
        list of (cell type, corr) tuples, sorted by decreasing correlation
    """
    results = []
    for cell_type_name, data in db_data.items():
        corr = scipy.stats.spearmanr(input_data, data[db_gene_ids])[0]
        results.append((cell_type_name, corr))
    # sort by decreasing correlation
    results.sort(key=lambda x: x[1], reverse=True)
    return results

def search(input_data, input_gene_names, db_data, db_gene_names, method='spearman'):
    """
    Finds the most similar cell types by the given method.

    Args:
        input_data (array): 1d array representing gene expression
        input_gene_names (array or list): sequence of gene names present in data
        db_data (dict): dict of cell type : 1d array of mean gene expression for that cell type
        db_gene_names (array or list): sequence of gene names present in db
        method (str): currently, only 'spearman' is supported.

    Returns:
        list of (cell type, similarity metric) sorted by decreasing similarity
    """
    data_gene_ids, db_gene_ids = gene_overlap_indices(input_gene_names, db_gene_names)
    input_data = input_data[data_gene_ids]
    if method == 'spearman':
        results = spearman_search(input_data, db_data, db_gene_ids)
    return results
