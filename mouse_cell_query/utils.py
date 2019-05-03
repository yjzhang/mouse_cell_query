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
    db_gene_names = [x.upper() for x in db_gene_names]
    input_gene_names = [x.upper() for x in input_gene_names]
    db_genes_set = set(db_gene_names)
    # map of gene names to IDs
    db_gene_ids_map = {}
    for i, gene in enumerate(db_gene_names):
        db_gene_ids_map[gene] = i
    db_gene_ids = []
    data_gene_ids = []
    for index, gene in enumerate(input_gene_names):
        if gene in db_genes_set:
            data_gene_ids.append(index)
            db_gene_ids.append(db_gene_ids_map[gene])
    db_gene_ids = np.array(db_gene_ids)
    data_gene_ids = np.array(data_gene_ids)
    return data_gene_ids, db_gene_ids

def spearman_search(input_data, db_data, data_gene_ids, db_gene_ids):
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
    input_data = input_data[data_gene_ids]
    for cell_type_name, data in db_data.items():
        corr = scipy.stats.spearmanr(input_data, data[db_gene_ids])[0]
        results.append((cell_type_name, corr))
    # sort by decreasing correlation
    results.sort(key=lambda x: x[1], reverse=True)
    return results

def spearman_nonzero_search(input_data, db_data, data_gene_ids, db_gene_ids):
    """
    only use nonzero elements in input_data
    """
    results = []
    input_data = input_data[data_gene_ids]
    nonzeros = np.nonzero(input_data)[0]
    input_data = input_data[nonzeros]
    for cell_type_name, data in db_data.items():
        corr = scipy.stats.spearmanr(input_data, data[db_gene_ids][nonzeros])[0]
        results.append((cell_type_name, corr))
    # sort by decreasing correlation
    results.sort(key=lambda x: x[1], reverse=True)
    return results

def poisson_search(input_data, db_data, data_gene_ids, db_gene_ids):
    """
    Search using Poisson distance
    """
    from uncurl_analysis import bulk_data
    input_data = input_data[data_gene_ids]
    results = []
    for cell_type_name, data in db_data.items():
        data = data/data.sum()
        data_subset = data[db_gene_ids]
        dist = bulk_data.log_prob_poisson(data_subset, input_data)
        results.append((cell_type_name, dist))
    # sort by decreasing correlation
    results.sort(key=lambda x: x[1], reverse=True)
    return results

def cosine_search(input_data, db_data, data_gene_ids, db_gene_ids):
    """
    Search using cosine similarity
    """
    from uncurl_analysis import bulk_data
    input_data = input_data[data_gene_ids]
    results = []
    for cell_type_name, data in db_data.items():
        data_subset = data[db_gene_ids]
        dist = bulk_data.cosine(data_subset, input_data)[0][0]
        results.append((cell_type_name, dist))
    # sort by decreasing cosine similarity
    results.sort(key=lambda x: x[1], reverse=True)
    return results


def hamming_search(input_data, db_data, data_gene_ids, db_gene_ids):
    """
    Search using Hamming distance on binarized data
    """
    import scipy
    input_data = input_data[data_gene_ids]
    input_data = (input_data != 0)
    results = []
    for cell_type_name, data in db_data.items():
        data_subset = data[db_gene_ids]
        data_subset = (data_subset != 0)
        dist = scipy.spatial.distance.hamming(input_data, data_subset)
        results.append((cell_type_name, dist))
    # sort by increasing distance
    results.sort(key=lambda x: x[1], reverse=False)
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
    if method == 'spearman':
        f = spearman_search
    elif method == 'poisson':
        f = poisson_search
    elif method == 'spearman_nonzero':
        f = spearman_nonzero_search
    elif method == 'cosine':
        f = cosine_search
    elif method == 'hamming':
        f = hamming_search
    results = f(input_data, db_data, data_gene_ids, db_gene_ids)
    return results
