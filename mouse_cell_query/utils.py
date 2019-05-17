# utils for querying...

import numpy as np
import scipy.stats

try:
    from functools import lru_cache
except ImportError:
    from backports.functools_lru_cache import lru_cache


def gene_overlap_indices(input_gene_names, db_gene_names):
    """
    This returns two arrays of integers: an array of indices in input_genes,
    and an array of indices in db_genes, such that they refer to the same genes.

    Args:
        input_gene_names: list or array of strings
        db_gene_names: list or array of strings
    """
    return tuple_gene_overlap_indices(tuple(input_gene_names), tuple(db_gene_names))

@lru_cache(maxsize=None)
def tuple_gene_overlap_indices(input_gene_names, db_gene_names):
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

def spearman_search(input_data, db_data):
    """
    Search using Spearman correlation. assumes that the input data has already been aligned
    by gene name.

    Args:
        input_data (array): 1d array of shape genes
        db_data: 1d array

    Returns: spearman correlation (float)
    """
    corr = scipy.stats.spearmanr(input_data, db_data)[0]
    return corr

def spearman_nonzero_search(input_data, db_data):
    """
    only use nonzero elements in input_data
    """
    nonzeros = np.nonzero(input_data)[0]
    input_data = input_data[nonzeros]
    corr = scipy.stats.spearmanr(input_data, db_data[nonzeros])[0]
    return corr

def kendall_tau(input_data, db_data):
    return scipy.stats.kendalltau(input_data, db_data)[0]

def poisson_search(input_data, db_data):
    """
    Search using Poisson distance
    """
    from uncurl_analysis import bulk_data
    data = db_data/db_data.sum()
    dist = bulk_data.log_prob_poisson(data, input_data)
    return dist

def cosine_search(input_data, db_data):
    """
    Search using cosine similarity
    """
    from uncurl_analysis import bulk_data
    dist = bulk_data.cosine(db_data, input_data)[0][0]
    return dist


def hamming_search(input_data, db_data):
    """
    Search using Hamming distance on binarized data
    """
    import scipy
    input_data = (input_data != 0)
    db_data = (db_data != 0)
    dist = scipy.spatial.distance.hamming(input_data, db_data)
    return dist

def random_result(input_data, db_data):
    """
    random similarity score between 0 and 1, just for testing purposes
    """
    import random
    return random.random()

def search(input_data, input_gene_names, db_data, db_gene_names=None, db_gene_data=None, method='spearman'):
    """
    Finds the most similar cell types by the given method.

    Args:
        input_data (array): 1d array representing gene expression
        input_gene_names (array or list): sequence of gene names present in data
        db_data (dict): dict of cell type : 1d array of mean gene expression for that cell type
        db_gene_names (array, list): sequence of gene names present in db
        db_gene_data (dict): dict of cell type : array of gene names.
        method (str): currently, only 'spearman' is supported.

    Returns:
        list of (cell type, similarity metric) sorted by decreasing similarity
    """
    if method == 'spearman':
        f = spearman_search
        reverse = True
    elif method == 'kendall':
        f = kendall_tau
        reverse = True
    elif method == 'poisson':
        f = poisson_search
        reverse = True
    elif method == 'spearman_nonzero':
        f = spearman_nonzero_search
        reverse = True
    elif method == 'cosine':
        f = cosine_search
        reverse = True
    elif method == 'hamming':
        f = hamming_search
        reverse = False
    elif method == 'random':
        f = random_result
        reverse = False
    if db_gene_names is not None:
        data_gene_ids, db_gene_ids = gene_overlap_indices(input_gene_names, db_gene_names)
        data_subset = input_data[data_gene_ids]
    results = []
    for cell_type_name, data in db_data.items():
        if db_gene_names is not None:
            db_data_subset = data[db_gene_ids]
        elif db_gene_data is not None:
            db_genes = db_gene_data[cell_type_name].astype(str)
            data_gene_ids, db_gene_ids = gene_overlap_indices(input_gene_names, db_genes)
            data_subset = input_data[data_gene_ids]
            db_data_subset = data[db_gene_ids]
        results.append((cell_type_name, f(data_subset, db_data_subset)))
    results.sort(key=lambda x: x[1], reverse=reverse)
    return results
