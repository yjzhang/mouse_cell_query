import os

import numpy as np
from uncurl_analysis import dense_matrix_h5

from .utils import search

PATH = os.path.dirname(__file__)
DATA_MEANS_PATH = os.path.join(PATH, 'data', 'cell_type_means.h5')
DATA_MEDIANS_PATH = os.path.join(PATH, 'data', 'cell_type_medians.h5')
GENES_PATH = os.path.join(PATH, 'data', 'gene_names.txt')

def search_db(input_array, input_gene_names, method='spearman', db='means'):
    """
    args:
        input_array (array): 1d array of shape genes
        input_gene_names (array): 1d array of gene names
        method (str): default: 'spearman'
        db (str): could be 'means' or 'medians'. default: 'means'
    """
    if db == 'means':
        data_dict = dense_matrix_h5.H5Dict(DATA_MEANS_PATH)
    elif db == 'medians':
        data_dict = dense_matrix_h5.H5Dict(DATA_MEDIANS_PATH)
    gene_names = np.loadtxt(GENES_PATH, dtype=str)
    return search(input_array, input_gene_names, data_dict, gene_names, method)

def get_cell_names(db='means'):
    """
    Returns a list of all cell names.
    """
    if db == 'means':
        data_dict = dense_matrix_h5.H5Dict(DATA_MEANS_PATH)
    elif db == 'medians':
        data_dict = dense_matrix_h5.H5Dict(DATA_MEDIANS_PATH)
    return data_dict.keys()

