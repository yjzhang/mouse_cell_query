import os

import numpy as np
from uncurl_analysis import dense_matrix_h5

from .utils import search

PATH = os.path.dirname(__file__)
DATA_MEANS_PATH = os.path.join(PATH, 'data', 'cell_type_means.h5')
DATA_MEDIANS_PATH = os.path.join(PATH, 'data', 'cell_type_medians.h5')
DATA_MCA_MEANS_PATH = os.path.join(PATH, 'data', 'cell_type_means_mca.h5')
GENES_PATH = os.path.join(PATH, 'data', 'gene_names.txt')
GENES_MCA_PATH = os.path.join(PATH, 'data', 'cell_type_genes_mca.h5')

def search_db(input_array, input_gene_names, method='spearman', db='tm_means'):
    """
    args:
        input_array (array): 1d array of shape genes
        input_gene_names (array): 1d array of gene names
        method (str): default: 'spearman'
        db (str): Default: 'tm_means' (Tabula Muris means). Also: 'tm_medians', 'mca_means'
    """
    gene_names = None
    gene_data = None
    if db == 'tm_means':
        data_dict = dense_matrix_h5.H5Dict(DATA_MEANS_PATH)
        gene_names = np.loadtxt(GENES_PATH, dtype=str)
    elif db == 'tm_medians':
        data_dict = dense_matrix_h5.H5Dict(DATA_MEDIANS_PATH)
        gene_names = np.loadtxt(GENES_PATH, dtype=str)
    elif db == 'mca_means':
        data_dict = dense_matrix_h5.H5Dict(DATA_MCA_MEANS_PATH)
        gene_data = dense_matrix_h5.H5Dict(GENES_MCA_PATH)
    else:
        raise Exception('db name is unknown')
    return search(input_array, input_gene_names, data_dict, db_gene_names=gene_names, db_gene_data=gene_data, method=method)

def get_cell_names(db='tm_means'):
    """
    Returns a list of all cell names.
    """
    if db == 'tm_means':
        data_dict = dense_matrix_h5.H5Dict(DATA_MEANS_PATH)
    elif db == 'tm_medians':
        data_dict = dense_matrix_h5.H5Dict(DATA_MEDIANS_PATH)
    elif db == 'mca_means' or db=='mca_medians':
        data_dict = dense_matrix_h5.H5Dict(DATA_MCA_MEANS_PATH)
    return data_dict.keys()

def get_gene_names(db='tm_means'):
    if db.startswith('tm'):
        gene_names = np.loadtxt(GENES_PATH, dtype=str)
        return gene_names
    elif db == 'mca_means' or db=='mca_medians':
        gene_data = dense_matrix_h5.H5Dict(GENES_MCA_PATH)
        return gene_data
