import os

import numpy as np
from uncurl_analysis import dense_matrix_h5

from .utils import search

PATH = os.path.dirname(__file__)
DATA_PATH = os.path.join(PATH, 'data', 'cell_type_means.h5')
GENES_PATH = os.path.join(PATH, 'data', 'gene_names.txt')

def search_db(input_array, input_gene_names, method='spearman'):
    """
    load gene names/data
    """
    data_dict = dense_matrix_h5.H5Dict(DATA_PATH)
    gene_names = np.loadtxt(GENES_PATH, dtype=str)
    return search(input_array, input_gene_names, data_dict, gene_names, method)
