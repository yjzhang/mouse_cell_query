import os

import numpy as np
from uncurl_analysis import dense_matrix_h5

from .utils import search

PATH = os.path.dirname(__file__)
DATA_MEANS_PATH = os.path.join(PATH, 'data', 'cell_type_means.h5')
DATA_FACS_MEANS_PATH = os.path.join(PATH, 'data', 'cell_type_means_tm_facs.h5')
DATA_MEDIANS_PATH = os.path.join(PATH, 'data', 'cell_type_medians.h5')
DATA_MCA_MEANS_PATH = os.path.join(PATH, 'data', 'cell_type_means_mca.h5')
DATA_MCA_COARSE_MEANS_PATH = os.path.join(PATH, 'data', 'cell_type_means_mca_coarse.h5')
DATA_ALLEN_CLUSTER_MEANS_PATH = os.path.join(PATH, 'data', 'allen_class_cluster_means.h5')
DATA_ALLEN_CLASS_MEANS_PATH = os.path.join(PATH, 'data', 'allen_class_means.h5')
DATA_ALLEN_SUBCLASS_MEANS_PATH = os.path.join(PATH, 'data', 'allen_class_subclass_means.h5')

GENES_PATH = os.path.join(PATH, 'data', 'gene_names.txt')
GENES_MCA_PATH = os.path.join(PATH, 'data', 'cell_type_genes_mca.h5')
GENES_MCA_COARSE_PATH = os.path.join(PATH, 'data', 'mca_coarse_gene_names.txt')
GENES_ALLEN_PATH = os.path.join(PATH, 'data', 'allen_gene_names.txt')

# TODO: build index with NMSlib

DB_TO_PATH = {
        'tm_means': DATA_MEANS_PATH,
        'tm_means_normalized': os.path.join(PATH, 'data', 'cell_type_means_tm_droplet_normalized.h5'),
        'tm_medians': DATA_MEDIANS_PATH,
        'mca_means': DATA_MCA_MEANS_PATH,
        'mca_coarse_means': DATA_MCA_COARSE_MEANS_PATH,
        'allen_cluster_means': DATA_ALLEN_CLUSTER_MEANS_PATH,
        'allen_class_means': DATA_ALLEN_CLASS_MEANS_PATH,
        'allen_subclass_means': DATA_ALLEN_SUBCLASS_MEANS_PATH,
        'tm_facs_means': DATA_FACS_MEANS_PATH,
        'tm_facs_means_normalized': os.path.join(PATH, 'data', 'cell_type_means_tm_facs_normalized.h5'),
}

DB_TO_GENE_PATH = {
        'tm_means': GENES_PATH,
        'tm_means_normalized': GENES_PATH,
        'tm_medians': GENES_PATH,
        'mca_means': GENES_MCA_PATH,
        'mca_coarse_means': GENES_MCA_COARSE_PATH,
        'allen_cluster_means': GENES_ALLEN_PATH,
        'allen_class_means': GENES_ALLEN_PATH,
        'allen_subclass_means': GENES_ALLEN_PATH,
        'tm_facs_means': os.path.join(PATH, 'data', 'tabula_muris_facs_genes.txt'),
        'tm_facs_means_normalized': os.path.join(PATH, 'data', 'tabula_muris_facs_genes.txt')
}

def get_dbs():
    """
    Returns a list of all valid db names.
    """
    return DB_TO_PATH.keys()

def search_db(input_array, input_gene_names, method='spearman', db='tm_means'):
    """
    Args:
        input_array (array): 1d array of shape genes
        input_gene_names (array): 1d array of gene names
        method (str): default: 'spearman'
        db (str): Default: 'tm_means' (Tabula Muris means). Also: 'tm_medians', 'mca_means', 'mca_coarse_means'

    Returns:
        list of tuples of (cell name, score), sorted by similarity to query
    """
    gene_names = None
    gene_data = None
    data_dict = dense_matrix_h5.H5Dict(DB_TO_PATH[db])
    gene_path = DB_TO_GENE_PATH[db]
    if gene_path.endswith('h5'):
        gene_data = dense_matrix_h5.H5Dict(gene_path)
    else:
        gene_names = np.loadtxt(gene_path, dtype=str)
    return search(input_array, input_gene_names, data_dict, db_gene_names=gene_names, db_gene_data=gene_data, method=method)

def get_cell_names(db='tm_means'):
    """
    Returns a list of all cell names.
    """
    data_dict = dense_matrix_h5.H5Dict(DB_TO_PATH[db])
    return data_dict.keys()

def get_gene_names(db='tm_means'):
    gene_path = DB_TO_GENE_PATH[db]
    if gene_path.endswith('h5'):
        gene_data = dense_matrix_h5.H5Dict(gene_path)
        return gene_data
    else:
        gene_names = np.loadtxt(gene_path, dtype=str)
        return gene_names
