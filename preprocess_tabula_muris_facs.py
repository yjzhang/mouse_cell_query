# TODO: preprocess the tabula muris data into a format that can be queried upon.
# data source: https://figshare.com/articles/Single-cell_RNA-seq_data_from_microfluidic_emulsion_v2_/5968960
# tabula_muris/droplet.zip

from collections import defaultdict
import os

import numpy as np
import pandas as pd
import scipy.io
from scipy import sparse

metadata_facs = pd.read_csv('tabula_muris/metadata_FACS.csv')
annotations_facs = pd.read_csv('tabula_muris/annotations_facs.csv')

current_channel = None
current_tissue = None
current_matrix = None
current_barcodes = None
tissue_values = defaultdict(lambda: [])
cell_type_values = defaultdict(lambda: [])
tissue_cell_types = defaultdict(lambda: [])
cell_type_ids = defaultdict(lambda: [])
genes = None
normalize = False
for i, row in annotations_facs.iterrows():
    tissue = row.tissue
    if tissue != current_tissue:
        dirname = tissue
        current_tissue = tissue
        print(tissue)
        current_matrix = pd.read_csv(os.path.join('tabula_muris', 'FACS', dirname + '-counts.csv'))
        if genes is None:
            genes = current_matrix.iloc[:,0].values.astype(str)
    cell_id = row.cell
    cell_type = row.cell_ontology_class
    index = row.cell
    values = current_matrix[index].values
    if normalize:
        values = values/values.sum()
    cell_type_values[cell_type].append(values)
    tissue_values[tissue].append(values)
    tissue_cell_types[tissue].append(cell_type)
    cell_type_ids[cell_type].append(cell_id)

np.savetxt('tabula_muris_facs_genes.txt', genes, fmt='%s')

# mapping of cell type name to cell ontology id
cell_name_to_ontology_id = {}
for index, row in annotations_facs.iterrows():
    cell_type = row.cell_ontology_class
    cell_id = row.cell_ontology_id
    cell_name_to_ontology_id[cell_type] = cell_id

# calculate cell type means
cell_type_means = {}
for cell_type, values in cell_type_values.items():
    values_mean = np.zeros(values[0].shape[0])
    for v in values:
        values_mean += v.flatten()
    values_mean /= len(values)
    cell_type_means[cell_type] = values_mean

# calculate cell type medians
cell_type_medians = {}
for cell_type, values in cell_type_values.items():
    v_matrix = np.vstack([x.flatten() for x in values])
    # do a median calculation here
    cell_type_medians[cell_type] = np.median(v_matrix, 0)

# dump as pickle
#with open('cell_type_means_dict_tm_facs.pkl', 'wb') as f:
#    pickle.dump(cell_type_means, f)

#with open('cell_type_ontology_ids_tm_facs.pkl', 'wb') as f:
#    pickle.dump(cell_name_to_ontology_id, f)

# convert to h5 using pytables
#from uncurl_analysis import dense_matrix_h5
#dense_matrix_h5.store_dict('cell_type_means_tm_facs.h5', cell_type_means)
#dense_matrix_h5.store_dict('cell_type_medians_tm_facs.h5', cell_type_medians)

# save the whole matrix of cell_type_values.
import subprocess
os.makedirs('cell_type_matrices_tm_facs', exist_ok=True)
for cell_type, values in cell_type_values.items():
    # this is of shape cells x genes
    v_matrix = np.vstack([x.flatten() for x in values])
    v_matrix = sparse.csc_matrix(v_matrix)
    scipy.io.mmwrite('cell_type_matrices_tm_facs/{0}.mtx'.format(cell_type), v_matrix)
    subprocess.call(['gzip', 'cell_type_matrices_tm_facs/{0}.mtx'.format(cell_type)])
    np.savetxt('cell_type_matrices_tm_facs/{0}_barcodes.txt'.format(cell_type), np.array(cell_type_ids[cell_type]), fmt='%s')

# TODO: save matrices of tissues
os.makedirs('tissue_type_matrices_tm_facs', exist_ok=True)
for tissue_type, values in tissue_values.items():
    # this is of shape cells x genes
    v_matrix = np.vstack([x.flatten() for x in values])
    v_matrix = sparse.csc_matrix(v_matrix)
    scipy.io.mmwrite('tissue_type_matrices_tm_facs/{0}.mtx'.format(tissue_type), v_matrix)
    subprocess.call(['gzip', 'tissue_type_matrices_tm_facs/{0}.mtx'.format(tissue_type)])
    np.savetxt('tissue_type_matrices_tm_facs/{0}_cell_types.txt'.format(tissue_type), tissue_cell_types[tissue_type], fmt='%s')


