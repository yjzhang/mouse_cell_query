# TODO: preprocess the tabula muris data into a format that can be queried upon.
# data source: https://figshare.com/articles/Single-cell_RNA-seq_data_from_microfluidic_emulsion_v2_/5968960
# tabula_muris/droplet.zip

from collections import defaultdict
import os

import numpy as np
import pandas as pd
import scipy.io
from scipy import sparse

metadata_droplet = pd.read_csv('tabula_muris/metadata_droplet.csv')
annotations_droplet = pd.read_csv('tabula_muris/annotations_droplet.csv')

current_channel = None
current_tissue = None
current_matrix = None
current_barcodes = None
tissue_values = defaultdict(lambda: [])
cell_type_values = defaultdict(lambda: [])
tissue_cell_types = defaultdict(lambda: [])
normalize = False
for index, row in annotations_droplet.iterrows():
    channel = row.channel
    tissue = row.tissue
    if channel != current_channel or tissue != current_tissue:
        dirname = tissue + '-' + channel
        current_channel = channel
        current_tissue = tissue
        print(tissue, channel)
        current_matrix = scipy.io.mmread(os.path.join('tabula_muris', 'droplet', dirname, 'matrix.mtx'))
        current_matrix = sparse.csc_matrix(current_matrix)
        # TODO: preprocess current_matrix?
        current_barcodes = pd.read_csv(os.path.join('tabula_muris', 'droplet', dirname, 'barcodes.tsv'), header=None)
    cell_barcode = row.cell[-16:]+'-1'
    index = np.where(current_barcodes[0] == cell_barcode)[0][0]
    cell_type = row.cell_ontology_class
    values = current_matrix[:, index]
    if normalize:
        values = values/values.sum()
    cell_type_values[cell_type].append(values)
    tissue_values[tissue].append(values)
    tissue_cell_types[tissue].append(cell_type)

# mapping of cell type name to cell ontology id
cell_name_to_ontology_id = {}
for index, row in annotations_droplet.iterrows():
    cell_type = row.cell_ontology_class
    cell_id = row.cell_ontology_id
    cell_name_to_ontology_id[cell_type] = cell_id

# calculate cell type means
cell_type_means = {}
for cell_type, values in cell_type_values.items():
    values_mean = np.zeros(values[0].shape[0])
    for v in values:
        values_mean += v.toarray().flatten()
    values_mean /= len(values)
    cell_type_means[cell_type] = values_mean

# calculate cell type medians
cell_type_medians = {}
for cell_type, values in cell_type_values.items():
    v_matrix = np.vstack([x.toarray().flatten() for x in values])
    # do a median calculation here
    cell_type_medians[cell_type] = np.median(v_matrix, 0)

# dump as pickle
import pickle
with open('cell_type_means_dict.pkl', 'wb') as f:
    pickle.dump(cell_type_means, f)

with open('cell_type_ontology_ids.pkl', 'wb') as f:
    pickle.dump(cell_name_to_ontology_id, f)

# convert to h5 using pytables
from uncurl_analysis import dense_matrix_h5
dense_matrix_h5.store_dict('cell_type_means_tm_droplet.h5', cell_type_means)
dense_matrix_h5.store_dict('cell_type_medians_tm_droplet.h5', cell_type_medians)

# save the whole matrix of cell_type_values.
import subprocess
os.makedirs('cell_type_matrices_tm_droplet', exist_ok=True)
for cell_type, values in cell_type_values.items():
    # this is of shape cells x genes
    v_matrix = np.vstack([x.toarray().flatten() for x in values])
    v_matrix = sparse.csc_matrix(v_matrix)
    scipy.io.mmwrite('cell_type_matrices_tm_droplet/{0}.mtx'.format(cell_type), v_matrix)
    subprocess.call(['gzip', 'cell_type_matrices_tm_droplet/{0}.mtx'.format(cell_type)])
# TODO: save by tissue type as well
os.makedirs('tissue_cell_type_matrices_tm_droplet', exist_ok=True)
for tissue, values in tissue_values.items():
    # this is of shape cells x genes
    v_matrix = np.vstack([x.toarray().flatten() for x in values])
    v_matrix = sparse.csc_matrix(v_matrix)
    scipy.io.mmwrite('tissue_cell_type_matrices_tm_droplet/{0}.mtx'.format(tissue), v_matrix)
    subprocess.call(['gzip', 'tissue_cell_type_matrices_tm_droplet/{0}.mtx'.format(tissue)])
    np.savetxt('tissue_cell_type_matrices_tm_droplet/{0}_cell_types.txt'.format(tissue), tissue_cell_types[tissue], fmt='%s')


# get gene names
# it appears that the gene names are constant across all subsets
gene_names = pd.read_table('tabula_muris/droplet/Bladder-10X_P4_3/genes.tsv', header=None)
genes = gene_names[0]
np.savetxt('gene_names.txt', genes, fmt='%s')
