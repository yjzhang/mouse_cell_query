# TODO: preprocess the tabula muris data into a format that can be queried upon.
# tabula_muris/droplet.zip

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
tissue_values = {}
cell_type_values = {}
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
        current_barcodes = pd.read_csv(os.path.join('tabula_muris', 'droplet', dirname, 'barcodes.tsv'), header=None)
    cell_barcode = row.cell[-16:]+'-1'
    index = np.where(current_barcodes[0] == cell_barcode)[0][0]
    cell_type = row.cell_ontology_class
    if cell_type not in cell_type_values:
        cell_type_values[cell_type] = []
    cell_type_values[cell_type].append(current_matrix[:, index])

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

# dump as pickle
import pickle
with open('cell_type_means_dict.pkl', 'wb') as f:
    pickle.dump(cell_type_means, f)

with open('cell_type_ontology_ids.pkl', 'wb') as f:
    pickle.dump(cell_name_to_ontology_id, f)

# convert to h5 using pytables
from uncurl_analysis import dense_matrix_h5
dense_matrix_h5.store_dict('cell_type_means.h5', cell_type_means)

# get gene names
# it appears that the gene names are constant across all subsets
gene_names = pd.read_table('tabula_muris/droplet/Bladder-10X_P4_3/genes.tsv', header=None)
genes = gene_names[0]
np.savetxt('gene_names.txt', genes, fmt='%s')
