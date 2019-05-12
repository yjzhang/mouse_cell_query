# data source: https://figshare.com/articles/MCA_DGE_Data/5435866

import os

import numpy as np
import pandas as pd
import scipy.io
from scipy import sparse

annotations = pd.read_csv('mca_microwell_seq/MCA_CellAssignments.csv')

current_batch = None
current_tissue = None
current_matrix = None
current_barcodes = None
current_genes = None
cell_type_values = {}
cell_type_genes = {}
# TODO: try to align gene names across dataset?
for index, row in annotations.iterrows():
    try:
        batch = row.Batch
        tissue = row.Tissue
        if batch != current_batch:
            dirname = ''.join(batch.split('_'))
            current_tissue = tissue
            current_batch = batch
            print(batch)
            try:
                current_matrix = pd.read_table(os.path.join('mca_microwell_seq', 'rmbatch_dge', '{0}_rm.batch_dge.txt.gz'.format(dirname)), sep=' ')
            except:
                dirname = batch.split('_')[0]
                current_matrix = pd.read_table(os.path.join('mca_microwell_seq', 'rmbatch_dge', '{0}_rm.batch_dge.txt.gz'.format(dirname)), sep=' ')
            current_barcodes = np.array([x.split('.')[1] for x in current_matrix.columns])
            current_genes = current_matrix.index.values.astype(str)
            current_matrix = current_matrix.values
        cell_barcode = row['Cell.Barcode']
        index = np.where(current_barcodes == cell_barcode)[0][0]
        cell_type = row['Annotation']
        if cell_type not in cell_type_values:
            cell_type = cell_type.replace('/', '-')
            cell_type_values[cell_type] = []
            cell_type_genes[cell_type] = current_genes
        cell_type_values[cell_type].append(current_matrix[:, index])
    except:
        continue

import pickle
with open('cell_type_values_mca_dict.pkl', 'wb') as f:
    pickle.dump(cell_type_values, f)

with open('cell_type_gene_names_mca_dict.pkl', 'wb') as f:
    pickle.dump(cell_type_genes, f)



# calculate cell type means
from uncurl_analysis import dense_matrix_h5
dense_matrix_h5.store_dict('cell_type_genes_mca.h5', cell_type_genes)
cell_type_means = {}
for cell_type, values in cell_type_values.items():
    print(cell_type)
    values_mean = np.zeros(values[0].shape[0])
    for v in values:
        try:
            values_mean += v.flatten()
        except:
            continue
    values_mean /= len(values)
    cell_type_means[cell_type] = values_mean
dense_matrix_h5.store_dict('cell_type_means_mca.h5', cell_type_means)

# calculate cell type medians
cell_type_medians = {}
for cell_type, values in cell_type_values.items():
    v_matrix = np.vstack([x.flatten() for x in values])
    # do a median calculation here
    cell_type_medians[cell_type] = np.median(v_matrix, 0)

# dump as pickle

# convert to h5 using pytables
dense_matrix_h5.store_dict('cell_type_medians_mca.h5', cell_type_medians)

# save the whole matrix of cell_type_values.
import subprocess
os.makedirs('cell_type_matrices_mca', exist_ok=True)
for cell_type, values in cell_type_values.items():
    # this is of shape cells x genes
    v_matrix = np.vstack([x.toarray().flatten() for x in values])
    v_matrix = sparse.csc_matrix(v_matrix)
    scipy.io.mmwrite('cell_type_matrices_mca/{0}.mtx'.format(cell_type), v_matrix)
    subprocess.call(['gzip', 'cell_type_matrices_mca/{0}.mtx'.format(cell_type)])

