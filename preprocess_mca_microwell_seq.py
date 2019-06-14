# data source: https://figshare.com/articles/MCA_DGE_Data/5435866

import os

import numpy as np
import pandas as pd
import scipy.io
from scipy import sparse

# TODO: make it so all datasets have the same genes?

class GeneNameTransform(object):

    def __init__(self, old_gene_names, new_gene_names):
        """
        Args:
            old_gene_names (array of str)
            new_gene_names (array of str)
        """
        self.names_set = set(old_gene_names)
        self.names_set.update(new_gene_names)
        self.names_list = list(self.names_set)
        old_indices = {x: i for i, x in enumerate(old_gene_names)}
        new_indices = {x: i for i, x in enumerate(new_gene_names)}
        self.old_to_unified = np.array([old_indices[x] if x in old_indices else -1 for x in self.names_list])
        self.old_zero_mask = np.array([i for i, x in enumerate(self.names_list) if x not in old_indices])
        self.new_to_unified = np.array([new_indices[x] if x in new_indices else -1 for x in self.names_list])
        self.new_zero_mask = np.array([i for i, x in enumerate(self.names_list) if x not in new_indices])

    def transform_old(self, old_data):
        """
        Transforms a dataset of the "old" gene names to the unified form.

        old_data: 1d array
        """
        # TODO: what if data is a sparse matrix... it might still work?
        if isinstance(old_data, list):
            if len(old_data) == 0:
                return old_data
            results = []
            for data in old_data:
                results.append(self.transform_old(data))
            return results
        else:
            # assume genes is first dimension
            if len(old_data.shape) == 2:
                new_data = old_data[self.old_to_unified, :]
                if len(self.old_zero_mask) > 0:
                    new_data[self.old_zero_mask, :] = 0
            else:
                new_data = old_data[self.old_to_unified]
                if len(self.old_zero_mask) > 0:
                    new_data[self.old_zero_mask] = 0
            return new_data


    def transform_new(self, new_data):
        """
        transforms a dataset of the "new" gene names to the unified form.

        new_data: 1d array
        """
        if isinstance(new_data, list):
            if len(new_data) == 0:
                return new_data
            results = []
            for data in new_data:
                results.append(self.transform_new(data))
            return results
        else:
            if len(new_data.shape) == 2:
                new_data = new_data[self.new_to_unified, :]
                if len(self.new_zero_mask) > 0:
                    new_data[self.new_zero_mask, :] = 0
            else:
                new_data = new_data[self.new_to_unified]
                if len(self.new_zero_mask) > 0:
                    new_data[self.new_zero_mask] = 0
            return new_data


annotations = pd.read_csv('mca_microwell_seq/MCA_CellAssignments.csv')

current_batch = None
current_tissue = None
current_table = None
current_matrix = None
current_barcodes = None
current_genes = None
prev_cell_type = None

cell_type_values = {}

unified_gene_list = None

for index, row in annotations.iterrows():
    batch = row.Batch
    tissue = row.Tissue
    cell_barcode = row['Cell.Barcode']
    cell_type = row['Annotation']
    print(index, cell_barcode, cell_type)
    if batch != current_batch:
        gene_mappers = {}
        dirname = ''.join(batch.split('_'))
        current_tissue = tissue
        current_batch = batch
        print(batch)
        try:
            current_table = pd.read_table(os.path.join('mca_microwell_seq', 'rmbatch_dge', '{0}_rm.batch_dge.txt.gz'.format(dirname)), sep=' ')
        except:
            dirname = batch.split('_')[0]
            current_table = pd.read_table(os.path.join('mca_microwell_seq', 'rmbatch_dge', '{0}_rm.batch_dge.txt.gz'.format(dirname)), sep=' ')
        current_barcodes = np.array([x.split('.')[1] for x in current_table.columns])
        current_genes = current_table.index.values.astype(str)
        current_matrix = current_table.values
        # transform current_matrix...
        if unified_gene_list is None or len(unified_gene_list) != len(current_genes) or (unified_gene_list != current_genes).any():
            if unified_gene_list is None:
                unified_gene_list = current_genes
            else:
                gene_name_mapper = GeneNameTransform(unified_gene_list, current_genes)
                cell_type_values = {k: gene_name_mapper.transform_old(v) for k, v in cell_type_values.items()}
                # transform new
                current_matrix = gene_name_mapper.transform_new(current_matrix)
                unified_gene_list = np.array(gene_name_mapper.names_list)
    try:
        index = np.where(current_barcodes == cell_barcode)[0][0]
    except:
        continue
    if cell_type not in cell_type_values:
        cell_type = cell_type.replace('/', '-')
        cell_type_values[cell_type] = []
    cell_type_values[cell_type].append(current_matrix[:, index])
    prev_cell_type = cell_type

np.savetxt('genes_mca.txt', unified_gene_list, fmt='%s')

import pickle
with open('cell_type_values_mca_dict.pkl', 'wb') as f:
    pickle.dump(cell_type_values, f)

# calculate cell type means
from uncurl_analysis import dense_matrix_h5
cell_type_means = {}
for cell_type, values in cell_type_values.items():
    print(cell_type)
    values_mean = np.zeros(values[0].shape[0])
    for v in values:
        try:
            values_mean += v.flatten()
        except:
            print('error in calculating means')
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
    v_matrix = np.vstack([x.flatten() for x in values])
    v_matrix = sparse.csc_matrix(v_matrix)
    scipy.io.mmwrite('cell_type_matrices_mca/{0}.mtx'.format(cell_type), v_matrix)
    subprocess.call(['gzip', 'cell_type_matrices_mca/{0}.mtx'.format(cell_type)])

