# data source: https://figshare.com/articles/MCA_DGE_Data/5435866

import gc
import os
import subprocess
import time

import numpy as np
import pandas as pd
import scipy.io
from scipy import sparse

# TODO: make it so all datasets have the same genes?

class GeneNameTransform(object):

    def __init__(self, old_gene_names, new_gene_names, fixed_new=False):
        """
        Args:
            old_gene_names (array of str)
            new_gene_names (array of str)
        """
        self.names_set = set(old_gene_names)
        self.names_to_add = list(set(new_gene_names).difference(old_gene_names))
        if len(self.names_to_add) == 0:
            print('No new genes to add')
        if fixed_new and len(self.names_set.difference(new_gene_names)) == 0:
            # if there are no genes that are in old_gene_names but not in new_gene_names
            print('using new gene names')
            self.names_list = np.array(new_gene_names)
        else:
            self.names_list = np.hstack([old_gene_names, self.names_to_add])
        old_indices = {x: i for i, x in enumerate(old_gene_names)}
        new_indices = {x: i for i, x in enumerate(new_gene_names)}
        self.old_to_unified = np.array([old_indices[x] if x in old_indices else -1 for x in self.names_list])
        self.old_zero_mask = np.array([i for i, x in enumerate(self.names_list) if x not in old_indices])
        self.new_to_unified = np.array([new_indices[x] if x in new_indices else -1 for x in self.names_list])
        self.new_zero_mask = np.array([i for i, x in enumerate(self.names_list) if x not in new_indices])

    def transform_old(self, old_data):
        """
        Transforms a dataset of the "old" gene names to the unified form.

        old_data: 1d array or 2d array. If 2d array, assume that genes are in the first dimension.
        """
        if len(self.names_to_add) == 0:
            return old_data
        # if data is a sparse matrix... it might still work?
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

        new_data: 1d or 2d array. if 2d array, assume that genes are in the 1st dimension.
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

def save_cell_types(cell_type_value, gene_names=None):
    os.makedirs('cell_type_matrices_mca', exist_ok=True)
    for cell_type, values in cell_type_values.items():
        # this is of shape cells x genes
        v_matrix = sparse.hstack(values).T
        v_matrix = sparse.csc_matrix(v_matrix)
        scipy.io.mmwrite('cell_type_matrices_mca/{0}.mtx'.format(cell_type), v_matrix)
        subprocess.call(['gzip', 'cell_type_matrices_mca/{0}.mtx'.format(cell_type)])
        if gene_names is not None:
            np.savetxt('cell_type_matrices_mca/{0}_genes.txt'.format(cell_type), gene_names, fmt='%s')


def save_cell_type(cell_type, values, gene_names=None):
    os.makedirs('cell_type_matrices_mca', exist_ok=True)
    v_matrix = sparse.hstack(values).T
    v_matrix = sparse.csc_matrix(v_matrix)
    scipy.io.mmwrite('cell_type_matrices_mca/{0}.mtx'.format(cell_type), v_matrix)
    subprocess.call(['gzip', '-f', 'cell_type_matrices_mca/{0}.mtx'.format(cell_type)])
    if gene_names is not None:
        np.savetxt('cell_type_matrices_mca/{0}_genes.txt'.format(cell_type), gene_names, fmt='%s')


annotations = pd.read_csv('mca_microwell_seq/MCA_CellAssignments.csv')

current_batch = None
current_tissue = None
current_table = None
current_matrix = None
current_barcodes = None
current_genes = None
prev_cell_type = None

visited_batches = {}
remaining_batches = {}
# TODO: keep track of all cell types that have files still remaining in the dataset...
# so that we can clean memory

cell_type_values = {}
cell_type_ids = {}

unified_gene_list = None

for cell_index, row in annotations.iterrows():
    # TODO: save barcodes???
    batch = row.Batch
    tissue = row.Tissue
    if batch != current_batch:
        # TODO: save cell types, re-make cell_type_values? but only after a certain point;
        # only after most of the cells have been processed. and then, save a copy of all of the gene name lists?
        print(cell_index)
        remaining_cell_types = set(annotations['Annotation'].values[cell_index+1:])
        cell_types_to_remove = []
        for cell_type, values in cell_type_values.items():
            print(cell_type, len(values))
            if cell_type not in remaining_cell_types:
                print('done with cell type ' + cell_type + ', writing out matrix...')
                save_cell_type(cell_type, values, unified_gene_list)
                cell_types_to_remove.append(cell_type)
        for cell_type in cell_types_to_remove:
            del cell_type_values[cell_type]
        del current_matrix
        gc.collect()
        gene_mappers = {}
        dirname = ''.join(batch.split('_'))
        if dirname == 'EmbryonicMesenchyme1':
            dirname = 'EmbryonicMesenchymeE14.5'
        elif dirname == 'FetalLiver1':
            dirname = 'FetalLiverE14.1'
        elif dirname == 'Male(fetal)Gonad1':
            dirname = 'FetalMaleGonad'
        elif dirname == 'NeonatalBrain1':
            # ugh wtf
            dirname = 'NeontalBrain1'
        elif dirname == 'NeonatalBrain2':
            dirname = 'NeontalBrain2'
        elif dirname == 'Placenta1':
            dirname = 'PlacentaE14.1'
        elif dirname == 'Placenta2':
            dirname = 'PlacentaE14.2'
        path = os.path.join('mca_microwell_seq', 'rmbatch_dge', '{0}_rm.batch_dge.txt.gz'.format(dirname))
        if not os.path.exists(path):
            dirname = batch.split('_')[0]
            path = os.path.join('mca_microwell_seq', 'rmbatch_dge', '{0}_rm.batch_dge.txt.gz'.format(dirname))
        if not os.path.exists(path):
            print('data not available:', dirname, path)
            continue
        current_tissue = tissue
        current_batch = batch
        print(batch)
        t0 = time.time()
        print('loading data...')
        current_table = pd.read_table(path, sep=' ')
        current_barcodes = np.array([x.split('.')[1] for x in current_table.columns])
        current_genes = current_table.index.values.astype(str)
        current_matrix = current_table.to_numpy(dtype=np.int16, copy=True)
        del current_table
        gc.collect()
        print('finished loading data: {0}'.format(time.time() - t0))
        # transform current_matrix...
        if unified_gene_list is None or len(unified_gene_list) != len(current_genes) or (unified_gene_list != current_genes).any():
            t0 = time.time()
            print('merging gene names...')
            # concatenate the matrices to save time
            cell_type_values = {k: [sparse.hstack(v)] for k, v in cell_type_values.items()}
            if unified_gene_list is None:
                unified_gene_list = current_genes
            else:
                gene_name_mapper = GeneNameTransform(unified_gene_list, current_genes)
                cell_type_values = {k: gene_name_mapper.transform_old(v) for k, v in cell_type_values.items()}
                # transform new
                current_matrix = gene_name_mapper.transform_new(current_matrix)
                unified_gene_list = np.array(gene_name_mapper.names_list)
            print('finished merging gene names: {0}'.format(time.time() - t0))
        current_matrix = sparse.csc_matrix(current_matrix)
        gc.collect()
    cell_barcode = row['Cell.Barcode']
    cell_type = row['Annotation']
    #try:
    #    index = np.where(current_barcodes == cell_barcode)[0][0]
    #except:
    #    print('Cell ID not found:', batch, cell_index, cell_barcode, cell_type)
    #    continue
    cell_type = cell_type.replace('/', '-')
    #if cell_type not in cell_type_values:
    #    cell_type_values[cell_type] = []
    #cell_type_values[cell_type].append(current_matrix[:, index])
    if cell_type not in cell_type_ids:
        cell_type_ids[cell_type] = []
    cell_type_ids[cell_type].append(row['Cell.name'])
    print(cell_index, cell_barcode, cell_type)
    prev_cell_type = cell_type


for cell_type, ids in cell_type_ids.items():
    print(cell_type, len(ids))
    print('done with cell type ' + cell_type + ', writing out cell ids...')
    np.savetxt('cell_type_matrices_mca/{0}_barcodes.txt'.format(cell_type), np.array(ids), fmt='%s')

for cell_type, values in cell_type_values.items():
    print(cell_type, len(values))
    print('done with cell type ' + cell_type + ', writing out matrix...')
    save_cell_type(cell_type, values, unified_gene_list)

np.savetxt('genes_mca.txt', unified_gene_list, fmt='%s')

"""
import pickle
with open('cell_type_values_mca_dict.pkl', 'wb') as f:
    pickle.dump(cell_type_values, f)

# calculate cell type means
from uncurl_analysis import dense_matrix_h5
cell_type_means = {}
for cell_type, values in cell_type_values.items():
    print(cell_type)
    v_matrix = np.hstack(values)
    values_mean = v_matrix.mean(1)
    cell_type_means[cell_type] = np.array(values_mean).flatten()
dense_matrix_h5.store_dict('cell_type_means_mca.h5', cell_type_means)

# calculate cell type medians
cell_type_medians = {}
for cell_type, values in cell_type_values.items():
    v_matrix = np.hstack(values).T.toarray()
    # do a median calculation here
    cell_type_medians[cell_type] = np.median(v_matrix, 0)

# dump as pickle

# convert to h5 using pytables
dense_matrix_h5.store_dict('cell_type_medians_mca.h5', cell_type_medians)

# save the whole matrix of cell_type_values.
os.makedirs('cell_type_matrices_mca', exist_ok=True)
for cell_type, values in cell_type_values.items():
    # this is of shape cells x genes
    v_matrix = sparse.hstack(values).T
    v_matrix = sparse.csc_matrix(v_matrix)
    scipy.io.mmwrite('cell_type_matrices_mca/{0}.mtx'.format(cell_type), v_matrix)
    subprocess.call(['gzip', 'cell_type_matrices_mca/{0}.mtx'.format(cell_type)])
"""

# TODO: unify gene names across datasets (again)
unified_gene_list = np.loadtxt('genes_mca.txt', dtype=str)
for filename in os.listdir('cell_type_matrices_mca'):
    if filename.endswith('.mtx.gz'):
        base_name = filename[:-7]
        gene_name = base_name + '_genes.txt'
        data_path = os.path.join('cell_type_matrices_mca', filename)
        genes_path = os.path.join('cell_type_matrices_mca', gene_name)
        if os.path.exists(genes_path):
            print(base_name)
            # TODO: load dataset
            data = scipy.io.mmread(data_path)
            data = data.T
            data = sparse.csc_matrix(data)
            current_genes = np.loadtxt(genes_path, dtype=str)
            gene_name_mapper = GeneNameTransform(current_genes, unified_gene_list, fixed_new=True)
            data = gene_name_mapper.transform_old(data)
            scipy.io.mmwrite(data_path[:-3], data)
            subprocess.call(['gzip', '-f', data_path])
            np.savetxt(genes_path, np.array(gene_name_mapper.names_list), fmt='%s')
