# loads the processed mca data, combines the same cell types from different tissues?
import numpy as np

from uncurl_analysis import dense_matrix_h5

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
            else:
                new_data = new_data[self.new_to_unified]
            new_data = new_data[self.new_to_unified]
            if len(self.new_zero_mask) > 0:
                new_data[self.new_zero_mask] = 0
            return new_data


gene_names_original = dense_matrix_h5.H5Dict('mouse_cell_query/data/cell_type_genes_mca.h5')
cell_type_means = dense_matrix_h5.H5Dict('mouse_cell_query/data/cell_type_means_mca.h5')

cell_types = cell_type_means.keys()

# sort cell types into coarser classes - merge tissues
cell_type_classes = set()
cell_type_classes_map = {}

for cell_type in cell_types:
    x = cell_type.split('(')[0]
    if x not in cell_type_classes:
        cell_type_classes.add(x)
    cell_type_classes_map[cell_type] = x

cell_type_classes_coarse = set()
cell_type_classes_map_coarse = {}

for cell_type in cell_types:
    cl0 = cell_type.split('(')[0].split('_')[0].strip()
    x = cl0.capitalize().strip()
    if cl0.endswith('s'):
        x = x[:-1]
    if x not in cell_type_classes_coarse:
        cell_type_classes_coarse.add(x)
    cell_type_classes_map_coarse[cell_type] = x

import os
import scipy.io
from scipy import sparse

base_dir = 'cell_type_matrices_mca'
new_cell_types_dict = {x: [] for x in cell_type_classes_coarse}
unified_gene_list = None
gene_name_mapper = None
for filename in os.listdir(base_dir):
    data = sparse.csc_matrix(scipy.io.mmread(os.path.join(base_dir, filename)).T)
    cell_type = filename.split('.')[0]
    gene_names = gene_names_original[cell_type].astype(str)
    coarse_cell_type = cell_type_classes_map_coarse[cell_type]
    print(cell_type, coarse_cell_type)
    # TODO: merge gene names
    if unified_gene_list is None or len(gene_names) != len(unified_gene_list) or (gene_names != unified_gene_list).any():
        if unified_gene_list is None:
            unified_gene_list = gene_names
        else:
            gene_name_mapper = GeneNameTransform(unified_gene_list, gene_names)
            # TODO: transform gene names of everything in set
            new_cell_types_dict = {k: gene_name_mapper.transform_old(v) for k, v in new_cell_types_dict.items()}
            # transform new
            data = gene_name_mapper.transform_new(data)
            unified_gene_list = np.array(gene_name_mapper.names_list)
    new_cell_types_dict[coarse_cell_type].append(data)

# calculate means
new_cell_types_dict = {k : sparse.hstack(v) for k, v in new_cell_types_dict.items()}
normalize = False
if normalize:
    new_cell_types_dict = {k : v/v.sum(1) for k, v in new_cell_types_dict.items()}
new_cell_types_means = {k : np.array(v.mean(1)).flatten() for k, v in new_cell_types_dict.items()}

dense_matrix_h5.store_dict('cell_type_means_mca_coarse_normalized.h5', new_cell_types_means)
np.savetxt('mca_coarse_gene_names.txt', unified_gene_list, fmt='%s')

import subprocess
os.makedirs('cell_type_matrices_mca_coarse', exist_ok=True)
for cell_type, v_matrix in new_cell_types_dict.items():
    # this is of shape cells x genes
    scipy.io.mmwrite('cell_type_matrices_mca_coarse/{0}.mtx'.format(cell_type), v_matrix)
    subprocess.call(['gzip', 'cell_type_matrices_mca_coarse/{0}.mtx'.format(cell_type)])

