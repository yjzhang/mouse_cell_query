# loads the processed mca data, combines the same cell types from different tissues?

import gc
import os

import numpy as np
import scipy.io
from scipy import sparse


# TODO: load cell types from folder
base_path = 'cell_type_matrices_mca'
cell_types = []
for filename in os.listdir(base_path):
    if filename.endswith('.mtx.gz'):
        base_name = filename[:-7]
        cell_types.append(base_name)

# sort cell types into coarser classes - merge tissues
cell_type_classes = set()
cell_type_classes_map = {}
coarse_to_fine_map = {}

for cell_type in cell_types:
    x = cell_type.split('(')[0]
    if x not in cell_type_classes:
        cell_type_classes.add(x)
        coarse_to_fine_map[x] = []
    cell_type_classes_map[cell_type] = x
    coarse_to_fine_map[x].append(cell_type)

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

print(cell_type_classes_coarse)

import subprocess
base_dir = 'cell_type_matrices_mca'
new_cell_types_dict = {x: [] for x in cell_type_classes_coarse}
os.makedirs('cell_type_matrices_mca_coarse', exist_ok=True)
for coarse_cell_type in cell_type_classes_coarse:
    print('Processing: ', coarse_cell_type)
    for filename in os.listdir(base_dir):
        if filename.endswith('.mtx.gz'):
            cell_type = filename.split('.')[0]
            if cell_type_classes_map_coarse[cell_type] == coarse_cell_type:
                data = sparse.csc_matrix(scipy.io.mmread(os.path.join(base_dir, filename)).T)
                print(cell_type, ' => ', coarse_cell_type)
                new_cell_types_dict[coarse_cell_type].append(data)
                gc.collect()
    v_matrix = sparse.vstack(new_cell_types_dict[coarse_cell_type])
    del new_cell_types_dict[coarse_cell_type]
    #new_cell_types_dict[coarse_cell_type] = v_matrix
    gc.collect()
    scipy.io.mmwrite('cell_type_matrices_mca_coarse/{0}.mtx'.format(coarse_cell_type), v_matrix)
    subprocess.call(['gzip', 'cell_type_matrices_mca_coarse/{0}.mtx'.format(coarse_cell_type)])

print('done processing')

# calculate means
#normalize = False
#if normalize:
#    new_cell_types_dict = {k : v/v.sum(1) for k, v in new_cell_types_dict.items()}
#new_cell_types_means = {k : np.array(v.mean(1)).flatten() for k, v in new_cell_types_dict.items()}

#dense_matrix_h5.store_dict('cell_type_means_mca_coarse_normalized.h5', new_cell_types_means)

