# loads the processed mca data, combines the same cell types from different tissues?
import numpy as np

from uncurl_analysis import dense_matrix_h5


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
for filename in os.listdir(base_dir):
    data = sparse.csc_matrix(scipy.io.mmread(os.path.join(base_dir, filename)).T)
    cell_type = filename.split('.')[0]
    coarse_cell_type = cell_type_classes_map_coarse[cell_type]
    print(cell_type, coarse_cell_type)
    new_cell_types_dict[coarse_cell_type].append(data)

# calculate means
new_cell_types_dict = {k : sparse.hstack(v) for k, v in new_cell_types_dict.items()}
normalize = False
if normalize:
    new_cell_types_dict = {k : v/v.sum(1) for k, v in new_cell_types_dict.items()}
new_cell_types_means = {k : np.array(v.mean(1)).flatten() for k, v in new_cell_types_dict.items()}

#dense_matrix_h5.store_dict('cell_type_means_mca_coarse_normalized.h5', new_cell_types_means)

import subprocess
os.makedirs('cell_type_matrices_mca_coarse', exist_ok=True)
for cell_type, v_matrix in new_cell_types_dict.items():
    # this is of shape cells x genes
    scipy.io.mmwrite('cell_type_matrices_mca_coarse/{0}.mtx'.format(cell_type), v_matrix)
    subprocess.call(['gzip', 'cell_type_matrices_mca_coarse/{0}.mtx'.format(cell_type)])

