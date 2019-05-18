import numpy as np
import pandas as pd
from scipy import sparse
import scipy.io

data = pd.read_csv('allen_brain_atlas/GSE115746_cells_exon_counts.csv.gz')
gene_names = data.iloc[:,0].values
cell_ids = data.columns.values

np.savetxt('allen_brain_atlas/cell_ids.txt', cell_ids)
np.savetxt('allen_brain_atlas/gene_names.txt', gene_names)

data_matrix = data.values[:,1:]
data_sparse = sparse.csc_matrix(data_matrix)
scipy.io.mmwrite('allen_brain_atlas/data_sparse.mtx', data_sparse)


# load annotations
from collections import defaultdict
annotations = pd.read_csv('allen_brain_atlas/Supplementary_Table_10_Full_Metadata.csv')
filtered_data = ['cell_id', 'class', 'subclass', 'cluster']
class_values = defaultdict(lambda: [])
class_subclass_values = defaultdict(lambda: [])
class_cluster_values = defaultdict(lambda: [])
for i, row in annotations.iterrows():
    cell_id = row['sample_name']
    cell_class = row['class']
    cell_subclass = row['subclass']
    cluster = row['cluster']
    if cell_class == 'Low Quality' or cell_class == 'No Class':
        continue
    try:
        cell_index = np.where(cell_ids == cell_id)[0][0]
        data_val = data_matrix[:, cell_index]
        class_values[cell_class].append(data_val)
        class_subclass_values[cell_class + ' - ' + cell_subclass].append(data_val)
        class_cluster_values[cell_class + ' - ' + cluster].append(data_val)
    except:
        continue

# calculate means
class_means = {k: np.sum(v, 0)/len(v) for k, v in class_values.items()}
class_subclass_means = {k: np.sum(v, 0)/len(v) for k, v in class_subclass_values.items()}
class_cluster_means = {k: np.sum(v, 0)/len(v) for k, v in class_cluster_values.items()}

from uncurl_analysis import dense_matrix_h5

dense_matrix_h5.save_dict('allen_class_means.h5', class_means)
dense_matrix_h5.save_dict('allen_class_subclass_means.h5', class_subclass_means)
dense_matrix_h5.save_dict('allen_class_cluster_means.h5', class_cluster_means)
