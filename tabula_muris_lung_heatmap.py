import pickle
with open('tm_droplet_lung_cellmesh_query_top_cells.pkl', 'rb') as f:
    label_cell_types = pickle.load(f)

# pick a subset of cell types
# include all cell types except for nan
cells_to_include = [
 'ciliated columnar cell of tracheobronchial tree',
 'type II pneumocyte',
 'lung endothelial cell',
 'stromal cell',
 'B cell',
 'T cell',
 'natural killer cell',
 'alveolar macrophage',
 'non-classical monocyte',
 'classical monocyte',
 'mast cell',
 'myeloid cell',
 'leukocyte',
]
# load the hand-created mappings file
cell_types_map = {}
cell_types_alternate_map = {}
with open('cell_ontology_to_cellmesh_tabula_muris.tsv') as f:
    l0 = f.readline()
    for line in f.readlines():
        line_data = [x.strip() for x in line.split('\t')]
        name = line_data[1]
        primary_cellmesh_name = line_data[2]
        alternate_cellmesh_names = line_data[3:]
        use_cellmesh = line_data[0]
        if use_cellmesh != 'n':
            cell_types_map[name] = primary_cellmesh_name
            cell_types_alternate_map[name] = alternate_cellmesh_names
cell_types_map['alveolar macrophage'] = 'Macrophages'
cell_types_map['stromal cell'] = 'Fibroblasts'
correct_cellmesh_results = []
for x in cells_to_include:
    correct_name = cell_types_map[x]
    if correct_name not in correct_cellmesh_results:
        correct_cellmesh_results.append(correct_name)
cellmesh_to_index = {c: i for i, c in enumerate(correct_cellmesh_results)}


# create a heatmap
import numpy as np
cell_type_results = np.zeros((len(cells_to_include), len(correct_cellmesh_results)))
for i, c in enumerate(cells_to_include):
    result_cells = label_cell_types[(c, 50, 'ratio', 'prob')]
    for i2, c2 in enumerate(result_cells):
        for x in correct_cellmesh_results:
            if c2 == x:
                index = cellmesh_to_index[x]
                if cell_type_results[i, index] == 0:
                    rr = 1.0/(i2+1)
                    cell_type_results[i, index] = rr
                break

"""
from sklearn.cluster.bicluster import SpectralCoclustering
spec = SpectralCoclustering(10)
spec.fit(cell_type_results + 1e-8)
row_labels = spec.row_labels_
column_labels = spec.column_labels_
row_order = np.argsort(row_labels)
col_order = np.argsort(column_labels)
cell_type_results_reordered = cell_type_results[row_order, :]
cell_type_results_reordered = cell_type_results_reordered[:, col_order]
# TODO: save reordered columns?
tm_labels = np.array(cells_to_include)[row_order]
cellmesh_labels = np.array(correct_cellmesh_results)[col_order]
"""

cell_type_results_reordered = cell_type_results
tm_labels = cells_to_include
cellmesh_labels = correct_cellmesh_results

import seaborn as sns
import matplotlib.pyplot as plt
plt.figure(figsize=(15,10))
plt.tight_layout()
plt.cla()
plt.clf()
plt.gcf().subplots_adjust(bottom=0.3)
plt.gcf().subplots_adjust(left=0.5)
sns.set(style='whitegrid', font_scale=1.5)
sns.heatmap(cell_type_results_reordered, xticklabels=cellmesh_labels, yticklabels=tm_labels)
plt.title('CellMesh results for lung cells')
plt.xlabel('CellMesh terms')
plt.ylabel('Tabula Muris cell type')
plt.savefig('heatmap_tm_droplet_lung.png', dpi=100)

