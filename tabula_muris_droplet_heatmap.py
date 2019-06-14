# code to generate the heatmaps
# requires output of tabula_muris_droplet_cellmesh_query.py

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

map_method_means = pd.read_csv('MAP_method_means.csv')
map_cell_types = pd.read_csv('MAP_cell_types.csv')
# load scQuery results
scQuery_results = []
with open('scQuery_top_gene_results.csv') as f:
    for i, line in enumerate(f.readlines()):
        line_data = line.split(',')[:2]
        if i > 0:
            line_data[1] = float(line_data[1])
        scQuery_results.append(line_data)
scQuery_results = pd.DataFrame(scQuery_results[1:], columns=['cell_type', 'mean_average_precision'])
scQuery_results['n_genes'] = 50
scQuery_results['query_method'] = 'scQuery'
scQuery_results['gene_method'] = 'ratio'

all_cell_types = sorted(list(set(map_cell_types.cell_type)))
map_cell_types = map_cell_types.append(scQuery_results)

import pickle
with open('tm_droplet_cellmesh_query_top_cells.pkl', 'rb') as f:
    label_cell_types = pickle.load(f)

# pick a subset of cell types
# pick cells within a group?
# immune cells: dendritic cell, B cell, natural killer cell, mast cell, T cell, macrophage, monocyte, basophil
cells_to_include_immune = [
        'natural killer cell',
        'T cell',
        'B cell',
        'mast cell',
        'basophil',
        'monocyte',
        'dendritic cell',
        'macrophage',
]

#cells_to_include = cells_to_include_immune
cells_to_include = all_cell_types

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
correct_cellmesh_results = list(set([cell_types_map[x] for x in cells_to_include]))
cellmesh_to_index = {c: i for i, c in enumerate(correct_cellmesh_results)}
actual_cellmesh_results = [label_cell_types[(c, 50, 'ratio', 'cellmesh')][0] for c in cells_to_include]


# create a heatmap
import numpy as np
cell_type_results = np.zeros((len(cells_to_include), len(correct_cellmesh_results)))
for i, c in enumerate(cells_to_include):
    result_cells = label_cell_types[(c, 50, 'ratio', 'cellmesh')]
    for i2, c2 in enumerate(result_cells):
        for x in correct_cellmesh_results:
            if c2 == x:
                index = cellmesh_to_index[x]
                if cell_type_results[i, index] == 0:
                    rr = 1.0/(i2+1)
                    cell_type_results[i, index] = rr

# TODO: do co-clustering on cell_type_results 
from sklearn.cluster.bicluster import SpectralCoclustering
spec = SpectralCoclustering(10)
spec.fit(cell_type_results + 1e-8)
row_labels = spec.row_labels_
column_labels = spec.column_labels_

row_order = np.argsort(row_labels)
col_order = np.argsort(column_labels)
cell_type_results_reordered = cell_type_results[row_order, :]
cell_type_results_reordered = cell_type_results_reordered[:, col_order]
tm_labels = np.array(cells_to_include)[row_order]
cellmesh_labels = np.array(correct_cellmesh_results)[col_order]


plt.figure(figsize=(51,51))
plt.tight_layout()
plt.cla()
plt.clf()
plt.gcf().subplots_adjust(bottom=0.2)
plt.gcf().subplots_adjust(left=0.2)
sns.set(style='whitegrid', font_scale=2.0)
sns.heatmap(cell_type_results_reordered, xticklabels=cellmesh_labels, yticklabels=tm_labels, linewidths=.5)
plt.title('CellMesh results for immune cells')
plt.xlabel('CellMesh terms')
plt.ylabel('Tabula Muris cell type')
plt.savefig('heatmap_tm_droplet_all_cell_types.png', dpi=100)

