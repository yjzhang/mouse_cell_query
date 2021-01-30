import pickle

import numpy as np
import scipy.io

labels = np.loadtxt('spinal_labels_1.txt', delimiter=';', dtype=str)
genes = np.loadtxt('genes.txt', dtype=str)

"""
data = scipy.io.mmread('spinal_cord.mtx.gz')
# TODO: diff exp
from uncurl_analysis import gene_extraction
import time
t0 =  time.time()
scores_t, pvals_t = gene_extraction.one_vs_rest_t(data, labels, eps=0.1, test='t')
print('diffexp time for t test:', time.time() - t0)

with open('tm_droplet_t_scores.pkl', 'wb') as f:
    pickle.dump(scores_t, f)
with open('tm_droplet_t_pvals.pkl', 'wb') as f:
    pickle.dump(pvals_t, f)

"""
############################################################################

# 3. for each cluster, run cellmesh and cellmarker

with open('tm_droplet_t_scores.pkl', 'rb') as f:
    scores_t = pickle.load(f)
with open('tm_droplet_t_pvals.pkl', 'rb') as f:
    pvals_t = pickle.load(f)
# get gene names - save ranked genes for each cell type
top_genes_ratio = {}
for cell, cell_genes in scores_t.items():
    top_genes_ratio[cell] = [genes[i[0]] for i in cell_genes]
import pandas as pd
top_genes_ratio = pd.DataFrame(top_genes_ratio)
top_genes_ratio.to_csv('split_seq_top_genes_ratio.csv')

import cellmesh
from cellmesh import prob_method
import cellmarker
from cellmarker import prob_method as cellmarker_prob_method
labels_set = set(labels)
label_results = {}
label_cell_types = {}
n_genes = [50]#, 200, 1000]
gene_methods = ['ratio']#, 't', 'u']
query_methods = ['cellmesh', 'prob']#, 'meshscan', 'meshscan_prob']#, 'aggregate', 'aggregate_2']
all_species = ['mouse']
# TODO: mesh subsets - None, nervous system, brain + spinal cord, neurons + neuroglia
cell_type_subsets = [None, ('D009420',), ('D013116', 'D001921'), ('D009457', 'D009474')]
all_mesh_cell_id_names = cellmesh.get_all_cell_id_names(include_cell_components=False)
all_mesh_terms = [x[1] for x in all_mesh_cell_id_names]
for label in labels_set:
    for n_gene in n_genes:
        for method in gene_methods:
            for query_method in query_methods:
                for species in all_species:
                    for cell_type_subset in cell_type_subsets:
                        if method == 'ratio':
                            top_genes = [genes[x[0]] for x in scores_t[label][:n_gene]]
                        elif method == 't':
                            top_genes = [genes[x[0]] for x in pvals_t[label][:n_gene]]
                        top_genes = [x.upper() for x in top_genes]
                        if query_method == 'cellmarker':
                            results = cellmarker.hypergeometric_test(top_genes, species=species)
                            top_cells = [x[0] for x in results]
                        elif query_method == 'cellmarker_prob':
                            results = cellmarker_prob_method.prob_test(top_genes, species=species)
                            top_cells = [x[0] for x in results]
                        elif query_method == 'cellmesh':
                            results = cellmesh.hypergeometric_test(top_genes, species=species, cell_type_subset=cell_type_subset)
                            top_cells = [x[1] for x in results]
                        elif query_method == 'prob':
                            results = prob_method.prob_test(top_genes, species=species, cell_type_subset=cell_type_subset)
                            top_cells = [x[1] for x in results]
                        label_cell_types[(label, n_gene, method, query_method, species, cell_type_subset)] = top_cells
                        print(label, n_gene, method, query_method, species, top_cells[:10])

with open('split_seq_cellmesh_query_top_cells_cell_subset.pkl', 'wb') as f:
    pickle.dump(label_cell_types, f)

cell_type_df = []
for k, v in label_cell_types.items():
    if k[4] != 'mouse' or k[1] != 50:
        continue
    for i, cell in enumerate(v[:5]):
        cell_type_df.append((k[0], k[3], k[5], cell, i))

cell_type_df = pd.DataFrame(cell_type_df, columns=['ground_truth', 'query_method', 'cell_subset', 'cell_type', 'rank'])

cell_type_df.to_csv('split_seq_50genes_mouse_top5_cells_cell_subset.csv', index=None)
