# TODO:
# 1. load all tabula muris data, with specific clusters
# 2. do diffexp using standard uncurl-app pipelines, take top k genes per cell type, 

# 1. load data
# start with droplet data

import numpy as np
import scipy.io
from scipy import sparse
import os

"""
path = 'cell_type_matrices_tm_facs'
all_matrices = []
all_labels = []
for f in os.listdir(path):
    file_path = os.path.join(path, f)
    label = f.split('.')[0]
    print(label)
    data = scipy.io.mmread(file_path).T
    data = sparse.csc_matrix(data)
    all_matrices.append(data)
    all_labels.extend([label]*data.shape[1])
    print('num cells: ', data.shape[1])

all_matrices = sparse.hstack(all_matrices)
all_labels = np.array(all_labels)

print(all_matrices.shape)
print(all_labels.shape)

# save all_matrices, all_labels
scipy.io.mmwrite('tm_facs_all_matrices.mtx', all_matrices)
np.savetxt('tm_facs_all_labels.txt', all_labels, fmt='%s')
"""

genes = np.loadtxt('tabula_muris_facs_genes.txt', dtype=str)

all_matrices = scipy.io.mmread('tm_facs_all_matrices.mtx')
all_labels = np.loadtxt('tm_facs_all_labels.txt', dtype=str, delimiter='##')

# 2. calculate diffexp for each cluster
from uncurl_analysis import gene_extraction
import time
t0 =  time.time()
scores_t, pvals_t = gene_extraction.one_vs_rest_t(all_matrices, all_labels, eps=0.1, test='t')
print('diffexp time for t test:', time.time() - t0)
t0 =  time.time()
scores_u, pvals_u = gene_extraction.one_vs_rest_t(all_matrices, all_labels, eps=0.1, test='u')
print('diffexp time for u test:', time.time() - t0)

# 3. for each cluster, run cellmesh and cellmarker
import cellmesh
import cellmarker
labels_set = set(all_labels)
label_results = {}
label_cell_types = {}
n_genes = [20, 50, 100, 200, 1000]
gene_methods = ['ratio', 't', 'u']
query_methods = ['cellmarker', 'cellmesh', 'cellmesh_tfidf']
for label in labels_set:
    for n_gene in n_genes:
        for method in gene_methods:
            for query_method in query_methods:
                if method == 'ratio':
                    top_genes = [genes[x[0]] for x in scores_t[label][:n_gene]]
                elif method == 't':
                    top_genes = [genes[x[0]] for x in pvals_t[label][:n_gene]]
                elif method == 'u':
                    top_genes = [genes[x[0]] for x in pvals_u[label][:n_gene]]
                top_genes = [x.upper() for x in top_genes]
                if query_method == 'cellmarker':
                    results = cellmarker.hypergeometric_test(top_genes)
                    top_cells = [x[0] for x in results]
                elif query_method == 'cellmesh':
                    results = cellmesh.hypergeometric_test(top_genes)
                    top_cells = [x[1] for x in results]
                elif query_method == 'cellmesh_tfidf':
                    results = cellmesh.normed_hypergeometric_test(top_genes)
                    top_cells = [x[1] for x in results]
                label_results[(label, n_gene, method, query_method)] = results
                label_cell_types[(label, n_gene, method, query_method)] = top_cells
                print(label, n_gene, method, query_method, top_cells[:10])

import pickle
with open('tm_facs_cellmesh_query_results.pkl', 'wb') as f:
    pickle.dump(label_results, f)
with open('tm_facs_cellmesh_query_top_cells.pkl', 'wb') as f:
    pickle.dump(label_cell_types, f)


with open('tm_facs_cellmesh_query_top_cells.pkl', 'rb') as f:
    label_cell_types = pickle.load(f)

# TODO: how to compare accuracy?
# can we use a precision-recall curve?
# load the hand-created mappings file
cell_types_map = {}
cell_types_alternate_map = {}
with open('cell_ontology_to_cellmesh_tabula_muris.csv') as f:
    l0 = f.readline()
    for line in f.readlines():
        line_data = line.split(',')
        name = line_data[1]
        primary_cellmesh_name = line_data[2]
        alternate_cellmesh_names = line_data[3:]
        use_cellmesh = line_data[0]
        if use_cellmesh != 'n':
            cell_types_map[name] = primary_cellmesh_name
            cell_types_alternate_map[name] = alternate_cellmesh_names

with open('tm_cell_onto_alternate_names.csv') as f:
    l0 = f.readline()
    for line in f.readlines():
        line_data = line.split(',')
        name = line_data[1]
        alternate_cell_type_names = line_data[2:]
        if name in cell_types_alternate_map:
            cell_types_alternate_map[name].extend(alternate_cell_type_names)


# TODO: calculate accuracy of each label list
label_accuracies = {}
label_extended_accuracies = {}
for key, value in label_cell_types.items():
    name = key[0]
    accuracies = []
    extended_accuracies = []
    if name not in cell_types_map:
        continue
    for v in value:
        if v == cell_types_map[name] or v.lower() in name.lower() or name.lower() in v.lower():
            accuracies.append(1)
            extended_accuracies.append(1)
        else:
            accuracies.append(0)
            if v in cell_types_alternate_map[name]:
                extended_accuracies.append(0.5)
            else:
                extended_accuracies.append(0)
    label_accuracies[key] = accuracies
    label_extended_accuracies[key] = extended_accuracies



def calculate_precision_recall_curve(accuracies, n_relevant=1, use_extended=False):
    """
    https://nlp.stanford.edu/IR-book/html/htmledition/evaluation-of-unranked-retrieval-sets-1.html
    https://nlp.stanford.edu/IR-book/html/htmledition/evaluation-of-ranked-retrieval-results-1.html
    """
    precision = []
    recall = []
    mean_avg_precision = 0
    num_correct = 0
    for i, a in enumerate(accuracies):
        if a == 1:
            num_correct += 1
        if use_extended and a > 0:
            num_correct += 1
        precision.append(float(min(num_correct, n_relevant))/(i+1))
        recall.append(float(min(num_correct, n_relevant))/n_relevant)
    prev_recall = 0
    for p, r in zip(precision, recall):
        mean_avg_precision += (r - prev_recall)*p
        prev_recall = r
    return precision, recall, mean_avg_precision

# calculate precision-recall curves for all keys
label_map = {}
for key, val in label_extended_accuracies.items():
    label_map[key] = calculate_precision_recall_curve(val, use_extended=True)[2]

# take mean over cell types
# average MAP for all methods, over all datasets
map_methods = {}
for key, val in label_map.items():
    method_key = key[1:]
    try:
        map_methods[method_key].append(val)
    except:
        map_methods[method_key] = [val]

map_method_means = {key: np.mean(val) for key, val in map_methods.items()}

with open('MAP_method_means.pkl', 'wb') as f:
    pickle.dump(map_method_means, f)

with open('MAP_method_means.pkl', 'rb') as f:
    map_method_means = pickle.load(f)

# convert map_method_means to a pandas dict
import pandas as pd

map_list = [key + (v,) for key, v in map_method_means.items()]
map_list.sort()
map_method_means = pd.DataFrame(map_list, columns=['n_genes',  'gene_method', 'query_method', 'mean_average_precision'])

map_cell_types_list = [key + (v,) for key, v in label_map.items()]
map_cell_types_list.sort()
map_cell_types = pd.DataFrame(map_cell_types_list, columns=['cell_type', 'n_genes',  'gene_method', 'query_method', 'mean_average_precision'])

map_method_means.to_csv('MAP_facs_method_means.csv', index=None)
map_cell_types.to_csv('MAP_facs_cell_types.csv', index=None)

##################################################################
