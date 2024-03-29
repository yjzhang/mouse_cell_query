# TODO:
# 1. load all mca data, with specific clusters
# 2. do diffexp using standard uncurl-app pipelines, take top k genes per cell type, 

# 1. load data
import random
import numpy as np
import scipy.io
from scipy import sparse
import os
import pickle

path = 'cell_type_matrices_mca'
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
scipy.io.mmwrite('mca_fine_all_matrices.mtx', all_matrices)
np.savetxt('mca_fine_all_labels.txt', all_labels, fmt='%s')

############################################################################################

genes = np.loadtxt('genes_mca.txt', dtype=str)
all_matrices = scipy.io.mmread('mca_fine_all_matrices.mtx')
all_labels = np.loadtxt('mca_fine_all_labels.txt', dtype=str, delimiter='##')

# 2. calculate diffexp for each cluster
from uncurl_analysis import gene_extraction
import time
t0 =  time.time()
scores_t, pvals_t = gene_extraction.one_vs_rest_t(all_matrices, all_labels, eps=0.1, test='t')
print('diffexp time for t test:', time.time() - t0)
t0 =  time.time()
#scores_u, pvals_u = gene_extraction.one_vs_rest_t(all_matrices, all_labels, eps=0.1, test='u')
print('diffexp time for u test:', time.time() - t0)

with open('mca_fine_t_scores.pkl', 'wb') as f:
    pickle.dump(scores_t, f)
with open('mca_fine_t_pvals.pkl', 'wb') as f:
    pickle.dump(pvals_t, f)
#with open('mca_fine_u_pvals.pkl', 'wb') as f:
#    pickle.dump(pvals_u, f)

###########################################################################################
# 3. for each cluster, run cellmesh and cellmarker

with open('mca_fine_t_scores.pkl', 'rb') as f:
    scores_t = pickle.load(f)
with open('mca_fine_t_pvals.pkl', 'rb') as f:
    pvals_t = pickle.load(f)
#with open('mca_fine_u_pvals.pkl', 'rb') as f:
#   pvals_u = pickle.load(f)
all_labels = np.loadtxt('mca_fine_all_labels.txt', dtype=str, delimiter='##')
genes = np.loadtxt('genes_mca.txt', dtype=str)

# get gene names - save ranked genes for each cell type
top_genes_ratio = {}
for cell, cell_genes in scores_t.items():
    top_genes_ratio[cell] = [genes[i[0]] for i in cell_genes]
import pandas as pd
top_genes_ratio = pd.DataFrame(top_genes_ratio)
top_genes_ratio.to_csv('mca_fine_top_genes_ratio.csv')

import cellmesh
from cellmesh import prob_method, gsva_ext_method
import cellmarker
from cellmarker import prob_method as cellmarker_prob_method
from mouse_cell_query import query_aggregation
labels_set = set(all_labels)
label_results = {}
label_cell_types = {}
n_genes = [50]
gene_methods = ['ratio']#, 't', 'u']
query_methods = ['cellmarker', 'cellmarker_prob','panglao', 'cellmesh', 'cellmesh_tfidf', 'prob', 'gsva', 'random_mesh', 'tissue_mesh_prob']#, 'aggregate', 'aggregate_2']
all_species = ['human', 'mouse', 'both']
all_mesh_cell_id_names = cellmesh.get_all_cell_id_names(include_cell_components=False)
all_mesh_terms = [x[1] for x in all_mesh_cell_id_names]
for label in labels_set:
    for n_gene in n_genes:
        for method in gene_methods:
            for query_method in query_methods:
                for species in all_species:
                    if method == 'ratio':
                        top_genes = [genes[x[0]] for x in scores_t[label][:n_gene]]
                    elif method == 't':
                        top_genes = [genes[x[0]] for x in pvals_t[label][:n_gene]]
                    elif method == 'u':
                        top_genes = [genes[x[0]] for x in pvals_u[label][:n_gene]]
                    top_genes = [x.upper() for x in top_genes]
                    if query_method == 'cellmarker':
                        results = cellmarker.hypergeometric_test(top_genes, species=species)
                        top_cells = [x[0] for x in results]
                    elif query_method == 'cellmarker_prob':
                        results = cellmarker_prob_method.prob_test(top_genes, species=species)
                        top_cells = [x[0] for x in results]
                    elif query_method == 'panglao':
                        results = cellmarker.hypergeometric_test(top_genes, species=species, db_dir=cellmarker.PANGLAO_DB_DIR)
                        top_cells = [x[0] for x in results]
                    elif query_method == 'cellmesh':
                        results = cellmesh.hypergeometric_test(top_genes, species=species)
                        top_cells = [x[1] for x in results]
                    elif query_method == 'cellmesh_tfidf':
                        results = cellmesh.normed_hypergeometric_test(top_genes, species=species)
                        top_cells = [x[1] for x in results]
                    elif query_method == 'aggregate':
                        results = query_aggregation.cellmarker_cellmesh_hypergeometric_test(top_genes)
                        top_cells = [x[1] for x in results[1:]]
                        results = results[1:]
                    elif query_method == 'aggregate_2':
                        results = query_aggregation.cellmarker_cellmesh_tfidf_hypergeometric_test(top_genes)
                        top_cells = [x[1] for x in results[1:]]
                        results = results[1:]
                    elif query_method == 'prob':
                        results = prob_method.prob_test(top_genes, species=species)
                        top_cells = [x[1] for x in results]
                    elif query_method == 'gsva':
                        results = gsva_ext_method.gsva_ext_test(top_genes, species=species)
                        top_cells = [x[1] for x in results]
                    elif query_method == 'tissue_mesh_prob':
                        results = prob_method.prob_test(top_genes, db_dir=cellmesh.ANATOMY_DB_DIR, species=species)
                        top_cells = [x[1] for x in results]
                    elif query_method == 'random_mesh':
                        # select a random MeSH term
                        top_cells = random.sample(all_mesh_terms, 10)
                        results = top_cells
                    #label_results[(label, n_gene, method, query_method, species)] = [x[:-1] for x in results]
                    label_cell_types[(label, n_gene, method, query_method, species)] = top_cells
                    print(label, n_gene, method, query_method, species, top_cells[:10])


with open('mca_fine_cellmesh_query_top_cells.pkl', 'wb') as f:
    pickle.dump(label_cell_types, f)


############################################################################################

import pandas as pd
with open('mca_fine_cellmesh_query_top_cells.pkl', 'rb') as f:
    label_cell_types = pickle.load(f)

# TODO: how to compare accuracy?
# TODO: convert MCA names to tabula muris
# can we use a precision-recall curve?
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

with open('tm_cell_onto_alternate_names.tsv') as f:
    l0 = f.readline()
    for line in f.readlines():
        line_data = [x.strip() for x in line.split('\t')]
        name = line_data[1]
        alternate_cell_type_names = line_data[2:]
        if name in cell_types_alternate_map:
            cell_types_alternate_map[name].extend(alternate_cell_type_names)

mca_cell_names_to_cellmarker = pd.read_table('mca_cell_names_to_cellmarker.tsv')
mca_cell_names_map = {}
import cellmesh
for i, row in mca_cell_names_to_cellmarker.iterrows():
    # TODO: match with cell_types_map and cell_types_alternate_map
    cell_name = row['mca_cell_name']
    if isinstance(row['tabula_muris'], str):
        alt_name = row['tabula_muris'].strip()
        mca_cell_names_map[cell_name] = cell_types_alternate_map[alt_name]
        mca_cell_names_map[cell_name].append(cell_types_map[alt_name])
    if isinstance(row['cellmarker_cellonto'], str):
        alt_name = row['cellmarker_cellonto'].strip()
        try:
            mca_cell_names_map[cell_name].append(alt_name)
        except:
            mca_cell_names_map[cell_name] = [alt_name]
        mca_cell_names_map[cell_name].extend([x[1] for x in cellmesh.cellonto_to_cellmesh(alt_name)])
    if isinstance(row['cellmesh'], str):
        alt_name = row['cellmesh'].strip()
        try:
            mca_cell_names_map[cell_name].append(alt_name)
        except:
            mca_cell_names_map[cell_name] = [alt_name]
        mca_cell_names_map[cell_name].extend([x[1] for x in cellmesh.cellonto_to_cellmesh(alt_name)])

# TODO: calculate accuracy of each label list
label_accuracies = {}
label_extended_accuracies = {}
for key, value in label_cell_types.items():
    name = key[0]
    name = name.split('(')[0].split('_')[0].strip()
    name = name.capitalize().strip()
    if name.endswith('s'):
        name = name[:-1]
    accuracies = []
    extended_accuracies = []
    if name not in mca_cell_names_map:
        continue
    for v in value:
        if v.lower() in name.lower() or name.lower() in v.lower():
            accuracies.append(1)
            extended_accuracies.append(1)
        else:
            accuracies.append(0)
            if v in mca_cell_names_map[name]:
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

with open('mca_fine_MAP_method_means.pkl', 'wb') as f:
    pickle.dump(map_method_means, f)

with open('mca_fine_MAP_method_means.pkl', 'rb') as f:
    map_method_means = pickle.load(f)

# convert map_method_means to a pandas dict
import pandas as pd

map_list = [key + (v,) for key, v in map_method_means.items()]
map_list.sort()
map_method_means = pd.DataFrame(map_list, columns=['n_genes',  'gene_method', 'query_method', 'species', 'mean_average_precision'])

map_cell_types_list = [key + (v,) for key, v in label_map.items()]
map_cell_types_list.sort()
map_cell_types = pd.DataFrame(map_cell_types_list, columns=['cell_type', 'n_genes',  'gene_method', 'query_method', 'species', 'mean_average_precision'])

map_method_means.to_csv('mca_fine_MAP_method_means.csv', index=None)
map_cell_types.to_csv('mca_fine_MAP_cell_types.csv', index=None)

#####################################################################################################################

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

map_method_means = pd.read_csv('mca_fine_MAP_method_means.csv')
map_cell_types = pd.read_csv('mca_fine_MAP_cell_types.csv')
# load scQuery results
scQuery_results = pd.read_csv('MAP_mca_scquery_cell_types.csv')
map_cell_types = map_cell_types.append(scQuery_results)
map_method_means = map_method_means.append({'n_genes': 50,'gene_method': 'ratio', 'query_method': 'scquery', 'species': 'mouse', 'mean_average_precision': scQuery_results['mean_average_precision'].mean()}, ignore_index=True)


# TODO: plot?
# plot cellmarker, cellmesh, cellmesh_tfidf. fix gene_method='ratio', group by query_method, plot over all n_genes.

# for each cell type, plot performance of all methods
all_cell_types = sorted(list(set(map_cell_types.cell_type)))
map_cell_types_subset = map_cell_types[map_cell_types.gene_method=='ratio']
sns.set(style='whitegrid', font_scale=1.0)
plot_y = int(np.ceil(len(all_cell_types)/10))
fig, axes = plt.subplots(10, plot_y, figsize=(120, 50))
for i, axes_1 in enumerate(axes):
    for j, ax in enumerate(axes_1):
        index = i*plot_y + j
        if index >= len(all_cell_types):
            break
        data_subset = map_cell_types_subset[map_cell_types_subset.cell_type==all_cell_types[index]]
        g = sns.categorical.barplot(x='query_method', y='mean_average_precision', hue='n_genes', data=data_subset, ax=ax)
        ax.set_title(all_cell_types[index])
        if j != 0:
            ax.set_ylabel('')
        if i != 6:
            ax.set_xlabel('')
        ax.set_ylim(0, 1)
        if i != j:
            ax.get_legend().remove()
plt.savefig('map_ratios_mca_fine_all_cell_types.png', dpi=100)



map_method_means_subset = map_method_means[map_method_means.gene_method=='ratio']

# plot performance of all methods
sns.set(style='whitegrid', font_scale=1.5)
fig, ax = plt.subplots(figsize=(18, 10))
g = sns.categorical.barplot(x='query_method', y='mean_average_precision', hue='n_genes', data=map_method_means_subset, ax=ax)
plt.ylim(0, 0.8)
plt.title('Cell Type Annotation Accuracy')
plt.savefig('map_ratios_mca_fine.png', dpi=100)

# analysis: which cell types did each of the methods perform poorly on?
best_methods_per_cell_type = map_cell_types.sort_values('mean_average_precision', ascending=False).groupby('cell_type').first()

cell_type_no_cellmarker = map_cell_types[map_cell_types.query_method != 'cellmarker']
best_methods_no_cellmarker = cell_type_no_cellmarker.sort_values('mean_average_precision', ascending=False).groupby('cell_type').first()

# add scQuery results
# TODO: calculate results with scquery
#new_map_method_means = map_method_means.append({'n_genes': 50,'gene_method': 'ratio', 'query_method': 'scQuery', 'mean_average_precision': 0.32136}, ignore_index=True)
new_mmm_subset = map_method_means[(map_method_means.gene_method=='ratio') & (map_method_means.n_genes==50)]

sns.set(style='whitegrid', font_scale=1.5)
fig, ax = plt.subplots(figsize=(18, 10))
g = sns.categorical.barplot(x='query_method', y='mean_average_precision', hue='n_genes', data=new_mmm_subset, ax=ax)
plt.ylim(0, 0.8)
plt.title('Cell Type Annotation Accuracy')
plt.savefig('map_ratios_mca_fine_with_scquery.png', dpi=100)

#sns.set(style='whitegrid', font_scale=1.5)
#fig, ax = plt.subplots(figsize=(8, 10))
#g = sns.categorical.barplot(x='query_method', y='mean_average_precision', hue='n_genes', data=new_mmm_subset, ax=ax)
#plt.ylim(0, 0.8)
#plt.title('Cell Type Annotation Accuracy')
#plt.savefig('map_ratios_mca_fine_with_scquery.png', dpi=100)

# convert MAP to top-1, top-3, and top-5 accuracy
map_cell_types['top_1'] = [1 if x > 0.9 else 0 for x in map_cell_types['mean_average_precision']]
map_cell_types['top_3'] = [1 if x > 0.3 else 0 for x in map_cell_types['mean_average_precision']]
map_cell_types['top_5'] = [1 if x >= 0.2 else 0 for x in map_cell_types['mean_average_precision']]
map_cell_types['top_10'] = [1 if x >= 0.1 else 0 for x in map_cell_types['mean_average_precision']]
# calculate mean top-1, top-3, and top-5 accuracy by method & gene count
map_cell_type_means = map_cell_types.groupby(['query_method', 'gene_method', 'n_genes', 'species']).mean().reset_index()
# plot top-1, top-3, and top-5 accuracy
sns.set(style='whitegrid', font_scale=1.5)
fig, axes = plt.subplots(1, 4, figsize=(55, 10))
categories = ['top_1', 'top_3', 'top_5', 'top_10'] 
titles = ['Top-1', 'Top-3', 'Top-5', 'Top-10'] 
for i, ax in enumerate(axes):
    g = sns.categorical.barplot(x='query_method', y=categories[i], hue='n_genes', data=map_cell_type_means[map_cell_type_means.gene_method=='ratio'], ax=ax)
    g.set_ylim(0, 1.0)
    g.set_title('Cell Type Annotation {0} Accuracy'.format(titles[i])) 
plt.savefig('top_1_accuracy_mca_fine.png', dpi=100)

fig, ax = plt.subplots(figsize=(18, 10))
g = sns.categorical.barplot(x='query_method', y='top_3', hue='n_genes', data=map_cell_type_means[map_cell_type_means.gene_method=='ratio'], ax=ax)
plt.ylim(0, 1.0)
plt.title('Cell Type Annotation Top-3 Accuracy')
plt.savefig('top_3_accuracy_mca_fine.png', dpi=100)

fig, ax = plt.subplots(figsize=(18, 10))
g = sns.categorical.barplot(x='query_method', y='top_5', hue='n_genes', data=map_cell_type_means[map_cell_type_means.gene_method=='ratio'], ax=ax)
plt.ylim(0, 1.0)
plt.title('Cell Type Annotation Top-5 Accuracy')
plt.savefig('top_5_accuracy_mca_fine.png', dpi=100)


# plot top 1 accuracy with 50 genes, color by species
fig, ax = plt.subplots(figsize=(18, 10))
g = sns.categorical.barplot(x='query_method', y='top_1', hue='species',
        data=map_cell_type_means[(map_cell_type_means.gene_method=='ratio') & (map_cell_type_means.n_genes==50)],
        ax=ax)
plt.ylim(0, 1.0)
plt.title('Cell Type Annotation Top-1 Accuracy')
plt.savefig('top_1_accuracy_mca_fine_50_genes.png', dpi=100)

# plot top 3 accuracy with 50 genes, color by species
fig, ax = plt.subplots(figsize=(18, 10))
g = sns.categorical.barplot(x='query_method', y='top_3', hue='species',
        data=map_cell_type_means[(map_cell_type_means.gene_method=='ratio') & (map_cell_type_means.n_genes==50)],
        ax=ax)
plt.ylim(0, 1.0)
plt.title('Cell Type Annotation Top-3 Accuracy')
plt.savefig('top_3_accuracy_mca_fine_50_genes.png', dpi=100)

# plot top 5 accuracy with 50 genes, color by species
fig, ax = plt.subplots(figsize=(18, 10))
g = sns.categorical.barplot(x='query_method', y='top_5', hue='species',
        data=map_cell_type_means[(map_cell_type_means.gene_method=='ratio') & (map_cell_type_means.n_genes==50)],
        ax=ax)
plt.ylim(0, 1.0)
plt.title('Cell Type Annotation Top-5 Accuracy')
plt.savefig('top_5_accuracy_mca_fine_50_genes.png', dpi=100)

