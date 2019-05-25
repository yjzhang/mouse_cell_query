# TODO: process data for cellatlassearch

# the data output should be a csv where the columns are the gene names (capitalized), and the header is "Cell 1", "Cell 2",...
# try two types of outputs: random cells from each cluster, and try cluster means

# 1. load data

import numpy as np
import scipy.io
from scipy import sparse
import os

genes = np.loadtxt('mouse_cell_query/data/gene_names.txt', dtype=str)
all_matrices = scipy.io.mmread('tm_droplet_all_matrices.mtx')
all_matrices = sparse.csc_matrix(all_matrices)
all_labels = np.loadtxt('tm_droplet_all_labels.txt', dtype=str, delimiter='##')

# 1. take label means
label_means = {}
for label in set(all_labels):
    label_means[label] = np.array(all_matrices[:, all_labels==label].mean(1)).flatten()

import pandas as pd
means_data = pd.DataFrame(label_means, index=genes)
# this data can be used for cellatlassearch?
means_data.to_csv('tm_means_data.csv')

# do diffexp, save top n genes

from uncurl_analysis import gene_extraction
scores_t, pvals_t = gene_extraction.one_vs_rest_t(all_matrices, all_labels, eps=0.1, test='t')

# extract top n genes
top_genes = {key: [genes[x[0]] for x in vals[:100]] for key, vals in scores_t.items()}
top_genes_table = pd.DataFrame(top_genes, index=np.arange(0, 100))
# try this for scQuery
top_genes_table.to_csv('tm_top_100_genes.csv')
