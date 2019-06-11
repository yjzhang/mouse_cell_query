# TODO: take a random subsample of n cells from each cell type in tabula muris droplet
# this will be used for testing CellAtlasSearch and scMatch, as well as possibly for uncurl's cell similarity search.

n_cells = 5

import random
import numpy as np
import scipy.io
from scipy import sparse
import os
import pickle

genes = np.loadtxt('mouse_cell_query/data/gene_names.txt', dtype=str)
# all_matrices is a genes x cells array
all_matrices = scipy.io.mmread('tm_droplet_all_matrices.mtx')
all_labels = np.loadtxt('tm_droplet_all_labels.txt', dtype=str, delimiter='##')

all_matrices = sparse.csc_matrix(all_matrices)
sub_data = []
sub_labels = []
for lab in sorted(list(set(all_labels))):
    data_indices = list(np.where(all_labels == lab)[0])
    n = min(n_cells, len(data_indices))
    samples = random.sample(data_indices, n)
    sub_data.append(all_matrices[:, samples])
    sub_labels.append(all_labels[samples])

sub_data = sparse.hstack(sub_data)
sub_labels = np.hstack(sub_labels)

# re-label with number
new_labels = []
from collections import Counter
lab_counts = Counter()
for lab in sub_labels:
    lab_counts[lab] += 1
    new_labels.append(lab + '_' + str(lab_counts[lab]))

# create a table
import pandas as pd
data = pd.DataFrame(sub_data.toarray())
data.columns = [x.replace(' ', '_') for x in new_labels]
data.index = [x.upper() for x in genes]
data.to_csv('tm_droplet_data_subset.csv')
data.to_csv('tm_droplet_data_subset.tsv', sep='\t')


##########################################################

columns_to_include = [int(x.split('_')[-1]) < 6 for x in data.columns]
data = data.loc[:, columns_to_include]

###########################################################

# TODO: save means
from uncurl_analysis import dense_matrix_h5
data_means = dense_matrix_h5.H5Dict('mouse_cell_query/data/cell_type_means.h5')
all_means = []
for cell_type in sorted(data_means.keys()):
    means = data_means[cell_type]
    all_means.append(means)
all_means = np.vstack(all_means.T)
data = pd.DataFrame(all_means)
data.columns = sorted(data_means.keys())
data.index = [x.upper() for x in genes]
data.to_csv('tm_droplet_data_means.csv')
