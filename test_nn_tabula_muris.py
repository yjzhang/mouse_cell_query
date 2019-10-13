# TODO: test FACS model on droplet-seq and test droplet-seq model on FACS

from mouse_cell_query import nn_query
model = nn_query.Classifier.load_from_file('mouse_cell_query/models/tm_combined_model_200_200.h5')

import numpy as np
import scipy.io
from scipy import sparse
import os

genes = np.loadtxt('mouse_cell_query/data/gene_names.txt', dtype=str)
# all_matrices is cells x genes
all_matrices = sparse.csr_matrix(scipy.io.mmread('tm_droplet_all_matrices.mtx.gz').T)
print('loaded matrix')

all_facs_labels = np.loadtxt('tm_facs_all_labels.txt', dtype=str, delimiter='##')
all_labels = np.loadtxt('tm_droplet_all_labels.txt', dtype=str, delimiter='##')


# map droplet_seq labels onto facs labels
# take subset of cells with overlapping labels
# run classification
results = model.predict(all_matrices, genes)
label_results = model.results_to_labels(results)
accuracy = float(sum(label_results==all_labels))/len(all_labels)
print('accuracy of combined model on droplet-seq:', accuracy)
