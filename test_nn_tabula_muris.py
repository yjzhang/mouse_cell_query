# TODO: test FACS model on droplet-seq and test droplet-seq model on FACS

import nn_query
model = nn_query.Classifier.load_from_file('tm_facs_model.h5')

import numpy as np
import scipy.io
from scipy import sparse
import os

genes = np.loadtxt('mouse_cell_query/data/gene_names.txt', dtype=str)
# all_matrices is cells x genes
all_matrices = sparse.csr_matrix(scipy.io.mmread('tm_droplet_all_matrices.mtx').T)
print('loaded matrix')

all_facs_labels = np.loadtxt('tm_facs_all_labels.txt', dtype=str, delimiter='##')
all_labels = np.loadtxt('tm_droplet_all_labels.txt', dtype=str, delimiter='##')

facs_labels_list = list(set(all_facs_labels))
facs_labels_list.sort()
droplet_labels_list = list(set(all_labels))
droplet_labels_list.sort()

# map droplet_seq labels onto facs labels
labels_map_facs = {label: i for i, label in enumerate(sorted(list(set(all_facs_labels))))}
# take subset of cells with overlapping labels
label_in_facs = [True if x in labels_map_facs else False for x in all_labels]
all_matrices_subset = all_matrices[label_in_facs, :]
all_labels_subset = all_labels[label_in_facs]
all_labels_ints = np.array([labels_map_facs[lab] for lab in all_labels_subset])
# run classification
results = model.predict(all_matrices_subset, genes)
accuracy = float(sum(results[label_in_facs].argmax(1)==all_labels_ints))/len(all_labels_ints)
print('accuracy of facs model on droplet-seq:', accuracy)
