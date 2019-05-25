# trains a neural network on tabula muris data.

import numpy as np
import scipy.io
from scipy import sparse
import os

genes = np.loadtxt('mouse_cell_query/data/gene_names.txt', dtype=str)
# all_matrices is cells x genes
all_matrices = sparse.csr_matrix(scipy.io.mmread('tm_droplet_all_matrices.mtx').T)
all_labels = np.loadtxt('tm_droplet_all_labels.txt', dtype=str, delimiter='##')
labels_map = {label: i for i, label in enumerate(sorted(list(set(all_labels))))}
all_labels_ints = np.array([labels_map[lab] for lab in all_labels])

n_labels = len(labels_map)
all_labels_one_hot = np.zeros((len(all_labels), n_labels))
for i, lab in enumerate(all_labels_ints):
    all_labels_one_hot[i, lab] = 1.0

all_matrices_normalized = all_matrices/np.array(all_matrices.sum(1))
all_matrices_normalized = sparse.csr_matrix(all_matrices_normalized)

# train/val split
from sklearn.model_selection import train_test_split
train_indices, test_indices = train_test_split(np.arange(0, len(all_labels)))
train_matrix = all_matrices_normalized[train_indices, :]
train_labels = all_labels_one_hot[train_indices, :]
test_matrix = all_matrices_normalized[test_indices, :]
test_labels = all_labels_one_hot[test_indices, :]

# train matrix
import nn_query
model = nn_query.Classifier(genes, len(labels_map))
model.fit(train_matrix, train_labels, n_epochs=2)
# test...
results = model.predict(test_matrix, genes)
accuracy = float(sum(results.argmax(1)==all_labels_ints[test_indices]))/len(test_indices)
print('accuracy:', accuracy)
model.save('tm_droplet_model.h5')

# load model, test on FACS data
model = nn_query.Classifier.load_from_file('tm_droplet_model.h5')
