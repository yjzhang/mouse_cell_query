import numpy as np
import scipy.io
import scipy.sparse
import mouse_cell_query

path = '/home/yjzhang/Grad_School/single_cell/uncurl_test_datasets/10x_pure_pooled/'

data = scipy.io.mmread(path + 'data_400_cells.mtx')
data = scipy.sparse.csc_matrix(data)
genes = np.loadtxt(path + 'gene_names_400.tsv', dtype=str)
labels = np.loadtxt(path + 'labs_400_cells.txt').astype(int)

mouse_cell_query.search_db(np.array(data[:, labels==0].mean(1)).flatten(), genes, method='spearman_nonzero')
mouse_cell_query.search_db(np.array(data[:, labels==0].mean(1)).flatten(), genes, method='poisson')
mouse_cell_query.search_db(np.array(data[:, labels==0].mean(1)).flatten(), genes, method='cosine')
