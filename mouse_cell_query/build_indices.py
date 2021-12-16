# TODO:
# 1. sub-set genes: find most over-expressed genes in each cluster/cell type, and genes with the highest variances overall (that are also highly expressed).
# 2. Sub-select data with these genes.
# 3. Build indices with nmslib - try to build spearman correlation indices. Actually that's not possible, so let's just build cosine similarity indices.
import os

import numpy as np
import scipy.io
from scipy import sparse

DATA_DIR = '../'
TM_DROPLET_PATH = os.path.join(DATA_DIR, 'tabula_muris/cell_type_matrices_tm_droplet')
TM_FACS_PATH = os.path.join(DATA_DIR, 'tabula_muris/cell_type_matrices_tm_facs')
TM_GENES_PATH = os.path.join(DATA_DIR, 'mouse_cell_query/data/gene_names.txt')
MCA_PATH = os.path.join(DATA_DIR, 'cell_type_matrices_mca_coarse')
MCA_GENES_PATH = os.path.join(DATA_DIR, 'mca_microwell_seq/genes_mca.txt')

# TODO: map cell type names together using that one table


def load_data(paths, genes_path, genes_subset_path=None, normalize=True, log=True):
    """
    Returns:
        data - cells x genes csr matrix
        genes - array of gene names
        all_labels - all cell labels
    """
    all_matrices = []
    all_labels = []
    for data_path in paths:
        for f in os.listdir(data_path):
            if f.endswith('.mtx.gz'):
                file_path = os.path.join(data_path, f)
                label = f.split('.')[0]
                print(label)
                data = scipy.io.mmread(file_path)
                data = sparse.csc_matrix(data)
                all_matrices.append(data)
                all_labels.extend([label]*data.shape[1])
                print('num cells: ', data.shape[1])
    all_matrices = sparse.hstack(all_matrices)
    all_labels = np.array(all_labels)

    n_genes, n_cells = all_matrices.shape

    print('data matrix shape:', all_matrices.shape)
    print('number of cells:', all_labels.shape)

    from uncurl import preprocessing
    if normalize:
        all_matrices = preprocessing.cell_normalize(all_matrices, multiply_means=False)
    if log:
        all_matrices = preprocessing.log1p(all_matrices)
    all_matrices = sparse.csr_matrix(all_matrices.T)
    print('matrices normalized')

    genes = np.loadtxt(genes_path, dtype=str)
    if genes_subset_path is not None:
        all_matrices, genes = load_genes_subset(all_matrices, genes, genes_subset_path)
    return all_matrices, genes, all_labels


def load_genes_subset(data, genes, genes_subset_path):
    """
    data is a cells x genes matrix
    """
    genes_subset = np.loadtxt(genes_subset_path, dtype=str)
    genes_dict = {g: i for i, g in enumerate(genes)}
    genes_subset_ids = np.array([genes_dict[g] for g in genes_subset])
    print('number of genes in subset:', len(genes_subset))
    data_subset = data[:, genes_subset_ids]
    print('data subset shape:', data_subset.shape)
    return data_subset, genes_subset


# 1. gene selection?
def get_genes_tm(all_matrices, genes, all_labels):
    """
    This function returns the top genes from Tabula Muris.
    """
    labels_set = set(all_labels)
    n_genes = len(genes)
    if len(genes) == all_matrices.shape[1]:
        all_matrices = all_matrices.T
    print(all_matrices.shape)
    print(all_labels.shape)

    from uncurl_analysis import gene_extraction
    scores_t, pvals_t = gene_extraction.one_vs_rest_t(all_matrices, all_labels, eps=0.1, test='t')
    # take the top 50 genes for each cell type
    cell_type_genes = {}
    all_top_genes = set()
    n_gene = 50
    for label in labels_set:
        top_genes = [genes[x[0]] for x in scores_t[label][:n_gene]]
        cell_type_genes[label] = top_genes
        all_top_genes.update(top_genes)
    # TODO: get genes with highest variance
    from uncurl import preprocessing
    means, vars = preprocessing.sparse_mean_var(all_matrices)
    nonzeros = (all_matrices > 0)
    genes_sum = np.array(nonzeros.sum(1)).flatten()
    # get the 2000 genes with the highest variance/mean that also have at least 10% nonzero
    mv = vars/(means + 1e-8)
    mv_indices = mv.argsort()[::-1]
    top_genes = mv_indices[:1000]
    top_genes = top_genes[genes_sum[top_genes] > (n_genes/10.)]
    top_genes = [genes[x] for x in top_genes]
    all_top_genes.update(top_genes)
    return all_top_genes


# train an NN
def train_nn_model(data_subset, genes_subset, all_labels, output_path='tm_combined_classifier.h5'):
    """
    Trains a NN model for classification.

    Args:
        data_subset - csr matrix of shape (cells, genes)
        genes_subset - np array of length (genes)
        all_labels - np array of strings of length (cells)
    """
    labels_set = set(all_labels)
    labels_map = {label: i for i, label in enumerate(sorted(list(labels_set)))}
    all_labels_ints = np.array([labels_map[lab] for lab in all_labels])
    all_labels_one_hot = np.zeros((len(all_labels), len(labels_set)))
    for i, lab in enumerate(all_labels_ints):
        all_labels_one_hot[i, lab] = 1.0

    from sklearn.model_selection import train_test_split
    train_indices, test_indices = train_test_split(np.arange(0, len(all_labels)))
    train_matrix = data_subset[train_indices, :]
    train_labels = all_labels_one_hot[train_indices, :]
  
    test_matrix = data_subset[test_indices, :]
    #test_labels = all_labels_one_hot[test_indices, :]
    test_labels_ints = all_labels_ints[test_indices]
 
    import nn_query
    layers = [500, 250]
    model = nn_query.Classifier(genes_subset, len(labels_set), layers=layers,
            class_names=np.array(sorted(list(labels_set))))
    model.fit(train_matrix, train_labels, n_epochs=30, batch_size=1000)
    results = model.predict(test_matrix, genes_subset)
    result_labels = results.argmax(1)
    accuracy = float(sum(results.argmax(1)==test_labels_ints))/len(test_indices)
    print('test accuracy of combined model on combined data:', accuracy)
    # calculate mean balanced accuracy
    from sklearn.metrics import balanced_accuracy_score
    print('Mean balanced accuracy:', balanced_accuracy_score(test_labels_ints, result_labels, adjusted=True))
    model.save(output_path)
    return model


def build_nmslib_index(data, genes, all_labels, output_path='tm_combined_classifier.h5'):
    pass


if __name__ == '__main__':
    #tm_droplet_top_genes = get_genes_tm(TM_DROPLET_PATH, TM_GENES_PATH)
    #np.savetxt('tm_droplet_top_genes.txt', np.array(list(tm_droplet_top_genes)), fmt='%s')
    #tm_facs_top_genes = get_genes_tm(TM_FACS_PATH, TM_GENES_PATH)
    #np.savetxt('tm_facs_top_genes.txt', np.array(list(tm_facs_top_genes)), fmt='%s')
    #tm_combined_top_genes = tm_droplet_top_genes.union(tm_facs_top_genes)
    #np.savetxt('tm_combined_top_genes.txt', np.array(list(tm_combined_top_genes)), fmt='%s')
    paths = [TM_DROPLET_PATH, TM_FACS_PATH]
    genes_path = TM_GENES_PATH
    genes_subset_path = 'tm_combined_top_genes.txt'

    data_tm, genes, labels_tm = load_data([TM_DROPLET_PATH], genes_path, normalize=False)
    tm_droplet_top_genes = get_genes_tm(data_tm, genes, labels_tm)
    data_facs, genes, labels_facs = load_data([TM_DROPLET_PATH], genes_path, normalize=False)
    tm_facs_top_genes = get_genes_tm(data_facs, genes, labels_facs)
    tm_combined_top_genes = tm_droplet_top_genes.union(tm_facs_top_genes)
    np.savetxt('tm_combined_top_genes.txt', np.array(list(tm_combined_top_genes)), fmt='%s')

    data_combined = sparse.vstack(data_tm, data_facs)
    labels_combined = np.hstack(labels_tm, labels_facs)
    data_subset, genes_subset = load_genes_subset(data_combined, labels_combined, genes_subset_path)
 
    output_path = 'tm_combined_classifier.h5'
    #model = train_nn_model(paths, genes_path, genes_subset_path)
