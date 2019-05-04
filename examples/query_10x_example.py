import numpy as np
import scipy.io
import scipy.sparse
import mouse_cell_query

path = '/home/yjzhang/Grad_School/single_cell/uncurl_test_datasets/10x_pure_pooled/'

data = scipy.io.mmread(path + 'data_400_cells.mtx')
data = scipy.sparse.csc_matrix(data)
genes = np.loadtxt(path + 'gene_names_400.tsv', dtype=str)
labels = np.loadtxt(path + 'labs_400_cells.txt').astype(int)
true_labels = {
        0:   'CD19+ b cells',
        1:   'CD14+ monocytes',
        2:   'CD34+',
        3:   'CD4+ t helper',
        4:   'CD56+ nk',
        5:   'CD8+ cytotoxic t',
        6:   'CD4+/CD45RO+ memory t',
        7:   'CD8+/CD45RA+ naive cytotoxic',
        8:   'CD4+/CD45RA+/CD25- naive t',
        9:   'CD4+/CD25 regulatory t'
}

# potentially correct labels
true_labels_cell_ontology = {
        0:   ['B cell'],
        1:   ['monocyte', 'classical monocyte', 'non-classical monocyte'],
        2:   ['hematopoietic precursor cell'],
        3:   ['T cell'],
        4:   ['natural killer cell'],
        5:   ['T cell'],
        6:   ['T cell'],
        7:   ['T cell'],
        8:   ['T cell', 'immature T cell'],
        9:   ['T cell']
}

# TODO: run a systematic experiment
import matplotlib.pyplot as plt
plt.figure(figsize=(12, 8))
methods = ['spearman', 'spearman_nonzero', 'cosine', 'poisson', 'random']
ticks = ['x', 'o', '*', '+', 's']
for method, tick in zip(methods, ticks):
    cell_label_results = {}
    for lab in sorted(list(set(labels))):
        results = mouse_cell_query.search_db(np.array(data[:, labels==lab].mean(1)).flatten(), genes, method=method)
        cell_label_results[lab] = results

    # get accuracy at top 1, top 5, top 10
    result_ranks = []
    for labi, results in cell_label_results.items():
        true_results = true_labels_cell_ontology[labi]
        rank = 0
        for i, xr in enumerate(results):
            if xr[0] in true_results:
                rank = i
                break
        result_ranks.append(rank)
    result_ranks = np.array(result_ranks)
    accuracies = np.zeros(20)
    for i in range(20):
        accuracies[i] = float(sum(result_ranks < i))/len(result_ranks)
    plt.plot(range(20), accuracies, '--' + tick, label=method)

plt.grid()
plt.title('Accuracy vs rank on 10x_400')
plt.xlabel('rank')
plt.ylabel('accuracy')
plt.xticks(range(0, 21, 2))
plt.legend()
plt.savefig('query_accuracy_10x_400.png')


# try a different dataset... Zeisel dataset?
path = '/home/yjzhang/Grad_School/single_cell/uncurl_test_datasets/zeisel/'
data = scipy.io.loadmat(path + 'Zeisel.mat')
