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
plt.figure(figsize=(12, 7))
#methods = ['spearman', 'spearman_nonzero', 'cosine', 'poisson', 'random']
methods = ['spearman', 'spearman_nonzero', 'cosine', 'poisson', 'random']
dbs = ['tm_means', 'mca_means']
ticks = ['x', 'o', '*', '+', 's']
ranks = 21
import time
timings = {}
for method, tick in zip(methods, ticks):
    print(method)
    for db in dbs:
        print(db)
        cell_label_results = {}
        timing_labels = []
        for lab in sorted(list(set(labels))):
            # time different methods
            t0 = time.time()
            results = mouse_cell_query.search_db(np.array(data[:, labels==lab].mean(1)).flatten(), genes, method=method, db=db)
            timing_labels.append(time.time() - t0)
            cell_label_results[lab] = results

        timings[(method, db)] = timing_labels
        print('mean time per query: ', np.mean(timing_labels))
        # get accuracy at top 1, top 5, top 10
        result_ranks = []
        for labi, results in cell_label_results.items():
            true_results = true_labels_cell_ontology[labi]
            rank = 0
            for i, xr in enumerate(results):
                has_rank = False
                if xr[0] in true_results:
                    has_rank = True
                    rank = i
                for tr in true_results:
                    if tr in xr[0]:
                        has_rank = True
                        rank = i
                if has_rank:
                    break
            result_ranks.append(rank)
        result_ranks = np.array(result_ranks)
        accuracies = []
        ranks_list = []
        prev_accuracy = -1
        prev_rank = 0
        # plot accuracies only when change happens?
        for i in range(ranks):
            ac = float(sum(result_ranks <= i))/len(result_ranks)
            if ac != prev_accuracy:
                #if prev_accuracy != -1:
                #    accuracies.append(prev_accuracy)
                #    ranks_list.append(i-1)
                accuracies.append(ac)
                ranks_list.append(i)
                prev_accuracy = ac
                prev_rank = i
        plt.plot(ranks_list, accuracies, '--' + tick, label=method+'_'+db)

plt.grid()
plt.title('Accuracy vs rank on 10x_400')
plt.xlabel('rank')
plt.ylabel('accuracy')
plt.xticks(range(0, ranks, int(ranks/10)))
plt.legend()
plt.savefig('query_accuracy_10x_400.png', dpi=200)
