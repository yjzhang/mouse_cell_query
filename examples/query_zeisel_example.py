import numpy as np
import scipy.io
import scipy.sparse
import mouse_cell_query

path = '/home/yjzhang/Grad_School/single_cell/uncurl_test_datasets/zeisel/'

data_mat = scipy.io.loadmat(path + 'Zeisel.mat')
data = data_mat['data']
genes = np.array([x.strip() for x in data_mat['genes']])
labels = np.array([x.strip() for x in data_mat['label_names']])
true_labels = [
        'astrocytes_ependymal',
        'endothelial-mural',
        'interneurons',
        'microglia',
        'oligodendrocytes',
        'pyramidal CA1',
        'pyramidal SS',
]

# potentially correct labels - TODO
true_labels_cell_ontology = {
        'astrocytes_ependymal': ['Astrocyte', 'Non-Neuronal - Astro'],
        'endothelial-mural': ['Endothelial', 'Non-Neuronal - Endo'],
        'interneurons': ['Neuron', 'GABAergic', 'Glutamatergic'],
        'microglia': ['Microglia', 'Non-Neuronal - Microglia'],
        'oligodendrocytes': ['Oligodendrocyte', 'Non-Neuronal - Oligo'],
        'pyramidal CA1': ['Pyramidal neuron cell', 'VISp'],
        'pyramidal SS': ['Pyramidal neuron cell', 'VISp'],
}

# TODO: run a systematic experiment
import matplotlib.pyplot as plt
plt.figure(figsize=(12, 7))
#methods = ['spearman', 'spearman_nonzero', 'cosine', 'poisson', 'random']
methods = ['spearman', 'kendall', 'cosine', 'poisson', 'random']
dbs = ['allen_cluster_means', 'mca_coarse_means']
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
            print(lab, results[:10])

        timings[(method, db)] = timing_labels
        print('mean time per query: ', np.mean(timing_labels))
        # get accuracy at top 1, top 5, top 10
        result_ranks = []
        for labi, results in cell_label_results.items():
            true_results = true_labels_cell_ontology[labi]
            rank = len(results)
            for i, xr in enumerate(results):
                has_rank = False
                if xr[0] in true_results:
                    has_rank = True
                    rank = i
                for tr in true_results:
                    if tr.lower() in xr[0].lower():
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
plt.title('Accuracy vs rank on zeisel 2015')
plt.xlabel('rank')
plt.ylabel('accuracy')
plt.xticks(range(0, ranks, int(ranks/10)))
plt.legend()
plt.savefig('query_accuracy_zeisel.png', dpi=200)
