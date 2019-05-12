import numpy as np
import scipy.io
import scipy.sparse
import mouse_cell_query

path = '/home/yjzhang/Grad_School/single_cell/uncurl_test_datasets/zeisel/'

data_mat = scipy.io.loadmat(path + 'Zeisel.mat')
data = data_mat['data']
genes = data_mat['genes']
labels = data_mat['label_names']
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
        'astrocytes_ependymal': [],
        'endothelial-mural': [],
        'interneurons': [],
        'microglia': [],
        'oligodendrocytes': [],
        'pyramidal CA1': [],
        'pyramidal SS': [],
}

# TODO: run a systematic experiment
import matplotlib.pyplot as plt
plt.figure(figsize=(12, 7))
#methods = ['spearman', 'spearman_nonzero', 'cosine', 'poisson', 'random']
methods = ['spearman', 'spearman_nonzero', 'cosine', 'poisson', 'random']
dbs = ['means', 'medians']
ticks = ['x', 'o', '*', '+', 's']
ranks = 21
for method, tick in zip(methods, ticks):
    for db in dbs:
        cell_label_results = {}
        for lab in sorted(list(set(labels))):
            results = mouse_cell_query.search_db(np.array(data[:, labels==lab].mean(1)).flatten(), genes, method=method, db=db)
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


# try a different dataset... Zeisel dataset?
path = '/home/yjzhang/Grad_School/single_cell/uncurl_test_datasets/zeisel/'
data = scipy.io.loadmat(path + 'Zeisel.mat')
