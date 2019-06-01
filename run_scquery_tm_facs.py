# TODO: use scquery_selenium to get cell type names

# 1. get top genes
import pickle
import numpy as np

genes = np.loadtxt('tabula_muris_facs_genes.txt', dtype=str)
with open('tm_facs_t_scores.pkl', 'rb') as f:
    scores_t = pickle.load(f)


# for all cell types...
import scquery_selenium
cell_types = scores_t.keys()
cell_type_genes = {}
cell_type_selenium_results = {}
for cell_type in cell_types:
    top_50_genes = [genes[x[0]] for x in scores_t[cell_type][:50]]
    cell_type_genes[cell_type] = top_50_genes

print('number of cell types:', len(cell_type_genes))
for cell_type, gene_set in cell_type_genes.items():
    cell_type_selenium_results[cell_type] = scquery_selenium.scquery_gene_set(gene_set)
    print(cell_type, cell_type_selenium_results[cell_type])

for cell_type, gene_set in cell_type_genes.items():
    if cell_type not in cell_type_selenium_results:
        cell_type_selenium_results[cell_type] = scquery_selenium.scquery_gene_set(gene_set)
        print(cell_type, cell_type_selenium_results[cell_type])


with open('tm_facs_scquery_results.pkl', 'wb') as f:
    pickle.dump(cell_type_selenium_results, f)

###########################################################################################
# TODO: calculate accuracies

with open('tm_facs_scquery_results.pkl', 'rb') as f:
    cell_type_selenium_results = pickle.load(f)

label_cell_types = {key: [v[0] for v in values] for key, values in cell_type_selenium_results.items()}
cell_types_map = {}
cell_types_alternate_map = {}
with open('cell_ontology_to_cellmesh_tabula_muris.tsv') as f:
    l0 = f.readline()
    for line in f.readlines():
        line_data = line.split('\t')
        name = line_data[1]
        primary_cellmesh_name = line_data[2]
        alternate_cellmesh_names = line_data[3:]
        use_cellmesh = line_data[0]
        if use_cellmesh != 'n':
            cell_types_map[name] = primary_cellmesh_name
            cell_types_alternate_map[name] = alternate_cellmesh_names

with open('tm_cell_onto_alternate_names.tsv') as f:
    l0 = f.readline()
    for line in f.readlines():
        line_data = line.split('\t')
        name = line_data[1]
        alternate_cell_type_names = line_data[2:]
        if name in cell_types_alternate_map:
            cell_types_alternate_map[name].extend(alternate_cell_type_names)

# append names to alternate_map
with open('tm_to_scquery.tsv') as f:
    l0 = f.readline()
    for line in f.readlines():
        line_data = line.split('\t')
        name = line_data[0].strip()
        if len(line_data) > 1 and name in cell_types_alternate_map:
            alternate_name = line_data[1].strip()
            cell_types_alternate_map[name].append(alternate_name)

# calculate accuracy
label_accuracies = {}
label_extended_accuracies = {}
for key, value in label_cell_types.items():
    name = key
    accuracies = []
    extended_accuracies = []
    if name not in cell_types_map:
        continue
    for v in value:
        if v == cell_types_map[name] or v.lower() in name.lower() or name.lower() in v.lower():
            accuracies.append(1)
            extended_accuracies.append(1)
        else:
            accuracies.append(0)
            if v in cell_types_alternate_map[name]:
                extended_accuracies.append(0.5)
            else:
                extended_accuracies.append(0)
    label_accuracies[key] = accuracies
    label_extended_accuracies[key] = extended_accuracies

def calculate_precision_recall_curve(accuracies, n_relevant=1, use_extended=False):
    """
    https://nlp.stanford.edu/IR-book/html/htmledition/evaluation-of-unranked-retrieval-sets-1.html
    https://nlp.stanford.edu/IR-book/html/htmledition/evaluation-of-ranked-retrieval-results-1.html
    """
    precision = []
    recall = []
    mean_avg_precision = 0
    num_correct = 0
    for i, a in enumerate(accuracies):
        if a == 1:
            num_correct += 1
        if use_extended and a > 0:
            num_correct += 1
        precision.append(float(min(num_correct, n_relevant))/(i+1))
        recall.append(float(min(num_correct, n_relevant))/n_relevant)
    prev_recall = 0
    for p, r in zip(precision, recall):
        mean_avg_precision += (r - prev_recall)*p
        prev_recall = r
    return precision, recall, mean_avg_precision

# calculate precision-recall curves for all keys
label_map = {}
for key, val in label_extended_accuracies.items():
    label_map[key] = calculate_precision_recall_curve(val, use_extended=True)[2]

# take mean over cell types
# average MAP for all methods, over all datasets
map_mean = np.mean([v for v in label_map.values()])

# convert map_method_means to a pandas dict
import pandas as pd

map_cell_types_list = [(key, 50, 'ratio', 'scquery', v) for key, v in label_map.items()]
map_cell_types_list.sort()
map_cell_types = pd.DataFrame(map_cell_types_list, columns=['cell_type', 'n_genes',  'gene_method', 'query_method', 'mean_average_precision'])

map_cell_types.to_csv('MAP_facs_scquery_cell_types.csv', index=None)



