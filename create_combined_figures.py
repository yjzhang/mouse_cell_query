# TODO: create figures that combine the data streams?

import pickle
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# datasets to use: tabula muris droplet, tabula muris facs, mca (microwell-seq)

# what to do for figures:

# figure 3b: Single cell comparison figure
# for each dataset: compare cellmarker + hypergeom, cellmesh + hypergeom, cellmesh + prob, and scquery.
# fix the number of genes: either top 50 or top 100

# 1. load datasets
tm_droplet_cell_types = pd.read_csv('results_no_chromosomes_no_cell_components_no_cell_lines/MAP_cell_types.csv')
tm_facs_cell_types = pd.read_csv('results_no_chromosomes_no_cell_components_no_cell_lines/MAP_facs_cell_types.csv')
mca_coarse_cell_types = pd.read_csv('mca_coarse_MAP_cell_types.csv')

# 2. add scquery results

# load droplet scquery results
scQuery_results = []
with open('scQuery_top_gene_results.csv') as f:
    for i, line in enumerate(f.readlines()):
        line_data = line.split(',')[:2]
        if i > 0:
            line_data[1] = float(line_data[1])
        scQuery_results.append(line_data)
scQuery_results = pd.DataFrame(scQuery_results[1:], columns=['cell_type', 'mean_average_precision'])
scQuery_results['n_genes'] = 50
scQuery_results['query_method'] = 'scquery'
scQuery_results['gene_method'] = 'ratio'

# load facs scquery results
facs_scQuery_results = pd.read_csv('MAP_facs_scquery_cell_types.csv')

# load mca scquery results
mca_scQuery_results = pd.read_csv('MAP_mca_scquery_cell_types.csv')

# combine scquery results with other results
tm_droplet_cell_types = tm_droplet_cell_types.append(scQuery_results)
tm_facs_cell_types = tm_facs_cell_types.append(facs_scQuery_results)
mca_coarse_cell_types = mca_coarse_cell_types.append(mca_scQuery_results)

# plot top-3 accuracy for ratio, 50 genes
# cellmarker, cellmesh, prob, scquery
tm_droplet_cell_types['top_3_accuracy'] = [1 if x > 0.3 else 0 for x in tm_droplet_cell_types['mean_average_precision']]
tm_droplet_cell_types['top_1_accuracy'] = [1 if x > 0.6 else 0 for x in tm_droplet_cell_types['mean_average_precision']]
tm_facs_cell_types['top_3_accuracy'] = [1 if x > 0.3 else 0 for x in tm_facs_cell_types['mean_average_precision']]
tm_facs_cell_types['top_1_accuracy'] = [1 if x > 0.6 else 0 for x in tm_facs_cell_types['mean_average_precision']]
mca_coarse_cell_types['top_3_accuracy'] = [1 if x > 0.3 else 0 for x in mca_coarse_cell_types['mean_average_precision']]
mca_coarse_cell_types['top_1_accuracy'] = [1 if x > 0.6 else 0 for x in mca_coarse_cell_types['mean_average_precision']]

# plot a grouped bar plot...
tm_droplet_cell_types_means = tm_droplet_cell_types.groupby(['query_method', 'gene_method', 'n_genes']).mean().reset_index()
tm_facs_cell_types_means = tm_facs_cell_types.groupby(['query_method', 'gene_method', 'n_genes']).mean().reset_index()
mca_coarse_cell_types_means = mca_coarse_cell_types.groupby(['query_method', 'gene_method', 'n_genes']).mean().reset_index()

tm_droplet_cell_types_means['dataset'] = 'tm_droplet'
tm_facs_cell_types_means['dataset'] = 'tm_facs'
mca_coarse_cell_types_means['dataset'] = 'mca'

all_datasets_cell_types_means = pd.concat([tm_droplet_cell_types_means, tm_facs_cell_types_means, mca_coarse_cell_types_means], join='inner')
all_datasets_cell_types_means.replace('cellmarker', 'cellmarker_hypergeom', inplace=True)
all_datasets_cell_types_means.replace('cellmesh', 'cellmesh_hypergeom', inplace=True)
all_datasets_cell_types_means.replace('prob', 'cellmesh_prob', inplace=True)
included_methods = ['cellmarker_hypergeom', 'cellmesh_hypergeom', 'cellmesh_prob', 'scquery']
all_datasets_cell_types_means = all_datasets_cell_types_means[all_datasets_cell_types_means.query_method.isin(included_methods)]


sns.set(style='whitegrid', font_scale=1.7)
fig, ax = plt.subplots(figsize=(18, 10))
g = sns.categorical.barplot(x='query_method', y='top_3_accuracy', hue='dataset',
        data=all_datasets_cell_types_means[(all_datasets_cell_types_means.gene_method=='ratio') & (all_datasets_cell_types_means.n_genes==50)],
        ax=ax)
plt.ylim(0, 1.0)
plt.title('Cell Type Annotation Top-3 Accuracy')
plt.savefig('top_3_accuracy_all_datasets_by_query_method.png', dpi=100)

sns.set(style='whitegrid', font_scale=1.7)
fig, ax = plt.subplots(figsize=(18, 10))
g = sns.categorical.barplot(x='query_method', y='top_1_accuracy', hue='dataset',
        data=all_datasets_cell_types_means[(all_datasets_cell_types_means.gene_method=='ratio') & (all_datasets_cell_types_means.n_genes==50)],
        ax=ax)
plt.ylim(0, 1.0)
plt.title('Cell Type Annotation Top-1 Accuracy')
plt.savefig('top_1_accuracy_all_datasets_by_query_method.png', dpi=100)

fig, ax = plt.subplots(figsize=(18, 10))
g = sns.categorical.barplot(x='dataset', y='top_3_accuracy', hue='query_method',
        data=all_datasets_cell_types_means[(all_datasets_cell_types_means.gene_method=='ratio') & (all_datasets_cell_types_means.n_genes==50)],
        ax=ax)
plt.ylim(0, 1.0)
plt.title('Cell Type Annotation Top-3 Accuracy')
plt.savefig('top_3_accuracy_all_datasets_by_dataset.png', dpi=100)

fig, ax = plt.subplots(figsize=(18, 10))
g = sns.categorical.barplot(x='dataset', y='top_1_accuracy', hue='query_method',
        data=all_datasets_cell_types_means[(all_datasets_cell_types_means.gene_method=='ratio') & (all_datasets_cell_types_means.n_genes==50)],
        ax=ax)
plt.ylim(0, 1.0)
plt.title('Cell Type Annotation Top-1 Accuracy')
plt.savefig('top_1_accuracy_all_datasets_by_dataset.png', dpi=100)


# create a line graph showing how performance changes with num genes

plt.cla()
sns.relplot(x='n_genes', y='top_3_accuracy', hue='query_method', style='query_method', col='dataset', kind='line',
        markers=True,
        data=all_datasets_cell_types_means[(all_datasets_cell_types_means.gene_method=='ratio')],
        linewidth=2)
plt.ylim(0, 1.0)
#plt.title('Cell Type Annotation Top-3 Accuracy')
plt.savefig('top_3_accuracy_droplet_lineplot.png', dpi=100)

plt.cla()
sns.relplot(x='n_genes', y='top_1_accuracy', hue='query_method', style='query_method', col='dataset', kind='line',
        markers=True,
        data=all_datasets_cell_types_means[(all_datasets_cell_types_means.gene_method=='ratio')],
        linewidth=2)
plt.ylim(0, 1.0)
#plt.title('Cell Type Annotation Top-3 Accuracy')
plt.savefig('top_1_accuracy_droplet_lineplot.png', dpi=100)





# figure 4: Detailed analysis pipeline for UNCURL-App on tabula muris droplet
# plot a heatmap here

