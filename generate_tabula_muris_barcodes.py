# TODO: preprocess the tabula muris data into a format that can be queried upon.
# data source: https://figshare.com/articles/Single-cell_RNA-seq_data_from_microfluidic_emulsion_v2_/5968960
# tabula_muris/droplet.zip

from collections import defaultdict

import numpy as np
import pandas as pd

metadata_droplet = pd.read_csv('tabula_muris/metadata_droplet.csv')
annotations_droplet = pd.read_csv('tabula_muris/annotations_droplet.csv')

# save cell barcodes for each cell type
cell_type_ids = defaultdict(lambda: [])
normalize = False
for index, row in annotations_droplet.iterrows():
    cell_id = row['cell']
    cell_type = row['cell_ontology_class']
    cell_type_ids[cell_type].append(cell_id)

for cell_type, ids in cell_type_ids.items():
    # this is of shape cells x genes
    np.savetxt('tabula_muris/cell_type_matrices_tm_droplet/{0}_barcodes.txt'.format(cell_type), np.array(cell_type_ids[cell_type]), fmt='%s')

########################################### facs data

annotations_droplet = pd.read_csv('tabula_muris/annotations_facs.csv')
# save cell barcodes for each cell type
cell_type_ids = defaultdict(lambda: [])
normalize = False
for index, row in annotations_droplet.iterrows():
    cell_id = row['cell']
    cell_type = row['cell_ontology_class']
    cell_type_ids[cell_type].append(cell_id)

for cell_type, ids in cell_type_ids.items():
    # this is of shape cells x genes
    np.savetxt('tabula_muris/cell_type_matrices_tm_facs/{0}_barcodes.txt'.format(cell_type), np.array(cell_type_ids[cell_type]), fmt='%s')


