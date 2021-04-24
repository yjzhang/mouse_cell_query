# TODO: quantify Seurat and singleR query accuracies for tabula muris

# 1. load files

def load_file_data(filename):
    results = {}
    with open(filename) as f:
        for l in f.readlines():
            line = [x.strip() for x in l.split('\t')]
            if len(line) <= 1:
                continue
            truth = line[0]
            results[truth] = line[1:]
    return results

def calculate_accuracies(cell_types_map, cell_types_alternate_map, results):
    """
    Returns two dicts, both dicts of cell type name to list of 1 or 0 indicating whether the label is correct.
    """
    label_accuracies = {}
    label_extended_accuracies = {}
    for key, values in results.items():
        name = key
        accuracies = []
        extended_accuracies = []
        if name not in cell_types_map:
            continue
        for v in values:
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
    return label_accuracies, label_extended_accuracies

def calc_topk_accuracy(label_accuracies, k=1):
    """
    Returns a single number, indicating the top 1 accuracy over all cell types.
    """
    acc_count = 0.0
    total_count = len(label_accuracies)
    for key, vals in label_accuracies.items():
        if sum(vals[:k]) >= 1:
            acc_count += 1
    return acc_count/total_count


droplet_to_facs_seurat = load_file_data('droplet_to_facs_predictions.txt')
droplet_to_facs_singler = load_file_data('droplet_to_facs_predictions_singleR.txt')
droplet_to_ref_singler = load_file_data('droplet_to_bulk_ref_predictions_singleR.txt')

# 2. load ground truth refs
cell_types_map = {}
cell_types_alternate_map = {}
with open('cell_ontology_to_cellmesh_tabula_muris.tsv') as f:
    l0 = f.readline()
    for line in f.readlines():
        line_data = [x.strip() for x in line.split('\t')]
        name = line_data[1]
        primary_cellmesh_name = line_data[2]
        alternate_cellmesh_names = line_data[3:]
        use_cellmesh = line_data[0]
        if use_cellmesh != 'n':
            cell_types_map[name] = primary_cellmesh_name
            cell_types_alternate_map[name] = alternate_cellmesh_names

with open('tabula_muris_to_panglao.tsv') as f:
    l0 = f.readline()
    for line in f.readlines():
        line_data = [x.strip() for x in line.split('\t')]
        name = line_data[1]
        primary_panglao_name = line_data[2]
        alternate_panglao_names = line_data[3:]
        use_cellmesh = line_data[0]
        if use_cellmesh != 'n':
            cell_types_alternate_map[name].extend(alternate_panglao_names)

with open('tm_cell_onto_alternate_names.tsv') as f:
    l0 = f.readline()
    for line in f.readlines():
        line_data = [x.strip() for x in line.split('\t')]
        name = line_data[1]
        alternate_cell_type_names = line_data[2:]
        if name in cell_types_alternate_map:
            cell_types_alternate_map[name].extend(alternate_cell_type_names)

# 2. calculate accuracy for each file
acc_seurat = calculate_accuracies(cell_types_map, cell_types_alternate_map, droplet_to_facs_seurat)
acc_singler = calculate_accuracies(cell_types_map, cell_types_alternate_map, droplet_to_facs_singler)
acc_singler_bulk = calculate_accuracies(cell_types_map, cell_types_alternate_map, droplet_to_ref_singler)

seurat_top1 = calc_topk_accuracy(acc_seurat[0])
singler_top1 = calc_topk_accuracy(acc_singler[0])
singler_bulk_top1 = calc_topk_accuracy(acc_singler_bulk[0])

seurat_top3 = calc_topk_accuracy(acc_seurat[0], 3)
singler_top3 = calc_topk_accuracy(acc_singler[0], 3)
singler_bulk_top3 = calc_topk_accuracy(acc_singler_bulk[0], 3)

seurat_top5 = calc_topk_accuracy(acc_seurat[0], 5)
singler_top5 = calc_topk_accuracy(acc_singler[0], 5)
singler_bulk_top5 = calc_topk_accuracy(acc_singler_bulk[0], 5)
