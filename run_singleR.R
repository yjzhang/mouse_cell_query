library(Seurat)
library(SingleR)

# 1. load facs and droplet data
facs_cell_ids <- readRDS('facs_cell_ids.rds')
droplet_cell_types <- readRDS('droplet_cell_types.rds')
facs_merged <- readRDS('facs_merged.rds')
facs_merged <- NormalizeData(facs_merged)
facs_merged <- as.SingleCellExperiment(facs_merged)
trained <- trainSingleR(facs_merged, facs_cell_ids, genes="de", sd.thresh=0.5, aggr.ref=TRUE)

# TODO: try using facs data, and also try using an existing reference.

all_predictions <- c()
for (cell in droplet_cell_types) {
    print(cell)
    query_data <- readRDS(paste('droplet_seurat', cell, '.rds'))
    query_data <- NormalizeData(query_data)
    query_data <- as.SingleCellExperiment(query_data)
    predic <- classifySingleR(test=query_data, trained=trained)
    p_table <- table(predic$labels)
    print('table:')
    print(p_table)
    data_table <- as.data.frame(p_table)
    data_table <- data_table[order(data_table$Freq),]
    print(data_table)
    all_predictions[[cell]] <- lapply(data_table[,'Var1'], as.character)
}

saveRDS(all_predictions, 'droplet_to_facs_predictions_singleR.rds')

for (cell in names(all_predictions)) {
    print(cell)
    results = rev(all_predictions[[cell]])
    outstr <- paste(c(cell, results), collapse="\t")
    write(outstr, 'droplet_to_facs_predictions_singleR.txt', append=TRUE)
    write("\n", 'droplet_to_facs_predictions_singleR.txt', append=TRUE)
}


# try prediction using an existing reference

gene_intersection <- function(data1, data2) {
    genes1 <- rownames(data1)
    genes2 <- rownames(data2)
    int <- intersect(genes1, genes2)
    m1 <- match(int, genes1)
    m2 <- match(int, genes2)
    return(list('m1'=m1, 'm2'=m2))
}

# 1. load ref from celldex
library(celldex)
ref <- MouseRNAseqData()

gene_subsets <- gene_intersection(ref, query_data)
m1 <- gene_subsets$m1
m2 <- gene_subsets$m2
ref <- ref[m1]
ref_trained <- trainSingleR(ref, ref$label.main, genes="de", sd.thresh=0.5, aggr.ref=FALSE)

ref_all_predictions <- c()
for (cell in droplet_cell_types) {
    print(cell)
    query_data <- readRDS(paste('droplet_seurat', cell, '.rds'))
    query_data <- query_data[m2]
    query_data <- NormalizeData(query_data)
    query_data <- as.SingleCellExperiment(query_data)
    predic <- classifySingleR(test=query_data, trained=ref_trained)
    p_table <- table(predic$labels)
    print('table:')
    print(p_table)
    data_table <- as.data.frame(p_table)
    data_table <- data_table[order(data_table$Freq),]
    print(data_table)
    ref_all_predictions[[cell]] <- lapply(data_table[,'Var1'], as.character)
}

saveRDS(ref_all_predictions, 'droplet_to_bulk_ref_predictions_singleR.rds')

ref_all_predictions <- readRDS('droplet_to_bulk_ref_predictions_singleR.rds')
for (cell in names(ref_all_predictions)) {
    print(cell)
    results = rev(ref_all_predictions[[cell]])
    outstr <- paste(c(cell, results), collapse="\t")
    write(outstr, 'droplet_to_bulk_ref_predictions_singleR.txt', append=TRUE)
    write("\n", 'droplet_to_bulk_ref_predictions_singleR.txt', append=TRUE)
}
