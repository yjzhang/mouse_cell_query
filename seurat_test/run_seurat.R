# TODO

# load MCA data

# load TM data, use as reference

library(dplyr)
library(Seurat)
library(patchwork)


# input dir for tm droplet cell types: ~/Grad_School/single_cell/mouse_cell_query/tabula_muris/cell_type_matrices
# input dir for mca cell types: /Grad_School/single_cell/mouse_cell_query/mca_microwell_seq/rmbatch_dge
# load all the files as mtx, load gene names and cell type annotations

# 1. load droplet data ######################################################
droplet_genes_path <- "tabula_muris/tabula_muris_droplet_genes.txt"
files_droplet_data <- list.files(path="tabula_muris/cell_type_matrices_tm_droplet/", pattern="*.mtx.gz", full.names=TRUE, recursive=FALSE)
droplet_matrices <- c()
droplet_seurat <- c()
# add cell types
droplet_cell_types <- c()
for (x in files_droplet_data) {
    # 1. get barcode name
    length <- nchar(x)
    base_path <- substr(x, 0, length-7)
    split_path <- unlist(strsplit(base_path, "/"))
    cell_name <- split_path[length(split_path)]
    print(base_path)
    print(cell_name)
    barcode_path <- paste(base_path, "_barcodes.txt", sep="")
    cell.data <- ReadMtx(x, barcode_path, droplet_genes_path, feature.column=1)
#    if (ncol(cell.data) < 30) {
#        next
#    }
    cell <- CreateSeuratObject(counts=cell.data)
    Idents(cell) <- base_path
#    cell <- NormalizeData(cell, verbose = FALSE)
#    cell <- FindVariableFeatures(cell, selection.method = "vst",
#        nfeatures = 2000, verbose = FALSE)
    droplet_seurat <- c(droplet_seurat, cell)
    droplet_cell_types <- c(cell_name, droplet_cell_types)
    saveRDS(cell, paste('droplet_seurat', cell_name, '.rds'))
}
droplet.data <- c()
for (i in 1:length(droplet_seurat)) {
    droplet.data[droplet_cell_types[i]] <- droplet_seurat[i]
}
droplet_merged <- merge(x=droplet_seurat[[1]], y=droplet_seurat[2:length(droplet_seurat)], add.cell.ids=droplet_cell_types)
droplet_merged <- FindVariableFeatures(droplet_merged, selection.method='vst', nfeatures=2000)

saveRDS(droplet_cell_types, 'droplet_cell_types.rds')
saveRDS(droplet.data, 'droplet_data.rds')
saveRDS(droplet_merged, 'droplet_merged.rds')

# 2. load facs data ########################################################

library(dplyr)
library(Seurat)
library(patchwork)

facs_genes_path <- "tabula_muris/tabula_muris_facs_genes.txt"
files_facs_data <- list.files(path="tabula_muris/cell_type_matrices_tm_facs/", pattern="*.mtx.gz", full.names=TRUE, recursive=FALSE)
facs_matrices <- c()
facs_seurat <- c()
# add cell types
facs_cell_types <- c()
for (x in files_facs_data) {
    # 1. get barcode name
    length <- nchar(x)
    base_path <- substr(x, 0, length-7)
    split_path <- unlist(strsplit(base_path, "/"))
    cell_name <- split_path[length(split_path)]
    print(base_path)
    print(cell_name)
    barcode_path <- paste(base_path, "_barcodes.txt", sep="")
    cell.data <- ReadMtx(x, barcode_path, facs_genes_path, feature.column=1)
#    if (ncol(cell.data) < 30) {
#        next
#    }
    cell <- CreateSeuratObject(counts=cell.data)
    Idents(cell) <- base_path
#    cell <- NormalizeData(cell, verbose = FALSE)
#    cell <- FindVariableFeatures(cell, selection.method = "vst",
#        nfeatures = 2000, verbose = FALSE)
    facs_seurat <- c(facs_seurat, cell)
    facs_cell_types <- c(facs_cell_types, cell_name)
}
#integrated_data <- c()
#for (i in 1:length(facs_seurat)) {
#    integrated_data[paste('facs', facs_cell_types[i], sep='_')] <- facs_seurat[i]
#}
#for (i in 1:length(droplet_seurat)) {
#    integrated_data[droplet_cell_types[i]] <- droplet_seurat[i]
#}
facs_cell_ids <- c()
for (i in 1:length(facs_seurat)) {
    facs_cell_ids <- c(facs_cell_ids, rep(facs_cell_types[i], ncol(facs_seurat[[i]])))
}
facs_merged <- merge(x=facs_seurat[[1]], y=facs_seurat[2:length(facs_seurat)], add.cell.ids=facs_cell_types)
#facs_merged <- FindVariableFeatures(facs_merged, selection.method='vst', nfeatures=2000)

saveRDS(facs_cell_types, 'facs_cell_types.rds')
saveRDS(facs_cell_ids, 'facs_cell_ids.rds')
saveRDS(facs_merged, 'facs_merged.rds')

# 3. Integrate facs with droplet data ########################################

library(dplyr)
library(Seurat)
library(patchwork)

# TODO: downsample cells - so that we only use like 10k cells per dataset? Or just downsample the droplet dataset?


#droplet_merged <- readRDS('droplet_merged.rds')
#facs_merged <- readRDS('facs_merged.rds')
#integrated.anchors <- FindIntegrationAnchors(object.list=c(droplet_merged, facs_merged), dims=1:20)
#data.integrated <- IntegrateData(anchorset=integrated.anchors, dims=1:20)
#saveRDS(data.integrated, 'data_integrated.rds')

# this causes memory issues - can we break up the data into multiple pieces and then do the integration?


# 4. Map integrated droplet data onto facs data #########################################

library(dplyr)
library(Seurat)
library(patchwork)

facs_cell_ids <- readRDS('facs_cell_ids.rds')
droplet_cell_types <- readRDS('droplet_cell_types.rds')
facs_merged <- readRDS('facs_merged.rds')
#droplet.data <- readRDS('droplet_data.rds')
facs_merged <- NormalizeData(facs_merged)
facs_merged <- FindVariableFeatures(facs_merged, selection.method='vst', nfeatures=2000)

all_predictions <- c()
for (cell in droplet_cell_types) {
    print(cell)
    query_data <- readRDS(paste('droplet_seurat', cell, '.rds'))
    query_data <- NormalizeData(query_data)
    result <- tryCatch({
        facs.anchors <- FindTransferAnchors(reference=facs_merged, query=query_data, dims=1:15)
        predictions <- TransferData(anchorset=facs.anchors, refdata=facs_cell_ids, dims=1:15, k.weight=14)

        p_table <- table(predictions$predicted.id) 
        p_table <- sort(p_table, decreasing=TRUE)
        p_table <- as.data.frame(p_table)
        all_predictions[[cell]] <- lapply(p_table[,'Var1'], as.character)
    },
    error=function(cond) {
        print('Error')
        print(cond)
    })
}

saveRDS(all_predictions, 'droplet_to_facs_predictions.rds')

for (cell in names(all_predictions)) {
    print(cell)
    results = all_predictions[[cell]]
    outstr <- paste(c(cell, results), collapse="\t")
    write(outstr, 'droplet_to_facs_predictions.txt', append=TRUE)
    write("\n", 'droplet_to_facs_predictions.txt', append=TRUE)
}
#droplet_seurat.anchors <- FindIntegrationAnchors(object.list = droplet_seurat, dims = 1:30);
#all_data.integrated <- IntegrateData(anchorset = all_data.anchors, dims = 1:30);
#saveRDS(all_data.integrated, 'all_data_integrated.rds')

