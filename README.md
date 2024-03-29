Querying various single-cell "atlases"

TODO:
- find datasets, download them
- download cluster data, calculate cluster means?
- test various methods for similarity search

## Datasets

### Tabula Muris

<https://figshare.com/projects/Tabula_Muris_Transcriptomic_characterization_of_20_organs_and_tissues_from_Mus_musculus_at_single_cell_resolution/27733>

This dataset contains cell type and tissue-annotated datasets, where the cell type annotation comes from Cell Ontology.

There are two single-cell sequencing datasets:
- 10x droplet seq dataset, with 422,803 droplets, 55,656 of which passed a QC cutoff of 500 genes and 1000 UMI
- Smart-Seq2 sequencing of FACS-sorted cells, with 53,760 cells, 44,879 of which passed a QC cutoff of at least 500 genes and 50,000 reads

There are x cells and y cell types.

### Microwell-seq

<https://figshare.com/articles/MCA_DGE_Data/5435866>

This dataset contains some batch-removed cells.

## Querying

the main problem with querying is that the input data and the database might have different gene sets. right now we just subset the genes present in both the db and query. Is there a way to do this better?

ideas:
- hamming distance with binarized data from input vs db?
- rank correlation using only nonzero elements in query?
- somehow combining hamming distance with rank correlation? average?
- extracting top genes from each cell type in tabula muris
