# src/

This folder contains the files used to run [on-disk PLIER](https://github.com/wgmao/DelayedPLIER).

## Pipeline
These files download data, preprocess it, and use it to train PLIER.
They are numbered in the order which they should be run.
If two files start with the same number, they can be run at the same time.

| File           | Description |
| -------------- | ----------- |
| 0_download_recount3.R  | Downloads all mouse samples from the recount3 compendium |
| 1_get_gene_lengths.R | Downloads the length of the genes present in the recount3 data for use in TPM normalizing the data |
| 1a_metadata_to_tsv.R | Converts the metadata from recount3 into a tsv for ease of use in python |
| 1b_remove_scrnaseq.py | Removes the single-cell RNAseq data from the dataset |
| 2_create_pathway_graph.py | Parses pathway files from the reactome database and converts them into a format usable by PLIER |
| 3_preprocess_expression.py | TPM normalizes, variance filters, and otherwise makes the recount expression data more manageable for PLIER |
| 4_convert_to_hdf5.R | On-disk PLIER expects the expression to live in an hdf5 file. This script converts the preprocessed tsv file and stores its data in an hdf5 file |
| 5_calculate_pcs.py | Calculates an initialization for PLIER using incremental PCA |
| 6_run_delayed_plier.R | Runs PLIER on the expression data |

## Libraries
These files contain useful functions used in the pipeline

| File           | Description |
| -------------- | ----------- |
| delayed_plier.R | Stores the functions used to run on-disk PLIER |
| transform.py | Contains a wrapper class that makes the output of PLIER easier to apply to expression data in Python |
| utils.py | Contains utility functions useful in the pipeline |