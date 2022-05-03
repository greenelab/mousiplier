"""This script adds marker genes for three brain regions to the pathway matrix """

import argparse
import os

import pandas as pd

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--marker_files',
                        nargs='*',
                        help='A series of files containing relevant genes',
                        default=['data/markers/cerebral_cortex.txt',
                                 'data/markers/midbrain.txt',
                                 'data/markers/striatum.txt'])
    parser.add_argument('--pathway_file',
                        help='The file containing a matrix mapping pathways to genes',
                        default='data/plier_pathways.tsv')
    args = parser.parse_args()

    # keep_default_na=False keeps us from clobbering the NA gene symbol
    pathway_df = pd.read_csv(args.pathway_file, sep='\t', index_col=0, keep_default_na=False)

    brain_df = None

    for file_path in args.marker_files:
        pathway = os.path.splitext(file_path)[0]
        pathway = os.path.basename(pathway)
        genes = pd.read_csv(file_path, header=None, index_col=0)
        genes[pathway] = 1.0
        genes.index.name = None

        if brain_df is None:
            brain_df = genes
        else:
            brain_df = brain_df.merge(genes, how='outer', left_index=True, right_index=True)

    original_genes = set(pathway_df.index)

    test_df = pd.DataFrame({'genes': ['fake_gene_1', 'fake_gene_2'], 'fake_pathway': [1., 1.]})
    test_df = test_df.set_index('genes')
    test2_df = pd.DataFrame({'genes': ['fake_gene_3', 'fake_gene_4'], 'fake_pathway2': [1., 1.]})
    test2_df = test2_df.set_index('genes')

    test_df.index.name = None

    pathway_df = pathway_df.merge(test_df, how='outer', left_index=True, right_index=True)
    pathway_df = pathway_df.fillna(0)
    pathway_df.index.name = None

    current_genes = set(pathway_df.index)

    pathway_df.to_csv('data/extended_plier_pathways.tsv', sep='\t')