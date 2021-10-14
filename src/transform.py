"""
This file implements a class that uses the weights from a delayedPLIER run to transform gene
expression data into LV space
"""

import random
from typing import List

import h5py
import numpy as np
import pandas as pd




class PlierTransform():
    def __init__(self, weight_file: str, genes: List[str], debug: bool=False):
        """
        Load the PLIER weights into a numpy array

        Arguments
        ---------
        weight_file: The path to the file output by PLIER storing the weights. By default it will
                     be called Z.hdf5
        genes: The genes used PLIER. These genes should be the exact list of genes in the loadings
               in the same order they are present in the loadings. The easiest way to get these
               is to parse the file produced by `3_preprocess_expression.py` (which is called
               no_scrna_tpm.tsv) if you ran the preprocessing scripts via snakemake
        debug: A flag that prints more information about the input data when set to True
        """
        with h5py.File(weight_file, 'r') as in_file:
            # TODO change count to counts when moving to recount data
            if debug:
                print('Reading loadings...')
            dataset = in_file.get('count')

            # h5py syntax for 'load the whole array from disk'
            loadings = dataset[()]
            # We want the loadings to be genes x LVs but PLIER gives LVs x genes
            loadings = loadings.T

            self.file = weight_file
            self.genes = genes

            assert len(genes) == loadings.shape[0]

            if debug:
                print('Loading values:')
                print(loadings)
                print('Loadings shape: {}'.format(loadings.shape))
                sparsity = (loadings.size - np.count_nonzero(loadings)) / loadings.size
                print('Loadings sparsity = {}'.format(sparsity))

            self.loadings = loadings

    def transform(self, expression: pd.DataFrame) -> pd.DataFrame:
        """
        Transform a samples x genes matrix into the PLIER latent space (a samples x LVs matrix)

        Arguments
        ---------
        expression: A dataframe containing the tpm-normalized expression data where the rows are
                    samples and the columns are genes in the same format as the PLIER genes

        Returns
        -------
        plier_expression: A dataframe where the rows are samples and the columns are latent
                          latent variables
        """
        reordered_expression = expression[self.genes]
        print(reordered_expression)

        # Ensure there aren't multiple columns with the same gene
        dup_columns = reordered_expression.columns.duplicated()
        reordered_expression = reordered_expression.loc[:,~dup_columns]

        # Ensure the same number of genes are present in the loadings and expression
        try:
            assert reordered_expression.shape[1] == self.loadings.shape[0]
        except AssertionError as e:
            print('Expression dims: {}'.format(reordered_expression.shape))
            print('Loading dims: {}'.format(self.loadings.shape))
            raise e

        expression_matrix = reordered_expression.values

        transformed_matrix = expression_matrix @ self.loadings

        transformed_df = pd.DataFrame(transformed_matrix, index=reordered_expression.index)

        return transformed_df

        # Convert genes?
        # Ensure data is TPM transformed?

    def __str__(self):
        """
        Creates a human readable string representation to work with `print`.
        """
        rep = 'PlierTransform object from {}:\n{}\n'.format(self.file, self.loadings)
        rep += 'First and last genes: {}'.format((self.genes[0], self.genes[-1]))
        return rep

if __name__ == '__main__':
    pathways = pd.read_csv('data/plier_pathways.tsv', sep='\t', index_col=0)
    genes = (list(pathways.index[:5485]))

    a = PlierTransform('DelayedPLIER/test_output/Z.hdf5', genes)
    print('Input loadings:')
    print(a)

    ones = np.ones((1337, 10000))
    example_genes = random.sample(list(pathways.index[:10000]), k=10000)
    example_expression = pd.DataFrame(ones, columns=example_genes)

    print('Input expression:')
    print(example_expression)

    transformed_df = a.transform(example_expression)
    print('Transformed df:')
    print(transformed_df)