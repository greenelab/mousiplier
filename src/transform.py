"""
This file implements a class that uses the weights from a delayedPLIER run to transform gene
expression data into LV space
"""

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
        genes: The genes used in the corresponding to the
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
                    samples and the columns are genes

        Returns
        -------
        plier_expression: A dataframe where the rows are samples and the columns are latent
                          latent variables
        """
        # Convert genes?
        # Match genes
        # Ensure data is TPM transformed?
        # Check dims
        # Multiply PLIER weights by expression values
        # Return result

    def __str__(self):
        """
        Creates a human readable string representation to work with `print`.
        """
        rep = 'PlierTransform object from {}:\n{}\n'.format(self.file, self.loadings)
        rep += 'First and last genes: {}'.format((self.genes[0], self.genes[-1]))
        return rep

pathways = pd.read_csv('data/plier_pathways.tsv', sep='\t', index_col=0)
genes = (list(pathways.index[:5485]))

a = PlierTransform('DelayedPLIER/test_output/Z.hdf5', genes, debug=True)
print(a)

ones = np.ones((1337, 5485))
example_expression = pd.DataFrame(np.ones, index=genes)