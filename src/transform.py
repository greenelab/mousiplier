"""
This file implements a class that uses the weights from a delayedPLIER run to transform gene
expression data into LV space
"""

import random
from typing import List

import numpy as np
import pandas as pd


class PlierTransform():
    def __init__(self, weight_file: str, lambda_file: str, debug: bool = False):
        """
        Load the PLIER weights into a numpy array

        Arguments
        ---------
        weight_file: The path to the file output by PLIER storing the weights. By default it will
                     be called Z.tsv
        lambda_file: The file containing the L2 norm used by PLIER for training
        debug: A flag that prints more information about the input data when set to True
        """
        lv_df = pd.read_csv(weight_file, sep='\t')
        with open(lambda_file) as in_file:
            self.l2 = float(in_file.readline().strip())
        self.lv_df = lv_df
        loadings = lv_df.to_numpy()
        self.loadings = loadings
        self.file = weight_file
        self.genes = list(lv_df.index)

        assert len(self.genes) == loadings.shape[0]

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
        expression: A dataframe containing the rpkm-normalized expression data where the rows are
                    samples and the columns are genes in the same format as the PLIER genes

        Returns
        -------
        plier_expression: A dataframe where the rows are samples and the columns are latent
                         variables
        """
        reordered_expression = expression[list(self.genes)]

        # Ensure there aren't multiple columns with the same gene
        dup_columns = reordered_expression.columns.duplicated()
        reordered_expression = reordered_expression.loc[:, ~dup_columns]

        # Ensure the same number of genes are present in the loadings and expression
        try:
            assert reordered_expression.shape[1] == self.loadings.shape[0]
        except AssertionError as e:
            print('Expression dims: {}'.format(reordered_expression.shape))
            print('Loading dims: {}'.format(self.loadings.shape))
            raise e

        expression_matrix = reordered_expression.values

        xTx = self.loadings.T @ self.loadings
        inv_term = np.linalg.inv(xTx + np.identity(self.loadings.shape[1]) * self.l2)
        transformed_matrix = expression_matrix @ self.loadings @ inv_term

        col_names = ['LV{}'.format(i+1) for i in range(transformed_matrix.shape[1])]

        transformed_df = pd.DataFrame(transformed_matrix, index=reordered_expression.index,
                                      columns=col_names)

        return transformed_df

    def __str__(self):
        """
        Creates a human readable string representation to work with `print`.
        """
        rep = 'PlierTransform object from {}:\n{}\n'.format(self.file, self.loadings)
        rep += 'First and last genes: {}'.format((self.genes[0], self.genes[-1]))
        return rep


if __name__ == '__main__':
    pathways = pd.read_csv('data/example_pathway_matrix.tsv', sep='\t', index_col=0)

    a = PlierTransform('output/Z.tsv', 'output/lambda.txt')
    print('Input loadings:')
    print(a)

    SAMPLE_COUNT = 1337

    ones = np.ones((SAMPLE_COUNT, len(pathways.index)))
    example_genes = pathways.index
    example_samples = ['sample{}'.format(i+1) for i in range(SAMPLE_COUNT)]
    example_expression = pd.DataFrame(ones, columns=example_genes, index=example_samples)

    print('Input expression:')
    print(example_expression)

    transformed_df = a.transform(example_expression)
    print('Transformed df:')
    print(transformed_df)
